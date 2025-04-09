module MO

export MOMetadatum, MO_field, MO_mask, MO_immersed_grid, adjusted_MO_tracers, initialize!
export MOMonthly, MODaily
export MOFieldTimeSeries, MORestoring, LinearlyTaperedPolarMask

using ClimaOcean
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: inpaint_mask!, NearestNeighborInpainting, download_progress, compute_native_date_range
using ClimaOcean.InitialConditions: three_dimensional_regrid!, interpolate!

using Oceananigans
using Oceananigans: location
using Oceananigans.Architectures: architecture, child_architecture
using Oceananigans.BoundaryConditions
using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: DistributedField, all_reduce, barrier!
using Oceananigans.Utils

using KernelAbstractions: @kernel, @index
using NCDatasets
using JLD2
using Downloads: download
using Dates
using Adapt
using Scratch

download_MO_cache::String = ""
function __init__()
    global download_MO_cache = @get_scratch!("MO")
end

include("MO_metadata.jl")
include("MO_mask.jl")

# Vertical coordinate
#const ECCO_z = [
#    -6128.75,
#    -5683.75,
#    -5250.25,
#    -4839.75,
#    -4452.25,
#    -4087.75,
#    -3746.25,
#    -3427.75,
#    -3132.25,
#    -2859.75,
#    -2610.25,
#    -2383.74,
#    -2180.13,
#    -1999.09,
#    -1839.64,
#    -1699.66,
#    -1575.64,
#    -1463.12,
#    -1357.68,
#    -1255.87,
#    -1155.72,
#    -1056.53,
#    -958.45,
#    -862.10,
#    -768.43,
#    -678.57,
#    -593.72,
#    -515.09,
#    -443.70,
#    -380.30,
#    -325.30,
#    -278.70,
#    -240.09,
#    -208.72,
#    -183.57,
#    -163.43,
#    -147.11,
#    -133.45,
#    -121.51,
#    -110.59,
#    -100.20,
#    -90.06,
#    -80.01,
#    -70.0,
#    -60.0,
#    -50.0,
#    -40.0,
#    -30.0,
#    -20.0,
#    -10.0,
#      0.0,
#]
const MO_z = [
    0.49,
    1.54,
    2.64,
    3.81,
    5.07,
    6.44,
    7.92   
    9.57,
    11.40, 
    13.46,
    15.81,
    18.49,
    21.59,
    25.21,   
    29.44,
    34.43,
    40.34,
    47.37,
    55.76,
    65.80,
    77.85,   
    92.32,
    109.72,
    130.66, 
    155.85,
    186.12,
    222.47,
    266.04,   
    318.12,
    380.21, 
    453.93,
    541.08,
    643.56,
    763.33,
    902.33,   
    1062.44, 
    1245.29,
    1452.25,
    1684.28,
    1941.89,
    2225.07,
    2533.33,  
    2865.70,
    3220.82, 
    3597.03,
    3992.48,
    4405.22,
    4833.29,
    5274.78,   
    5727.91,
]

empty_MO_field(variable_name::Symbol; kw...) = empty_MO_field(Metadatum(variable_name, dataset=MOMonthly()); kw...)

function empty_MO_field(metadata::MOMetadata;
                          architecture = CPU(), 
                          horizontal_halo = (7, 7))

    Nx, Ny, Nz, _ = size(metadata)
    loc = location(metadata)
    longitude = (0, 360)
    latitude = (-90, 90)
    TX, TY = (Periodic, Bounded)

    if variable_is_three_dimensional(metadata)
        TZ = Bounded
        LZ = Center
        z = MO_z
        halo = (horizontal_halo..., 3)
        sz = (Nx, Ny, Nz)
    else # the variable is two-dimensional
        TZ = Flat
        LZ = Nothing
        z = nothing
        halo = horizontal_halo
        sz = (Nx, Ny)
    end

    grid = LatitudeLongitudeGrid(architecture, Float32; halo, longitude, latitude, z,
                                 size = sz,
                                 topology = (TX, TY, TZ))

    return Field{loc...}(grid)
end

# Only temperature and salinity need a thorough inpainting because of stability,
# other variables can do with only a couple of passes. Sea ice variables 
# cannot be inpainted because zeros in the data are physical, not missing values.
function default_inpainting(metadata::MOMetadata)
    if metadata.name in [:temperature, :salinity]
        return NearestNeighborInpainting(Inf)
    elseif metadata.name in [:sea_ice_fraction, :sea_ice_thickness]
        return nothing
    else
        return NearestNeighborInpainting(5)
    end
end

"""
    MO_field(metadata::MOMetadata;
               architecture = CPU(),
               inpainting = nothing,
               mask = nothing,
               horizontal_halo = (7, 7),
               cache_inpainted_data = false)

Return a `Field` on `architecture` described by `MOMetadata` with
`horizontal_halo` size.
If not `nothing`, the `inpainting` method is used to fill the cells
within the specified `mask`. `mask` is set to `MO_mask` for non-nothing
`inpainting`.
"""
function MO_field(metadata::MOMetadata;
                    architecture = CPU(),
                    inpainting = default_inpainting(metadata),
                    mask = nothing,
                    horizontal_halo = (7, 7),
                    cache_inpainted_data = true)
                    
    field = empty_MO_field(metadata; architecture, horizontal_halo)
    inpainted_path = inpainted_metadata_path(metadata)

    if !isnothing(inpainting) && isfile(inpainted_path)
        file = jldopen(inpainted_path, "r")
        maxiter = file["inpainting_maxiter"]

        # read data if generated with the same inpainting
        if maxiter == inpainting.maxiter
            data = file["data"]
            close(file)
            copyto!(parent(field), data)
            return field
        end

        close(file)
    end

    download_dataset(metadata)
    path = metadata_path(metadata)
    ds = Dataset(path)
    shortname = short_name(metadata)

    if variable_is_three_dimensional(metadata)
        data = ds[shortname][:, :, :, 1]
        data = reverse(data, dims=3)
    else
        data = ds[shortname][:, :, 1]
    end        

    close(ds)
    
    # Convert data from Union(FT, missing} to FT
    FT = eltype(field)
    data[ismissing.(data)] .= 1e10 # Artificially large number!
    data = if location(field)[2] == Face # ?
        new_data = zeros(FT, size(field))
        new_data[:, 1:end-1, :] .= data
        new_data    
    else
        data = Array{FT}(data)
    end
    
    # MO4 data is on a -180, 180 longitude grid as opposed to MO2 data that
    # is on a 0, 360 longitude grid. To make the data consistent, we shift MO4
    # data by 180 degrees in longitude
    if metadata.dataset isa MOMonthly 
        Nx = size(data, 1)
        if variable_is_three_dimensional(metadata)
            shift = (Nx ÷ 2, 0, 0)
        else
            shift = (Nx ÷ 2, 0)
        end
        data = circshift(data, shift)
    end

    set!(field, data)
    fill_halo_regions!(field)

    if !isnothing(inpainting)
        # Respect user-supplied mask, but otherwise build default MO mask.
        if isnothing(mask)
            mask = MO_mask(metadata, architecture; data_field=field)
        end

        # Make sure all values are extended properly
        name = string(metadata.name)
        date = string(metadata.dates)
        dataset = summary(metadata.dataset)
        @info string("Inpainting ", dataset, " ", name, " data from ", date, "...")
        start_time = time_ns()
        
        inpaint_mask!(field, mask; inpainting)
        fill_halo_regions!(field)

        elapsed = 1e-9 * (time_ns() - start_time)
        @info string(" ... (", prettytime(elapsed), ")")
    
        # We cache the inpainted data to avoid recomputing it
        @root if cache_inpainted_data
            file = jldopen(inpainted_path, "w+")
            file["data"] = on_architecture(CPU(), parent(field))
            file["inpainting_maxiter"] = inpainting.maxiter
            close(file)
        end
    end

    return field
end

# Fallback
MO_field(var_name::Symbol; kw...) = MO_field(MOMetadata(var_name); kw...)

function inpainted_metadata_filename(metadata::MOMetadata)
    original_filename = metadata_filename(metadata)
    without_extension = original_filename[1:end-3]
    return without_extension * "_inpainted.jld2"
end

inpainted_metadata_path(metadata::MOMetadata) = joinpath(metadata.dir, inpainted_metadata_filename(metadata))

function set!(field::Field, MO_metadata::MOMetadatum; kw...)

    # Fields initialized from MO
    grid = field.grid
    arch = child_architecture(grid)
    mask = MO_mask(MO_metadata, arch)

    f = MO_field(MO_metadata; mask,
                   architecture = arch,
                   kw...)

    interpolate!(field, f)

    return field
end

include("MO_restoring.jl")

end # Module 

