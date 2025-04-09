using CFTime
using Dates
using Downloads

using Oceananigans.DistributedComputations

using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: Metadata, metadata_path, download_progress, AnyDateTime

import Dates: year, month, day
import Oceananigans.Fields: set!
import Base

import Oceananigans.Fields: set!, location
import ClimaOcean.DataWrangling: all_dates, metadata_filename, download_dataset, default_download_directory

struct MultiYearIFS end
struct HourlyIFS end

const IFSMetadata{D} = Metadata{<:Union{<:MultiYearIFS, <:HourlyIFS}, D}
const IFSMetadatum   = Metadatum{<:Union{<:MultiYearIFS, <:HourlyIFS}}

# ===============================================================================
# STEVE:HARDCODED
const varname::String = "seaSurfaceTemperature"
const YYYYMMDD::String = "20250318"
const HH::String = "00"
const year::Int = 2025
const month::Int = 3
const day::Int = 18
const hour::Int = 0
const datadir::String = "/fsx/climaocean/data"
const model::String = "ECMWFHRes"
const forecast_hour = nothing
const nlon::Int = 1440
const nlat::Int = 721
const forecast_days::Int = 3
# ===============================================================================

default_download_directory(::Union{<:MultiYearIFS, <:HourlyIFS}) = download_IFS_cache

# STEVE:HARDCODED
Base.size(data::IFSMetadata) = (nlon, nlat, length(data.dates))
Base.size(::IFSMetadatum)    = (nlon, nlat, 1)

# IFS is a spatially 2D dataset
variable_is_three_dimensional(data::IFSMetadata) = false

# The whole range of dates in the different dataset datasets
# NOTE! rivers and icebergs have a different frequency! (typical IFS data is three-hourly while rivers and icebergs are daily)
function all_dates(::HourlyIFS, name)
    if name == :river_freshwater_flux || name == :iceberg_freshwater_flux
        return DateTime(1990, 1, 1) : Day(1) : DateTime(1990, 12, 31)
    else
        return DateTime(year, month, day, hour, 0, 0) : Hour(1) : DateTime(year, month, day+forecast_days, hour, 0, 0)
	# Hourly 1-90 (3.75 days)
	# 3-hourly 93-145
	# 6-hourly 145-240
    end
end

# Valid for all IFS datasets
function IFS_time_indices(dataset, dates, name)
    all_IFS_dates = all_dates(dataset, name)
    indices = Int[]

    for date in dates
        index = findfirst(x -> x == date, all_IFS_dates)
        !isnothing(index) && push!(indices, index)
    end

    return indices
end

# File name generation specific to each Dataset dataset
# Note that `HourlyIFS` has only one file associated, so we can define
# the filename directly for the whole `Metadata` object, independent of the `dates`
function metadata_filename(metadata::Metadata{<:HourlyIFS}) # No difference 

#   shortname::String = "seaSurfaceTemperature",
#   YYYYMMDD::String = "20250318",
#   HH::String = "00",
#   datadir::String = "/fsx/climaocean/data",
#   model::String = "ECMWFHRes",
#   forecast_hour = nothing

    shortname = short_name(metadata)
 
    infile_prefix="$datadir/$model.$YYYYMMDD.$HH.full_forecast"
    infile = "$infile_prefix.$shortname.nc"

    return infile
end


# Convenience functions
short_name(data::IFSMetadata) = IFS_short_names[data.name]
location(::IFSMetadata) = (Center, Center, Center)

# A list of all variables provided in the IFS dataset:
IFS_variable_names = (:river_freshwater_flux,
                        :rain_freshwater_flux,
                        :snow_freshwater_flux,
                        :iceberg_freshwater_flux,
                        :specific_humidity,
                        :sea_level_pressure,
                        :downwelling_longwave_radiation,
                        :downwelling_shortwave_radiation,
                        :temperature,
                        :eastward_velocity,
                        :northward_velocity)


# var_array=( 'airTemperature' 'meanSeaLevelPressure' 'precipitationRate' 'windVelocity10MeterEastward' 'windVelocity10MeterNorthward' )
IFS_short_names = Dict(
    :river_freshwater_flux           => "friver",   # Freshwater fluxes from rivers
    :rain_freshwater_flux            => "prra",     # Freshwater flux from rainfall
    :snow_freshwater_flux            => "prsn",     # Freshwater flux from snowfall
    :iceberg_freshwater_flux         => "licalvf",  # Freshwater flux from calving icebergs
    :specific_humidity               => "huss",     # Surface specific humidity
    :sea_level_pressure              => "meanSeaLevelPressure",      # Sea level pressure
    :downwelling_longwave_radiation  => "rlds",     # Downwelling longwave radiation
    :downwelling_shortwave_radiation => "rsds",     # Downwelling shortwave radiation
    :temperature                     => "airTemperature",      # Near-surface air temperature
    :eastward_velocity               => "windVelocity10MeterEastward",      # Eastward near-surface wind
    :northward_velocity              => "windVelocity10MeterNorthward",      # Northward near-surface wind
)


IFS_repeat_year_urls = Dict(
    :shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                            "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :river_freshwater_flux => "https://www.dropbox.com/scl/fi/21ggl4p74k4zvbf04nb67/" *
                              "RYF.friver.1990_1991.nc?rlkey=ny2qcjkk1cfijmwyqxsfm68fz&dl=0",

    :rain_freshwater_flux => "https://www.dropbox.com/scl/fi/5icl1gbd7f5hvyn656kjq/" *
                             "RYF.prra.1990_1991.nc?rlkey=iifyjm4ppwyd8ztcek4dtx0k8&dl=0",

    :snow_freshwater_flux => "https://www.dropbox.com/scl/fi/1r4ajjzb3643z93ads4x4/" *
                             "RYF.prsn.1990_1991.nc?rlkey=auyqpwn060cvy4w01a2yskfah&dl=0",

    :iceberg_freshwater_flux => "https://www.dropbox.com/scl/fi/44nc5y27ohvif7lkvpyv0/" *
                                "RYF.licalvf.1990_1991.nc?rlkey=w7rqu48y2baw1efmgrnmym0jk&dl=0",

    :specific_humidity => "https://www.dropbox.com/scl/fi/66z6ymfr4ghkynizydc29/" *
                          "RYF.huss.1990_1991.nc?rlkey=107yq04aew8lrmfyorj68v4td&dl=0",

    :sea_level_pressure => nothing,

    :downwelling_longwave_radiation  => "https://www.dropbox.com/scl/fi/y6r62szkirrivua5nqq61/" *
                                        "RYF.rlds.1990_1991.nc?rlkey=wt9yq3cyrvs2rbowoirf4nkum&dl=0",

    :downwelling_shortwave_radiation => "https://www.dropbox.com/scl/fi/z6fkvmd9oe3ycmaxta131/" *
                                        "RYF.rsds.1990_1991.nc?rlkey=r7q6zcbj6a4fxsq0f8th7c4tc&dl=0",

    :temperature => nothing,

    :eastward_velocity => nothing,

    :northward_velocity => nothing,
)

metadata_url(metadata::Metadata{<:HourlyIFS}) = IFS_repeat_year_urls[metadata.name]  

function download_dataset(metadata::IFSMetadata)

    @root for metadatum in metadata

        fileurl  = metadata_url(metadatum)

	#STEVE: if the file url is 'nothing', then we're going to use the ECMWF IFS HRES forecast instead for that field
        if fileurl == nothing
            #STEVE: assume we've already downloaded the data to the desired storage location
	    continue
	end

        filepath = metadata_path(metadatum)

        if !isfile(filepath)
            Downloads.download(fileurl, filepath; progress=download_progress)
        end
    end
 

    return nothing
end
