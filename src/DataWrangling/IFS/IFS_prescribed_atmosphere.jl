const AA = Oceananigans.Architectures.AbstractArchitecture

IFSPrescribedAtmosphere(arch::Distributed, FT = Float32; kw...) =
    IFSPrescribedAtmosphere(child_architecture(arch); kw...)


"""
    IFSPrescribedAtmosphere([architecture = CPU(), FT = Float32];
                              dataset = HourlyIFS(),
                              start_date = first_date(dataset, :temperature),
                              end_date = last_date(dataset, :temperature),
                              backend = IFSNetCDFBackend(10),
                              time_indexing = Cyclical(),
                              surface_layer_height = 10,  # meters
                              include_rivers_and_icebergs = false,
                              other_kw...)

Return a [`PrescribedAtmosphere`](@ref) representing IFS forecast data.
The atmospheric data will be held in `IFSFieldTimeSeries` objects containing.
For a detailed description of the keyword arguments, see the [`IFSFieldTimeSeries`](@ref) constructor.
"""
function IFSPrescribedAtmosphere(architecture = CPU(), FT = Float32;
                                   dataset = HourlyIFS(),
                                   start_date = first_date(dataset, :temperature),
                                   end_date = last_date(dataset, :temperature),
                                   backend = IFSNetCDFBackend(10),
                                   time_indexing = Cyclical(),
                                   surface_layer_height = 10,  # meters
                                   include_rivers_and_icebergs = false,
                                   other_kw...)

    kw = (; time_indexing, backend, start_date, end_date, dataset)
    kw = merge(kw, other_kw) 

    # var_array=( 'airTemperature' 'meanSeaLevelPressure' 'precipitationRate' 'windVelocity10MeterEastward' 'windVelocity10MeterNorthward' )
    # (Defined via 'short_name()')
    ua  = IFSFieldTimeSeries(:eastward_velocity, architecture, FT;               kw...)
    va  = IFSFieldTimeSeries(:northward_velocity, architecture, FT;              kw...)
    Ta  = IFSFieldTimeSeries(:temperature, architecture, FT;                     kw...)
    qa  = IFSFieldTimeSeries(:specific_humidity, architecture, FT;               kw...)
    pa  = IFSFieldTimeSeries(:sea_level_pressure, architecture, FT;              kw...)
    Fra = IFSFieldTimeSeries(:rain_freshwater_flux, architecture, FT;            kw...)
    Fsn = IFSFieldTimeSeries(:snow_freshwater_flux, architecture, FT;            kw...)
    Ql  = IFSFieldTimeSeries(:downwelling_longwave_radiation, architecture, FT;  kw...)
    Qs  = IFSFieldTimeSeries(:downwelling_shortwave_radiation, architecture, FT; kw...)

    freshwater_flux = (rain = Fra,
                       snow = Fsn)

    # Remember that rivers and icebergs are on a different grid and have
    # a different frequency than the rest of the IFS data. We use `PrescribedAtmospheres`
    # "auxiliary_freshwater_flux" feature to represent them.
    if include_rivers_and_icebergs    
        Fri = IFSFieldTimeSeries(:river_freshwater_flux, architecture;   kw...)
        Fic = IFSFieldTimeSeries(:iceberg_freshwater_flux, architecture; kw...)
        auxiliary_freshwater_flux = (rivers = Fri, icebergs = Fic)
    else
        auxiliary_freshwater_flux = nothing
    end

    times = ua.times
    grid  = ua.grid

    velocities = (u = ua,
                  v = va)

    tracers = (T = Ta,
               q = qa)

    pressure = pa

    downwelling_radiation = TwoBandDownwellingRadiation(shortwave=Qs, longwave=Ql)

    FT = eltype(ua)
    surface_layer_height = convert(FT, surface_layer_height)

    atmosphere = PrescribedAtmosphere(grid, times;
                                      velocities,
                                      freshwater_flux,
                                      auxiliary_freshwater_flux,
                                      tracers,
                                      downwelling_radiation,
                                      surface_layer_height,
                                      pressure)

    return atmosphere
end
