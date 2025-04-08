module NetcdfReader

export read_netcdf
export read_IFSforecast

using NCDatasets

function read_netcdf(infile::String, varname::String)

    println("read_netcdf:: reading infile = $infile using varname = $varname")

    ds = NCDataset(infile,"r")
    #
    # Show listing of the netcdf file
    display(ds)

    # Load specific variable
    var = ds[varname]

    # load all data
    data = var[:,:,:]

    # load all data ignoring attributes like scale_factor, add_offset, _FillValue and time units
    #data2 = var.var[:,:,:];

    close(ds)

    return data

end

function read_IFSforecast(;
		          varname::String = "seaSurfaceTemperature", 
		          YYYYMMDD::String = "20250318", 
			  HH::String = "00", 
			  datadir::String = "/fsx/climaocean/data", 
			  model::String = "ECMWFHRes", 
			  forecast_hour = nothing
			  )

    println("read_IFSforecast: Loading variable: $varname...")

    if forecast_hour == nothing
        infile_prefix="$datadir/$model.$YYYYMMDD.$HH.full_forecast"
    else
        infile_prefix="$datadir/$model.$YYYYMMDD.$HH.f$forecast_hour"
    end

    infile = "$infile_prefix.$varname.nc"
    data = read_netcdf(infile, varname)

    return data

end


end
