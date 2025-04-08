module NetcdfReader

export read_netcdf
using NCDatasets

function read_netcdf(infile, varname)

    println("read_netcdf:: reading infile = $infile using varname = $varname")

    # Show listing of the netcdf file
#   NCDataset(infile,"r")

    ds = NCDataset(infile,"r")
    display(ds)
#   var_ds = ds[varname]

    # load a subset
    #subdata = v[10:30,30:5:end]

    # load all data
#   data = var_ds[:,:]

    # load all data ignoring attributes like scale_factor, add_offset, _FillValue and time units
    #data2 = v.var[:,:];

    # load an attribute
    #unit = v.attrib["units"]
    #
    #
    close(ds)

    return ds

end


end
