include("./NetcdfReader.jl")
using .NetcdfReader

varname="seaSurfaceTemperature"
YYYYMMDD="20250318"
HH="00"
datadir="/fsx/climaocean/data"
model="ECMWFHRes"
forecast_hour="000"

infile_prefix="$datadir/$model.$YYYYMMDD.$HH.f$forecast_hour"
infile="$infile_prefix.$varname.nc"
ds = read_netcdf(infile, varname)

var_array = [ "airTemperature", "airDensity", "meanSeaLevelPressure", "precipitationRate", "seaSurfaceTemperature", "windVelocity10MeterEastward", "windVelocity10MeterNorthward" ]

for (i, varname) in enumerate(var_array)

    println("Loading variable: $varname...")
    local infile = "$infile_prefix.$varname.nc"
    local ds = read_netcdf(infile, varname)

    ds
end
