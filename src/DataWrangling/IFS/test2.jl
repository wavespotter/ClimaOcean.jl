include("./NetcdfReader.jl")
using .NetcdfReader

varname="seaSurfaceTemperature"
YYYYMMDD="20250318"
HH="00"
datadir="/fsx/climaocean/data"
model="ECMWFHRes"
forecast_hour="000"

data = read_IFSforecast(varname=varname, YYYYMMDD=YYYYMMDD, HH=HH, datadir=datadir, model=model, forecast_hour=forecast_hour)

display(data)
