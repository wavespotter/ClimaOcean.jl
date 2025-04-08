using CFTime
using Dates
using ClimaOcean.DataWrangling
using ClimaOcean.DataWrangling: netrc_downloader, metadata_path, AnyDateTime
using Oceananigans.DistributedComputations

import Dates: year, month, day
using Downloads

import Oceananigans.Fields: set!, location
import Base
import ClimaOcean.DataWrangling: all_dates, metadata_filename, download_dataset, default_download_directory

struct MOMonthly end
struct MODaily end

const MOMetadata{D} = Metadata{<:Union{<:MOMonthly, <:MODaily}, D}
const MOMetadatum   = Metadatum{<:Union{<:MOMonthly, <:MODaily}}

"""
    MOMetadatum(name; 
                  date = first_date(MO4Monthly()), 
                  dir = download_MO_cache)

an alias to construct a [`Metadatum`](@ref) of [`MO4Montly`](@ref)
"""
function MOMetadatum(name; 
                       date = first_date(MO4Monthly()), 
                       dir = download_MO_cache)
  
    return Metadatum(name; date, dir, dataset=MO4Monthly())
end

default_download_directory(::Union{<:MOMonthly, <:MODaily}) = download_MO_cache

datasetstr(md::MOMetadata) = string(md.dataset)

datestr(md::MOMetadata) = string(first(md.dates), "--", last(md.dates))
datestr(md::MOMetadatum) = string(md.dates)

Base.summary(md::MOMetadata) = string("MOMetadata{", datasetstr(md), "} of ",
                                        md.name, " for ", datestr(md))

Base.size(data::Metadata{<:MODaily})   = (4320, 2041, 50, length(data.dates))
Base.size(data::Metadata{<:MOMonthly}) = (4320, 2041, 50, length(data.dates))

Base.size(::Metadatum{<:MODaily})   = (4320, 2041, 50, 1)
Base.size(::Metadatum{<:MOMonthly}) = (4320, 2041, 50, 1)

# The whole range of dates in the different dataset datasets
all_dates(::MOMonthly, name) = DateTime(2023, 5, 1) : Month(1) : DateTime(2025, 5, 1)
all_dates(::MODaily, name)   = DateTime(2023, 5, 1) : Day(1)   : DateTime(2025, 5, 1)

# Fallback, actually, we do not really need the name for MO since all
# variables have the same frequency and the same time-range, differently from JRA55
all_dates(dataset::Union{<:MOMonthly, <:MODaily}) = all_dates(dataset, :temperature)

# File name generation specific to each Dataset dataset
function metadata_filename(metadata::Metadatum{<:MOMonthly})
    shortname = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    return shortname * "_" * yearstr * "_" * monthstr * ".nc"
end

function metadata_filename(metadata::Metadatum{<:Union{MO2Daily, MO2Monthly}})
    shortname   = short_name(metadata)
    yearstr  = string(Dates.year(metadata.dates))
    monthstr = string(Dates.month(metadata.dates), pad=2)
    postfix = variable_is_three_dimensional(metadata) ? ".4320x2041x50." : ".4320x2041."

    if metadata.dataset isa MOMonthly
        return shortname * postfix * yearstr * monthstr * ".nc"
    elseif metadata.dataset isa MODaily
        daystr = string(Dates.day(metadata.dates), pad=2)
        return shortname * postfix * yearstr * monthstr * daystr * ".nc"
    end
end

# Convenience functions
short_name(data::Metadata{<:MODaily})   = MO_short_names[data.name]
short_name(data::Metadata{<:MOMonthly}) = MO_short_names[data.name]

location(data::MOMetadata) = MO_location[data.name]

variable_is_three_dimensional(data::MOMetadata) =
    data.name == :temperature ||
    data.name == :salinity ||
    data.name == :u_velocity ||
    data.name == :v_velocity

#MO4_short_names = Dict(
#    :temperature           => "THETA",
#    :salinity              => "SALT",
#    :u_velocity            => "EVEL",
#    :v_velocity            => "NVEL",
#    :free_surface          => "SSH",
#    :sea_ice_thickness     => "SIheff",
#    :sea_ice_concentration => "SIarea",
#    :net_heat_flux         => "oceQnet"
#)

MO_short_names = Dict(
    :temperature           => "THETA",
    :salinity              => "SALT",
    :u_velocity            => "UVEL",
    :v_velocity            => "VVEL",
    :free_surface          => "SSH",
    :sea_ice_thickness     => "SIheff",
    :sea_ice_concentration => "SIarea",
    :net_heat_flux         => "oceQnet"
)

MO_location = Dict(
    :temperature           => (Center, Center, Center),
    :salinity              => (Center, Center, Center),
    :free_surface          => (Center, Center, Nothing),
    :sea_ice_thickness     => (Center, Center, Nothing),
    :sea_ice_concentration => (Center, Center, Nothing),
    :net_heat_flux         => (Center, Center, Nothing),
    :u_velocity            => (Face,   Center, Center),
    :v_velocity            => (Center, Face,   Center),
)


