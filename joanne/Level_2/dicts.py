import datetime
import glob

import xarray as xr
import pandas as pd

import joanne

# %%

nc_meta = {
    #         'trajectory' : {'cf_role' : 'trajectory_id'},
    "time": {
        "standard_name": "time",
        "long_name": "Time of recorded measurement",
        "units": "seconds since 1970-01-01 00:00:00 UTC",
        "calendar": "gregorian",
        "axis": "T",
    },
    "height": {
        "standard_name": "geopotential_height",
        "long_name": "Geopotential Height obtained by integrating upwards the atmospheric thickness estimated using the hypsometric equation",
        "units": "m",
        "axis": "Z",
        "positive": "up",
    },
    "latitude": {
        "standard_name": "latitude",
        "long_name": "North Latitude",
        "units": "degrees_north",
        "axis": "X",
    },
    "longitude": {
        "standard_name": "longitude",
        "long_name": "East Longitude",
        "units": "degrees_east",
        "axis": "Y",
    },
    "pressure": {
        "standard_name": "air_pressure",
        "long_name": "Atmospheric Pressure",
        "units": "hPa",
        "coordinates": "time longitude latitude height",
    },
    "temperature": {
        "standard_name": "air_temperature",
        "long_name": "Dry Bulb Temperature",
        "units": "degree_Celsius",
        "coordinates": "time longitude latitude height",
    },
    "relative_humidity": {
        "standard_name": "relative_humidity",
        "long_name": "Relative Humidity",
        "units": "%",
        "coordinates": "time longitude latitude height",
    },
    "wind_speed": {
        "standard_name": "wind_speed",
        "long_name": "Wind Speed",
        "units": "m s-1",
        "coordinates": "time longitude latitude height",
    },
    "wind_direction": {
        "standard_name": "wind_from_direction",
        "long_name": "Wind Direction",
        "units": "degrees",
        "coordinates": "time longitude latitude height",
    },
}

list_of_flight_attrs = [
    "True Heading (deg)",
    "True Air Speed (m/s)",
    "Ground Track (deg)",
    "Ground Speed (m/s)",
    "Longitude (deg)",
    "Latitude (deg)",
    "MSL Altitude (m)",
    "Geopotential Altitude (m)",
    "Software Notes",
    "Format Notes",
]


def get_flight_attrs(a_filepath, list_of_flight_attrs=list_of_flight_attrs):

    flight_attrs = {}

    file_path = a_filepath  # a_filepaths[i][0]
    f = open(file_path, "r")

    lines = f.readlines()

    # checking which line number in the file has the string we are looking for: "Launch Obs Done?"
    # this line number changes for different files
    for attr in list_of_flight_attrs:
        for j in range(len(lines)):
            if attr in lines[j]:
                line_id = j
                break

        if attr == "True Air Speed (m/s)":
            attr = "True-Air-Speed-(ms-1)"
        elif attr == "Ground Speed (m/s)":
            attr = "Ground-Speed-(ms-1)"
        elif attr == "Software Notes":
            attr = "AVAPS-Software-Notes"
        elif attr == "Format Notes":
            attr = "AVAPS-Format-Notes"

        attr = attr.replace(" ", "-")

        if "AVAPS" in attr:
            flight_attrs[attr] = lines[line_id].split("= ")[1]
        else:
            flight_attrs[attr] = float(lines[line_id].split("= ")[1])

    f.close()

    return flight_attrs


def get_global_attrs(Platform, file_time, sonde_ds):
    ###----- Global Attributes -----###

    nc_global_attrs = {
        "Title": "Sounding data containing temperature, pressure, humidity,"
        " latitude, longitude, wind direction, wind speed, and time",
        "Campaign": "EUREC4A-ATOMIC (Jan-Feb, 2020)",
        "Platform": Platform,
        "Instrument": "Vaisala RD41",
        "Launch-date": str(pd.to_datetime(file_time).date()),
        "Launch-time-(UTC)": str(sonde_ds.launch_time.values),
        "Sonde-Serial-ID": sonde_ds.SondeId,
        "ASPEN-Version": sonde_ds.AspenVersion,
        "Processing-Time": sonde_ds.ProcessingTime,
        "Mission-PI": "Mission PI",
        "Author": "Geet George",
        "Author-Email": "geet.george@mpimet.mpg.de",
        "JOANNE-version": joanne.__version__,
        "Conventions": "CF-1.7",
        "featureType": "trajectory",
        "Creation-Time": str(datetime.datetime.utcnow()) + " UTC",
    }

    return nc_global_attrs
