import datetime
import glob

import xarray as xr
import pandas as pd

import joanne

# %%

list_of_vars = ["alt", "time", "wspd", "wdir", "ta", "p", "rh", "lat", "lon"]

nc_meta = {
    "time": {
        "standard_name": "time",
        "long_name": "Time of recorded measurement",
        # "units": "seconds since 2020-01-01 00:00:00 UTC",
        # "calendar": "gregorian",
        "axis": "T",
    },
    "alt": {
        "standard_name": "geopotential_height",
        "long_name": "Geopotential Height",  # obtained by integrating upwards the atmospheric thickness estimated using the hypsometric equation",
        "units": "m",
        "axis": "Z",
        "positive": "up",
    },
    "lat": {
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees_north",
        "axis": "Y",
    },
    "lon": {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
        "axis": "X",
    },
    "p": {
        "standard_name": "air_pressure",
        "long_name": "atmospheric pressure",
        "units": "Pa",
        "coordinates": "time lon lat alt",
    },
    "ta": {
        "standard_name": "air_temperature",
        "long_name": "air temperature",
        "units": "K",
        "coordinates": "time lon lat alt",
    },
    "rh": {
        "standard_name": "relative_humidity",
        "long_name": "Relative Humidity",
        "units": "",
        "coordinates": "time lon lat alt",
    },
    "wspd": {
        "standard_name": "wind_speed",
        "long_name": "Wind Speed",
        "units": "m s-1",
        "coordinates": "time lon lat alt",
    },
    "wdir": {
        "standard_name": "wind_from_direction",
        "long_name": "Wind Direction",
        "units": "degrees",
        "coordinates": "time lon lat alt",
    },
    # "sonde_id": {
    #     "descripion": "unique sonde ID in the format PLATFORM_FLIGHT-ID_sSONDE-NUMBER-FOR-THE-FLIGHT",
    #     "long_name": "sonde identifier",
    #     "cf_role": "trajectory_id",
    # },
}

# nc_dims = {"alt":["time"], "time", "wspd", "wdir", "ta", "p", "rh", "lat", "lon"}

list_of_flight_attrs = [
    "True Heading (deg)",
    "True Air Speed (m/s)",
    "Ground Track (deg)",
    "Ground Speed (m/s)",
    "Longitude (deg)",
    "Latitude (deg)",
    "MSL Altitude (m)",
    "Geopotential Altitude (m)",
    # "Software Notes",
    # "Format Notes",
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
        "title": "EUREC4A JOANNE Level-2",
        "Conventions": "CF-1.8",
        "campaign_id": "EUREC4A",
        "project_id": "JOANNE",
        "platform_id": Platform,
        "instrument_id": "Vaisala RD-41",
        "product_id": "Level-2",
        "AVAPS-Software-version": "Version 4.1.2",
        "ASPEN-version": sonde_ds.AspenVersion,
        "JOANNE-version": joanne.__version__,
        "launch-date": str(pd.to_datetime(file_time).date()),
        "launch-time-(UTC)": str(sonde_ds.launch_time.values),
        "sonde-serial-ID": sonde_ds.SondeId,
        # "sounding_id": sonde_id,
        "processing-time": sonde_ds.ProcessingTime,
        # "Mission-PI": "Mission PI",
        "author": "Geet George",
        "author_email": "geet.george@mpimet.mpg.de",
        "featureType": "trajectory",
        "Creation-Time": str(datetime.datetime.utcnow()) + " UTC",
    }

    return nc_global_attrs
