import datetime
import subprocess
from joanne import pkg_ver

# git_dir = "/Users/geet/Documents/EUREC4A/JOANNE/"

# try:
#     git_module_version = (
#         subprocess.check_output(["git", "describe"], cwd=git_dir)
#         .strip()
#         .decode("utf-8")
#     )
# except:
#     git_module_version = "--"

list_of_vars = [
    "launch_time",
    "height",
    "latitude",
    "longitude",
    "pressure",
    "temperature",
    "relative_humidity",
    "wind_speed",
    "wind_direction",
    "u_wind",
    "v_wind",
    "potential_temperature",
    "specific_humidity",
    "precipitable_water",
    "static_stability",
    "low_height_flag",
    "cloud_flag",
    "Platform",
    "flight_height",
    "flight_lat",
    "flight_lon",
]

nc_attrs = {
    "launch_time": {
        "standard_name": "time",
        "long_name": "Time of dropsonde launch",
        "units": "seconds since 1970-01-01 00:00:00 UTC",
        "calendar": "gregorian",
        "axis": "T",
    },
    "height": {
        "standard_name": "geopotential_height",
        "long_name": "Geopotential Height",
        "description": "Height obtained by integrating upwards the atmospheric thickness estimated from the hypsometric equation",
        "units": "gpm",
        "axis": "Z",
        "positive": "up",
    },
    "latitude": {
        "standard_name": "latitude",
        "long_name": "North Latitude",
        "units": "degree",
        "axis": "X",
    },
    "longitude": {
        "standard_name": "longitude",
        "long_name": "East Longitude",
        "units": "degree",
        "axis": "Y",
    },
    "pressure": {
        "standard_name": "air_pressure",
        "long_name": "Atmospheric Pressure",
        "units": "hPa",
        "coordinates": "launch_time longitude latitude height",
    },
    "temperature": {
        "standard_name": "air_temperature",
        "long_name": "Dry Bulb Temperature",
        "units": "degree_Celsius",
        "coordinates": "launch_time longitude latitude height",
    },
    "potential_temperature": {
        "standard_name": "potential_temperature",
        "long_name": "potential temperature",
        "units": "K",
        "coordinates": "launch_time longitude latitude height",
    },
    "relative_humidity": {
        "standard_name": "relative_humidity",
        "long_name": "Relative Humidity",
        "units": "%",
        "coordinates": "launch_time longitude latitude height",
    },
    "specific_humidity": {
        "standard_name": "specific_humidity",
        "long_name": "Specific humidity",
        "units": "kg kg-1",
        "coordinates": "launch_time longitude latitude height",
    },
    "wind_speed": {
        "standard_name": "wind_speed",
        "long_name": "Wind Speed",
        "units": "m s-1",
        "coordinates": "launch_time longitude latitude height",
    },
    "u_wind": {
        "standard_name": "eastward_wind",
        "long_name": "u-component of the wind",
        "units": "m s-1",
        "coordinates": "launch_time longitude latitude height",
    },
    "v_wind": {
        "standard_name": "northward_wind",
        "long_name": "v-component of the wind",
        "units": "m s-1",
        "coordinates": "launch_time longitude latitude height",
    },
    "wind_direction": {
        "standard_name": "wind_from_direction",
        "long_name": "Wind Direction",
        "units": "m/s",
        "coordinates": "launch_time longitude latitude height",
    },
    "precipitable_water": {
        "standard_name": "precipitable_water",
        "long_name": "integrated water vapour in the measured column",
        "units": "kg m-2",
        "coordinates": "launch_time",
    },
    "static_stability": {
        "standard_name": "static_stability",
        "long_name": "static stability",
        "description": "gradient of potential temperature along the pressure grid",
        "units": " K hPa-1",
        "coordinates": "launch_time longitude latitude height",
    },
    "low_height_flag": {
        "long_name": "flag to indicate if flight height at launch was low",
        "flag_values": "1, 0",
        "flag_meanings": "flight height below 4 km flight height at or above 4 km",
        "valid_range": "0, 1",
    },
    "cloud_flag": {
        "long_name": "flag to indicate presence of cloud",
        "flag_values": "1, 0",
        "flag_meanings": "cloud no_cloud",
        "valid_range": "0, 1",
    },
    "Platform": {
        "standard_name": "platform",
        "long_name": "platform from which the sounding was made",
        "coordinates": "launch_time",
    },
    "flight_height": {
        "standard_name": "height",
        "long_name": "height of the aircraft when the dropsonde was launched",
        "units": "m",
        "coordinates": "launch_time",
    },
    "flight_lat": {
        "standard_name": "latitude",
        "long_name": "north latitude of the aircraft when the dropsonde was launched",
        "units": "degree north",
        "coordinates": "launch_time",
    },
    "flight_lon": {
        "standard_name": "longitude",
        "long_name": "east longitude of the aircraft when the dropsonde was launched",
        "units": "degree east",
        "coordinates": "launch_time",
    },
}

nc_dims = {
    "launch_time": ["sounding"],
    "height": ["obs"],
    "latitude": ["sounding", "obs"],
    "longitude": ["sounding", "obs"],
    "pressure": ["sounding", "obs"],
    "temperature": ["sounding", "obs"],
    "relative_humidity": ["sounding", "obs"],
    "wind_speed": ["sounding", "obs"],
    "wind_direction": ["sounding", "obs"],
    "u_wind": ["sounding", "obs"],
    "v_wind": ["sounding", "obs"],
    "potential_temperature": ["sounding", "obs"],
    "specific_humidity": ["sounding", "obs"],
    "precipitable_water": ["sounding"],
    "static_stability": ["sounding", "obs"],
    "low_height_flag": ["sounding"],
    "cloud_flag": ["sounding", "obs"],
    "Platform": ["sounding"],
    "flight_height": ["sounding"],
    "flight_lat": ["sounding"],
    "flight_lon": ["sounding"],
}

nc_global_attrs = {
    "Title": "Gridded, sounding data from JOANNE Level-2",
    "Campaign": "EUREC4A-ATOMIC (Jan-Feb, 2020)",
    "Instrument": "Vaisala RD41 (AVAPS receiver aboard aircraft)",
    "Data Processing for Level-2": 'AvapsEditorVersion "BatchAspen V3.4.3"',
    "Author": "Geet George (MPI-M, Hamburg); geet.george@mpimet.mpg.de",
    "version": pkg_ver.pkg_version,
    "Conventions": "CF-1.7",
    "featureType": "trajectory",
    "Creation Time": str(datetime.datetime.utcnow()) + " UTC",
}
