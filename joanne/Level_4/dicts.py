import datetime
import joanne

list_of_vars = [
    "sounding",
    "circle",
    # "launch_time",
    "alt",
    # "latitude",
    # "longitude",
    # "pressure",
    # "temperature",
    # "relative_humidity",
    # "wind_speed",
    # "wind_direction",
    # "u_wind",
    # "v_wind",
    # "potential_temperature",
    # "specific_humidity",
    # "precipitable_water",
    # "static_stability",
    # "low_height_flag",
    # "cloud_flag",
    "platform",
    "flight_height",
    # "flight_lat",
    # "flight_lon",
    "circle_lon",
    "circle_lat",
    "circle_diameter",
    "circle_time",
    # "dx",
    # "dy",
    "u",
    "dudx",
    "dudy",
    # "sondes_regressed",
    "segment_id",
    "v",
    "dvdx",
    "dvdy",
    "q",
    "dqdx",
    "dqdy",
    "ta",
    "dtadx",
    "dtady",
    "p",
    "dpdx",
    "dpdy",
    "D",
    "vor",
    # "density",
    # "mean_density",
    "W",
    "omega",
    # "h_adv_q",
    # "h_adv_ta",
    # "h_adv_p",
    # "h_adv_u",
    # "h_adv_v",
]

nc_attrs = {
    "sounding": {
        "standard_name": "sounding",
        "long_name": "Sonde number",
        "units": "",
        "axis": "T",
    },
    "circle": {
        "standard_name": "time",
        "long_name": "Circle number",
        "units": "",
        "axis": "T",
    },
    "launch_time": {
        "standard_name": "time",
        "long_name": "Time of dropsonde launch",
        "units": "seconds since 1970-01-01 00:00:00 UTC",
        "calendar": "gregorian",
        "axis": "T",
    },
    "alt": {
        "standard_name": "geopotential_height",
        "long_name": "Geopotential Height",
        "description": "Height obtained by integrating upwards the atmospheric thickness estimated from the hypsometric equation",
        "units": "m",
        "axis": "Z",
        "positive": "up",
    },
    "latitude": {
        "standard_name": "latitude",
        "long_name": "North Latitude",
        "units": "degree_north",
        "axis": "Y",
    },
    "longitude": {
        "standard_name": "longitude",
        "long_name": "East Longitude",
        "units": "degree_east",
        "axis": "X",
    },
    "p": {
        "standard_name": "air_pressure",
        "long_name": "Atmospheric Pressure",
        "units": "hPa",
        "coordinates": "launch_time longitude latitude height",
    },
    "ta": {
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
        "units": "degree",
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
        "units": "",
        "flag_meanings": "flight height below 4 km flight height at or above 4 km",
        "valid_range": "0, 1",
    },
    "cloud_flag": {
        "long_name": "flag to indicate presence of cloud",
        "flag_values": "1, 0",
        "units": "",
        "flag_meanings": "cloud no_cloud",
        "valid_range": "0, 1",
    },
    "platform": {
        "standard_name": "platform",
        "long_name": "platform of the flown circle",
        "coordinates": "circle",
        "units": "",
    },
    "flight_height": {
        "standard_name": "alt",
        "long_name": "mean height of the aircraft during the circle",
        "units": "m",
        "coordinates": "circle",
    },
    "segment_id": {
        "description": "unique segment ID in the format PLATFORM_FLIGHT-ID_cCIRCLE-NUMBER-FOR-THE-FLIGHT",
        "long_name": "segment (circle) identifier",
        "cf_role": "trajectory_id",
        "units": "",
    },
    "flight_lat": {
        "standard_name": "latitude",
        "long_name": "north latitude of the aircraft when the dropsonde was launched",
        "units": "degree_north",
        "coordinates": "launch_time",
    },
    "flight_lon": {
        "standard_name": "longitude",
        "long_name": "east longitude of the aircraft when the dropsonde was launched",
        "units": "degree_east",
        "coordinates": "launch_time",
    },
    "circle_lon": {  # mean lon for circles at all levels fitted as least square fit to all sondes
        "standard_name": "circle_center_longitude",
        "long_name": "east longitude of fitted circle for all regressed sondes in circle",
        "units": "degree_east",
        "coordinates": "circle",
    },
    "circle_lat": {  # mean y for circles at all levels fitted as least square fit to all sondes
        "standard_name": "circle_center_latitude",
        "long_name": "north latitude of fitted circle for all regressed sondes in circle",
        "units": "degree_north",
        "coordinates": "circle",
    },
    "circle_diameter": {  # mean diameter for circles at all levels fitted as least square fit to all sondes
        "standard_name": "circle_diameter",
        "long_name": "mean diameter of circle",
        "units": "m",
        "coordinates": "circle",
    },
    "circle_time": {
        "standard_name": "circle_time",
        "long_name": "mean time of circle",
        "units": "seconds since 1970-01-01 00:00:00 UTC",
        "calendar": "gregorian",
        "coordinates": "circle",
    },
    "dx": {
        "standard_name": "delta_x",
        "long_name": "zonal distance of sonde from xc of circle",
        "units": "m",
        "coordinates": "launch_time height",
    },
    "dy": {
        "standard_name": "delta_y",
        "long_name": "meridional distance of sonde from yc of circle",
        "units": "m",
        "coordinates": "launch_time height",
    },
    "u": {
        "standard_name": "eastward_wind",
        "long_name": "intercept value from regressed eastward wind in circle",
        "units": "m s-1",
        "coordinates": "circle height",
    },
    "dudx": {
        "long_name": "zonal gradient of eastward wind",
        "units": "s-1",
        "coordinates": "circle height",
    },
    "dudy": {
        "long_name": "meridional gradient of eastward wind",
        "units": "s-1",
        "coordinates": "circle height",
    },
    "sondes_regressed": {
        "long_name": "number of sondes regressed",
        "units": "",
        "coordinates": "circle height",
    },
    "v": {
        "standard_name": "northward_wind",
        "long_name": "intercept value from regressed northward wind in circle",
        "units": "m s-1",
        "coordinates": "circle height",
    },
    "dvdx": {
        "long_name": "zonal gradient of northward wind",
        "units": "s-1",
        "coordinates": "circle height",
    },
    "dvdy": {
        "long_name": "meridional gradient of northward wind",
        "units": "s-1",
        "coordinates": "circle height",
    },
    "q": {
        "standard_name": "specific_humidity",
        "long_name": "intercept value from regressed specific humidity in circle",
        "units": "kg kg-1",
        "coordinates": "circle height",
    },
    "dqdx": {
        "long_name": "zonal gradient of specific humidity",
        "units": "m-1",
        "coordinates": "circle height",
    },
    "dqdy": {
        "long_name": "meridional gradient of specific humidity",
        "units": "m-1",
        "coordinates": "circle height",
    },
    "ta": {
        "standard_name": "temperature",
        "long_name": "intercept value from regressed temperature in circle",
        "units": "degree_Celsius",
        "coordinates": "circle height",
    },
    "dtadx": {
        "long_name": "zonal gradient of temperature",
        "units": "degree_Celsius m-1",
        "coordinates": "circle height",
    },
    "dtady": {
        "long_name": "meridional gradient of temperature",
        "units": "degree_Celsius m-1",
        "coordinates": "circle height",
    },
    "p": {
        "standard_name": "pressure",
        "long_name": "intercept value from regressed pressure in circle",
        "units": "hPa",
        "coordinates": "circle height",
    },
    "dpdx": {
        "long_name": "zonal gradient of pressure",
        "units": "hPa m-1",
        "coordinates": "circle height",
    },
    "dpdy": {
        "long_name": "meridional gradient of pressure",
        "units": "hPa m-1",
        "coordinates": "circle height",
    },
    "D": {
        "standard_name": "divergence",
        "long_name": "horizontal mass divergence",
        "units": "s-1",
        "coordinates": "circle height",
    },
    "vor": {
        "standard_name": "vorticity",
        "long_name": "horizontal mass vorticity",
        "units": "s-1",
        "coordinates": "circle height",
    },
    "density": {
        "standard_name": "air_density",
        "long_name": "air density",
        "units": "kg m-3",
        "coordinates": "launch_time height",
    },
    "mean_density": {
        "standard_name": "air_density",
        "long_name": "mean air density across all sondes in circle",
        "units": "kg m-3",
        "coordinates": "circle height",
    },
    "W": {
        "standard_name": "vertical_velocity",
        "long_name": "large-scale atmospheric vertical velocity",
        "units": "m s-1",
        "coordinates": "circle height",
    },
    "omega": {
        "standard_name": "pressure_velocity",
        "long_name": "large-scale atmospheric pressure velocity",
        "units": "hPa h-1",
        "coordinates": "circle height",
    },
    "h_adv_q": {
        "standard_name": "q_advection",
        "long_name": "horizontal advection of specific humidity",
        "units": "kg kg-1 s-1",
        "coordinates": "circle height",
    },
    "h_adv_ta": {
        "standard_name": "T_advection",
        "long_name": "horizontal advection of temperature",
        "units": "degree_Celsius s-1",
        "coordinates": "circle height",
    },
    "h_adv_p": {
        "standard_name": "p_advection",
        "long_name": "horizontal advection of pressure",
        "units": "hPa s-1",
        "coordinates": "circle height",
    },
    "h_adv_u": {
        "standard_name": "u_advection",
        "long_name": "horizontal advection of eastward wind",
        "units": "m s-1 s-1",
        "coordinates": "circle height",
    },
    "h_adv_v": {
        "standard_name": "v_advection",
        "long_name": "horizontal advection of northward wind",
        "units": "m s-1 s-1",
        "coordinates": "circle height",
    },
}

nc_dims = {
    "sounding": ["sounding"],
    "circle": ["circle"],
    "segment_id": ["circle"],
    "launch_time": ["circle", "sounding"],
    "alt": ["alt"],
    "latitude": ["circle", "sounding", "alt"],
    "longitude": ["circle", "sounding", "alt"],
    "pressure": ["circle", "sounding", "alt"],
    "temperature": ["circle", "sounding", "alt"],
    "relative_humidity": ["circle", "sounding", "alt"],
    "wind_speed": ["circle", "sounding", "alt"],
    "wind_direction": ["circle", "sounding", "alt"],
    "u_wind": ["circle", "sounding", "alt"],
    "v_wind": ["circle", "sounding", "alt"],
    "potential_temperature": ["circle", "sounding", "alt"],
    "specific_humidity": ["circle", "sounding", "alt"],
    "precipitable_water": ["circle", "sounding"],
    "static_stability": ["circle", "sounding", "alt"],
    "low_height_flag": ["circle", "sounding"],
    "cloud_flag": ["circle", "sounding", "alt"],
    "platform": ["circle"],
    "flight_height": ["circle"],
    "flight_lat": ["circle", "sounding"],
    "flight_lon": ["circle", "sounding"],
    "circle_lon": ["circle"],
    "circle_lat": ["circle"],
    "circle_diameter": ["circle"],
    "circle_time": ["circle"],
    "dx": ["circle", "sounding", "alt"],
    "dy": ["circle", "sounding", "alt"],
    "u": ["circle", "alt"],
    "dudx": ["circle", "alt"],
    "dudy": ["circle", "alt"],
    "sondes_regressed": ["circle", "alt"],
    "v": ["circle", "alt"],
    "dvdx": ["circle", "alt"],
    "dvdy": ["circle", "alt"],
    "q": ["circle", "alt"],
    "dqdx": ["circle", "alt"],
    "dqdy": ["circle", "alt"],
    "ta": ["circle", "alt"],
    "dtadx": ["circle", "alt"],
    "dtady": ["circle", "alt"],
    "p": ["circle", "alt"],
    "dpdx": ["circle", "alt"],
    "dpdy": ["circle", "alt"],
    "D": ["circle", "alt"],
    "vor": ["circle", "alt"],
    "density": ["sounding", "circle", "alt"],
    "mean_density": ["circle", "alt"],
    "W": ["circle", "alt"],
    "omega": ["circle", "alt"],
    "h_adv_q": ["circle", "alt"],
    "h_adv_ta": ["circle", "alt"],
    "h_adv_p": ["circle", "alt"],
    "h_adv_u": ["circle", "alt"],
    "h_adv_v": ["circle", "alt"],
}

nc_global_attrs = {
    "title": "EUREC4A JOANNE Level-4",
    "Conventions": "CF-1.8",
    "campaign_id": "EUREC4A",
    "project_id": "JOANNE",
    "instrument_id": "Vaisala RD-41",
    "product_id": "Level-4",
    "AVAPS-Software-version": "Version 4.1.2",
    "ASPEN-version": "BatchAspen v3.4.3",
    "JOANNE-version": joanne.__version__,
    "author": "Geet George",
    "author_email": "geet.george@mpimet.mpg.de",
    "featureType": "trajectory",
    "creation_time": str(datetime.datetime.utcnow()) + " UTC",
    "Reference-Study": "https://doi.org/10.1175/JAS-D-18-0141.1",
}
