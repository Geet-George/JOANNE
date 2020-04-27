# %%
import datetime
import glob

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import metpy.interpolate as mpinterp
import numpy as np
import requests
import xarray as xr
from metpy import constants as mpconsts
from metpy.future import precipitable_water
from metpy.units import units
from tqdm import tqdm
import subprocess
# %%

try:
    git_module_version = subprocess.check_output(["git", "describe"]).strip().decode("utf-8")
except:
    git_module_version = "--"

# %%

### Defining functions


def retrieve_all_files(directory, file_ext="*.nc"):
    """
    Input : 

        directory : string
                    directory where the files are stored
        file_ext : string
                   extension of files needed; default is '*.nc'
    Output :

        list_of_files : list
                        list containing the file paths to all NC files in specified directory
    """

    list_of_files = sorted(glob.glob(directory + file_ext))

    return list_of_files


def remove_non_mono_incr_alt(lv2dataset):

    """
    This function removes the indices in the 
    geopotential height ('alt') array that are not monotonically 
    increasing and return a list of the remaining indices
    """
    ds_ka_alt = lv2dataset.height

    mono_ind = []
    g = 0
    while g < len(ds_ka_alt) - 1:
        n = 1
        prv = ds_ka_alt[g]
        nxt = ds_ka_alt[g + n]
        if prv < nxt:
            g += 1
            continue
        else:
            while prv >= nxt:
                mono_ind.append(g + n)
                if g + n == len(ds_ka_alt) - 1:
                    break
                else:
                    n += 1
                    nxt = ds_ka_alt[g + n]
            g += n

    x = list(set(range(len(ds_ka_alt))) - set(mono_ind))

    return x


def strictly_increasing(L):
    """
    This function checks if the provided array is strictly increasing
    """
    return all(x < y for x, y in zip(L, L[1:]))


def interp_along_height(dataset, height_limit=10000, vertical_spacing=10):
    """
    Input :
    
        dataset : Dataset with variables along "height' dimension
        height_limit : altitude up to which interpolation is to be carried out;
        default = 10 km
        vertical_spacing : vertical spacing to which values are to be interpolated;
        default = 10 m

    Output :

        new_interpolated_ds : New dataset with given dataset's variables 
        interpolated at given vertical_spacing up to given height_limit


    Function to interpolate all values along the height dimension of a netCDF dataset
    to a specified vertical spacing (10 m default) upto a given height level (10 km default)
    Given dataset must have data variables along the height dimension
    """
    new_index = np.arange(0, height_limit + vertical_spacing, vertical_spacing)

    new_interpolated_ds = dataset.interp(height=new_index)

    return new_interpolated_ds


def calc_q_from_rh(dataset):
    """
    Input :

        dataset : Dataset 

    Output :

        q : Specific humidity values
    
    Function to estimate specific humidity from the relative humidity,
    temperature and pressure in the given dataset. This function uses MetPy's
    functions to get q:

    (i) mpcalc.dewpoint_from_relative_humidity()
    (ii) mpcalc.specific_humidity_from_dewpoint()
                        
    """
    dp = mpcalc.dewpoint_from_relative_humidity(
        dataset.temperature.values * units.degC, dataset.relative_humidity.values / 100,
    ).magnitude

    q = mpcalc.specific_humidity_from_dewpoint(
        dp * units.degC, dataset.pressure.values * units.hPa
    ).magnitude

    return q


def calc_theta_from_T(dataset):
    """
    Input :

        dataset : Dataset

    Output :

        theta : Potential temperature values 

    Function to estimate potential temperature from the
    temperature and pressure in the given dataset. This function uses MetPy's
    functions to get theta:

    (i) mpcalc.potential_temperature()
    
    """
    theta = mpcalc.potential_temperature(
        dataset.pressure.values * units.hPa, dataset.temperature.values * units.degC
    ).magnitude

    return theta


def calc_T_from_theta(dataset):
    """
    Input :

        dataset : Dataset

    Output :

        T : Temperature values 

    Function to estimate temperature from potential temperature and pressure,
    in the given dataset. This function uses MetPy's
    functions to get T:

    (i) mpcalc.temperature_from_potential_temperature()
    
    """
    T = (
        mpcalc.temperature_from_potential_temperature(
            dataset.pressure.values * units.hPa,
            dataset.potential_temperature.values * units.kelvin,
        ).magnitude
        - 273.15
    )

    return T


def calc_rh_from_q(dataset, T=None):
    """
    Input :

        dataset : Dataset

        T : Temperature values estimated from interpolated potential_temperature;
            if not specified, function will calculate this from given dataset using
            calc_T_from_theta()

    Output :

        rh : Relative humidity values 

    Function to estimate relative humidity from specific humidity, temperature
    and pressure in the given dataset. This function uses MetPy's
    functions to get rh:

    (i) mpcalc.relative_humidity_from_specific_humidity()
    
    """
    if T is None:
        T = calc_T_from_theta(dataset)

    rh = mpcalc.relative_humidity_from_specific_humidity(
        dataset.specific_humidity.values,
        T * units.degC,
        dataset.pressure.values * units.hPa,
    ).magnitude

    return rh


def add_wind_components_to_dataset(dataset):
    """
    Input :

        dataset : xarray dataset
    
    Output :

        dataset : xarray dataset
                  Original dataset, with added variables 'u_wind' and 'v_wind' calculated
                  from wind_speed and wind_direction of the given dataset
    
    Function to compute u and v components of wind, from wind speed and direction in the given dataset,
    and add them as variables to the dataset.                   
    """
    u_wind, v_wind = mpcalc.wind_components(
        dataset.wind_speed.values * units["m/s"],
        dataset.wind_direction.values * units.deg,
    )

    dataset["u_wind"] = (dataset.pressure.dims, u_wind.magnitude)
    dataset["v_wind"] = (dataset.pressure.dims, v_wind.magnitude)

    return dataset


def adding_q_and_theta_to_dataset(dataset):

    """
    Input :
        
        dataset : xarray dataset 
    
    Output :

        dataset : xarray dataset
                  Original dataset with added variables of 
                  'specific_humidity' and 'potential_temperature'
                 
    Function to add variables of 'specific_humidity' and 'potential_temperature' in original
    dataset using functions 
    (i) calc_q 
    (ii) calc_theta
    
    """
    if "potential_temperature" not in list(dataset.data_vars):

        theta = calc_theta_from_T(dataset)
        dataset["potential_temperature"] = (dataset.pressure.dims, theta)

    if "specific_humidity" not in list(dataset.data_vars):

        q = calc_q_from_rh(dataset)
        dataset["specific_humidity"] = (dataset.pressure.dims, q)

    return dataset


def adding_precipitable_water_to_dataset(dataset, altitude_limit=None):
    """
    Input :
        dataset : xarray dataset

    Output :
        dataset : xarray dataset
                  Original dataset with added variable of precipitable_water

    Function to add variable 'precipitable_water' to given dataset, with no dimension,
    using MetPy functions :

    (i) mpcalc.precipitable_water()
    (ii) mpcalc.dewpoint_from_relative_humidity()
    """
    dp = mpcalc.dewpoint_from_relative_humidity(
        dataset.temperature.values * units.degC, dataset.relative_humidity.values / 100
    ).magnitude

    pw = precipitable_water(
        dataset.pressure.values * units.hPa, dp * units.degC, top=altitude_limit
    ).magnitude

    dataset["precipitable_water"] = pw

    return dataset


def adding_static_stability_to_dataset(dataset, method="gradient"):
    """
    Input :
        dataset : xarray dataset

    Output :
        dataset : xarray dataset
                  Original dataset with added variable of static_stability

    Function to add variable 'static_stability' to given dataset, along height dimension,
    using gradient of theta with p or using the MetPy functions of mpcalc.static_stability(). 
    The former is the default method, the latter can be selected by giving keyword argument as
    method = 'B92', which stands for Bluestein(1992).
    """
    if method == "gradient":
        pot = calc_theta_from_T(dataset)
        pres = dataset.pressure.values
        d_pot = pot[:-1] - pot[1:]
        d_pres = pres[1:] - pres[:-1]
        ss = d_pot / d_pres

    if method == "B92":
        ss = mpcalc.static_stability(
            dataset.pressure.values * units.hPa, dataset.temperature.values * units.degC
        ).magnitude

    static_stability = np.full(len(dataset.temperature), np.nan)

    static_stability[0] = 0
    static_stability[1:] = ss

    dataset["static_stability"] = (dataset.temperature.dims, static_stability)

    return dataset


def substitute_T_and_RH_for_interpolated_dataset(dataset):

    """
    Input :

        dataset : Dataset interpolated along height

    Output :

        dataset : Original dataset with new T and RH

    Function to remove interoplated values of T and RH in the original dataset and 
    replace with new values of T and RH,
    calculated from values of interpolated theta and q, respetively
    """
    T = calc_T_from_theta(dataset)
    rh = calc_rh_from_q(dataset, T=T)

    dataset["temperature"] = (dataset.pressure.dims, T)
    dataset["relative_humidity"] = (dataset.pressure.dims, rh * 100)

    return dataset


def add_cloud_flag(dataset):
    """
    Function under construction
    """
    cloud_flag = np.full(len(dataset.obs), 0)

    dataset["cloud_flag"] = (dataset.pressure.dims, cloud_flag)

    return dataset


def pressure_interpolation(
    pressures, altitudes, output_altitudes, convergence_error=0.05
):
    """
    Interpolates pressure on altitude grid
    
    The pressure is interpolated logarithmically.
    
    Input
    -----
    pressure : array
        pressures in hPa
    altitudes : array
        altitudes in m belonging to pressure values
    output_altitudes : array
        altitudes (m) on which the pressure should
        be interpolated to
    convergence_error : float
        Error that needs to be reached to stop
        convergence iteration
    
    Output
    -------
    pressure_interpolated : array
        array of interpolated pressure values
        on altitudes
    """

    pressure_interpolated = np.empty(len(output_altitudes))
    pressure_interpolated[:] = np.nan

    # Exclude heights outside of the intersection of measurements heights
    # and output_altitudes
    altitudes_above_measurements = output_altitudes > max(altitudes)
    range_of_alt_max = (
        np.min(
            np.where(
                altitudes_above_measurements
                | (output_altitudes == output_altitudes[-1])
            )
        )
        - 1
    )

    altitudes_below_measurements = output_altitudes < min(altitudes)
    range_of_alt_min = (
        np.max(
            np.where(
                altitudes_below_measurements | (output_altitudes == output_altitudes[0])
            )
        )
        + 1
    )

    for i in range(range_of_alt_min, range_of_alt_max):

        target_h = output_altitudes[i]

        lower_idx = np.nanmax(np.where(altitudes < target_h))
        upper_idx = np.nanmin(np.where(altitudes > target_h))

        p1 = pressures[lower_idx]  # pressure at lower altitude
        p2 = pressures[upper_idx]  # pressure at higher altitude
        a1 = altitudes[lower_idx]  # lower altitude
        a2 = altitudes[upper_idx]  # higher altitude

        xp = np.array([p1, p2])
        arr = np.array([a1, a2])

        err = 10

        if a2 - a1 < 100:
            while err > convergence_error:
                x = np.mean([p1, p2])
                ah = mpinterp.log_interpolate_1d(x, xp, arr, fill_value=np.nan)
                if ah > target_h:
                    p2 = x
                else:
                    p1 = x
                err = abs(ah - target_h)
            pressure_interpolated[i] = x

    return pressure_interpolated


def add_log_interp_pressure_to_dataset(
    dataset, interp_dataset=None, height_limit=10000, vertical_spacing=10
):
    """
    Input : 

        dataset : dataset

        interp_dataset : dataset
                         interpolated values of dataset;
                         if not specified, function will calculate this 
                         from given dataset using interp_along_height()

    Output :

        interp_dataset : dataset
                         returns modified interp_dataset with original linearly 
                         interpolated pressure variable replaced with logarithmically
                         interpolated pressure variable                      
    """

    if interp_dataset is None:
        interp_dataset = interp_along_height(
            dataset, height_limit=height_limit, vertical_spacing=vertical_spacing
        )

    pressure = pressure_interpolation(
        dataset.pressure.values, dataset.height.values, interp_dataset.height.values
    )

    interp_dataset["pressure"] = (dataset.pressure.dims, pressure)

    return interp_dataset


def add_platform_details_as_var(dataset):
    """
    Input :
        dataset : xarray dataset
                  dataset to which launch_time is to be included as a variable
    Output :
        dataset : xarray dataset
                  modified dataset, now with some platform details as variables, 
                  with no dimension attached to it
    """

    dataset["launch_time"] = (
        np.datetime64(dataset.attrs["Launch time (UTC)"]).astype("float") / 1e9
    )
    dataset["Platform"] = dataset.attrs["Platform"]
    dataset["flight_height"] = dataset.attrs["Geopotential Altitude (m)"]
    dataset["flight_lat"] = dataset.attrs["Latitude (deg)"]
    dataset["flight_lon"] = dataset.attrs["Longitude (deg)"]

    if dataset.attrs["Geopotential Altitude (m)"] < 4000:
        low_height_flag = 1
    else:
        low_height_flag = 0

    dataset["low_height_flag"] = low_height_flag

    return dataset


def ready_to_interpolate(file_path):

    """
    Input :

        file_path : string
                    path to NC file containing Level-2 data
    Output :

        dataset_to_interpolate : xarray dataset
                                 dataset ready for interpolation

    Function that takes in the path to Level-2 NC file and makes it ready for interpolation,
    by swapping dimension from 'obs' to 'height', and adding 'specific_humidity',
    'potential_temperature','wind_components','precipitable_water',
    and platform details variables to the dataset.                                
    """

    dataset_to_interpolate = xr.open_dataset(file_path).swap_dims({"obs": "height"})
    dataset_to_interpolate = adding_q_and_theta_to_dataset(dataset_to_interpolate)
    dataset_to_interpolate = add_wind_components_to_dataset(dataset_to_interpolate)
    dataset_to_interpolate = add_platform_details_as_var(dataset_to_interpolate)
    dataset_to_interpolate = adding_precipitable_water_to_dataset(
        dataset_to_interpolate
    )

    return dataset_to_interpolate


def interpolate_for_level_3(
    file_path_OR_dataset,
    height_limit=10000,
    vertical_spacing=10,
    pressure_log_interp=True,
):

    """
    Input :

        file_path_OR_dataset : string or  dataset
                               if file path to Level-2 NC file is provided as string, 
                               dataset will be created using the ready_to_interpolate() function,
                               if dataset is provided, it will be used directly

    Output :

        interpolated_dataset : xarray dataset
                               interpolated dataset

    Function to interpolate a dataset with Level-2 data, in the format 
    for Level-3 gridding
    """

    if type(file_path_OR_dataset) is str:
        dataset = ready_to_interpolate(file_path_OR_dataset)
    else:
        dataset = file_path_OR_dataset

    interpolated_dataset = interp_along_height(
        dataset, height_limit=height_limit, vertical_spacing=vertical_spacing
    )

    if pressure_log_interp is True:
        interpolated_dataset = add_log_interp_pressure_to_dataset(
            dataset, interpolated_dataset
        )

    interpolated_dataset = substitute_T_and_RH_for_interpolated_dataset(
        interpolated_dataset
    )

    interpolated_dataset = add_cloud_flag(interpolated_dataset)
    interpolated_dataset = adding_static_stability_to_dataset(interpolated_dataset)

    return interpolated_dataset


def concatenate_soundings(list_of_interpolated_dataset):
    """
    Input : 
        list_of_interpolated_dataset : list
                                       list containing individual, interpolated sounding profiles
                                       as xarray datasets
    Output :
        concatenated_dataset : xarray dataset
                               dataset with all soundings in list_of_soundings concatenated to a new
                               dimension called 'soundings', and swap 'height' dimension with 'obs', 
                               making 'height' a variable
    """
    concatenated_dataset = xr.concat(list_of_interpolated_dataset, dim="sounding")
    concatenated_dataset = concatenated_dataset.drop("obs").swap_dims({"height": "obs"})

    return concatenated_dataset


def lv3_structure_from_lv2(
    directory_OR_list_of_files,
    height_limit=10000,
    vertical_spacing=10,
    pressure_log_interp=True,
):
    """
    Input :
        directory_OR_list_of_files : string or list
                                     if directory where NC files are stored is provided as a string,
                                     a list of file paths for all NC files in the directory is created,
                                     otherwise a list of file paths needed to be gridded can also be 
                                     provided directly
    Output :
        dataset : xarray dataset
                  dataset with Level-3 structure
                  
    Function to create Level-3 gridded dataset from Level-2 files
    """

    if type(directory_OR_list_of_files) is str:
        list_of_files = retrieve_all_files(directory_OR_list_of_files, file_ext="*.nc")
    else:
        list_of_files = directory_OR_list_of_files

    interp_list = [None] * len(list_of_files)

    for id_, file_path in enumerate(tqdm(list_of_files)):
        interp_list[id_] = interpolate_for_level_3(
            file_path,
            height_limit=height_limit,
            vertical_spacing=vertical_spacing,
            pressure_log_interp=pressure_log_interp,
        )

    dataset = concatenate_soundings(interp_list)

    return dataset


# %%


# def main():
lv2_data_directory = (
    "/Users/geet/Documents/EUREC4A/JOANNE/Data/Level_2/"  # code_testing_data/"
)
lv3_data_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Data/Level_3/Test_data/"

lv3_dataset = lv3_structure_from_lv2(lv2_data_directory)
# lv3_dataset.to_netcdf(lv3_data_directory + "Level_3A.nc")

# if __name__ == "__main__":
#     lv3_dataset = main()
# %%

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
        #                       'valid_range' : [-90.  90.],
        "axis": "X",
    },
    "longitude": {
        "standard_name": "longitude",
        "long_name": "East Longitude",
        "units": "degree",
        #                       'valid_range' : [-180.  180.],
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

nc_data = {}

for var in nc_attrs.keys():
    nc_data[var] = lv3_dataset[var].values

nc_global_attrs = {
    "Title": "Gridded, sounding data from JOANNE Level-2",
    "Campaign": "EUREC4A-ATOMIC (Jan-Feb, 2020)",
    "Instrument": "Vaisala RD41 (AVAPS receiver aboard aircraft)",
    "Data Processing for Level-2": 'AvapsEditorVersion "BatchAspen V3.4.3"',
    "Author": "Geet George (MPI-M, Hamburg); geet.george@mpimet.mpg.de",
    "version": "0.1.0-alpha",
    "Conventions": "CF-1.7",
    "featureType": "trajectory",
    "Creation Time": str(datetime.datetime.utcnow()) + " UTC",
}


def create_variable(ds, var, **kwargs):
    """Insert the data into a variable in an :class:`xr.Dataset`"""
    data = nc_data[var]  # must be of type array
    attrs = nc_attrs[var].copy()
    dims = nc_dims[var]

    v = xr.Variable(dims, data, attrs=attrs)
    ds[var] = v

    return var


# %%

obs = lv3_dataset.obs.values
sounding = lv3_dataset.sounding.values

to_save_ds = xr.Dataset(coords={"obs": obs, "sounding": sounding})

for var in list_of_vars:
    create_variable(to_save_ds, var)

file_name = "EUREC4A" + "_Dropsonde-RD41_" + "Level_3A" + ".nc"

save_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Data/Level_3/"

comp = dict(zlib=True, complevel=4, fletcher32=True, _FillValue=np.finfo("float32").max)

encoding = {var: comp for var in to_save_ds.data_vars if var != "Platform"}

for key in nc_global_attrs.keys():
    to_save_ds.attrs[key] = nc_global_attrs[key]

to_save_ds.to_netcdf(
    save_directory + file_name, mode="w", format="NETCDF4", encoding=encoding
)
# %%

# df = pd.DataFrame(nc_attrs.values(), index=nc_attrs.keys())

# %%
