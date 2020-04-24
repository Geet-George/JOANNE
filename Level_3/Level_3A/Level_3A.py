# %%
import glob

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import metpy.interpolate as mpinterp
import numpy as np
import requests
import xarray as xr
from metpy.units import units
from metpy import constants as mpconsts
from tqdm import tqdm

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

    pw = mpcalc.precipitable_water(
        dp * units.degC, dataset.pressure.values * units.hPa, top=altitude_limit
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
    dataset["relative_humidity"] = (dataset.pressure.dims, rh)

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

    dataset["launch_time"] = np.datetime64(dataset.attrs["Launch time (UTC)"])
    dataset["Platform"] = dataset.attrs["Platform"]
    dataset["flight_height"] = dataset.attrs["Geopotential Altitude (m)"]
    dataset["flight_lat"] = dataset.attrs["Latitude (deg)"]
    dataset["flight_lon"] = dataset.attrs["Longitude (deg)"]

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
    by swapping dimension from 'obs' to 'height', and adding the 'specific_humidity' and
    'potential_temperature' variables to the dataset.                                
    """

    dataset_to_interpolate = xr.open_dataset(file_path).swap_dims({"obs": "height"})
    dataset_to_interpolate = adding_q_and_theta_to_dataset(dataset_to_interpolate)

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
    for Level-3 gridding, following these steps :

        (i) Variables 'specific_humidity','potential_temperature' and 'precipitable_water' are added to the dataset

        (ii) All variables along 'height' dimension in dataset are linearly interpolated along the height dimension, at specified height 
        intervals (default 10 m) and up to specified altitude (default 10 km) 

        (iii) Pressure values are interpolated using a logarithmic interpolation scheme
        and these values replace the linearly interpolated pressure values. 
        
        Caveat : The difference between these different interpolations is in the order of 0.005 hPa, which is lower than 
        the measurement uncertainty of the pressure sensor itself of the RD-41 itself. There is a significant computational 
        cost involved, however, with carrying out this added logarithmic interpolation of pressure, since it employs an 
        iterative numerical method. If linearly interpolated values of pressure suffice for the user, they can choose
        to specify keyword pressure_interpolation = 'linear' and this will skip the logarithmic interpolation, making 
        the function perform faster.

        (iv) The temperature and moisture variables are to be interpolated with values of theta and q, respectively. Thus,
        after interpolation, T and RH variables are recomputed from the interpolated values of theta and q. The new values for
        T and RH will replace the previously interpolated T and RH variables. Although, T and RH are the originally measured 
        properties by the dropsonde sensors, for interpolation q and theta are preferred, as these variables are conserved.
    """

    if type(file_path_OR_dataset) is str:
        dataset = ready_to_interpolate(file_path_OR_dataset)
    else:
        dataset = file_path_OR_dataset

    dataset = add_platform_details_as_var(dataset)
    dataset = adding_precipitable_water_to_dataset(dataset)

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

    for id_, file_path in enumerate(list_of_files):
        interp_list[id_] = interpolate_for_level_3(
            file_path,
            height_limit=height_limit,
            vertical_spacing=vertical_spacing,
            pressure_log_interp=pressure_log_interp,
        )

    dataset = concatenate_soundings(interp_list)

    return dataset


# %%


def main():
    lv2_data_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Data/code_testing_data/"
    lv3_data_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Data/Level_3/Test_data/"

    lv3_dataset = lv3_structure_from_lv2(lv2_data_directory)
    lv3_dataset.to_netcdf(lv3_data_directory + "Level_3A.nc")


if __name__ == "__main__":
    main()

# %%

lv2_data_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Level_2/Test_data/"
list_of_files = retrieve_all_files(lv2_data_directory)

dataset = xr.open_dataset(list_of_files[0])


# %%
