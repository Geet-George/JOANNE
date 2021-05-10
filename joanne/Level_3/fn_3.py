# %%
import datetime
import glob
import os.path
import subprocess
import warnings
from importlib import reload

import joanne
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import metpy.interpolate as mpinterp
import numpy as np
import pandas as pd
import requests
import xarray as xr
from eurec4a_snd.interpolate import postprocessing as pp
from joanne.Level_3 import dicts
from metpy import constants as mpconsts

# from metpy.future import precipitable_water
from metpy.units import units
from tqdm import tqdm

reload(dicts)

warnings.filterwarnings("ignore", message="Mean of empty slice")

#  %%

### Interpolation values (Thanks to Hauke Schulz, MPIM (eurec4a_snd))

interpolation_grid = np.arange(0, 10010, 10)  # meters
interpolation_bins = np.arange(-5, 10015, 10).astype(
    "int"
)  # Bins len(interpolation_grid)+1; (a,b]; (meters)
max_gap_fill = (
    50  # Maximum data gap size that should be filled by interpolation (meters)
)

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
    ds_ka_alt = lv2dataset.alt

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


def interp_along_height(
    dataset, height_limit=10000, vertical_spacing=10, max_gap=50, method="bin"
):
    """
    Input :
    
        dataset : Dataset with variables along 'alt' dimension
        height_limit : altitude up to which interpolation is to be carried out;
        default = 10 km
        vertical_spacing : vertical spacing to which values are to be interpolated;
        default = 10 m
        max_gap : no interpolation if gap between two datapoints is > max_gap; 
        default = 50 m

    Output :

        new_interpolated_ds : New dataset with given dataset's variables 
        interpolated at given vertical_spacing up to given height_limit


    Function to interpolate all values along the alt dimension of a netCDF dataset
    to a specified vertical spacing (10 m default) upto a given height level (10 km default)
    Given dataset must have data variables along the alt dimension
    """

    if method == "linear_interpolate":

        new_index = np.arange(0, height_limit + vertical_spacing, vertical_spacing)

        new_interpolated_ds = dataset.interp(alt=interpolation_grid)
        # new_interpolated_ds = new_interpolated_ds.interpolate_na(
        #     "alt", max_gap=max_gap, use_coordinate=True
        # )
    elif method == "bin":

        new_interpolated_ds = dataset.groupby_bins(
            "alt",
            interpolation_bins,
            labels=interpolation_grid,
            restore_coord_dims=True,
        ).mean()
        # for some reason, the groupby does not bin lat,lon and time since they are coordinates
        dataset["time"] = (
            ["alt"],
            dataset.time.values.astype(float),
        )
        # adding them as extra variables
        for coords in ["lat", "lon", "time"]:
            new_interpolated_ds[coords] = (
                dataset[coords]
                .groupby_bins(
                    "alt",
                    interpolation_bins,
                    labels=interpolation_grid,
                    restore_coord_dims=False,
                )
                .mean()
            )
        new_interpolated_ds = new_interpolated_ds.transpose()
        new_interpolated_ds = new_interpolated_ds.rename({"alt_bins": "alt"})

        new_interpolated_ds = new_interpolated_ds.interpolate_na(
            "alt", max_gap=max_gap_fill, use_coordinate=True
        )

        new_interpolated_ds["time"] = (
            ["alt"],
            pd.DatetimeIndex(new_interpolated_ds.time.values),
        )
        new_interpolated_ds.encoding["time"] = {
            "units": "seconds since 2020-01-01",
            "dtype": "int32",
            "_FillValue": np.iinfo("int32").max,
        }
        new_interpolated_ds = new_interpolated_ds.rename({"time": "interpolated_time"})

    return new_interpolated_ds


def calc_q_from_rh(ds):
    """
    Input :

        ds : Dataset 

    Output :

        q : Specific humidity values
    
    Function to estimate specific humidity from the relative humidity,
    temperature and pressure in the given dataset. This function uses MetPy's
    functions to get q:

    (i) mpcalc.dewpoint_from_relative_humidity()
    (ii) mpcalc.specific_humidity_from_dewpoint()
                        
    """
    # dp = mpcalc.dewpoint_from_relative_humidity(
    #     dataset.T.values * units.degC, dataset.rh.values / 100,
    # ).magnitude

    # q = mpcalc.specific_humidity_from_dewpoint(
    #     dp * units.degC, dataset.p.values * units.hPa
    # ).magnitude

    e_s = pp.calc_saturation_pressure(ds.ta.values)
    w_s = mpcalc.mixing_ratio(e_s * units.Pa, ds.p.values * units.Pa).magnitude
    w = ds.rh.values * w_s
    q = w / (1 + w)

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
        dataset.p.values * units.Pa, dataset.ta.values * units.kelvin
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
    ta = mpcalc.temperature_from_potential_temperature(
        dataset.p.values * units.Pa, dataset.theta.values * units.kelvin,
    ).magnitude

    return ta


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

    # rh = mpcalc.relative_humidity_from_specific_humidity(
    #     dataset.q.values, T * units.degC, dataset.p.values * units.hPa,
    # ).magnitude

    w = dataset.q / (1 - dataset.q)
    e_s = pp.calc_saturation_pressure(dataset.ta.values)
    w_s = mpcalc.mixing_ratio(e_s * units.Pa, dataset.p.values * units.Pa).magnitude

    rh = w / w_s

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
    u, v = mpcalc.wind_components(
        dataset.wspd.values * units["m/s"], dataset.wdir.values * units.deg,
    )

    dataset["u"] = (dataset.p.dims, u.magnitude)
    dataset["v"] = (dataset.p.dims, v.magnitude)

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
    if "theta" not in list(dataset.data_vars):

        theta = calc_theta_from_T(dataset)
        dataset["theta"] = (dataset.p.dims, theta)

    if "q" not in list(dataset.data_vars):

        q = calc_q_from_rh(dataset)
        dataset["q"] = (dataset.p.dims, q)

    return dataset


# def adding_precipitable_water_to_dataset(dataset, altitude_limit=None):
#     """
#     Input :
#         dataset : xarray dataset

#     Output :
#         dataset : xarray dataset
#                   Original dataset with added variable of precipitable_water

#     Function to add variable 'precipitable_water' to given dataset, with no dimension,
#     using MetPy functions :

#     (i) mpcalc.precipitable_water()
#     (ii) mpcalc.dewpoint_from_relative_humidity()
#     """
#     dp = mpcalc.dewpoint_from_relative_humidity(
#         (dataset.ta.values - 273.15) * units.degC, dataset.rh.values
#     ).magnitude

#     pw = precipitable_water(
#         dataset.p.values * units.Pa, dp * units.degC, top=altitude_limit
#     ).magnitude

#     dataset["PW"] = pw

#     return dataset


# def adding_static_stability_to_dataset(dataset, method="gradient"):
#     """
#     Input :
#         dataset : xarray dataset

#     Output :
#         dataset : xarray dataset
#                   Original dataset with added variable of static_stability

#     Function to add variable 'static_stability' to given dataset, along height dimension,
#     using gradient of theta with p or using the MetPy functions of mpcalc.static_stability().
#     The former is the default method, the latter can be selected by giving keyword argument as
#     method = 'B92', which stands for Bluestein(1992).
#     """
#     if method == "gradient":
#         pot = calc_theta_from_T(dataset)
#         pres = dataset.pressure.values
#         d_pot = pot[:-1] - pot[1:]
#         d_pres = pres[1:] - pres[:-1]
#         ss = d_pot / d_pres

#     if method == "B92":
#         ss = mpcalc.static_stability(
#             dataset.pressure.values * units.hPa, dataset.temperature.values * units.degC
#         ).magnitude

#     static_stability = np.full(len(dataset.temperature), np.nan)

#     static_stability[0] = 0
#     static_stability[1:] = ss

#     dataset["static_stability"] = (dataset.temperature.dims, static_stability)

#     return dataset


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

    dataset["ta"] = (dataset.p.dims, T)
    dataset["rh"] = (dataset.p.dims, rh)

    return dataset


def compute_wdir_from_u_and_v(u, v):

    """
    Input :

        u : u (zonal) component of wind
        v : v (meridional) component of wind

    Output :

        w_dir : wind direction from north in degree
        w_spd : magnitude of wind speed in m/s
    """

    w_dir = 90 - np.arctan2(-v, -u) * (180 / np.pi)

    mask = w_dir <= 0
    if np.any(mask):
        w_dir[mask] += 360

    w_spd = np.sqrt(u ** 2 + v ** 2)

    return w_dir, w_spd


def substitute_wdir_for_interpolated_dataset(dataset):

    """
    Input :

        dataset : Dataset interpolated along height

    Output :

        dataset : Original dataset with new wind direction and wind speed

    Function to remove interoplated values of wdir in the original dataset and 
    replace with new values of wdir computed from u and v
    """

    w_dir, w_spd = compute_wdir_from_u_and_v(dataset.u, dataset.v)

    dataset["w_dir"] = (dataset.p.dims, w_dir)
    dataset["w_spd"] = (dataset.p.dims, w_spd)

    return dataset


def add_cloud_flag(dataset):
    """
    Function under construction
    """
    cloud_flag = np.full(len(dataset.alt), 0)

    dataset["cloud_flag"] = (dataset.p.dims, cloud_flag)

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

    p = pressure_interpolation(
        dataset.p.values, dataset.alt.values, interp_dataset.alt.values
    )

    interp_dataset["p"] = (dataset.p.dims, p)

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

    dataset["launch_time"] = np.datetime64(dataset.attrs["launch_time_(UTC)"])
    dataset["platform_id"] = dataset.attrs["platform_id"]
    dataset["flight_altitude"] = dataset.attrs["aircraft_geopotential_altitude_(m)"]
    dataset["flight_lat"] = dataset.attrs["aircraft_latitude_(deg_N)"]
    dataset["flight_lon"] = dataset.attrs["aircraft_longitude_(deg_E)"]

    if dataset.attrs["aircraft_geopotential_altitude_(m)"] < 4000:
        low_height_flag = np.int8(1)
    else:
        low_height_flag = np.int8(0)

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
    'potential_temperature','wind_components',
    and platform details variables to the dataset.                                
    """

    dataset_to_interpolate = xr.open_dataset(file_path).swap_dims({"time": "alt"})
    dataset_to_interpolate = adding_q_and_theta_to_dataset(dataset_to_interpolate)
    dataset_to_interpolate = add_wind_components_to_dataset(dataset_to_interpolate)

    # dataset_to_interpolate = adding_precipitable_water_to_dataset(
    # dataset_to_interpolate
    # )

    return dataset_to_interpolate


def get_N_and_m_values(interp_dataset, original_dataset, bin_length=10):
    """
    Input :

        dataset : Dataset

    Output :

        N_ptu : number of PTU values in a bin
        m_ptu : method for binning PTU values
        N_gps : number of GPS values in a bin
        m_gps : method for binning GPS values

    Function to estimate number of observations in bin and the method for retrieving 
    data in the bin, i.e. either no data, interpolation or averaging    
    """
    interp_dataset["N_p"] = xr.DataArray(
        original_dataset.p.groupby_bins(
            "alt",
            interpolation_bins,
            labels=interpolation_grid,
            restore_coord_dims=True,
        )
        .count()
        .values,
        dims=["alt"],
        coords={"alt": interp_dataset.alt.values},
    )
    interp_dataset["N_ta"] = xr.DataArray(
        original_dataset.ta.groupby_bins(
            "alt",
            interpolation_bins,
            labels=interpolation_grid,
            restore_coord_dims=True,
        )
        .count()
        .values,
        dims=["alt"],
        coords={"alt": interp_dataset.alt.values},
    )
    interp_dataset["N_rh"] = xr.DataArray(
        original_dataset.rh.groupby_bins(
            "alt",
            interpolation_bins,
            labels=interpolation_grid,
            restore_coord_dims=True,
        )
        .count()
        .values,
        dims=["alt"],
        coords={"alt": interp_dataset.alt.values},
    )
    interp_dataset["N_gps"] = xr.DataArray(
        original_dataset.u.groupby_bins(
            "alt",
            interpolation_bins,
            labels=interpolation_grid,
            restore_coord_dims=True,
        )
        .count()
        .values,
        dims=["alt"],
        coords={"alt": interp_dataset.alt.values},
    )

    m_p = interp_dataset["N_p"].values.astype(int)
    m_ta = interp_dataset["N_ta"].values.astype(int)
    m_rh = interp_dataset["N_rh"].values.astype(int)
    m_gps = interp_dataset["N_gps"].values.astype(int)

    m_p[(m_p == np.isnan)] = np.int8(0)
    m_p[(m_p == 1)] = np.int8(1)
    m_p[(m_p > 1)] = np.int8(2)

    m_ta[(m_ta == np.isnan)] = np.int8(0)
    m_ta[(m_ta == 1)] = np.int8(1)
    m_ta[(m_ta > 1)] = np.int8(2)

    m_rh[(m_rh == np.isnan)] = np.int8(0)
    m_rh[(m_rh == 1)] = np.int8(1)
    m_rh[(m_rh > 1)] = np.int8(2)

    m_gps[(m_gps == np.isnan)] = np.int8(0)
    m_gps[(m_gps == 1)] = np.int8(1)
    m_gps[(m_gps > 1)] = np.int8(2)

    interp_dataset["m_p"] = xr.DataArray(
        m_p, dims=["alt"], coords={"alt": interp_dataset.alt.values},
    ).astype("int8")

    interp_dataset["m_ta"] = xr.DataArray(
        m_ta, dims=["alt"], coords={"alt": interp_dataset.alt.values},
    ).astype("int8")

    interp_dataset["m_rh"] = xr.DataArray(
        m_rh, dims=["alt"], coords={"alt": interp_dataset.alt.values},
    ).astype("int8")

    interp_dataset["m_gps"] = xr.DataArray(
        m_gps, dims=["alt"], coords={"alt": interp_dataset.alt.values}
    ).astype("int8")

    return interp_dataset


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

    interpolated_dataset = get_N_and_m_values(
        interpolated_dataset, dataset, bin_length=vertical_spacing
    )

    if pressure_log_interp is True:
        interpolated_dataset = add_log_interp_pressure_to_dataset(
            dataset, interpolated_dataset
        )

    interpolated_dataset = substitute_T_and_RH_for_interpolated_dataset(
        interpolated_dataset
    )

    interpolated_dataset = substitute_wdir_for_interpolated_dataset(
        interpolated_dataset
    )

    dataset = add_platform_details_as_var(dataset)

    for var in [
        "platform_id",
        "flight_altitude",
        "flight_lat",
        "flight_lon",
        "launch_time",
        "low_height_flag",
        "sonde_id",
    ]:
        interpolated_dataset[var] = dataset[var]

    # interpolated_dataset = add_cloud_flag(interpolated_dataset)
    # interpolated_dataset = adding_static_stability_to_dataset(interpolated_dataset)

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
    concatenated_dataset = concatenated_dataset.swap_dims({"sounding": "sonde_id"})
    # concatenated_dataset = concatenated_dataset.drop(
    #     "time"
    # )  # .swap_dims({"height": "obs"})

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

    save_directory = "/Users/geet/Documents/JOANNE/Data/Level_3/Interim_files/"

    for id_, file_path in enumerate(tqdm(list_of_files)):

        file_name = (
            "EUREC4A_JOANNE_Dropsonde-RD41_"
            + str(file_path[file_path.find("RD41_") + 3 : file_path.find("RD41_") + 19])
            + "Level_3_v"
            + str(joanne.__version__)
            + ".nc"
        )

        if os.path.exists(save_directory + file_name):

            interp_list[id_] = xr.open_dataset(save_directory + file_name)

        else:

            interp_list[id_] = interpolate_for_level_3(
                file_path,
                height_limit=height_limit,
                vertical_spacing=vertical_spacing,
                pressure_log_interp=pressure_log_interp,
            )

            interp_list[id_].to_netcdf(save_directory + file_name)

    concat_list = []
    for i in interp_list:

        if "ta" in i.var():
            concat_list.append(i)

    dataset = concatenate_soundings(concat_list)

    return dataset


def create_variable(ds, var, data, dims=dicts.nc_dims, attrs=dicts.nc_attrs, **kwargs):
    """Insert the data into a variable in an :class:`xr.Dataset`"""
    data = data[var]  # must be of type array
    attrs = attrs[var].copy()
    dims = dims[var]

    v = xr.Variable(dims, data, attrs=attrs)
    ds[var] = v

    return var


# %%
