# %%
import glob

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import metpy.interpolate as mpinterp
import numpy as np
import requests
import xarray as xr
from metpy.units import units
from tqdm import tqdm

# %%

### Defining functions


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
            calc_T_from_theta

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
        
        dataset : Dataset 
    
    Output :

        dataset : Original dataset with added variables of 'specific_humidity' and 'potential_temperature'
                 
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

def add_log_interp_pressure_to_dataset(dataset,interp_dataset=None) :
    """
    Input : 

        dataset : dataset

        interp_dataset : 
# %%

len(pressure_interpolation(t1.pressure.values, t1.height.values, interp_ds.height.values))

# %%


def main():
    lv2_data_directory = "JOANNE/Data/Test_data/"
    all_lv2_nc_files = sorted(glob.glob(lv2_data_directory + "*.nc"))

    for i in range(len(all_lv2_nc_files)):
        t1 = xr.open_dataset(all_lv2_nc_files[i]).swap_dims({"obs": "height"})
        t1 = adding_q_and_theta_to_dataset(t1)
        interp_ds = interp_along_height(t1)

        interp_ds = substitute_T_and_RH_for_interpolated_dataset(interp_ds)

    return interp_ds


if __name__ == "__main__":
    main()
    print("Kya re?")

# %%
