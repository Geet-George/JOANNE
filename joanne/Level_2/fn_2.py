# %%
import datetime
import glob
import sys
import warnings

from importlib import reload

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import seaborn as sb
import xarray as xr

# from seaborn import distplot
from tqdm import tqdm

import joanne
from joanne.Level_2 import dicts

reload(dicts)
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice encounter")
warnings.filterwarnings(
    "ignore", message="Attempted to set non-positive bottom ylim on a log-scaled axis"
)
# %%
def get_all_sondes_list(Platform):

    directory = "/Users/geet/Documents/JOANNE/Data/Level_1/" + Platform + "/"
    # directory where all sonde files are present

    a_dir = "/Users/geet/Documents/JOANNE/Data/Level_0/" + Platform + "/All_A_files/"
    # directory where all the A files are present

    logs_directory = "/Users/geet/Documents/JOANNE/Data/Level_2/logs_and_stats/"
    # directory to store logs and stats

    sonde_paths = sorted(glob.glob(directory + "*QC.nc"))
    # paths to the individual sonde files

    file_time_str = [None] * len(sonde_paths)
    file_time = [None] * len(sonde_paths)
    # list to store extracted sonde time from file name as string and as time

    a_files = [None] * len(sonde_paths)
    # list to store file names for the log files starting with A.

    sonde_ds = [None] * len(sonde_paths)
    # list to store individual datasets of all sondes from PQC files

    for i in range(len(sonde_paths)):
        file_time_str[i] = sonde_paths[i][-20:-5]
        file_time[i] = np.datetime64(
            pd.to_datetime(file_time_str[i], format="%Y%m%d_%H%M%S"), "s"
        )
        a_files[i] = "A" + file_time_str[i]
        sonde_ds[i] = xr.open_dataset(sonde_paths[i])

    return sonde_ds, directory, a_dir, logs_directory, a_files, file_time, sonde_paths


def get_var_count_sums(list_nc):

    s_time = [None] * len(list_nc)
    s_t = [None] * len(list_nc)
    s_rh = [None] * len(list_nc)
    s_p = [None] * len(list_nc)
    s_z = [None] * len(list_nc)
    s_u = [None] * len(list_nc)
    s_v = [None] * len(list_nc)
    # creating lists to store non-NaN index sums for all parameters

    list_of_variables = ["s_time", "s_t", "s_rh", "s_p", "s_z", "s_u", "s_v"]
    # list of parameter names as strings

    # sorting the non-NaN index sums from list_nc to the respective parameters' lists
    for j, var in enumerate(list_of_variables):
        for i in range(len(list_nc)):
            eval(var)[i] = list_nc[i][j]

    return list_of_variables, s_time, s_t, s_rh, s_p, s_z, s_u, s_v


def get_ld_flag_from_a_files(a_dir, a_files, logs_directory, Platform, logs=True):

    a_filepaths = []
    # list to store individual file paths for all A files

    for i in a_files:
        a_filepaths.append(sorted(glob.glob(a_dir + i + "*")))

    ld_FLAG = np.full(len(a_files), np.nan)
    # array to store ld_FLAG values

    # create and start writing a log file which will store sonde info about sondes with failed launch detection
    file = open(
        f"{logs_directory}no_launch_detect_logs_{Platform}_v{joanne.__version__}.txt",
        "w",
    )

    g = 0
    # counter of failed sondes

    for id_, i in enumerate(a_filepaths):

        try:
            if len(i) == 0:
                # if the file does not exist, i will be None, as no file would have been found and added to a_filepaths
                raise Exception
            else:
                file_path = i[0]
                f = open(file_path, "r")

        except Exception:
            print(
                f"{a_files[id_]} : This file does not exist, at least not in the given directory"
            )
            continue

        else:

            lines = f.readlines()

            # checking which line number in the file has the string we are looking for: "Launch Obs Done?"
            # this line number changes for different files
            for i, line in enumerate(lines):
                if "Launch Obs Done?" in line:
                    line_id = i
                    break

            # after line number is obtained, check if the value is an integer, or an exception
            try:
                if int(lines[line_id][25]) not in [0, 1]:
                    raise ValueError
                else:
                    a = int(lines[line_id][25])

            except ValueError:
                print("An exception flew by : Value is neither 0 nor 1")
                continue

            else:
                if a == 0:  # if value is 0, then the launch detection failed
                    ld_FLAG[id_] = False
                    g += 1
                    for line in lines:
                        if "Sonde ID/Type/Rev" in line:
                            # storing the sonde ID information and relevant details to the log file we created
                            file.write(line)
                        if "START Time:" in line:
                            # storing the sonde start time to the log file we created
                            file.write(line)
                            # line breaker in our log file as a break between two file records
                            file.write("------------------------------------------\n")
                            break
                else:
                    ld_FLAG[id_] = True

    file.write(f"In total, there were {g} sondes that didn't detect a launch.\n")
    # writing summary of failed sondes to the log file

    return ld_FLAG


def init_status_ds(
    list_of_variables, s_time, s_t, s_rh, s_p, s_z, s_u, s_v, ld_FLAG, file_time
):

    # adding ld_FLAG to the list of variables, for easy addition to the dataset
    if "ld_FLAG" not in list_of_variables:
        list_of_variables.append("ld_FLAG")

    data_vars = {}
    # dictionary to store variables to be added to the dataset

    # populating the dictionary with the variable names and respective data
    for var in list_of_variables:
        data_vars[var] = (["time"], eval(var))

    # Creating the dataset
    status_ds = xr.Dataset(data_vars, coords={"time": file_time})

    return status_ds


# # Plotting the distribution of ratio of parameter's non-NaN measurement counts to counts of time

# f, ax = plt.subplots(2, 1, figsize=(10, 14))

# # looping over all parameters
# for var in list_of_variables[1:-1]:

#     # calculating ratio
#     rat = status_ds[var] / status_ds.s_time

#     # plotting distribution on linear Y-axis
#     sb.distplot(
#         rat,
#         hist=False,
#         ax=ax[0],
#         kde_kws={"linewidth": 3, "alpha": 0.5},  # "marker":'o',
#         label=var[2:],
#     )

#     # plotting distribution on log Y-axis
#     sb.distplot(
#         rat,
#         hist=False,
#         ax=ax[1],
#         kde_kws={"linewidth": 0.5, "marker": "o", "alpha": 0.5},
#         label=var[2:],
#     )

#     # setting the second plot's Y-axis to log scale
#     ax[1].set_yscale("log")
#     plt.legend()
#     ax[1].set_ylim(0.01, 100)
#     ax[1].set_xlabel(r"$\dfrac{Count\ of\ Parameter}{Count\ of\ Time}$", fontsize=14)
#     [ax[m].set_ylabel("Distribution / %", fontsize=14) for m in [0, 1]]
#     [
#         [ax[m].spines[n].set_visible(False) for m in [0, 1]]
#         for n in ["left", "right", "top"]
#     ]
#     plt.suptitle(
#         f"Distribution for all {int(len(status_ds.time))} sondes launched from {Platform}",
#         fontsize=20,
#         verticalalignment="bottom",
#     )

# sb.despine(offset=10)
# plt.savefig(
#     f"{logs_directory}Count_of_measurements_{Platform}_v{joanne.__version__}.png",
#     dpi=300,
# )

# Time values are recorded every 0.25 seconds. Although, the PTU and GPS sensors have a measurement frequency of 2 Hz and 4 Hz, respectively, the distribution of measurements have a slightly different story. Based on the distribution, we know that the ideal case is for all parameters (except u,v) to have measurements at every other time record, and for u,v to have measurements at every time record. Since, the time records also include values during initialisation as well as during a little before and after the launch, when no signal can be sent back to the AVAPS PC, the actual ratio will always be lower than the ideal estimate of 1 (for u,v) and 0.5 (for the remaining parameters).

# The true distribution shows that peaks start to flatten around 0.8 and 0.4 for u,v and other parameters, respectively. Thus, sondes with ratios lower than these values are taken as not having a complete profile, and termed as 'ugly' sondes. These ugly sondes still have data, but because they are expected to have more NaN fields than most sondes, they are kept for more QC, NaN-filling later, depending upon the extent of the dearth of measurements in that sonde.

# P.S. It is still unclear to me why GPS values are only at every other time record, while u,v are always being estimated.


def add_ind_flags_to_statusds(status_ds, list_of_variables):

    # Determining ind_flags and adding to the dataset

    rat = [None] * len(list_of_variables[1:-1])
    rat_id = [None] * len(list_of_variables[1:-1])

    # looping over all parameters except time
    for i, var in enumerate(list_of_variables[1:-1]):

        # estimating the ratio of the parameter counts to time counts
        rat[i] = status_ds[var].values / status_ds.s_time.values

        if var == "s_u" or var == "s_v":
            thresh = 0.8

        else:
            thresh = 0.4

        # assigning flag values
        rat_id[i] = np.where(rat[i] > thresh, "good", "ugly")

        for j in range(len(rat[i])):
            if rat[i][j] == 0:
                rat_id[i][j] = "bad"

    ind_flag_vars = ["t_flag", "rh_flag", "p_flag", "z_flag", "u_flag", "v_flag"]
    # list of all ind_flags

    # adding the flags to the dataset
    for i, j in zip(ind_flag_vars, rat_id):
        status_ds[i] = (["time"], j)

    return status_ds, ind_flag_vars


def add_srf_flags_to_statusds(status_ds, sonde_paths):
    # Determining srf_flags and adding to the dataset

    srf_p_flag = np.full(len(sonde_paths), np.nan)
    srf_z_flag = np.full(len(sonde_paths), np.nan)
    srf_t_flag = np.full(len(sonde_paths), np.nan)
    srf_rh_flag = np.full(len(sonde_paths), np.nan)
    rms_palt_gpsalt = np.full(len(sonde_paths), np.nan)
    # creating arrays to store the srf_flags

    # looping over all sondes
    for id_, i in enumerate(sonde_paths):

        # estimating all srf_flags
        sonde = xr.open_dataset(i)
        srf_p_flag[id_] = pres_bounds(sonde)
        srf_z_flag[id_] = gps_bounds(sonde)
        srf_t_flag[id_] = tdry_bounds(sonde)
        srf_rh_flag[id_] = rh_bounds(sonde)
        rms_palt_gpsalt[id_] = palt_gpsalt_rms_check(sonde)

    srf_flag_vars = [
        "srf_p_flag",
        "srf_z_flag",
        "srf_t_flag",
        "srf_rh_flag",
        "rms_palt_gpsalt",
    ]
    # list of all srf_flags

    # adding the flags to the dataset
    for var in srf_flag_vars:
        status_ds[var] = (["time"], eval(var))

    return status_ds, srf_flag_vars


def get_the_ind_FLAG_to_statusds(status_ds, ind_flag_vars):
    # Determining ind_FLAG

    ind_FLAG = [None] * len(status_ds.time)

    ind_all_good = np.where(
        [
            all(
                item == "good"
                for item in status_ds[ind_flag_vars].sel(time=i).to_array().values
            )
            for i in status_ds.time
        ]
    )[0]

    ind_all_bad = np.where(
        [
            all(
                item == "bad"
                for item in status_ds[ind_flag_vars].sel(time=i).to_array().values
            )
            for i in status_ds.time
        ]
    )[0]

    ind_any_ugly_or_bad = np.where(
        [
            any(
                (item == "ugly") or (item == "bad")
                for item in status_ds[ind_flag_vars].sel(time=i).to_array().values
            )
            for i in status_ds.time
        ]
    )[0]

    for i in ind_all_good:
        ind_FLAG[i] = "GOOD"

    for i in ind_any_ugly_or_bad:
        ind_FLAG[i] = "UGLY"

    for i in ind_all_bad:
        ind_FLAG[i] = "BAD"

    status_ds["ind_FLAG"] = (["time"], ind_FLAG)

    return status_ds, ind_FLAG


def get_the_srf_FLAG_to_statusds(status_ds, srf_flag_vars):

    # Determining srf_FLAG

    srf_FLAG = [None] * len(status_ds.time)

    srf_all_good = set(
        np.where(
            [
                all(
                    item == 1
                    for item in status_ds[srf_flag_vars].sel(time=i).to_array().values
                )
                for i in status_ds.time
            ]
        )[0]
    )

    srf_all_bad = set(
        np.where(
            [
                all(
                    item == 0
                    for item in status_ds[srf_flag_vars].sel(time=i).to_array().values
                )
                for i in status_ds.time
            ]
        )[0]
    )

    any_bad = set(
        np.where(
            (
                [
                    any(
                        item == 0
                        for item in status_ds[srf_flag_vars]
                        .sel(time=i)
                        .to_array()
                        .values
                    )
                    for i in status_ds.time
                ]
            )
        )[0]
    )

    srf_any_bad = any_bad - any_bad.intersection(srf_all_bad)

    for i in srf_all_good:
        srf_FLAG[i] = "GOOD"

    for i in srf_any_bad:
        srf_FLAG[i] = "UGLY"

    for i in srf_all_bad:
        srf_FLAG[i] = "BAD"

    status_ds["srf_FLAG"] = (["time"], srf_FLAG)

    return status_ds, srf_FLAG


def get_the_FLAG(status_ds, ind_FLAG, srf_FLAG):

    # Determining sonde FLAG

    FLAG = [None] * len(status_ds.time)

    no_launch = np.where(status_ds.ld_FLAG != 1)[0]

    for i in range(len(FLAG)):

        if i in no_launch:
            FLAG[i] = "BAD"
        else:

            if srf_FLAG[i] == "BAD" and ind_FLAG[i] == "BAD":
                FLAG[i] = "BAD"
            elif srf_FLAG[i] == "GOOD" and ind_FLAG[i] == "GOOD":
                FLAG[i] = "GOOD"
            else:
                FLAG[i] = "UGLY"

    status_ds["FLAG"] = (["time"], FLAG)

    return status_ds


def get_total_non_nan_indices(sonde):

    """
    Retrieving the non-NaN indices for all parameters.
    
    'c' terms are complete arrays of indices that have non-NaN values
    's' terms are the counts of indices in the respective 'c' terms.
    
    
    Input : 
        sonde_path : Path to the sonde PQC file as a string
    Output :
        s_var : where, var is one of [time,t,rh,p,z,u,v]
        
    """

    #     import xarray as xr
    import numpy as np

    #     sonde = xr.open_dataset(sonde_path)

    c_time = ~np.isnan(sonde.time).values
    c_t = ~np.isnan(sonde.tdry).values
    c_rh = ~np.isnan(sonde.rh).values
    c_p = ~np.isnan(sonde.pres).values
    c_z = ~np.isnan(sonde.gpsalt).values
    c_u = ~np.isnan(sonde.u_wind).values
    c_v = ~np.isnan(sonde.v_wind).values

    s_time = c_time.sum()
    s_t = c_t.sum()
    s_rh = c_rh.sum()
    s_p = c_p.sum()
    s_z = c_z.sum()
    s_u = c_u.sum()
    s_v = c_v.sum()

    return s_time, s_t, s_rh, s_p, s_z, s_u, s_v


def pres_bounds(sonde, u_lim=1020, l_lim=1000):

    """
    Checking if maximum pressure of sonde is within bounds: 1000 hPa - 1020 hPa. 
    Value higher than bound is unrealistic, 
    and value lower than bound means sonde did not measure the bottommost levels of the atmosphere.
    
    This flag does not include any GPS values. Even if there were no pressure values above 1000 hPa, 
    there may still be GPS measurements in the lowest levels. 
    Such sondes can still be useful for wind and wind-derived products.
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if max pressure value is within bounds
               False, if max pressure value is out of bounds
    """

    if (sonde.pres.max() < 1000) | (sonde.pres.max() > 1020):
        return False
    else:
        return True


def gps_bounds(sonde, limit=30):

    """
    Checking if maximum GPS altitude of sonde is within bounds: <= limit (default assigned as 30 m)
    Value higher than bound has no near-surface measurements
    
    This flag does not include any pressure values. Even if there were no GPS values below 30 m,
    there may still be PTU measurements in the lowest levels. 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if min GPS altitude value is within bounds
               False, if min GPS altitude value is out of bounds
        
    """

    if sonde.gpsalt.min() > limit:
        return False
    else:
        return True


def tdry_bounds(sonde, u_limit=30, srf_limit=20):
    """
    Checking if tdry (air temperature) is within bounds:
    
    1. Maximum air temperature recorded should not be greater than the upper limit (u_limit), 
        set to a default value of 30 deg C.
    2. Mean air temperature in the bottom 100 m (by gpsalt) should not be lesser than srf_limit, 
        set to a default of 20 deg C.
    
    If any of the above limits is violated, the tdry for the sonde is considered out of bounds,
    and marked as 'False'. The sonde is also marked 'False', if there are no measurements in the 
    bottom 100 m (by gpsalt). 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if tdry value is within bounds
               False, if tdry is out of bounds or if measurements are not available
    """
    if sonde.tdry.max() >= u_limit:
        return False
    elif sonde.tdry.where(sonde.gpsalt < 100, drop=True).sum() == 0:
        return False
    elif sonde.tdry.where(sonde.gpsalt < 100, drop=True).mean() < srf_limit:
        return False
    else:
        return True


def rh_bounds(sonde, srf_limit=50):
    """
    Checking if rh (relative humidity) is within bounds:
    
    1. Mean RH in the bottom 100 m (by gpsalt) should not be lesser than srf_limit, 
        set to a default of 50 %.
    
    If the above limit is violated, the rh for the sonde is considered out of bounds,
    and marked as 'False'. The sonde is also marked 'False', if there are no measurements in the 
    bottom 100 m (by gpsalt). 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if rh value is within bounds
               False, if rh is out of bounds or if measurements are not available
    """
    if sonde.rh.where(sonde.gpsalt < 100, drop=True).sum() == 0:
        return False
    elif sonde.rh.where(sonde.gpsalt < 100, drop=True).mean() < srf_limit:
        return False
    else:
        return True


def palt_gpsalt_rms_check(sonde, rms_limit=100):
    """
    This function estimates the root mean square (RMS) difference between geopotential altitude (palt) 
    and the GPS altitude (gpsalt),for values below 4 km, and based on a limit (rms_limit; 
    set to a default value of 100 m), is flagged accordingly.
    
    If the estimated RMS difference is below the limit, then the sonde is flagged as 'True' for this test.
    
    If the estimated RMS difference is greater than the limit, or if there are no values of either palt or gpsalt
    overlapping in the lower 4 km, then the sonde is flagged as 'False' for this test. The lack of overlap could be 
    because either there are no palt values or no gpsalt values or both.
    
    If the above limit is violated, the rh for the sonde is considered out of bounds,
    and marked as 'False'. The sonde is also marked 'False', if there are no measurements in the 
    bottom 100 m (by gpsalt). 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if RMS below limit
               False, if RMS above limit or no overlapping values
    """

    x = (
        (
            sonde.alt.where(sonde.alt < 4000, drop=True)
            - sonde.gpsalt.where(sonde.alt < 4000, drop=True)
        )
        ** 2
    ).values

    nn = ~np.isnan(x)

    if nn.sum() == 0:
        return False  # all x are NaNs
    else:
        zdiff_rms = np.sqrt(np.nanmean(x))
        if zdiff_rms < rms_limit:
            return True
        else:
            return False


# Function to check if geopotential altitude estimation by ASPEN failed (for whatever reason)


def check_alt_values(sonde):
    """
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if at least 50 values for 'alt' are available
               False, if less than 50 values for 'alt' are available
    """


# Function to check if sonde failed due to no detection of launch


def check_launch_detect(sonde_path):
    """
    Alternative method to check automatic launch detection of the sonde
    
    Input : path to PQC file
    
    Output : If launch detected, True; else, False
    
    Function to check if the dropsonde detected no launch, and thus failed. This was a common
    cause of dropsonde failure during EUREC4A. Vaisala's best guess is that for some reason,
    the IR sensor near the parachute did not detect the parachute coming out, thus not 
    detecting a launch and thus, not switching from low-power transmission to high-power 
    transmission. This caused AVAPS to lose the sonde's signal very soon after 
    launch (~350 hPa).
    """

    xrdataset = xr.open_dataset(sonde_path)

    if str(xrdataset.reference_time.values)[-21:-2] == "T00:00:00.000000000":
        launch_detect_flag = (
            False  # this means the sonde failed and launch was not detected
        )
        print(
            xrdataset.attrs["SoundingDescription"][21:30],
            xrdataset.reference_time.values,
            xrdataset.launch_time.values,
        )
    else:
        launch_detect_flag = True  # this means launch was detected

    return launch_detect_flag


def create_variable(ds, vname, data, **kwargs):
    """Insert the data into a variable in an :class:`xr.Dataset`"""
    attrs = dicts.nc_meta[vname].copy()
    dims = ["time"]  # nc_dims[vname]

    v = xr.Variable(dims, np.asarray(data), attrs=attrs)
    ds[vname] = v

    return vname


def get_status_ds_for_platform(Platform):

    (
        sonde_ds,
        directory,
        a_dir,
        logs_directory,
        a_files,
        file_time,
        sonde_paths,
    ) = get_all_sondes_list(Platform)

    # Retrieving all non NaN index sums in to a list for all sondes
    list_nc = list(map(get_total_non_nan_indices, sonde_ds))

    list_of_variables, s_time, s_t, s_rh, s_p, s_z, s_u, s_v = get_var_count_sums(
        list_nc
    )

    ld_FLAG = get_ld_flag_from_a_files(a_dir, a_files, logs_directory, Platform)

    status_ds = init_status_ds(
        list_of_variables, s_time, s_t, s_rh, s_p, s_z, s_u, s_v, ld_FLAG, file_time
    )

    status_ds, ind_flag_vars = add_ind_flags_to_statusds(status_ds, list_of_variables)
    status_ds, srf_flag_vars = add_srf_flags_to_statusds(status_ds, sonde_paths)
    status_ds, ind_FLAG = get_the_ind_FLAG_to_statusds(status_ds, ind_flag_vars)
    status_ds, srf_FLAG = get_the_srf_FLAG_to_statusds(status_ds, srf_flag_vars)
    status_ds = get_the_FLAG(status_ds, ind_FLAG, srf_FLAG)

    status_ds.to_netcdf(
        f"{logs_directory}Status_of_sondes_{Platform}_v{joanne.__version__}.nc"
    )

    return status_ds


# %%
