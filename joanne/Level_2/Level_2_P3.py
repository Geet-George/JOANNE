# %%
from importlib import reload
import datetime
import glob
import sys
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
import xarray as xr
from seaborn import distplot
from tqdm import tqdm

import joanne

from joanne.Level_2 import fn_2 as f2

reload(f2)

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice encounter")
warnings.filterwarnings(
    "ignore", message="Attempted to set non-positive bottom ylim on a log-scaled axis"
)

# %%
###------ Platform Name ------###
Platform_Name = "P3"

# %%
directory = (
    "/Users/geet/Documents/EUREC4A/Dropsondes/" + Platform_Name + "/All_dropsondes/"
)
# directory where all sonde files are present

a_dir = (
    "/Users/geet/Documents/EUREC4A/Dropsondes_raw/" + Platform_Name + "/All_A_files/"
)
# directory where all the A files are present

logs_directory = "/Users/geet/Documents/JOANNE/joanne/Level_2/logs_and_stats/"
# directory to store logs and stats

sonde_paths = sorted(glob.glob(directory + "*_PQC.nc"))
# paths to the individual sonde files

file_time_str = [None] * len(sonde_paths)
file_time = [None] * len(sonde_paths)
# list to store extracted sonde time from file name as string and as time

a_files = [None] * len(sonde_paths)
# list to store file names for the log files starting with A.

sonde_ds = [None] * len(sonde_paths)
# list to store individual datasets of all sondes from PQC files

for i in range(len(sonde_paths)):
    file_time_str[i] = sonde_paths[i][-22:-7]
    file_time[i] = np.datetime64(
        pd.to_datetime(file_time_str[i], format="%Y%m%d_%H%M%S"), "s"
    )
    a_files[i] = "A" + file_time_str[i]
    sonde_ds[i] = xr.open_dataset(sonde_paths[i])


# %%
# Retrieving all non NaN index sums in to a list for all sondes
list_nc = list(map(f2.get_total_non_nan_indices, sonde_ds))


# %%
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


# %%
a_filepaths = []
# list to store individual file paths for all A files

for i in a_files:
    a_filepaths.append(sorted(glob.glob(a_dir + i + "*")))

ld_FLAG = np.full(len(a_files), np.nan)
# array to store ld_FLAG values

# create and start writing a log file which will store sonde info about sondes with failed launch detection
file = open(
    f"{logs_directory}no_launch_detect_logs_{Platform_Name}_v{joanne.__version__}.txt",
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

# %%
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

# %%
# Plotting the distribution of ratio of parameter's non-NaN measurement counts to counts of time

f, ax = plt.subplots(2, 1, figsize=(10, 14))

# looping over all parameters
for var in list_of_variables[1:-1]:

    # calculating ratio
    rat = status_ds[var] / status_ds.s_time

    # plotting distribution on linear Y-axis
    sb.distplot(
        rat,
        hist=False,
        ax=ax[0],
        kde_kws={"linewidth": 3, "alpha": 0.5},  # "marker":'o',
        label=var[2:],
    )

    # plotting distribution on log Y-axis
    sb.distplot(
        rat,
        hist=False,
        ax=ax[1],
        kde_kws={"linewidth": 0.5, "marker": "o", "alpha": 0.5},
        label=var[2:],
    )

    # setting the second plot's Y-axis to log scale
    ax[1].set_yscale("log")
    plt.legend()
    ax[1].set_ylim(0.01, 100)
    ax[1].set_xlabel(r"$\dfrac{Count\ of\ Parameter}{Count\ of\ Time}$", fontsize=14)
    [ax[m].set_ylabel("Distribution / %", fontsize=14) for m in [0, 1]]
    [
        [ax[m].spines[n].set_visible(False) for m in [0, 1]]
        for n in ["left", "right", "top"]
    ]
    plt.suptitle(
        f"Distribution for all {int(len(status_ds.time))} sondes launched from {Platform_Name}",
        fontsize=20,
        verticalalignment="bottom",
    )

sb.despine(offset=10)
plt.savefig(
    f"{logs_directory}Count_of_measurements_{Platform_Name}_v{joanne.__version__}.png",
    dpi=300,
)
# %% [markdown]
# Time values are recorded every 0.25 seconds. Although, the PTU and GPS sensors have a measurement frequency of 2 Hz and 4 Hz, respectively, the distribution of measurements have a slightly different story. Based on the distribution, we know that the ideal case is for all parameters (except u,v) to have measurements at every other time record, and for u,v to have measurements at every time record. Since, the time records also include values during initialisation as well as during a little before and after the launch, when no signal can be sent back to the AVAPS PC, the actual ratio will always be lower than the ideal estimate of 1 (for u,v) and 0.5 (for the remaining parameters).
#
# The true distribution shows that peaks start to flatten around 0.8 and 0.4 for u,v and other parameters, respectively. Thus, sondes with ratios lower than these values are taken as not having a complete profile, and termed as 'ugly' sondes. These ugly sondes still have data, but because they are expected to have more NaN fields than most sondes, they are kept for more QC, NaN-filling later, depending upon the extent of the dearth of measurements in that sonde.
#
# P.S. It is still unclear to me why GPS values are only at every other time record, while u,v are always being estimated.

# %%
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


# %%
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


# %%
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


# %%
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
                    for item in status_ds[srf_flag_vars].sel(time=i).to_array().values
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


# %%
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


# %%
good_ind = len(status_ds.where(status_ds.ind_FLAG == "GOOD", drop=True).time)
ugly_ind = len(status_ds.where(status_ds.ind_FLAG == "UGLY", drop=True).time)
bad_ind = len(status_ds.where(status_ds.ind_FLAG == "BAD", drop=True).time)

file.write("----------------------------------------------\n")
file.write(f"As per the ind_FLAG tests,\n")
file.write(f"{good_ind} are good sondes,\n")
file.write(f"{bad_ind} are bad sondes\nand {ugly_ind} are ugly sondes.\n")


# %%
good_srf = len(status_ds.where(status_ds.srf_FLAG == "GOOD", drop=True).time)
ugly_srf = len(status_ds.where(status_ds.srf_FLAG == "UGLY", drop=True).time)
bad_srf = len(status_ds.where(status_ds.srf_FLAG == "BAD", drop=True).time)

file.write("----------------------------------------------\n")
file.write(f"As per the srf_FLAG tests,\n")
file.write(f"{good_srf} are good sondes,\n")
file.write(f"{bad_srf} are bad sondes\nand {ugly_srf} are ugly sondes.\n")


# %%
good_sondes = len(status_ds.where(status_ds.FLAG == "GOOD", drop=True).time)
ugly_sondes = len(status_ds.where(status_ds.FLAG == "UGLY", drop=True).time)
bad_sondes = len(status_ds.where(status_ds.FLAG == "BAD", drop=True).time)

file.write("----------------------------------------------\n")
file.write(f"There are a total of {len(status_ds.time)} sondes\n")
file.write(f"out of which {good_sondes} are good sondes,\n")
file.write(
    f"{bad_sondes} are bad sondes\nand {ugly_sondes} are ugly sondes that can be salvaged with some effort.\n"
)

file.close()

status_ds.to_netcdf(
    f"{logs_directory}Status_of_sondes_{Platform_Name}_v{joanne.__version__}.nc"
)
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
# %%
flight_attrs = [None] * len(a_filepaths)

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

# mission_pi = []

for i in tqdm(range(len(sonde_ds))):

    if status_ds.FLAG[i] == "GOOD":

        ht_indices = ~np.isnan(sonde_ds[i].alt)
        # retrieving non-NaN indices of geopotential height (sonde_ds[i].alt)
        # only time values at these indices will be used in Level-2 trajectory data;
        # this means that only alternate u,v values are included in the Level-2 data
        # PTU has 2 Hz measurement frequency, while GPS has a 4 Hz measurement frequency

        ###----- Dimensions -----###

        obs = np.arange(1, ht_indices.sum() + 1, 1)
        # creating the observations dimension of the NC file

        ###----- Variables -----###

        height = np.float32(sonde_ds[i].alt[ht_indices].values)
        # Variable array: geopotential height

        time = sonde_ds[i].time[ht_indices].astype("float").values / 1e9
        # Variable array: time

        wind_speed = np.float32(sonde_ds[i].wspd[ht_indices].values)
        # Variable array: wind speed

        wind_direction = np.float32(sonde_ds[i].wdir[ht_indices].values)
        # Variable array: wind direction

        temperature = np.float32(sonde_ds[i].tdry[ht_indices].values)
        # Variable array: temperature

        pressure = np.float32(sonde_ds[i].pres[ht_indices].values)
        # Variable array: pressure

        relative_humidity = np.float32(sonde_ds[i].rh[ht_indices].values)
        # Variable array: relative humidity

        latitude = np.float32(sonde_ds[i].lat[ht_indices].values)
        # Variable array: latitude

        longitude = np.float32(sonde_ds[i].lon[ht_indices].values)
        # Variable array: longitude

        ###----- Global Attributes -----###

        nc_global_attrs = {
            "Title": "Sounding data containing temperature, pressure, humidity,"
            " latitude, longitude, wind direction, wind speed, and time",
            "Campaign": "EUREC4A-ATOMIC (Jan-Feb, 2020)",
            "Platform": Platform_Name,
            "Instrument": "Vaisala RD41",
            "Launch-date": str(pd.to_datetime(file_time[i]).date()),
            "Launch-time-(UTC)": str(sonde_ds[i].launch_time.values),
            "Sonde-Serial-ID": sonde_ds[i].SondeId,
            "ASPEN-Version": sonde_ds[i].AspenVersion,
            "Processing-Time": sonde_ds[i].ProcessingTime,
            "Mission-PI": "Mission PI",
            "Author": "Geet George",
            "Author-Email": "geet.george@mpimet.mpg.de",
            "JOANNE-version": joanne.__version__,
            "Conventions": "CF-1.7",
            "featureType": "trajectory",
            "Creation-Time": str(datetime.datetime.utcnow()) + " UTC",
        }

        ###--------- Retrieving flight parameters during launch --------###

        flight_attrs[i] = {}

        file_path = a_filepaths[i][0]
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
                flight_attrs[i][attr] = lines[line_id].split("= ")[1]
            else:
                flight_attrs[i][attr] = float(lines[line_id].split("= ")[1])

        f.close()

        ###--------- Creating and populating dataset --------###

        to_save_ds = xr.Dataset(coords={"time": obs})

        for var in nc_meta.keys():
            # v = var
            create_variable(to_save_ds, var, eval(var))

        file_name = (
            "EUREC4A_JOANNE_"
            + str(Platform_Name)
            + "_Dropsonde-RD41_"
            + file_time_str[i]
            + "_v"
            + str(joanne.__version__)
            + ".nc"
        )
        save_directory = "/Users/geet/Documents/JOANNE/Data/Level_2/"

        comp = dict(
            zlib=True, complevel=4, fletcher32=True, _FillValue=np.finfo("float32").max
        )

        encoding = {var: comp for var in to_save_ds.data_vars}

        for key in nc_global_attrs.keys():
            to_save_ds.attrs[key] = nc_global_attrs[key]

        for key in flight_attrs[i]:
            to_save_ds.attrs[key] = flight_attrs[i][key]

        ###--------- Saving dataset to NetCDF file --------###

        to_save_ds.to_netcdf(
            save_directory + file_name, mode="w", format="NETCDF4", encoding=encoding
        )
# %%
