# %%
import datetime
import glob
import sys
import warnings
from importlib import reload

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
import xarray as xr
from seaborn import distplot
from tqdm import tqdm

import joanne
from joanne.Level_2 import fn_2 as f2
from joanne.Level_2 import dicts

reload(f2)
reload(dicts)

# %%

varname_L1 = ["height", "time", "wspd", "wdir", "tdry", "pres", "rh", "lat", "lon"]
varname_L2 = ["alt", "time", "wspd", "wdir", "ta", "p", "rh", "lat", "lon"]

for Platform in ["HALO", "P3"]:

    (
        sonde_ds,
        directory,
        a_dir,
        logs_directory,
        a_files,
        file_time,
        sonde_paths,
    ) = f2.get_all_sondes_list(Platform)

    # look for status file with same major and minor version-bit
    # (patch number and modifiers can be different)

    status_filename = glob.glob(
        f"{logs_directory}Status_of_sondes_{Platform}_v{joanne.__version__[:3]}*.nc"
    )

    status_ds = xr.open_dataset(status_filename[0])

    launch_time = [None] * len(sonde_ds)

    for i in range(len(sonde_ds)):
        launch_time[i] = min(sonde_ds[i].time.values)

    sonde_id = [None] * len(launch_time)
    platform = [None] * len(launch_time)
    flight_id = [None] * len(launch_time)

    months = list(pd.DatetimeIndex(launch_time).month.astype(str).str.zfill(2))
    days = list(pd.DatetimeIndex(launch_time).day.astype(str).str.zfill(2))

    flight_id = [months[x] + days[x] for x in range(len(months))]

    for i in range(len(flight_id)):
        if i == 0:
            cntr = 1

        elif flight_id[i] == flight_id[i - 1]:
            cntr += 1
        elif flight_id[i] != flight_id[i - 1]:
            cntr = 1

        sonde_id[i] = Platform + "-" + flight_id[i] + "_s" + str(cntr).zfill(2)
    # platform[i] = "HALO"

    a_filepaths = []

    for i in a_files:
        a_filepaths.append(sorted(glob.glob(a_dir + i + "*")))

    file_time_str = [None] * len(sonde_paths)

    for i in range(len(sonde_paths)):
        # file_time_str[i] = sonde_paths[i][-22:-7]
        lt_pd = pd.to_datetime(min(sonde_ds[i].time.values))

        file_time_str[
            i
        ] = f"{(lt_pd.year)}{(lt_pd.month):02}{(lt_pd.day):02}_{(lt_pd.hour):02}{(lt_pd.minute):02}{(lt_pd.second):02}"

    for i in tqdm(range(len(sonde_ds))):

        if status_ds.FLAG[i] == "GOOD":  # or status_ds.FLAG[i] == "UGLY":

            # ht_indices = ~np.isnan(sonde_ds[i].alt)
            ht_indices = (
                ~np.isnan(sonde_ds[i].alt)
                & ~np.isnan(sonde_ds[i].lat)
                & ~np.isnan(sonde_ds[i].lon)
            )
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

            time = sonde_ds[i].time[ht_indices].values  # .astype("float").values / 1e9
            # Variable array: time

            variables = {}

            variables["time"] = time
            variables["alt"] = height

            ###--------- Unit Conversions --------###

            variables["rh"] = np.float32(sonde_ds[i]["rh"][ht_indices].values / 100)
            variables["lat"] = np.float32(sonde_ds[i]["lat"][ht_indices].values)
            variables["lon"] = np.float32(sonde_ds[i]["lon"][ht_indices].values)
            variables["p"] = np.float32(sonde_ds[i]["pres"][ht_indices].values * 100)
            variables["ta"] = np.float32(
                sonde_ds[i]["tdry"][ht_indices].values + 273.15
            )
            # variables["sonde_id"] = sonde_id[i]

            for var1, var2 in zip(varname_L1, varname_L2):
                if var2 not in variables.keys():
                    variables[var2] = np.float32(sonde_ds[i][var1][ht_indices].values)

            ###--------- Creating and populating dataset --------###

            to_save_ds = xr.Dataset(coords={"time": obs})

            for var in dicts.nc_meta.keys():
                # v = var
                f2.create_variable(to_save_ds, var, variables[var])

            attrs = {
                "descripion": "unique sonde ID in the format PLATFORM_FLIGHT-ID_sSONDE-NUMBER-FOR-THE-FLIGHT",
                "long_name": "sonde identifier",
                "cf_role": "trajectory_id",
            }
            sonde_id_var = xr.Variable([], sonde_id[i], attrs=attrs)
            to_save_ds["sonde_id"] = sonde_id_var

            # to_save_ds["sonde_id"] = ([],sonde_id[i],attrs={"descripion":"unique sonde ID in the format PLATFORM_FLIGHT-ID_sSONDE-NUMBER-FOR-THE-FLIGHT"})#,"long_name": "sonde identifier","cf_role": "trajectory_id",})

            file_name = (
                "EUREC4A_JOANNE_"
                + str(Platform)
                + "_Dropsonde-RD41_"
                + file_time_str[i]
                + "_v"
                + str(joanne.__version__)
                + ".nc"
            )
            save_directory = "/Users/geet/Documents/JOANNE/Data/Level_2/"

            comp = dict(
                zlib=True,
                complevel=4,
                fletcher32=True,
                _FillValue=np.finfo("float32").max,
            )

            encoding = {var: comp for var in to_save_ds.data_vars if var != "sonde_id"}
            encoding["time"] = {"units": "seconds since 2020-01-01", "dtype": "float"}

            nc_global_attrs = dicts.get_global_attrs(
                Platform, file_time[i], sonde_ds[i]
            )

            for key in nc_global_attrs.keys():
                to_save_ds.attrs[key] = nc_global_attrs[key]

            flight_attrs = dicts.get_flight_attrs(a_filepaths[i][0])

            for key in flight_attrs:
                to_save_ds.attrs[key] = flight_attrs[key]

            ###--------- Saving dataset to NetCDF file --------###

            to_save_ds.to_netcdf(
                save_directory + file_name,
                mode="w",
                format="NETCDF4",
                encoding=encoding,
            )
    # %%
