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

    status_ds = xr.open_dataset(
        f"{logs_directory}Status_of_sondes_{Platform}_v0.5.1+1.g0c163c5.dirty.nc"
    )

    a_filepaths = []

    for i in a_files:
        a_filepaths.append(sorted(glob.glob(a_dir + i + "*")))

    file_time_str = [None] * len(sonde_paths)

    for i in range(len(sonde_paths)):
        file_time_str[i] = sonde_paths[i][-22:-7]

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

            variables = {}

            variables["time"] = time
            variables["alt"] = height

            for var1, var2 in zip(varname_L1, varname_L2):
                if var2 not in variables.keys():
                    variables[var2] = np.float32(sonde_ds[i][var1][ht_indices].values)

            ###--------- Creating and populating dataset --------###

            to_save_ds = xr.Dataset(coords={"time": obs})

            for var in dicts.nc_meta.keys():
                # v = var
                f2.create_variable(to_save_ds, var, variables[var])

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

            encoding = {var: comp for var in to_save_ds.data_vars}

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
