# %%
import datetime
# import glob
import subprocess
import warnings

import numpy as np
# import requests
import xarray as xr

import fn_3a as f3a
import dicts

warnings.filterwarnings(
    "ignore", module="metpy.calc.thermo", message="invalid value encountered"
)
# %%
try:
    git_module_version = (
        subprocess.check_output(["git", "describe"]).strip().decode("utf-8")
    )
except:
    git_module_version = "--"

# %%
lv2_data_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Data/code_testing_data/"  # Level_2/"  # code_testing_data/"

lv3_dataset = f3a.lv3_structure_from_lv2(lv2_data_directory)

# %%

nc_data = {}

for var in dicts.list_of_vars:
    nc_data[var] = lv3_dataset[var].values

# %%

obs = lv3_dataset.obs.values
sounding = lv3_dataset.sounding.values

to_save_ds = xr.Dataset(coords={"obs": obs, "sounding": sounding})

for var in dicts.list_of_vars:
    f3a.create_variable(
        to_save_ds, var, data=nc_data, dims=dicts.nc_dims, attrs=dicts.nc_attrs
    )

file_name = (
    "EUREC4A" + "_Dropsonde-RD41_" + "Level_3A_" + str(git_module_version) + ".nc"
)

save_directory = "/Users/geet/Documents/EUREC4A/JOANNE/Data/Level_3/Test_data/"

comp = dict(zlib=True, complevel=4, fletcher32=True, _FillValue=np.finfo("float32").max)

encoding = {var: comp for var in to_save_ds.data_vars if var != "Platform"}

for key in dicts.nc_global_attrs.keys():
    to_save_ds.attrs[key] = dicts.nc_global_attrs[key]

to_save_ds.to_netcdf(
    save_directory + file_name, mode="w", format="NETCDF4", encoding=encoding
)