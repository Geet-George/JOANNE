# %%
import datetime
import os
import subprocess
import warnings
from importlib import reload

import numpy as np
import xarray as xr

# import dicts
from joanne.Level_3 import fn_3 as f3
from joanne.Level_3 import dicts as dicts
import joanne

warnings.filterwarnings(
    "ignore", module="metpy.calc.thermo", message="invalid value encountered"
)

reload(dicts)
reload(f3)
reload(joanne)
# %%
lv2_data_directory = "/Users/geet/Documents/JOANNE/Data/code_testing_data/"  # Level_2/"  # code_testing_data/"

lv3_dataset = f3.lv3_structure_from_lv2(lv2_data_directory)

# %%
nc_data = {}

for var in dicts.list_of_vars:
    if lv3_dataset[var].values.dtype == "float64":
        nc_data[var] = np.float32(lv3_dataset[var].values)
    else:
        nc_data[var] = lv3_dataset[var].values
# %%

obs = np.arange(0, len(lv3_dataset.alt), 1)
sounding = lv3_dataset.sounding.values

to_save_ds = xr.Dataset(coords={"alt": obs, "sounding": sounding})

for var in dicts.list_of_vars:
    f3.create_variable(
        to_save_ds, var, data=nc_data, dims=dicts.nc_dims, attrs=dicts.nc_attrs
    )

file_name = (
    "EUREC4A_JOANNE_Dropsonde-RD41_" + "Level_3_v" + str(joanne.__version__) + ".nc"
)

save_directory = (
    "/Users/geet/Documents/JOANNE/Data/Level_3/Test_data/"  # Test_data/" #Level_3/"
)

comp = dict(zlib=True, complevel=4, fletcher32=True, _FillValue=np.finfo("float32").max)

encoding = {
    var: comp for var in to_save_ds.data_vars if var not in ["platform", "sonde_id"]
}

for key in dicts.nc_global_attrs.keys():
    to_save_ds.attrs[key] = dicts.nc_global_attrs[key]

to_save_ds.to_netcdf(
    save_directory + file_name, mode="w", format="NETCDF4", encoding=encoding
)
# %%


# %%
