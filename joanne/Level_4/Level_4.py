# %%
import glob
import time
from datetime import date
from importlib import reload

import numpy as np
import xarray as xr
import yaml

import joanne
from tqdm import tqdm

from joanne.Level_4 import rgr_fn as rf
from joanne.Level_4 import ready_ds_for_regression as prep
from joanne.Level_4 import dicts

reload(prep)
reload(rf)
reload(dicts)
# %%
circles = prep.get_circles()
prep.get_xy_coords_for_circles(circles)

all_cir = xr.concat(circles, dim="circle")

### TEMPORARY FILE TO GET PREPPED XY COORDINATES OF CIRCLES ###
# all_cir = xr.open_dataset("all_cir_v0.8.1+7.g5.nc")

# %%
for par in tqdm(["u", "v", "q", "ta", "p"]):

    mean_var_name = par
    var_dx_name = "d" + par + "dx"
    var_dy_name = "d" + par + "dy"

    all_cir[mean_var_name], all_cir[var_dx_name], all_cir[var_dy_name] = rf.fit2d_xr(
        all_cir.dx, all_cir.dy, all_cir[par], "launch_time"
    )

    # all_cir[mean_var_name] = (["circle", "alt"], intercept)
    # all_cir[var_dx_name] = (["circle", "alt"], dpardx)
    # all_cir[var_dy_name] = (["circle", "alt"], dpardy)

lv4_dataset = rf.get_circle_products(all_cir)

lv4_dataset = prep.reswap_launchtime_sounding(lv4_dataset)

# %%
# start = time.perf_counter()

# list_of_parameters = [
#     "u",
#     "v",
#     "q",
#     "ta",
#     "p",
# ]

# rf.get_circle_products(circles, list_of_parameters)

# circles = prep.reswap_launchtime_sounding(circles)

# finish = time.perf_counter()

# print(f"Finished in {round(finish - start,2)} seconds ...")

# %%

# lv4_dataset = xr.concat(circles, dim="circle")

nc_data = {}

for var in dicts.list_of_vars:
    if (var != "platform") and (var != "segment_id"):
        nc_data[var] = np.float32(lv4_dataset[var].values)

    if (var == "platform") or (var == "segment_id"):
        nc_data[var] = lv4_dataset[var].values

    if (var == "launch_time") or (var == "circle_time"):
        nc_data[var] = np.float32(lv4_dataset[var].astype("float").values / 1e9)

# %%

alt = lv4_dataset.alt.values
sounding = lv4_dataset.sounding.values
circle = lv4_dataset.circle.values

to_save_ds = xr.Dataset(coords={"alt": alt, "sounding": sounding, "circle": circle})

for var in dicts.list_of_vars:
    prep.create_variable(
        to_save_ds, var, data=nc_data, dims=dicts.nc_dims, attrs=dicts.nc_attrs
    )

file_name = (
    "EUREC4A_JOANNE_Dropsonde-RD41_" + "Level_4_v" + str(joanne.__version__) + ".nc"
)

save_directory = "/Users/geet/Documents/JOANNE/Data/Level_4/"  # Test_data/" #Level_3/"

comp = dict(zlib=True, complevel=4, fletcher32=True, _FillValue=np.finfo("float32").max)

encoding = {}

encoding = {
    var: comp for var in to_save_ds.data_vars if var not in ["platform", "segment_id"]
}

for key in dicts.nc_global_attrs.keys():
    to_save_ds.attrs[key] = dicts.nc_global_attrs[key]

to_save_ds.to_netcdf(
    save_directory + file_name, mode="w", format="NETCDF4", encoding=encoding
)

# %%
