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

# %%

all_cir = xr.concat(circles, dim="circle")

# %%
for par in tqdm(["u", "v", "q", "ta", "p"]):

    mean_var_name = par
    var_dx_name = "d" + par + "dx"
    var_dy_name = "d" + par + "dy"

    (
        all_cir[mean_var_name],
        all_cir[var_dx_name],
        all_cir[var_dy_name],
        all_cir[par + "_sounding"],
    ) = rf.fit2d_xr(all_cir.dx, all_cir.dy, all_cir[par], ["sounding"], ["sounding"],)

lv4_dataset = rf.get_circle_products(all_cir)
# %%

lv4_dataset = lv4_dataset.drop("sounding")
# important to remove launch_time as dim and duplicate sounding variable

# %%
# nc_data = {}

# for var in dicts.list_of_vars:
#     if var not in ["platform", "segment_id", "sonde_id", "circle_time"]:
#         nc_data[var] = np.float32(lv4_dataset[var].values)

#     elif var in ["platform", "segment_id", "sonde_id"]:
#         nc_data[var] = lv4_dataset[var].values

#     elif var == "circle_time":
#         nc_data[var] = lv4_dataset[var].values.astype(float) / 1e9

#     else:
#         nc_data[var] = lv4_dataset[var].values
nc_data = {}

for var in dicts.list_of_vars:
    if lv4_dataset[var].values.dtype == "float64":
        nc_data[var] = np.float32(lv4_dataset[var].values)
    else:
        nc_data[var] = lv4_dataset[var].values
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

save_directory = "/Users/geet/Documents/JOANNE/Data/Level_4/"

comp = dict(zlib=True, complevel=4, fletcher32=True, _FillValue=np.finfo("float32").max)

encoding = {}

encoding = {
    var: comp
    for var in to_save_ds.data_vars
    if var not in ["platform", "segment_id", "sonde_id"]
}

encoding["circle_time"] = {"units": "seconds since 2020-01-01"}

# to_save_ds.encoding["circle_time"] = {
#     "units": "seconds since 2020-01-01",
# "dtype" :'datetime64'
# }

for key in dicts.nc_global_attrs.keys():
    to_save_ds.attrs[key] = dicts.nc_global_attrs[key]

to_save_ds.to_netcdf(
    save_directory + file_name, mode="w", format="NETCDF4", encoding=encoding
)

# %%
