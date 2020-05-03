# %%
import glob
import time
from datetime import date
from importlib import reload

import numpy as np
import xarray as xr
import yaml

import rgr_fn as rf
import ready_ds_for_regression as prep
import dicts

reload(prep)
reload(rf)
reload(dicts)
# %%
circles = prep.get_circles()
prep.get_xy_coords_for_circles(circles)

# %%
start = time.perf_counter()

list_of_parameters = [
    "u_wind",
    "v_wind",
    "specific_humidity",
    "temperature",
    "pressure",
]

rf.get_circle_products(circles, list_of_parameters)

circles = prep.reswap_launchtime_sounding(circles)

finish = time.perf_counter()

print(f"Finished in {round(finish - start,2)} seconds ...")
# %%


# %%

lv3b_dataset = xr.concat(circles,dim='circle')

# %%

nc_data = {}

for var in dicts.list_of_vars:
    nc_data[var] = lv3_dataset[var].values

# %%
for i in flight_date:
    for j in range(1, 7):
        circle.append((f"{i.astype(object).month:02d}{i.astype(object).day:02d}C{j}"))

ds = xr.Dataset(
    {
        "dudx": (["circle", "height"], m_u),  # {'ad':'dfw','qd':'df'}),
        "dudy": (["circle", "height"], c_u),
        "dvdx": (["circle", "height"], m_v),
        "dvdy": (["circle", "height"], c_v),
        "dTdx": (["circle", "height"], m_T),
        "dTdy": (["circle", "height"], c_T),
        "dqdx": (["circle", "height"], m_q),
        "dqdy": (["circle", "height"], c_q),
        "dqudx": (["circle", "height"], m_qu),
        "dqudy": (["circle", "height"], c_qu),
        "dqvdx": (["circle", "height"], m_qv),
        "dqvdy": (["circle", "height"], c_qv),
        "dTudx": (["circle", "height"], m_Tu),
        "dTudy": (["circle", "height"], c_Tu),
        "dTvdx": (["circle", "height"], m_Tv),
        "dTvdy": (["circle", "height"], c_Tv),
        "pressure": (["circle", "height"], pres),
        "divergence": (["circle", "height"], D),
        "vorticity": (["circle", "height"], vor),
        "u": (["circle", "height"], u),
        "v": (["circle", "height"], v),
        "w": (["circle", "height"], w_vel),
        "omega": (["circle", "height"], p_vel),
        "density": (["circle", "height"], mean_den),
        "regressed_sondes": (["circle", "height"], Ns),
    },
    coords={"height": circles[0].height.values, "circle": circle},
    attrs={
        "Author": "Geet George (MPI-M, Hamburg); geet.george@mpimet.mpg.de",
        "Instrument": "Vaisala RD41 (AVAPS receiver aboard HALO)",
        "Data Processing": 'AvapsEditorVersion "BatchAspen V3.4.3"',
        "Campaign": "EUREC4A (Jan-Feb,2020)",
        "Reference Study": "Bony, S. & Stevens, B., 2019, “Measuring Area-Averaged Vertical Motions with Dropsondes.” Journal of the Atmospheric Sciences 76 (3): 767–83",
        "Creation Time": str(datetime.datetime.utcnow()) + " UTC",
    },
)

file_name = "circle_products"
# file_name = str(datetime.datetime.utcnow().date()) + '_' + str(datetime.datetime.utcnow().hour)
# file_name = str(circle['launch_time'].mean().values)[0:13]
ds.to_netcdf(file_name + ".nc")
