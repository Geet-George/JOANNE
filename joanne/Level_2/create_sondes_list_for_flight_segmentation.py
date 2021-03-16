# %%

import xarray as xr
import numpy as np
from joanne.Level_2 import fn_2 as f2
import pandas as pd
import yaml
from datetime import datetime
from packaging import version
import glob
import joanne

# %%

dic = []
status_dir = "../../Data/Level_2/logs_and_stats/"

for Platform in ["P3", "HALO"]:
    # look for status file with same major and minor version-bit
    # (patch number and modifiers can be different)

    status_filename = glob.glob(
        f"{status_dir}Status_of_sondes_{Platform}_v{joanne.__version__[:3]}*.nc"
    )

    vers = [None] * len(status_filename)

    for n, i in enumerate(status_filename):
        vers[n] = version.parse(i)

    status_ds = xr.open_dataset(str(max(vers)))

    for i in range(len(status_ds.launch_time)):
        dic.append(
            {
                "flag": str(status_ds.FLAG.isel(launch_time=i).values),
                "launch_time": pd.to_datetime(
                    status_ds.launch_time.isel(launch_time=i).values
                ).to_pydatetime(),
                "platform": str(status_ds.platform.isel(launch_time=i).values),
                "sonde_id": str(status_ds.sonde_id.isel(launch_time=i).values),
            }
        )
# %%
filename = "sondes.yaml"
with open(filename, "w") as outfile:
    yaml.dump(dic, outfile)
    #   default_flow_style=False, sort_keys=False)

# %%
# all_P3_sondes, _, _, _, _, file_time, _, = f2.get_all_sondes_list("P3")

# # %%
# launch_time = [None] * len(all_P3_sondes)

# for i in range(len(all_P3_sondes)):
#     launch_time[i] = min(all_P3_sondes[i].time.values)

# # For all Level-1 QC files, the launch_time is indeed correctly obtained regardless of the A-file status. For sondes with the ld_FLAG as 0 (sondes which did not detect a launch) also have "correct" launch times, despite the A-file showing no launch detected. The filename however will not be the same as the launch_time for the BAD sondes. In any case of conflict, the launch_time of a sonde must be taken as that indicated by the variable launch_time in the file, and not that indicated in the filename.

# # %%

# sonde_id = [None] * len(launch_time)
# platform = [None] * len(launch_time)
# flight_id = [None] * len(launch_time)

# months = list(
#     pd.DatetimeIndex(launch_time).sort_values().month.astype(str).str.zfill(2)
# )
# days = list(pd.DatetimeIndex(launch_time).sort_values().day.astype(str).str.zfill(2))

# flight_id = [months[x] + days[x] for x in range(len(months))]

# for i in range(len(flight_id)):
#     if i == 0:
#         cntr = 1

#     elif flight_id[i] == flight_id[i - 1]:
#         cntr += 1
#     elif flight_id[i] != flight_id[i - 1]:
#         cntr = 1

#     sonde_id[i] = "P3-" + flight_id[i] + "_s" + str(cntr).zfill(2)
#     platform[i] = "P3"

# status_ds["sonde_id"] = (["time"], sonde_id)
# status_ds["platform"] = (["time"], platform)

# # %%
# #### --- UNCOMMENT FROM HERE TO CREATE YAML AND NC FOR LIST FOR FLIGHT SEGMENTATION --- ####
# status_ds["launch_time"] = (["time"], launch_time)
# # WITH ORDER OF launch_time SORTED, the 'time' and 'launch_time' are no longer the same
# status_ds = status_ds.swap_dims({"time": "launch_time"}).reset_coords()

# # %%

# dropsondeslist_for_flightsegmentation = xr.concat(status_ds.FLAG, status_ds.sonde_id)
# dropsondeslist_for_flightsegmentation.to_netcdf("ds_list.nc")

# # %%
# print(f"launch time is {min(all_HALO_sondes[44].time.values)} ")
# print(f"reference time is {all_HALO_sondes[44].reference_time.values}")
# print(f"file time is {status_ds.time[44].values}")

# # %%
# print(f"launch time is {min(all_HALO_sondes[43].time.values)} ")
# print(f"reference time is {all_HALO_sondes[43].reference_time.values}")
# print(f"file time is {status_ds.time[43].values}")

# %%

# dic = []

# for i in range(len(status_ds.launch_time)):
#     dic.append(
#         {
#             "flag": str(status_ds.FLAG.isel(launch_time=i).values),
#             "launch_time": pd.to_datetime(
#                 status_ds.launch_time.isel(launch_time=i).values
#             ).to_pydatetime(),
#             "platform": str(status_ds.platform.isel(launch_time=i).values),
#             "sonde_id": str(status_ds.sonde_id.isel(launch_time=i).values),
#         }
#     )
# # %%
# filename = "sondes.yaml"
# with open(filename, "w") as outfile:
#     yaml.dump(dic, outfile)
#     #   default_flow_style=False, sort_keys=False)

# %%
