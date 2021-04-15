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

for Platform in [
    "HALO",
    "P3",
]:
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
