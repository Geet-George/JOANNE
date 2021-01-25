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

reload(f2)

warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice encounter")
warnings.filterwarnings(
    "ignore", message="Attempted to set non-positive bottom ylim on a log-scaled axis"
)

# %%
###------ Platform Name ------###

for Platform in ["P3", "HALO"]:

    print(f"{Platform} running now...")

    directory = "/Users/geet/Documents/JOANNE/Data/Level_1/" + Platform + "/"
    # directory where all sonde files are present

    a_dir = "/Users/geet/Documents/JOANNE/Data/Level_0/" + Platform + "/All_A_files/"
    # directory where all the A files are present

    logs_directory = "/Users/geet/Documents/JOANNE/Data/Level_2/logs_and_stats/"
    # directory to store logs and stats

    status_ds = f2.get_status_ds_for_platform(Platform)

    file = open(
        f"{logs_directory}no_launch_detect_logs_{Platform}_v{joanne.__version__}.txt",
        "a",
    )

    status_dict = {}

    for x in ["ind", "srf", "FLAG"]:

        for y in ["GOOD", "UGLY", "BAD"]:

            if x != "FLAG":
                status_dict[f"{y}_{x}"] = len(
                    status_ds.where(status_ds[f"{x}_FLAG"] == y, drop=True).time
                )
            else:
                status_dict[f"{y}_{x}"] = len(
                    status_ds.where(status_ds[x] == y, drop=True).time
                )

    file.write("----------------------------------------------\n")
    file.write(
        f"As per the ind_FLAG tests,\n{status_dict['GOOD_ind']} are good sondes,\n"
    )
    file.write(
        f"{status_dict['BAD_ind']} are bad sondes\nand {status_dict['UGLY_ind']} are ugly sondes.\n"
    )
    file.write("----------------------------------------------\n")
    file.write(
        f"As per the srf_FLAG tests,\n{status_dict['GOOD_srf']} are good sondes,\n"
    )
    file.write(
        f"{status_dict['BAD_srf']} are bad sondes\nand {status_dict['UGLY_srf']} are ugly sondes.\n"
    )
    file.write("----------------------------------------------\n")
    file.write(f"There are a total of {len(status_ds.time)} sondes\n")
    file.write(f"out of which {status_dict['GOOD_FLAG']} are good sondes,\n")
    file.write(
        f"{status_dict['BAD_FLAG']} are bad sondes\nand {status_dict['UGLY_FLAG']} are ugly sondes that can be salvaged with some effort.\n"
    )

    file.close()

