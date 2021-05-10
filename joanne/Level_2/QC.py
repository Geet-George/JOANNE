# %%
import datetime
import glob
import sys
import warnings
from importlib import reload
import argparse

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
parser = argparse.ArgumentParser(
    description="This script takes .QC or .PQC files processed by ASPEN and runs the JOANNE QC tests on them. It outputs a status file that provides the success of the provided soundings in the JOANNE QC tests. There is an additional option to print out the sonde IDs which failed to detect a launch, along with a summary of the soundings' performance in the QC test by setting the argument 'logs' to True"
)

parser.add_argument(
    "-l",
    "--logs",
    help="Set this as True if a summary of launch detect failures is to be printed along with a summary of soundings' performance in JOANNE QC tests. This is set as False by default.",
    type=str,
)

args = parser.parse_args()

# %%

logs_directory = "/Users/geet/Documents/JOANNE/Data/Level_2/logs_and_stats/"
# directory to store logs and stats
data_directory = "/Users/geet/Documents/JOANNE/Data/"
save_directory = "/Users/geet/Documents/JOANNE/Data/QC/"
# %%


def run_qc(
    logs_directory=logs_directory,
    data_directory=data_directory,
    save_directory=save_directory,
):
    ###------ Platform Name ------###

    for Platform in ["P3", "HALO"]:

        print(f"{Platform} running now...")

        directory = f"{data_directory}Level_1/{Platform}/"
        # directory where all sonde files are present

        a_dir = f"{data_directory}Level_0/{Platform}/All_A_files/"
        # directory where all the A files are present

        status_ds = f2.get_status_ds_for_platform(Platform, save_dir=save_directory)

    return print("JOANNE QC tests finished running")


def create_QC_summary_logs(
    logs_directory=logs_directory, data_directory=data_directory,save_directory=save_directory,
):

    for Platform in ["P3", "HALO"]:

        file = open(
            f"{logs_directory}summary_logs_{Platform}_v{joanne.__version__}.txt", "w",
        )

        status_ds = f2.get_status_ds_for_platform(Platform, save_dir=save_directory)

        status_dict = {}

        for x in ["sat", "low", "qc_flag"]:

            for y in ["GOOD", "UGLY", "BAD"]:

                if x != "qc_flag":
                    status_dict[f"{y}_{x}"] = len(
                        status_ds.where(
                            status_ds[f"{x}_test"] == y, drop=True
                        ).launch_time
                    )
                else:
                    status_dict[f"{y}_{x}"] = len(
                        status_ds.where(status_ds[x] == y, drop=True).launch_time
                    )

        file.write("----------------------------------------------\n")
        file.write(
            f"As per the sat_test tests,\n{status_dict['GOOD_sat']} are good sondes,\n"
        )
        file.write(
            f"{status_dict['BAD_sat']} are bad sondes\nand {status_dict['UGLY_sat']} are ugly sondes.\n"
        )
        file.write("----------------------------------------------\n")
        file.write(
            f"As per the low_test tests,\n{status_dict['GOOD_low']} are good sondes,\n"
        )
        file.write(
            f"{status_dict['BAD_low']} are bad sondes\nand {status_dict['UGLY_low']} are ugly sondes.\n"
        )
        file.write("----------------------------------------------\n")
        file.write(f"There are a total of {len(status_ds.launch_time)} sondes\n")
        file.write(f"out of which {status_dict['GOOD_qc_flag']} are good sondes,\n")
        file.write(
            f"{status_dict['BAD_qc_flag']} are bad sondes\nand {status_dict['UGLY_qc_flag']} are ugly sondes that can be salvaged with some effort.\n"
        )

        file.close()

    return print("JOANNE QC summary logs finished printing")


if __name__ == "__main__":
    run_qc()
    if args.logs:
        create_QC_summary_logs()
# %%
