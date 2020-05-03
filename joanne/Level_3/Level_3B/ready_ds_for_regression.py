import yaml
import glob
import xarray as xr
import numpy as np

# from datetime import date
import datetime
from pylab import cos

# %%

yaml_directory = "/Users/geet/Documents/JOANNE/joanne/flight_segments/"
lv3_directory = "/Users/geet/Documents/JOANNE/Data/Level_3/"

lv3a_filename = "EUREC4A_Dropsonde-RD41_Level_3A.nc"


def get_level3a_dataset(lv3_directory=lv3_directory, lv3a_filename=lv3a_filename):
    return xr.open_dataset(lv3_directory + lv3a_filename)


def dim_ready_ds(ds_lv3a=get_level3a_dataset(), Platform="HALO"):

    dims_to_drop = ["obs", "sounding"]

    all_sondes = (
        ds_lv3a.where(ds_lv3a.Platform == Platform, drop=True)
        .swap_dims({"sounding": "launch_time"})
        .swap_dims({"obs": "height"})
        .drop(dims_to_drop)
    )

    return all_sondes


def get_circle_times_from_yaml(yaml_directory=yaml_directory):
    allyamlfiles = sorted(glob.glob(yaml_directory + "*.yaml"))

    circle_times = []
    flight_date = []

    for i in allyamlfiles:
        with open(i) as source:
            flightinfo = yaml.load(source, Loader=yaml.SafeLoader)

        circle_times.append(
            [
                (c["start"], c["end"])
                for c in flightinfo["segments"]
                if c["kind"] == "circle"
            ]
        )

        flight_date.append(np.datetime64(date.strftime(flightinfo["date"], "%Y-%m-%d")))

    return circle_times, flight_date


def get_circles(
    lv3_directory=lv3_directory,
    lv3a_filename=lv3a_filename,
    Platform="HALO",
    yaml_directory=yaml_directory,
):

    ds_lv3a = get_level3a_dataset(lv3_directory, lv3a_filename)

    all_sondes = dim_ready_ds(ds_lv3a, Platform)

    circle_times, flight_date = get_circle_times_from_yaml(yaml_directory)

    circles = []

    for i in range(len(flight_date)):
        for j in range(len(circle_times[i])):
            circles.append(
                all_sondes.sel(
                    launch_time=slice(
                        circle_times[i][j][0] - datetime.timedelta(minutes=2),
                        circle_times[i][j][1] + datetime.timedelta(minutes=2),
                    )
                )
            )

    return circles


def reswap_launchtime_sounding(circles):

    swapped_circles = []

    for circle in circles:
        circle["sounding"] = (
            ["launch_time"],
            np.arange(1, len(circle.launch_time) + 1, 1),
        )
        swapped_circles.append(circle.swap_dims({"launch_time": "sounding"}))

    return swapped_circles


# %%


def get_xy_coords_for_circles(circles):

    for i in range(len(circles)):

        x_coor = (
            circles[i]["longitude"]
            * 111.320
            * cos(np.radians(circles[i]["latitude"]))
            * 1000
        )
        y_coor = circles[i]["latitude"] * 110.54 * 1000

        # converting from lat, lon to coordinates in metre from (0,0).

        xc = [None] * len(x_coor.T)
        yc = [None] * len(y_coor.T)

        xc = np.mean(x_coor, axis=0)
        yc = np.mean(y_coor, axis=0)

        delta_x = x_coor - xc  # *111*1000 # difference of sonde long from mean long
        delta_y = y_coor - yc  # *111*1000 # difference of sonde lat from mean lat

        circles[i]["xc"] = (["height"], xc)
        circles[i]["yc"] = (["height"], yc)
        circles[i]["dx"] = (["launch_time", "height"], delta_x)
        circles[i]["dy"] = (["launch_time", "height"], delta_y)

    return print("Circles ready for regression")
