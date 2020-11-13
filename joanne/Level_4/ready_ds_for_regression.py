# %%
import yaml
import glob
import xarray as xr
import numpy as np

from datetime import date
import datetime
from pylab import cos

from joanne.Level_4 import dicts
import circle_fit as cf

# %%

yaml_directory = "/Users/geet/Documents/JOANNE/joanne/flight_segments/"
lv3_directory = "/Users/geet/Documents/JOANNE/Data/Level_3/"

lv3_filename = "EUREC4A_JOANNE_Dropsonde-RD41_Level_3_v0.8.1+0.g60c3587.dirty.nc"


def get_level3_dataset(lv3_directory=lv3_directory, lv3_filename=lv3_filename):
    return xr.open_dataset(lv3_directory + lv3_filename)


def get_circle_times_from_yaml(yaml_directory=yaml_directory):
    allyamlfiles = sorted(glob.glob(yaml_directory + "*P3*.yaml"))

    circle_times = []
    flight_date = []
    platform_name = []

    for i in allyamlfiles:
        with open(i) as source:
            flightinfo = yaml.load(source, Loader=yaml.SafeLoader)

        circle_times.append(
            [
                (c["start"], c["end"])
                for c in flightinfo["segments"]
                if "circle" in c["kinds"]
                if len(c["dropsondes"]["GOOD"]) >= 6
            ]
        )

        if "HALO" in i:
            platform_name.append("HALO")
        elif "P3" in i:
            platform_name.append("P3")
        else:
            platform_name.append("")

        flight_date.append(np.datetime64(date.strftime(flightinfo["date"], "%Y-%m-%d")))

    return circle_times, flight_date, platform_name


def dim_ready_ds(ds_lv3=get_level3_dataset()):

    dims_to_drop = ["sounding"]

    all_sondes = (
        ds_lv3.swap_dims({"sounding": "launch_time"})
        # .swap_dims({"obs": "alt"})
        .drop(dims_to_drop)
    )

    return all_sondes


def get_circles(
    lv3_directory=lv3_directory,
    lv3_filename=lv3_filename,
    # platform="HALO",
    yaml_directory=yaml_directory,
):

    ds_lv3 = get_level3_dataset(lv3_directory, lv3_filename)

    all_sondes = dim_ready_ds(ds_lv3)

    circle_times, flight_date, platform_name = get_circle_times_from_yaml(
        yaml_directory
    )

    circles = []

    for i in range(len(flight_date)):
        for j in range(len(circle_times[i])):
            circles.append(
                all_sondes.where(
                    all_sondes.platform == platform_name[i], drop=True
                ).sel(
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

        x_coor = circles[i]["lon"] * 111.320 * cos(np.radians(circles[i]["lat"])) * 1000
        y_coor = circles[i]["lat"] * 110.54 * 1000

        # converting from lat, lon to coordinates in metre from (0,0).
        c_xc = np.full(np.size(x_coor, 1), np.nan)
        c_yc = np.full(np.size(x_coor, 1), np.nan)
        c_r = np.full(np.size(x_coor, 1), np.nan)

        for j in range(np.size(x_coor, 1)):
            a = ~np.isnan(x_coor.values[:, j])
            if a.sum() > 4:
                c_xc[j], c_yc[j], c_r[j], _ = cf.least_squares_circle(
                    list(zip(x_coor.values[:, j], y_coor.values[:, j]))
                )

        circle_x = np.mean(c_xc)
        circle_y = np.mean(c_yc)
        circle_radius = np.mean(c_r)

        xc = [None] * len(x_coor.T)
        yc = [None] * len(y_coor.T)

        xc = np.mean(x_coor, axis=0)
        yc = np.mean(y_coor, axis=0)

        delta_x = x_coor - xc  # *111*1000 # difference of sonde long from mean long
        delta_y = y_coor - yc  # *111*1000 # difference of sonde lat from mean lat

        circles[i]["platform"] = circles[i].platform.values[0]
        circles[i]["flight_height"] = circles[i].flight_height.mean().values
        circles[i]["circle_time"] = circles[i].launch_time.mean().values
        circles[i]["circle_x"] = np.nanmean(c_xc)
        circles[i]["circle_y"] = np.nanmean(c_yc)
        circles[i]["circle_radius"] = np.nanmean(c_r)
        circles[i]["dx"] = (["launch_time", "alt"], delta_x)
        circles[i]["dy"] = (["launch_time", "alt"], delta_y)

    return print("Circles ready for regression")


def create_variable(ds, var, data, dims=dicts.nc_dims, attrs=dicts.nc_attrs, **kwargs):
    """Insert the data into a variable in an :class:`xr.Dataset`"""
    data = data[var]  # must be of type array
    attrs = attrs[var].copy()
    dims = dims[var]

    v = xr.Variable(dims, data, attrs=attrs)
    ds[var] = v

    return var
