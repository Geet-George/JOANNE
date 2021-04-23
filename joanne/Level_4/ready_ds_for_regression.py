# %%
import yaml
import glob
import xarray as xr
import numpy as np
from packaging import version

from datetime import date
import datetime
from pylab import cos
import joanne
from joanne.Level_4 import dicts
import circle_fit as cf

# %%

yaml_directory = "/Users/geet/Documents/JOANNE/joanne/flight_segments/"
lv3_directory = "/Users/geet/Documents/JOANNE/Data/Level_3/"

lv3_files = sorted(
    glob.glob(lv3_directory + f"EUREC4A_JOANNE_Dropsonde-RD41_Level_3_v*.nc")
)

vers = [None] * len(lv3_files)

for n, i in enumerate(lv3_files):
    vers[n] = version.parse(i)

lv3_filename = str(max(vers))


def get_level3_dataset(lv3_directory=lv3_directory, lv3_filename=lv3_filename):
    return xr.open_dataset(lv3_filename)


def get_circle_times_from_yaml(yaml_directory=yaml_directory):
    allyamlfiles = sorted(glob.glob(yaml_directory + "*.yaml"))

    circle_times = []
    sonde_ids = []
    flight_date = []
    platform_name = []
    segment_id = []

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

        sonde_ids.append(
            [
                c["dropsondes"]["GOOD"]
                for c in flightinfo["segments"]
                if "circle" in c["kinds"]
                if len(c["dropsondes"]["GOOD"]) >= 6
            ]
        )

        segment_id.append(
            [
                (c["segment_id"])
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

    return sonde_ids, circle_times, flight_date, platform_name, segment_id


def dim_ready_ds(ds_lv3=get_level3_dataset()):

    dims_to_drop = ["sounding"]

    all_sondes = (
        ds_lv3.swap_dims({"sounding": "sonde_id"})
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

    ds_fn = get_level3_dataset(lv3_directory, lv3_filename)

    (
        sonde_ids,
        circle_times,
        flight_date,
        platform_name,
        segment_id,
    ) = get_circle_times_from_yaml(yaml_directory)

    circles = []

    for i in range(len(flight_date)):
        for j in range(len(circle_times[i])):
            if len(sonde_ids[i]) != 0:
                circles.append(ds_fn.sel(sonde_id=sonde_ids[i][j]))
                # .swap_dims(
                # {"sonde_id": "sonde_id"}
                # )
                # )

            circles[-1]["segment_id"] = segment_id[i][j]

            circles[-1] = circles[-1].pad(
                sonde_id=(0, 13 - int(len(circles[-1].sonde_id))), mode="constant"
            )

            circles[-1]["sounding"] = (["sonde_id"], np.arange(0, 13, 1, dtype="int"))

            circles[-1] = circles[-1].swap_dims({"sonde_id": "sounding"})

    return circles


def reswap_launchtime_sounding(circle):

    # swapped_circles = []

    # for circle in circles:
    circle["sounding"] = (
        ["launch_time"],
        np.arange(1, len(circle.launch_time) + 1, 1),
    )
    circle = circle.swap_dims({"launch_time": "sounding"})

    return circle


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
                    [
                        (k, l)
                        for k, l in zip(x_coor.values[:, j], y_coor.values[:, j])
                        if ~np.isnan(k)
                    ]
                )

        circle_y = np.nanmean(c_yc) / (110.54 * 1000)
        circle_x = np.nanmean(c_xc) / (111.320 * cos(np.radians(circle_y)) * 1000)

        circle_diameter = np.nanmean(c_r) * 2

        xc = [None] * len(x_coor.T)
        yc = [None] * len(y_coor.T)

        xc = np.mean(x_coor, axis=0)
        yc = np.mean(y_coor, axis=0)

        delta_x = x_coor - xc  # *111*1000 # difference of sonde long from mean long
        delta_y = y_coor - yc  # *111*1000 # difference of sonde lat from mean lat

        circles[i]["platform_id"] = circles[i].platform_id.values[0]
        circles[i]["flight_altitude"] = circles[i].flight_altitude.mean().values
        circles[i]["circle_time"] = (
            circles[i].launch_time.mean().values.astype("datetime64")
        )
        # circles[i].encoding["circle_time"] = {
        #     "units": "seconds since 2020-01-01",
        #     "dtype": "datetime64[ns]",
        # }
        circles[i]["circle_lon"] = circle_x
        circles[i]["circle_lat"] = circle_y
        circles[i]["circle_diameter"] = circle_diameter
        circles[i]["dx"] = (["sounding", "alt"], delta_x)
        circles[i]["dy"] = (["sounding", "alt"], delta_y)

    return print("Circles ready for regression")


def create_variable(ds, var, data, dims=dicts.nc_dims, attrs=dicts.nc_attrs, **kwargs):
    """Insert the data into a variable in an :class:`xr.Dataset`"""
    data = data[var]  # must be of type array
    attrs = attrs[var].copy()
    dims = dims[var]

    v = xr.Variable(dims, data, attrs=attrs)
    ds[var] = v

    return var
