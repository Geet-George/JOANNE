# %%
import xarray as xr
import numpy as np
import pandas as pd
from pylab import *
from datetime import date

from scipy.interpolate import griddata

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.dates as mdates
import matplotlib.colors

import yaml
import glob

import metpy
import metpy.calc as mpcalc
from metpy.units import units

# from metpy.future import precipitable_water

from tqdm import tqdm_notebook as tqdm
from tqdm import trange

from pandas import DataFrame
from sklearn import linear_model


# %%
allyamlfiles = sorted(glob.glob("*.yaml"))
all_sondes = xr.open_dataset("all_sondes.nc")
all_sondes = all_sondes.where(all_sondes.Platform == "HALO", drop=True)

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

# %%
all_circles = xr.open_dataset("circle_products.nc")
all_sondes = xr.open_dataset("all_sondes.nc")

# %%
allyamlfiles = sorted(glob.glob("*.yaml"))
all_sondes = xr.open_dataset("all_sondes.nc")
all_sondes = all_sondes.where(all_sondes.Platform == "HALO", drop=True)

circle_times = []
flight_date = []

for i in allyamlfiles[0:1]:
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

# %%
for i in range(len(circles)):

    x_coor = circles[i]["lon"] * 111.320 * cos(np.radians(circles[i]["lat"])) * 1000
    y_coor = circles[i]["lat"] * 110.54 * 1000

    # converting from lat, lon to coordinates in metre from (0,0).

    xc = [None] * len(x_coor.T)
    yc = [None] * len(y_coor.T)

    xc = np.mean(x_coor, axis=0)
    yc = np.mean(y_coor, axis=0)

    delta_x = x_coor - xc  # *111*1000 # difference of sonde long from mean long
    delta_y = y_coor - yc  # *111*1000 # difference of sonde lat from mean lat

    u_o = np.mean(circles[i]["u_wind"], axis=0)  # mean u velocity
    v_o = np.mean(circles[i]["v_wind"], axis=0)  # mean v velocity

    delta_u = circles[i]["u_wind"] - u_o  # difference of sonde u from mean u
    delta_v = circles[i]["v_wind"] - v_o  # difference of sonde v from mean v

    y_u = delta_u / delta_y  # yu (ratio; delta u & delta y: refer notes)
    y_v = delta_v / delta_y  # yv (ratio; delta v & delta y: refer notes)

    x = delta_x / delta_y  # x (ratio; delta x & delta y: refer notes)

    circles[i]["dx"] = (["launch_time", "alt"], delta_x)
    circles[i]["dy"] = (["launch_time", "alt"], delta_y)
    circles[i]["qu"] = (["launch_time", "alt"], circles[i]["q"] * circles[i]["u_wind"])
    circles[i]["qv"] = (["launch_time", "alt"], circles[i]["q"] * circles[i]["v_wind"])
    circles[i]["Tu"] = (
        ["launch_time", "alt"],
        (circles[i]["tdry"] + 273.15) * circles[i]["u_wind"],
    )
    circles[i]["Tv"] = (
        ["launch_time", "alt"],
        (circles[i]["tdry"] + 273.15) * circles[i]["v_wind"],
    )

# %%
m_u = np.full((len(circles), len(circles[0].alt)), np.nan)
m_v = np.full((len(circles), len(circles[0].alt)), np.nan)
c_u = np.full((len(circles), len(circles[0].alt)), np.nan)
c_v = np.full((len(circles), len(circles[0].alt)), np.nan)
m_q = np.full((len(circles), len(circles[0].alt)), np.nan)
c_q = np.full((len(circles), len(circles[0].alt)), np.nan)
m_T = np.full((len(circles), len(circles[0].alt)), np.nan)
c_T = np.full((len(circles), len(circles[0].alt)), np.nan)

m_qu = np.full((len(circles), len(circles[0].alt)), np.nan)
c_qu = np.full((len(circles), len(circles[0].alt)), np.nan)
m_Tu = np.full((len(circles), len(circles[0].alt)), np.nan)
c_Tu = np.full((len(circles), len(circles[0].alt)), np.nan)

m_qv = np.full((len(circles), len(circles[0].alt)), np.nan)
c_qv = np.full((len(circles), len(circles[0].alt)), np.nan)
m_Tv = np.full((len(circles), len(circles[0].alt)), np.nan)
c_Tv = np.full((len(circles), len(circles[0].alt)), np.nan)

D = np.full((len(circles), len(circles[0].alt)), np.nan)
vor = np.full((len(circles), len(circles[0].alt)), np.nan)
den_m = np.full((len(circles)), [None])
mean_den = np.full((len(circles), len(circles[0].alt)), np.nan)
w_vel = np.full((len(circles), len(circles[0].alt)), np.nan)
p_vel = np.full((len(circles), len(circles[0].alt)), np.nan)

Ns = np.full((len(circles), len(circles[0].alt)), np.nan)


for i in tqdm(range(len(circles))):
    # for i in tqdm(range(1)) :

    #     # loop for linear regression at every level
    for k in range(0, len(circles[i].alt)):

        df = {"dx": circles[i]["dx"].isel(alt=k), "dy": circles[i]["dy"].isel(alt=k)}
        df_q = {"q": circles[i]["q"].isel(alt=k)}

        ddf = DataFrame(df, columns=["dx", "dy"])
        ddf_q = DataFrame(df_q, columns=["q"])

        id_ = np.where((isnan(ddf["dx"]) == False) & (isnan(ddf_q["q"]) == False))[0]

        Ns[i][k] = int(len(id_))

        if id_.size >= 6:
            X = ddf.values[id_]
            Y_u = circles[i]["u_wind"].isel(alt=k).isel(launch_time=id_)
            Y_v = circles[i]["v_wind"].isel(alt=k).isel(launch_time=id_)
            Y_q = circles[i]["q"].isel(alt=k).isel(launch_time=id_)
            Y_T = circles[i]["tdry"].isel(alt=k).isel(launch_time=id_) + 273.15

            Y_qu = Y_u * Y_q
            Y_qv = Y_v * Y_q

            Y_Tu = Y_u * Y_T
            Y_Tv = Y_v * Y_T

            regr_u = linear_model.LinearRegression()
            regr_u.fit(X, Y_u)

            regr_v = linear_model.LinearRegression()
            regr_v.fit(X, Y_v)

            regr_q = linear_model.LinearRegression()
            regr_q.fit(X, Y_q)

            regr_T = linear_model.LinearRegression()
            regr_T.fit(X, Y_T)

            regr_qu = linear_model.LinearRegression()
            regr_qu.fit(X, Y_qu)

            regr_qv = linear_model.LinearRegression()
            regr_qv.fit(X, Y_qv)

            regr_Tu = linear_model.LinearRegression()
            regr_Tu.fit(X, Y_Tu)

            regr_Tv = linear_model.LinearRegression()
            regr_Tv.fit(X, Y_Tv)

            # #         mean_u = regr.intercept_
            m_u[i][k], c_u[i][k] = regr_u.coef_
            m_v[i][k], c_v[i][k] = regr_v.coef_
            m_q[i][k], c_q[i][k] = regr_q.coef_
            m_T[i][k], c_T[i][k] = regr_T.coef_

            m_qu[i][k], c_qu[i][k] = regr_qu.coef_
            m_Tu[i][k], c_Tu[i][k] = regr_Tu.coef_

            m_qv[i][k], c_qv[i][k] = regr_qv.coef_
            m_Tv[i][k], c_Tv[i][k] = regr_Tv.coef_

    #         pres[i] = np.full((len(circles[i].alt)),np.nan)
    #         D[i] = np.full((len(circles[i].alt)),np.nan)
    #         vor[i] = np.full((len(circles[i].alt)),np.nan)

    #         pres[i] = circles[i]['pres'].mean(dim='launch_time')
    D[i] = m_u[i] + c_v[i]
    vor[i] = m_v[i] - c_u[i]

    den_m[i] = mpcalc.density(
        circles[i]["pres"].values * units.hPa,
        circles[i]["tdry"].values * units.degC,
        circles[i]["mr"].values / 1000,
    ).magnitude
    mean_den[i] = np.nanmean(den_m[i], axis=0)

    # Calculating the vertical velocity and pressure velocity

    #         nan_ids = np.isnan(D[i])

    #         w_vel[i] = np.nan * len(circles[i].alt)
    w_vel[i][0] = 0

    #         for m in range(1,len(circles[i].alt)) :
    #             if m-1 is in nan_ids :
    #                 w_vel[i][m] = np.nansum(D[i][0:m+1]*10)*(-1)

    #         for i in range(len(D)) :

    nan_ids = np.where(np.isnan(D[i]) == True)[0]

    w_vel[i][0] = 0
    last = 0

    for m in range(1, len(circles[i].alt)):

        if m in nan_ids:
            w_vel[i][m] = np.nan
        else:
            w_vel[i][m] = w_vel[i][last] + D[i][m] * 10 * (m - last)
            last = m

    for n in range(1, 910):

        p_vel[i][n] = -mean_den[i][n] * 9.81 * w_vel[i][n] * 60 * 60 / 100

# %%

# %%
u = [None] * len(circles)
v = [None] * len(circles)
pres = [None] * len(circles)

for i in range(len(circles)):
    u[i] = circles[i].u_wind.mean(dim="launch_time").values
    v[i] = circles[i].v_wind.mean(dim="launch_time").values
    pres[i] = circles[i].pres.mean(dim="launch_time").values


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
    coords={"height": circles[0].alt.values, "circle": circle},
    attrs={
        "Author": "Geet George (MPI-M, Hamburg); geet.george@mpimet.mpg.de",
        "Instrument": "Vaisala RD94 (AVAPS receiver aboard HALO)",
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
