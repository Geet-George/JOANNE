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
plt.plot(all_circles.sel(circle="0205C6").omega)


# %%
circles[i].w_wind  # .mean(dim='launch_time')


# %%
i = 22
for j in circles[i].launch_time:
    plt.plot(
        circles[i].w_wind.sel(launch_time=j),
        circles[i].alt,
        c="grey",
        alpha=0.5,
        linewidth=0.5,
    )

plt.plot(all_circles.isel(circle=i).w, all_circles.height, c="r", linewidth=2)


# %%
pre = xr.open_dataset("pre_eurec4a_soundings_1989_to_2019_v0.1.nc")
rh = mpcalc.relative_humidity_from_dewpoint(
    pre.temperature.values * units.degC, pre.dewpoint.values * units.degC
).magnitude

pre["rh"] = (["date_and_time", "height"], rh * 100)
pre["speed"] = (["date_and_time", "height"], pre.speed * 0.5144)
pre["u_wind"] = (["date_and_time", "height"], pre.u_wind * 0.5144)
pre["v_wind"] = (["date_and_time", "height"], pre.v_wind * 0.5144)

eur = xr.open_dataset("all_sondes.nc")


# %%
pre_var = ("rh", "temperature", "speed", "direction", "u_wind", "v_wind", "pressure")
eur_var = ("rh", "tdry", "wspd", "wdir", "u_wind", "v_wind", "pres")
xlabels = (
    "RH\n/ %",
    "T\n/ $\degree$C",
    "Wind Spd\n/ m s$^{-1}$",
    "Wind Dir\n/ $\degree E$",
    "u\n/ m s$^{-1}$",
    "v\n/ m s$^{-1}$",
    "Pressure\n/ hPa",
)

f, ax = plt.subplots(1, len(pre_var), sharey=True, figsize=(12, 4))

for g, i in enumerate(pre_var):
    mean = pre[i].mean(dim="date_and_time")
    std = pre[i].std(dim="date_and_time")
    ax[g].plot(mean, pre.height / 1000, linewidth=2, color="grey", label="TBPB")

    ax[g].fill_betweenx(
        pre.height / 1000, mean - std, mean + std, color="grey", alpha=0.3
    )

    ax[g].set_xlabel(xlabels[g], fontsize=14)

for g, i in enumerate(eur_var):
    mean = eur[i].mean(dim="launch_time")
    std = eur[i].std(dim="launch_time")
    ax[g].plot(mean, eur.alt / 1000, linewidth=2, color="steelblue", label="EUREC4A")

    ax[g].fill_betweenx(
        eur.alt / 1000, mean - std, mean + std, color="steelblue", alpha=0.3
    )
    ax[g].set_ylim(0, 9)

    if g == len(eur_var) - 1:
        ax[g].legend(loc=1, bbox_to_anchor=(1.7, 0.95))

    [ax[g].spines[m].set_visible(False) for m in ["right", "top", "left"]]

    if g == 0:
        ax[g].set_ylabel("Altitude / km", fontsize=14)

plt.savefig("EUREC4A_vs_TBPB_9km.jpg", dpi=300)


# %%
pre_var = ("rh", "temperature", "speed", "direction", "u_wind", "v_wind", "pressure")
eur_var = ("rh", "tdry", "wspd", "wdir", "u_wind", "v_wind", "pres")
xlabels = (
    "RH\n/ %",
    "T\n/ $\degree$C",
    "Wind Spd\n/ m s$^{-1}$",
    "Wind Dir\n/ $\degree E$",
    "u\n/ m s$^{-1}$",
    "v\n/ m s$^{-1}$",
    "Pressure\n/ hPa",
)

f, ax = plt.subplots(1, len(pre_var), sharey=True, figsize=(12, 4))

for g, i in enumerate(pre_var):
    mean = pre[i].sel(height=slice(0, 4000)).mean(dim="date_and_time")
    std = pre[i].sel(height=slice(0, 4000)).std(dim="date_and_time")
    ax[g].plot(
        mean,
        pre.sel(height=slice(0, 4000)).height / 1000,
        linewidth=2,
        color="grey",
        label="TBPB",
    )

    ax[g].fill_betweenx(
        pre.sel(height=slice(0, 4000)).height / 1000,
        mean - std,
        mean + std,
        color="grey",
        alpha=0.3,
    )

    ax[g].set_xlabel(xlabels[g], fontsize=14)

for g, i in enumerate(eur_var):
    mean = eur[i].sel(alt=slice(0, 4000)).mean(dim="launch_time")
    std = eur[i].sel(alt=slice(0, 4000)).std(dim="launch_time")
    ax[g].plot(
        mean,
        eur.sel(alt=slice(0, 4000)).alt / 1000,
        linewidth=2,
        color="steelblue",
        label="EUREC4A",
    )

    ax[g].fill_betweenx(
        eur.sel(alt=slice(0, 4000)).alt / 1000,
        mean - std,
        mean + std,
        color="steelblue",
        alpha=0.3,
    )
    ax[g].set_ylim(0, 4)

    if g == len(eur_var) - 1:
        ax[g].legend(loc=1, bbox_to_anchor=(1.7, 0.95))

    [ax[g].spines[m].set_visible(False) for m in ["right", "top", "left"]]

    if g == 0:
        ax[g].set_ylabel("Altitude / km", fontsize=14)

plt.savefig("EUREC4A_vs_TBPB_4km.jpg", dpi=300)


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
circles[0].q


# %%
nan_counts = []
launch_times = []

for j in range(len(circles)):
    for i in range(len(circles[j].launch_time)):
        nan_counts.append(np.isnan(circles[j].u_wind.isel(launch_time=i)).sum().values)
        launch_times.append(circles[j].launch_time.isel(launch_time=i).values)


# %%
plt.figure(figsize=(12, 5))
plt.plot(launch_times, nan_counts)
plt.xlabel("Time / UTC")
plt.ylabel("NaN Counts for theta")
plt.xticks(plt.gca().get_xticks()[1::2])
plt.savefig("NaN Counts for profiles.pdf")
# nan_counts


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

    #         p_vel[i] = np.full_like(D[i],np.nan) # Pressure velocity in hPa/h

    for n in range(1, 910):

        p_vel[i][n] = -mean_den[i][n] * 9.81 * w_vel[i][n] * 60 * 60 / 100


# %%
# for i in range(len(D)) :

#     nan_ids = np.where(np.isnan(D[i]) == True)[0]

#     w_vel[i][0] = 0
#     last = 0

#     for m in range(1,len(circles[i].alt)) :

#         if m in nan_ids :
#             w_vel[i][m] = np.nan
#         #     counter.append(m)
#         else :
#             w_vel[i][m] = w_vel[i][last] + D[i][m]*10*(m-last)
#             last = m


# %%
last


# %%
plt.scatter(w_vel[2], circles[0].alt)
plt.ylim(0, 200)


# %%
m = 4
# np.nansum(D[i][0:m+1]*10)
D[i][0 : m + 1]


# %%
m = circles[30].u_wind.mean(dim="launch_time")
std = circles[30].u_wind.std(dim="launch_time")

plt.plot(m, circles[30].alt)
plt.fill_betweenx(circles[30].alt, m - std, m + std, alpha=0.5)


# %%
m = circles[30].u_wind.mean(dim="launch_time")
std = circles[30].u_wind.std(dim="launch_time")

plt.plot(m, circles[30].alt)
plt.fill_betweenx(circles[30].alt, m - std, m + std, alpha=0.5)


# %%
u = [None] * len(circles)
v = [None] * len(circles)
pres = [None] * len(circles)

for i in range(len(circles)):
    u[i] = circles[i].u_wind.mean(dim="launch_time").values
    v[i] = circles[i].v_wind.mean(dim="launch_time").values
    pres[i] = circles[i].pres.mean(dim="launch_time").values


# %%
circles[0]


# %%
circle = []

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


# %%
ds


# %%
f = plt.figure(figsize=(10, 10))

for i in range(6, 12):
    p = plt.plot(ds.isel(circle=i).omega, ds.height, label=i)
    plt.plot(
        ds.isel(circle=i).regressed_sondes,
        ds.height,
        linestyle="--",
        c=p[0].get_color(),
    )
    #     plt.ylim(0,200)
    plt.legend()


# %%
p[0].get_color()


# %%
# for i in range(len(circles)) :

#         x_coor = circles[i]['lon']*111.320*cos(np.radians(circles[i]['lat']))*1000
#         y_coor = circles[i]['lat']*110.54*1000

#         # converting from lat, lon to coordinates in metre from (0,0).

#         xc = [None] * len(x_coor.T)
#         yc = [None] * len(y_coor.T)

#         xc = np.mean(x_coor,axis=0)
#         yc = np.mean(y_coor,axis=0)

#         delta_x = (x_coor - xc)#*111*1000 # difference of sonde long from mean long
#         delta_y = (y_coor - yc)#*111*1000 # difference of sonde lat from mean lat

#         u_o = np.mean(circles[i]['u_wind'],axis=0) # mean u velocity
#         v_o = np.mean(circles[i]['v_wind'],axis=0) # mean v velocity

#         delta_u = circles[i]['u_wind'] - u_o # difference of sonde u from mean u
#         delta_v = circles[i]['v_wind'] - v_o # difference of sonde v from mean v

#         y_u = delta_u/delta_y # yu (ratio; delta u & delta y: refer notes)
#         y_v = delta_v/delta_y # yv (ratio; delta v & delta y: refer notes)

#         x = delta_x/delta_y # x (ratio; delta x & delta y: refer notes)

#         circles[i]['dx'] = (['launch_time','alt'],delta_x)
#         circles[i]['dy'] = (['launch_time','alt'],delta_y)
#         circles[i]['qu'] = (['launch_time','alt'],
#                             circles[i]['q']*circles[i]['u_wind'])
#         circles[i]['qv'] = (['launch_time','alt'],
#                             circles[i]['q']*circles[i]['v_wind'])
#         circles[i]['Tu'] = (['launch_time','alt'],
#                             (circles[i]['tdry']+ 273.15)*circles[i]['u_wind'] )
#         circles[i]['Tv'] = (['launch_time','alt'],
#                             (circles[i]['tdry']+ 273.15)*circles[i]['v_wind'] )

# m_u = np.full((len(circles),len(circles[0].alt)),np.nan)
# m_v = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_u = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_v = np.full((len(circles),len(circles[0].alt)),np.nan)
# m_q = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_q = np.full((len(circles),len(circles[0].alt)),np.nan)
# m_T = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_T = np.full((len(circles),len(circles[0].alt)),np.nan)

# m_qu = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_qu = np.full((len(circles),len(circles[0].alt)),np.nan)
# m_Tu = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_Tu = np.full((len(circles),len(circles[0].alt)),np.nan)

# m_qv = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_qv = np.full((len(circles),len(circles[0].alt)),np.nan)
# m_Tv = np.full((len(circles),len(circles[0].alt)),np.nan)
# c_Tv = np.full((len(circles),len(circles[0].alt)),np.nan)

# D = np.full((len(circles),len(circles[0].alt)),np.nan)
# vor = np.full((len(circles),len(circles[0].alt)),np.nan)
# den_m = np.full((len(circles)),[None])
# mean_den = np.full((len(circles),len(circles[0].alt)),np.nan)
# w_vel = np.full((len(circles),len(circles[0].alt)),np.nan)
# p_vel = np.full((len(circles),len(circles[0].alt)),np.nan)

# # for i in trange(len(circles)) :
# for i in trange(4) :

# #     for j in range(len(circle_times[i])) :

# #         m_u[i] = np.full((len(circles[i].alt)),np.nan)
# #         m_v[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_u[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_v[i] = np.full((len(circles[i].alt)),np.nan)
# #         m_q[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_q[i] = np.full((len(circles[i].alt)),np.nan)
# #         m_T[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_T[i] = np.full((len(circles[i].alt)),np.nan)

# #         m_qu[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_qu[i] = np.full((len(circles[i].alt)),np.nan)
# #         m_Tu[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_Tu[i] = np.full((len(circles[i].alt)),np.nan)

# #         m_qv[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_qv[i] = np.full((len(circles[i].alt)),np.nan)
# #         m_Tv[i] = np.full((len(circles[i].alt)),np.nan)
# #         c_Tv[i] = np.full((len(circles[i].alt)),np.nan)

#         #     # loop for linear regression at every level
#         for k in range(0,len(circles[i].alt)) :

#             df = {'dx':delta_x.isel(alt=k),
#                  'dy':delta_y.isel(alt=k)}
#             df_q = {'q':circles[i]['q'].isel(alt=k)}

#             ddf = DataFrame(df,columns=['dx','dy'])
#             ddf_q = DataFrame(df_q,columns=['q'])

#             id = np.where((isnan(ddf['dx'])==False) & (isnan(ddf_q['q'])==False))[0]

#             if id.size != 0 :
#                 X = ddf.values[id]
#                 Y_u = circles[i]['u_wind'].isel(alt=k).isel(launch_time=id)
#                 Y_v = circles[i]['v_wind'].isel(alt=k).isel(launch_time=id)
#                 Y_q = circles[i]['q'].isel(alt=k).isel(launch_time=id)
#                 Y_T = circles[i]['tdry'].isel(alt=k).isel(launch_time=id) + 273.15

#                 Y_qu = Y_u*Y_q
#                 Y_qv = Y_v*Y_q

#                 Y_Tu = Y_u*Y_T
#                 Y_Tv = Y_v*Y_T

#                 regr_u = linear_model.LinearRegression()
#                 regr_u.fit(X,Y_u)

#                 regr_v = linear_model.LinearRegression()
#                 regr_v.fit(X,Y_v)

#                 regr_q = linear_model.LinearRegression()
#                 regr_q.fit(X,Y_q)

#                 regr_T = linear_model.LinearRegression()
#                 regr_T.fit(X,Y_T)

#                 regr_qu = linear_model.LinearRegression()
#                 regr_qu.fit(X,Y_qu)

#                 regr_qv = linear_model.LinearRegression()
#                 regr_qv.fit(X,Y_qv)

#                 regr_Tu = linear_model.LinearRegression()
#                 regr_Tu.fit(X,Y_Tu)

#                 regr_Tv = linear_model.LinearRegression()
#                 regr_Tv.fit(X,Y_Tv)

#         # #         mean_u = regr.intercept_
#                 m_u[i][k],c_u[i][k]= regr_u.coef_
#                 m_v[i][k],c_v[i][k]= regr_v.coef_
#                 m_q[i][k],c_q[i][k]= regr_q.coef_
#                 m_T[i][k],c_T[i][k]= regr_T.coef_

#                 m_qu[i][k],c_qu[i][k] = regr_qu.coef_
#                 m_Tu[i][k],c_Tu[i][k] = regr_Tu.coef_

#                 m_qv[i][k],c_qv[i][k] = regr_qv.coef_
#                 m_Tv[i][k],c_Tv[i][k] = regr_Tv.coef_

# #         pres[i] = np.full((len(circles[i].alt)),np.nan)
# #         D[i] = np.full((len(circles[i].alt)),np.nan)
# #         vor[i] = np.full((len(circles[i].alt)),np.nan)

# #         pres[i] = circles[i]['pres'].mean(dim='launch_time')
#         D[i] = (m_u[i] + c_v[i])
#         vor[i] = (m_v[i] - c_u[i])

#         den_m[i] = mpcalc.density(circles[i]['pres'].values*units.hPa,
#                          circles[i]['tdry'].values*units.degC,
#                          circles[i]['mr'].values/1000).magnitude
#         mean_den[i] = np.nanmean(den_m[i],axis=0)

#         # Calculating the vertical velocity and pressure velocity

# #         w_vel[i] = np.nan * len(circles[i].alt)
#         w_vel[i][0] = 0

#         for m in range(1,len(circles[i].alt)) :
#             w_vel[i][m] = np.nansum(D[i][0:m+1]*10)*(-1)

# #         p_vel[i] = np.full_like(D[i],np.nan) # Pressure velocity in hPa/h

#         for n in range(1,910) :

#             p_vel[i][n] = -mean_den[i][n]*9.81*w_vel[i][n]*60*60/100


# %%
# Pressure Gradient

upto = 200

lat_points = np.full(size(circles), [None])
lon_points = np.full(size(circles), [None])
p_values = np.full(size(circles), [None])
grid_x = np.full(size(circles), [None])
grid_y = np.full(size(circles), [None])
grid = np.full(size(circles), [None])

mean_time = np.full(size(circles), [None])
da = []

g = 0
for i in range(len(circles)):
    for j in range(len(circle_times[i])):

        lat_points[g] = circles[i, j].sel(alt=slice(0, upto)).lat.mean(dim="alt").values
        lon_points[g] = (
            circles[i, j].sel(alt=slice(0, upto)).lon.mean(dim="alt").values + 360
        )
        l_points = list(set(zip(lon_points[g], lat_points[g])))
        p_values[g] = circles[i, j].sel(alt=slice(0, upto)).pres.mean(dim="alt").values
        mean_time[g] = circles[i, j].launch_time.mean().values

        grid_x[g], grid_y[g] = np.mgrid[
            floor(np.min(lon_points[g])) : ceil(np.max(lon_points[g])) : 0.2,
            floor(np.min(lat_points[g])) : ceil(np.max(lat_points[g])) : 0.2,
        ]

        grid[g] = griddata(
            l_points, p_values[g], (grid_x[g], grid_y[g]), method="linear"
        )

        try:
            da.append(
                xr.DataArray(
                    grid[g],
                    coords={
                        "lon": (grid_x[g][:, 0]),
                        "lat": (grid_y[g][0]),
                        "circle": (circles[i, j].launch_time.mean().values),
                    },
                    dims={"lon", "lat"},
                )
            )
        except Exception:
            print(i, j)
        #         else :
        #             mean_time.append(circles[i,j].launch_time.mean().values)

        g += 1


# %%
ds1 = xr.concat(da[:], dim="circle")
ds1.to_netcdf("spatially_gridded_pressure.nc")


# %%
# ds = xr.Dataset({'lat':(['circle','sonde_no'],lat_points),
#                  'lon':(['circle','sonde_no'],lon_points),
#                  'pressure':(['circle','sonde_no'],p_values)},
#                   coords={'circle':(mean_time)})
# ds.to_netcdf('a.nc')


# %%


# %%
# ds1 = xr.Dataset({'pressure':(['circle','lon','lat'],grid)},
#                 coords={'circle':(mean_time),
#                         'lon':(grid_x),
#                         'lat':(grid_y)})


# %%
u = np.gradient(grid)[0]
v = np.gradient(grid)[1]

plt.quiver(grid_x, grid_y, u, v)
plt.scatter(lon_points, lat_points, s=30, c="r")


# %%
n = 10
plt.imshow(
    grid[n],
    interpolation="none",
    extent=[
        np.min(lon_points[n]),
        np.max(lon_points[n]),
        np.min(lat_points[n]),
        np.max(lat_points[n]),
    ],
)
plt.colorbar(orientation="horizontal")


# %%
shape(grid[0])


# %%
ncds = all_sondes.sel(launch_time="2020-01-24")
ncds = ncds.where(ncds.Platform == "HALO", drop=True)
ncds = ncds.where(
    ncds.reference_time != np.datetime64("2020-01-24T00:00:00.000000000"), drop=True
).isel(launch_time=slice(0, 72))
ncds


# %%
plt.plot(ncds.lon[0:72].isel(alt=10), ncds.lat[0:72].isel(alt=10), marker="o")


# %%
sondes_per_reg = 18
total_sondes = len(ncds.launch_time)

iterations = total_sondes - (sondes_per_reg - 1)


# %%
shape(grid[k])


# %%
upto = 200  # height up to which the average is to be taken

grid_x = [None] * iterations
grid_y = [None] * iterations
grid = [None] * iterations

u = [None] * iterations
v = [None] * iterations
mean_time = [None] * iterations

for k in tqdm(range(0, iterations)):

    lat_points = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(0, upto))
        .lat.mean(dim="alt")
        .values
    )
    lon_points = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(0, upto))
        .lon.mean(dim="alt")
        .values
        + 360
    )

    l_points = list(set(zip(lon_points, lat_points)))

    p_values = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(0, upto))
        .pres.mean(dim="alt")
        .values
    )

    grid_x[k], grid_y[k] = np.mgrid[
        floor(np.min(lon_points)) : ceil(np.max(lon_points)) : 0.2,
        floor(np.min(lat_points)) : ceil(np.max(lat_points)) : 0.2,
    ]

    grid[k] = griddata(l_points, p_values, (grid_x[k], grid_y[k]), method="linear")

    u[k] = np.gradient(grid[k])[0]
    v[k] = np.gradient(grid[k])[1]

    f, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12, 6))

    norm = matplotlib.colors.Normalize(vmin=0.1, vmax=2, clip=False)
    cmap = matplotlib.colors.ListedColormap(
        plt.cm.cividis(np.linspace(0.25, 1, 10)), "name"
    )

    vmin1 = 0.1
    vmax1 = 0.5
    q = ax[0].quiver(
        grid_x[k],
        grid_y[k],
        u[k],
        v[k],
        np.hypot(u[k], v[k]),
        units="width",
        pivot="mid",
        width=0.01,  # vmin=vmin1,vmax=vmax1,
        scale=1 / 0.1,
        cmap=cmap,
        norm=norm,
    )
    #     ax[0].quiverkey(q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
    #                    coordinates='figure')

    cbar1 = plt.colorbar(q, ax=ax[0], norm=norm)
    cbar1.set_label("Pressure Gradient (hPa / $\degree$)")
    #     cbar1 = q.clim(0.1,2)

    ax[0].scatter(lon_points, lat_points, s=45, c="r", label="Dropsonde considered")
    ax[0].legend(loc=1)
    ax[0].set_ylabel("Latitude ($\degree$N)")
    ax[0].set_xlabel("Longitude ($\degree$E)")

    vmin2 = 996
    vmax2 = 1002.5
    levels2 = np.linspace(vmin2, vmax2 + 0.5, 50)

    #     CS = ax[1].contourf(grid_x[k],grid_y[k],
    #                   grid[k],levels=levels2,vmin=vmin2,vmax=vmax2,cmap='Blues'
    #                   )
    CS = ax[1].imshow(
        grid[k].T,
        interpolation="none",
        origin="lower",
        extent=[
            floor(np.min(lon_points)),
            ceil(np.max(lon_points)),
            floor(np.min(lat_points)),
            ceil(np.max(lat_points)),
        ],
        vmin=vmin2,
        vmax=vmax2,
    )
    # plt.colorbar(orientation='horizontal')
    ax[1].set_xlabel("Longitude ($\degree$E)")
    cbar2 = plt.colorbar(CS, ax=ax[1], ticks=[vmin2, vmax2])
    cbar2.set_label("Pressure (hPa)")

    mean_time[k] = np.array(
        ncds.isel(launch_time=slice(k, k + sondes_per_reg)).launch_time.mean().values,
        dtype="datetime64[m]",
    )

    plt.suptitle(f"Mean time of circle: {mean_time[k]}", fontsize=16)
    plt.savefig(f"{k}.jpeg")
    plt.close()


# %%
# upto = 200 # height up to which the average is to be taken

# grid_x = [None] * iterations
# grid_y = [None] * iterations
# # grid = [None] * iterations

# grid_u = [None] * iterations
# grid_v = [None] * iterations

# for k in tqdm(range(0,iterations)) :

#     lat_points = ncds.isel(launch_time=slice(k,k+sondes_per_reg)).\
#                     sel(alt=slice(0,upto)).lat.mean(dim='alt').values
#     lon_points = ncds.isel(launch_time=slice(k,k+sondes_per_reg)).\
#                     sel(alt=slice(0,upto)).lon.mean(dim='alt').values + 360

#     l_points = list(set(zip(lon_points,lat_points)))

#     u_values = ncds.isel(launch_time=slice(k,k+sondes_per_reg)).\
#                     sel(alt=slice(0,upto)).u_wind.mean(dim='alt').values
#     v_values = ncds.isel(launch_time=slice(k,k+sondes_per_reg)).\
#                     sel(alt=slice(0,upto)).v_wind.mean(dim='alt').values

#     grid_x[k], grid_y[k] = np.mgrid[floor(np.min(lon_points)):ceil(np.max(lon_points)):0.2,
#                           floor(np.min(lat_points)):ceil(np.max(lat_points)):0.2]

#     grid_u[k] = griddata(l_points,u_values, (grid_x[k],grid_y[k]), method = 'linear')
#     grid_v[k] = griddata(l_points,v_values, (grid_x[k],grid_y[k]), method = 'linear')

# #     u[k] = np.gradient(grid[k])[0]
# #     v[k] = np.gradient(grid[k])[1]

#     f,ax = plt.subplots(1,1,sharex=True,sharey=True,figsize=(12,12))

#     norm = matplotlib.colors.Normalize(vmin=0.1,vmax=2,clip=False)
#     cmap = matplotlib.colors.ListedColormap(plt.cm.cividis(np.linspace(0.25,1,10)), "name")

#     vmin1 = 1; vmax1 = 10
#     q = ax.quiver(grid_x[k],grid_y[k],
#                      grid_u[k],grid_v[k],np.hypot(grid_u[k],grid_v[k]),
#                      units='width', pivot='mid', width=0.01,#vmin=vmin1,vmax=vmax1,
#                scale=1 / 0.1,cmap=cmap,norm=norm)
# #     ax[0].quiverkey(q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
# #                    coordinates='figure')

#     cbar1 = plt.colorbar(q,ax=ax,norm=norm)
#     cbar1.set_label('Pressure Gradient (hPa / $\degree$)')
# #     cbar1 = q.clim(0.1,2)

# #     ax[0].scatter(lon_points,lat_points,s=45,c='r',label='Dropsonde considered')
# #     ax[0].legend()
# #     ax[0].set_ylabel('Latitude ($\degree$N)')
# #     ax[0].set_xlabel('Longitude ($\degree$E)')

# #     vmin2=996; vmax2=1002.5
# #     levels2 = np.linspace(vmin2, vmax2+0.5,50)

# #     CS = ax[1].contourf(grid_x[k],grid_y[k],
# #                   grid[k],levels=levels2,vmin=vmin2,vmax=vmax2,cmap='Blues'
# #                   )
# #     ax[1].set_xlabel('Longitude ($\degree$E)')
# #     cbar2 = plt.colorbar(CS,ax=ax[1],ticks=[vmin2,vmax2])
# #     cbar2.set_label('Pressure (hPa)')

#     mean_time = np.array(ncds.isel(launch_time=slice(k,
#                     k+sondes_per_reg)).launch_time.mean().values, dtype='datetime64[m]')

#     plt.suptitle(f'Mean time of circle: {mean_time}')
#     plt.savefig(f'{k}_uv.jpeg')
#     plt.close()


# %%
ncds.mr


# %%
upto = 800  # height up to which the average is to be taken
fro = 600  # height from which the average is to be taken

grid_x = [None] * iterations
grid_y = [None] * iterations
grid = [None] * iterations

u = [None] * iterations
v = [None] * iterations

for k in tqdm(range(0, iterations)):

    lat_points = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(fro, upto))
        .lat.mean(dim="alt")
        .values
    )
    lon_points = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(fro, upto))
        .lon.mean(dim="alt")
        .values
        + 360
    )

    l_points = list(set(zip(lon_points, lat_points)))

    mr_values = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(fro, upto))
        .mr.mean(dim="alt")
        .values
    )

    grid_x[k], grid_y[k] = np.mgrid[
        floor(np.min(lon_points)) : ceil(np.max(lon_points)) : 0.2,
        floor(np.min(lat_points)) : ceil(np.max(lat_points)) : 0.2,
    ]

    grid[k] = griddata(l_points, mr_values, (grid_x[k], grid_y[k]), method="linear")

    u[k] = np.gradient(grid[k])[0]
    v[k] = np.gradient(grid[k])[1]

    f, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12, 6))

    vmin1 = 3e-01
    vmax1 = 6e-01

    norm = matplotlib.colors.Normalize(vmin=vmin1, vmax=vmax1, clip=False)
    cmap = matplotlib.colors.ListedColormap(
        plt.cm.Blues(np.linspace(0.25, 1, 10)), "name"
    )

    q = ax[0].quiver(
        grid_x[k],
        grid_y[k],
        u[k],
        v[k],
        np.hypot(u[k], v[k]),
        units="width",
        pivot="mid",
        width=0.01,  # vmin=vmin1,vmax=vmax1,
        cmap=cmap,
        norm=norm,
    )
    #     ax[0].quiverkey(q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
    #                    coordinates='figure')

    cbar1 = plt.colorbar(q, ax=ax[0], norm=norm)
    cbar1.set_label("Moisture Gradient ( g/kg / $\degree$)")
    #     cbar1 = q.clim(0.1,2)

    ax[0].scatter(lon_points, lat_points, s=45, c="r", label="Dropsondes considered")
    ax[0].legend(loc=1)
    ax[0].set_ylabel("Latitude ($\degree$N)")
    ax[0].set_xlabel("Longitude ($\degree$E)")

    #     vmin2=0.015; vmax2=0.017
    #     levels2 = np.linspace(vmin2, vmax2+0.5,50)

    vmin2 = 13
    vmax2 = 15.5
    levels2 = np.linspace(vmin2, vmax2 + 0.5, 50)

    #     CS = ax[1].contourf(grid_x[k],grid_y[k],
    #                   grid[k],levels=levels2,vmin=vmin2,vmax=vmax2,cmap='Blues'
    #                   )
    CS = ax[1].imshow(
        grid[k].T,
        interpolation="none",
        origin="lower",
        extent=[
            floor(np.min(lon_points)),
            ceil(np.max(lon_points)),
            floor(np.min(lat_points)),
            ceil(np.max(lat_points)),
        ],
        vmin=vmin2,
        vmax=vmax2,
    )
    # plt.colorbar(orientation='horizontal')
    ax[1].set_xlabel("Longitude ($\degree$E)")
    cbar2 = plt.colorbar(CS, ax=ax[1], ticks=[vmin2, vmax2])
    cbar2.set_label("Mixing Ratio (g/kg)")

    mean_time = np.array(
        ncds.isel(launch_time=slice(k, k + sondes_per_reg)).launch_time.mean().values,
        dtype="datetime64[m]",
    )

    plt.suptitle(f"Mean time of circle: {mean_time}")
    plt.savefig(f"{k}_mr_LCL.jpeg")
    plt.close()


# %%
grid[k]


# %%
upto = 150  # height up to which the average is to be taken
fro = 0  # height from which the average is to be taken

grid_x = [None] * iterations
grid_y = [None] * iterations
grid = [None] * iterations

u = [None] * iterations
v = [None] * iterations

for k in tqdm(range(0, iterations)):

    lat_points = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(fro, upto))
        .lat.mean(dim="alt")
        .values
    )
    lon_points = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(fro, upto))
        .lon.mean(dim="alt")
        .values
        + 360
    )

    l_points = list(set(zip(lon_points, lat_points)))

    q_values = (
        ncds.isel(launch_time=slice(k, k + sondes_per_reg))
        .sel(alt=slice(fro, upto))
        .q.mean(dim="alt")
        .values
    )

    grid_x[k], grid_y[k] = np.mgrid[
        floor(np.min(lon_points)) : ceil(np.max(lon_points)) : 0.2,
        floor(np.min(lat_points)) : ceil(np.max(lat_points)) : 0.2,
    ]

    grid[k] = griddata(l_points, q_values, (grid_x[k], grid_y[k]), method="linear")

    u[k] = np.gradient(grid[k])[0]
    v[k] = np.gradient(grid[k])[1]

    f, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12, 6))

    vmin1 = 0
    vmax1 = 6e-04

    norm = matplotlib.colors.Normalize(vmin=vmin1, vmax=vmax1, clip=False)
    cmap = matplotlib.colors.ListedColormap(
        plt.cm.Blues(np.linspace(0.25, 1, 10)), "name"
    )

    q = ax[0].quiver(
        grid_x[k],
        grid_y[k],
        u[k],
        v[k],
        np.hypot(u[k], v[k]),
        units="width",
        pivot="mid",
        width=0.01,  # vmin=vmin1,vmax=vmax1,
        cmap=cmap,
        norm=norm,
    )
    #     ax[0].quiverkey(q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
    #                    coordinates='figure')

    cbar1 = plt.colorbar(q, ax=ax[0], norm=norm)
    cbar1.set_label("Moisture Gradient ( / $\degree$)")
    #     cbar1 = q.clim(0.1,2)

    ax[0].scatter(lon_points, lat_points, s=45, c="r", label="Dropsondes considered")
    ax[0].legend()
    ax[0].set_ylabel("Latitude ($\degree$N)")
    ax[0].set_xlabel("Longitude ($\degree$E)")

    #     vmin2=0.015; vmax2=0.017
    #     levels2 = np.linspace(vmin2, vmax2+0.5,50)

    CS = ax[1].contourf(
        grid_x[k], grid_y[k], grid[k], cmap="Blues"
    )  # ,levels=levels2,vmin=vmin2,vmax=vmax2
    ax[1].set_xlabel("Longitude ($\degree$E)")
    cbar2 = plt.colorbar(CS, ax=ax[1])  # ,ticks=[vmin2,vmax2]
    cbar2.set_label("Specific Humidity (hPa)")

    mean_time = np.array(
        ncds.isel(launch_time=slice(k, k + sondes_per_reg)).launch_time.mean().values,
        dtype="datetime64[m]",
    )

    plt.suptitle(f"Mean time of circle: {mean_time}")
    plt.savefig(f"{k}_q_srf.jpeg")
    plt.close()


# %%
