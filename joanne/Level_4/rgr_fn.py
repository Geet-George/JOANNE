# Module to store functions for regression

from sklearn import linear_model
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units


def run_regression(circle, parameter):

    """
    Input :
        circle : xarray dataset
                 dataset with sondes of a circle, with dx and dy calculated
                
        parameter : string
                    the parameter on which regression is to be carried out
    
    Output :

        mean_parameter : mean of parameter (intercept)
        m_parameter, c_parameter    : coefficients of regression

    """
    id_u = ~np.isnan(circle.u.values)
    id_v = ~np.isnan(circle.v.values)
    id_q = ~np.isnan(circle.q.values)
    id_x = ~np.isnan(circle.dx.values)
    id_y = ~np.isnan(circle.dy.values)
    id_t = ~np.isnan(circle.T.values)
    id_p = ~np.isnan(circle.p.values)

    id_quxv = np.logical_and(np.logical_and(id_q, id_u), np.logical_and(id_x, id_v))
    id_ = np.logical_and(np.logical_and(id_t, id_p), id_quxv)

    mean_parameter = np.full(len(circle.height), np.nan)
    m_parameter = np.full(len(circle.height), np.nan)
    c_parameter = np.full(len(circle.height), np.nan)

    Ns = np.full(
        len(circle.height), np.nan
    )  # number of sondes available for regression

    for k in range(len(circle.height)):
        Ns[k] = id_[:, k].sum()
        if Ns[k] >= 6:
            X_dx = circle["dx"].isel(height=k).isel(launch_time=id_[:, k]).values
            X_dy = circle["dy"].isel(height=k).isel(launch_time=id_[:, k]).values

            X = list(zip(X_dx, X_dy))

            Y_parameter = (
                circle[parameter].isel(height=k).isel(launch_time=id_[:, k]).values
            )

            regr = linear_model.LinearRegression()
            regr.fit(X, Y_parameter)

            mean_parameter[k] = regr.intercept_
            m_parameter[k], c_parameter[k] = regr.coef_

    # if "sondes_regressed" not in circle:
    #     circle["sondes_regressed"] = (["height"], Ns)

    return (mean_parameter, m_parameter, c_parameter, Ns)


# rename_dict = {
#     "u_wind": "u",
#     "v_wind": "v",
#     "specific_humidity": "q",
#     "temperature": "T",
#     "pressure": "p",
# }


def regress_for_all_parameters(circle, list_of_parameters):

    for par in list_of_parameters:
        (par_mean, par_m, par_c, Ns) = run_regression(circle, par)

        circle[par] = (["height"], par_mean)
        circle["d" + par + "dx"] = (["height"], par_m)
        circle["d" + par + "dy"] = (["height"], par_c)

        if "sondes_regressed" not in list(circle.data_vars):
            circle["sondes_regressed"] = (["height"], Ns)

    return print(f"Finished all parameters ...")


def regress_for_all_circles(circles, list_of_parameters):

    # mean = [None] * len(circles)
    # m = [None] * len(circles)
    # c = [None] * len(circles)

    map_iterators = map(
        regress_for_all_parameters,
        circles,
        [list_of_parameters for _ in range(len(circles))],
    )

    for id_, i in enumerate(map_iterators):
        i
        print(f"{id_+1}/{len(circles)} circles completed ...")

    return print(f"Regressed over all circles ...")


def get_div_and_vor(circles):

    for circle in circles:

        D = circle.dudx.values + circle.dvdy.values

        vor = circle.dvdx.values - circle.dudy.values

        circle["D"] = (["height"], D)
        circle["vor"] = (["height"], vor)

    return print("Finished estimating divergence and vorticity for all circles....")


def get_density_vertical_velocity_and_omega(circles):

    for circle in circles:
        den_m = [None] * len(circle.launch_time)

        for sounding in range(len(circle.launch_time)):
            mr = mpcalc.mixing_ratio_from_specific_humidity(
                circle.isel(launch_time=sounding).q.values
            )
            den_m[sounding] = mpcalc.density(
                circle.isel(launch_time=sounding).p.values * units.hPa,
                circle.isel(launch_time=sounding).T.values * units.degC,
                mr,
            ).magnitude

        circle["density"] = (["launch_time", "height"], den_m)
        circle["mean_density"] = (["height"], np.nanmean(den_m, axis=0))

        D = circle.D.values
        mean_den = circle.mean_density

        nan_ids = np.where(np.isnan(D) == True)[0]

        w_vel = np.full(len(circle.height), np.nan)
        p_vel = np.full(len(circle.height), np.nan)

        w_vel[0] = 0
        last = 0

        for m in range(1, len(circle.height)):

            if m in nan_ids:
                w_vel[m] = np.nan
            else:
                w_vel[m] = w_vel[last] - D[m] * 10 * (m - last)
                last = m

        for n in range(1, len(circle.height)):

            p_vel[n] = -mean_den[n] * 9.81 * w_vel[n] * 60 * 60 / 100

        circle["W"] = (["height"], w_vel)
        circle["omega"] = (["height"], p_vel)

    return print("Finished estimating density, W and omega ...")


def get_advection(circles, list_of_parameters=["u", "v", "q", "T", "p"]):

    for id_, circle in enumerate(circles):
        adv_dicts = {}
        for var in list_of_parameters:
            adv_dicts[f"h_adv_{var}"] = -(circle.u * eval(f"circle.d{var}dx")) - (
                circle.v * eval(f"circle.d{var}dy")
            )
            circle[f"h_adv_{var}"] = (["height"], adv_dicts[f"h_adv_{var}"])

    return print("Finished estimating advection terms ...")


def get_circle_products(circles, list_of_parameters):

    # for id_,circle in enumerate(circles) :

    regress_for_all_circles(circles, list_of_parameters)

    get_div_and_vor(circles)

    get_density_vertical_velocity_and_omega(circles)

    get_advection(circles)

    return print(f"All circle products retrieved!")
