def get_total_non_nan_indices(sonde):

    """
    Retrieving the non-NaN indices for all parameters.
    
    'c' terms are complete arrays of indices that have non-NaN values
    's' terms are the counts of indices in the respective 'c' terms.
    
    
    Input : 
        sonde_path : Path to the sonde PQC file as a string
    Output :
        s_var : where, var is one of [time,t,rh,p,z,u,v]
        
    """

    #     import xarray as xr
    import numpy as np

    #     sonde = xr.open_dataset(sonde_path)

    c_time = ~np.isnan(sonde.time).values
    c_t = ~np.isnan(sonde.tdry).values
    c_rh = ~np.isnan(sonde.rh).values
    c_p = ~np.isnan(sonde.pres).values
    c_z = ~np.isnan(sonde.gpsalt).values
    c_u = ~np.isnan(sonde.u_wind).values
    c_v = ~np.isnan(sonde.v_wind).values

    s_time = c_time.sum()
    s_t = c_t.sum()
    s_rh = c_rh.sum()
    s_p = c_p.sum()
    s_z = c_z.sum()
    s_u = c_u.sum()
    s_v = c_v.sum()

    return s_time, s_t, s_rh, s_p, s_z, s_u, s_v


def pres_bounds(sonde, u_lim=1020, l_lim=1000):

    """
    Checking if maximum pressure of sonde is within bounds: 1000 hPa - 1020 hPa. 
    Value higher than bound is unrealistic, 
    and value lower than bound means sonde did not measure the bottommost levels of the atmosphere.
    
    This flag does not include any GPS values. Even if there were no pressure values above 1000 hPa, 
    there may still be GPS measurements in the lowest levels. 
    Such sondes can still be useful for wind and wind-derived products.
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if max pressure value is within bounds
               False, if max pressure value is out of bounds
    """

    if (sonde.pres.max() < 1000) | (sonde.pres.max() > 1020):
        return False
    else:
        return True


def gps_bounds(sonde, limit=30):

    """
    Checking if maximum GPS altitude of sonde is within bounds: <= limit (default assigned as 30 m)
    Value higher than bound has no near-surface measurements
    
    This flag does not include any pressure values. Even if there were no GPS values below 30 m,
    there may still be PTU measurements in the lowest levels. 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if min GPS altitude value is within bounds
               False, if min GPS altitude value is out of bounds
        
    """

    if sonde.gpsalt.min() > limit:
        return False
    else:
        return True


def tdry_bounds(sonde, u_limit=30, srf_limit=20):
    """
    Checking if tdry (air temperature) is within bounds:
    
    1. Maximum air temperature recorded should not be greater than the upper limit (u_limit), 
        set to a default value of 30 deg C.
    2. Mean air temperature in the bottom 100 m (by gpsalt) should not be lesser than srf_limit, 
        set to a default of 20 deg C.
    
    If any of the above limits is violated, the tdry for the sonde is considered out of bounds,
    and marked as 'False'. The sonde is also marked 'False', if there are no measurements in the 
    bottom 100 m (by gpsalt). 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if tdry value is within bounds
               False, if tdry is out of bounds or if measurements are not available
    """
    if sonde.tdry.max() >= u_limit:
        return False
    elif sonde.tdry.where(sonde.gpsalt < 100, drop=True).sum() == 0:
        return False
    elif sonde.tdry.where(sonde.gpsalt < 100, drop=True).mean() < srf_limit:
        return False
    else:
        return True


def rh_bounds(sonde, srf_limit=50):
    """
    Checking if rh (relative humidity) is within bounds:
    
    1. Mean RH in the bottom 100 m (by gpsalt) should not be lesser than srf_limit, 
        set to a default of 50 %.
    
    If the above limit is violated, the rh for the sonde is considered out of bounds,
    and marked as 'False'. The sonde is also marked 'False', if there are no measurements in the 
    bottom 100 m (by gpsalt). 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if rh value is within bounds
               False, if rh is out of bounds or if measurements are not available
    """
    if sonde.rh.where(sonde.gpsalt < 100, drop=True).sum() == 0:
        return False
    elif sonde.rh.where(sonde.gpsalt < 100, drop=True).mean() < srf_limit:
        return False
    else:
        return True


def palt_gpsalt_rms_check(sonde, rms_limit=100):
    """
    This function estimates the root mean square (RMS) difference between geopotential altitude (palt) 
    and the GPS altitude (gpsalt),for values below 4 km, and based on a limit (rms_limit; 
    set to a default value of 100 m), is flagged accordingly.
    
    If the estimated RMS difference is below the limit, then the sonde is flagged as 'True' for this test.
    
    If the estimated RMS difference is greater than the limit, or if there are no values of either palt or gpsalt
    overlapping in the lower 4 km, then the sonde is flagged as 'False' for this test. The lack of overlap could be 
    because either there are no palt values or no gpsalt values or both.
    
    If the above limit is violated, the rh for the sonde is considered out of bounds,
    and marked as 'False'. The sonde is also marked 'False', if there are no measurements in the 
    bottom 100 m (by gpsalt). 
    
    Input : 
        sonde : Opened xarray dataset of ASPEN-processed PQC dropsonde file
    Output :
        bool : True, if RMS below limit
               False, if RMS above limit or no overlapping values
    """

    x = (
        (
            sonde.alt.where(sonde.alt < 4000, drop=True)
            - sonde.gpsalt.where(sonde.alt < 4000, drop=True)
        )
        ** 2
    ).values

    nn = ~np.isnan(x)

    if nn.sum() == 0:
        return False  # all x are NaNs
    else:
        zdiff_rms = np.sqrt(np.nanmean(x))
        if zdiff_rms < rms_limit:
            return True
        else:
            return False


# Function to check if sonde failed due to no detection of launch


def check_launch_detect(sonde_path):
    """
    Alternative method to check automatic launch detection of the sonde
    
    Input : path to PQC file
    
    Output : If launch detected, True; else, False
    
    Function to check if the dropsonde detected no launch, and thus failed. This was a common
    cause of dropsonde failure during EUREC4A. Vaisala's best guess is that for some reason,
    the IR sensor near the parachute did not detect the parachute coming out, thus not 
    detecting a launch and thus, not switching from low-power transmission to high-power 
    transmission. This caused AVAPS to lose the sonde's signal very soon after 
    launch (~350 hPa).
    """

    xrdataset = xr.open_dataset(sonde_path)

    if str(xrdataset.reference_time.values)[-21:-2] == "T00:00:00.000000000":
        launch_detect_flag = (
            False  # this means the sonde failed and launch was not detected
        )
        print(
            xrdataset.attrs["SoundingDescription"][21:30],
            xrdataset.reference_time.values,
            xrdataset.launch_time.values,
        )
    else:
        launch_detect_flag = True  # this means launch was detected

    return launch_detect_flag


def create_variable(ds, vname, data, **kwargs):
    """Insert the data into a variable in an :class:`xr.Dataset`"""
    attrs = nc_meta[vname].copy()
    dims = ["time"]  # nc_dims[vname]

    v = xr.Variable(dims, np.asarray(data), attrs=attrs)
    ds[vname] = v

    return vname