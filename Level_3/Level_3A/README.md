# Level - 3A

<div style="text-align: justify">

Level-3A is a dataset where all Level-2 soundings are gridded to a uniform vertical resolution of 10 m, along the geopotential height dimension, up to an altitude of 10 km. For Level-3A, the P3 and HALO dropsondes are also integrated into a single dataset.

<div id="TOC">
<b>Contents :</b>
    <ul>
        <li>
            <a href="#file_str">File Structure</a>
        </li>
        <li>
            <a href="#variables">Added Variables</a>
        </li>
<!--             <ul>
                <li>
                    <a href="#launch_time">Launch Time</a>
                </li>
                <li>
                    <a href="#potential_temperature">Potential Temperature</a>
                </li>
                <li>
                    <a href="#specific_humidity">Specific Humidity</a>
                </li>
                <li>
                    <a href="#precipitable_water">Precipitable Water</a>
                </li>
                <li>
                    <a href="#static_stability">Static Stability</a>
                </li>
                <li>
                    <a href="#platform">Platform</a>
                </li>
            </ul> -->
        <li>
            <a href="#gridding">Gridding</a>
        </li>
        <li>
            <a href="#interpolation">Interpolation</a>
        </li>
        <li>
            <a href="#flags">Flags</a>
        </li>
<!--             <ul>
                <li>
                    <a href="#low_height_flag">Low Flight Height Flag</a>
                </li>
                <li>
                    <a href="#cloud_flag">Cloud Flag</a>
                </li>
            </ul> -->
    </ul>
</div>

<div id="file_str"><h2><a href="#TOC"> File Structure</a></h2></div>

<table>
<tbody>
<tr>
<td><b>OBJECT</b></td>
<td><b>NAME</b></td>
<td><b>UNITS</b></td>
<td><b>DIMENSIONS</b></td>
</tr>
<tr>
<td><b>Dimensions</b></td>
<td>obs</td>
<td>&nbsp;</td>
<td>obs</td>
</tr>
<tr>
<td>&nbsp;</td>
<td>sounding</td>
<td></td>
<td>sounding</td>
</tr>
<tr>
<td><b>Coordinates</b></td>
<td>launch time</td>
<td>seconds since<br>1970-01-01 00:00:00 UTC</td>
<td>sounding</td>
</tr>
<tr>
<td>&nbsp;</td>
<td>height </td>
<td>m</td>
<td>obs</td>
</tr>
<tr>
<td>&nbsp;</td>
<td>latitude </td>
<td>&#176N</td>
<td>obs, sounding</td> 
</tr>
<tr>
<td>&nbsp;</td>
<td>longitude </td>
<td>&#176E</td>
<td>obs, sounding</td> 
</tr>
<tr>
<td><b>Variables</b></td>
<td>pressure </td>
<td>hPa</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>temperature </td>
<td>&#x2103;</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>relative humidity </td>
<td>%</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>wind speed </td>
<td>m / s</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>wind direction </td>
<td>&#176</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>GPS altitude<sup>*</sup></td>
<td>m</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>vertical velocity<sup>*</sup></td>
<td>m / s</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>horizontal wind (u)</td>
<td>m / s</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>horizontal wind (v)</td>
<td>m / s</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>potential temperature</td>
<td>K</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>specific humidity</td>
<td>kg / kg</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>precipitable water</td>
<td>kg / m<sup>2</sup></td>
<td>sounding</td>
</tr>
<tr>
<td></td>
<td>static stability</td>
<td>K / hPa</td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>low height flag</td>
<td></td>
<td>sounding</td>
</tr>
<tr>
<td></td>
<td>cloud flag</td>
<td></td>
<td>obs, sounding</td>
</tr>
<tr>
<td></td>
<td>platform</td>
<td></td>
<td>sounding</td>
</tr>
<tr>
<td></td>
<td>flight height</td>
<td>m</td>
<td>sounding</td>
</tr>
<tr>
<td></td>
<td>flight latitude</td>
<td>&#176N</td>
<td>sounding</td>
</tr>
<tr>
<td></td>
<td>flight longitude</td>
<td>&#176E</td>
<td>sounding</td>
</tr>
</tbody>
</table>

<sup>*</sup> - variables haven't been included in the data in the current version, as they need to go through further quality checks.

<div id="variables"><h2><a href="#TOC">Added Variables</a></h2></div>

<!-- <div id="launch_time"><a href="#TOC"> -->
<h4>Launch Time</h4>
<!-- </a></h4></div> -->

Level-3A data are of the trajectory type with a single timestamp associated with each sounding, i.e. the launch time. The soundings in Level-3 do not have time at all levels, because for the original sounding, time is the independent variable, and interpolating that does not make sense. If it is essential for the user to obtain time for the relevant levels here, the data are still available in the Level-1 files, and thus, can still be retrieved. It would not be recommended though, to use "time of recording" in conjunction with interpolated variables.

<!-- <div id="potential_temperature"><a href="#TOC"> -->
<h4>Potential Temperature</h4>

The values of potential temperature are estimated from the sounding profile on their respective, raw vertical grid, before interpolating them on to a common grid.

<!-- <div id="specific_humidity"><a href="#TOC"> -->
<h4>Specific Humidity</h4>

For the estimation of saturated vapour pressure, method by <a href="#hardy1998">Hardy (1998)</a> is used. Specific humidity is estimated from the sounding profile on its raw vertical grid, before interpolation.

<!-- <div id="precipitable_water"><a href="#TOC"> -->
<h4>Precipitable Water</h4>
    
Precipitable water (PW) is computed before the interpolation of the soundings is carried out. So, values of PW is from measurements on the raw grid. Since, there will only be a single value per sounding, it makes sense to stay true to the raw sounding, in this case.

The value of PW is obtained by integrating from surface up to the top of the measured column and thus, will naturally depend on the height of the atmospheric column in consideration. This means that PW cannot be compared across soundings without considering flight altitude. Generally since the moisture in the upper layers of the atmosphere is very low, this does not cause a major difference in PW, even if the flight altitude varies by a couple of km, above ~6-7 km. However, the difference between flight altitudes during EUREC<sup>4</sup>A sonde launches was significant. See more in the discussion of the <a href="#low_height_flag">`low_height_flag`</a> variable.

<!-- <div id="static_stability"><a href="#TOC"> -->
<h4>Static Stability</h4>

Static stability is considered here as the gradient of potential temperature with respect to pressure. This value is estimated after the interpolation of variables to the common grid.

<!-- <div id="platform"><a href="#TOC"> -->
<h4>Platform</h4>

Although all soundings are in a single file in Level-3, they can still be separated into HALO and P3 sondes, using this variable, which specifies the platform from which the dropsonde was launched. The values of the variable are strings, and have two possible values - `"HALO"` and `"P3"`.

<div id="gridding"><h2><a href="#TOC">Gridding</a></h2></div>

The primary objective behind the Level-3A product is gridding all soundings on a common, vertical grid, thus making it easier to use the soundings for different analyses in bulk. The vertical grid spacing for the dataset is kept at 10 m, up to an altitude of 10 km. 

In the case of a regular drop, i.e. if there are no issues like a fast fall, or a failed parachute, the average descent rate of the dropsondes is ~21 m/s at 12 km altitude and ~11 m/s near to the surface<a href="#vaisala_datasheet"><sup>1</sup></a>. The PTU (pressure, temperature and humidity(U)) sensors have a measurement frequency of 2 Hz, while the GPS has a 4 Hz measurement frequency. This would translate to a vertical resolution of roughly 9-10 m at HALO's flight altitude, and 5-6 m close to the surface for the PTU values. For the more frequently measured wind values, the resolution will accordingly be higher. However, this vertical spacing will vary depending on the actual vertical wind velocity.

Taking 7-8 m to be the typical vertical resolution of PTU measurements, at a grid of 10 m spacing, this would be a slight downscaling of the raw measurements. For the wind values, however, this would mean neglecting almost half the measured values. An option is to keep the wind values at 5 m grid, and the others at a 10 m grid, although this would mean that it becomes slightly inconvenient to work with the dataset as a whole. This is left to the users' feedback for now. 

<div id="interpolation"><h2><a href="#TOC">Interpolation</a></h2></div>
The interpolation to the common grid is carried out through the following steps: 

(i) Variables `specific_humidity`, `potential_temperature`, `u_wind`, `v_wind`, `precipitable_water` and `static stability` are added to the dataset.

(ii) All variables along `height` coordinates in dataset are linearly interpolated along the `height` dimension, at specified height intervals (default 10 m) and up to specified altitude (default 10 km) 

(iii) Pressure values are interpolated using a logarithmic interpolation scheme and these values replace the linearly interpolated pressure values. (Note : The difference between these different interpolations is in the order of 0.005 hPa, which is lower than the measurement uncertainty of the pressure sensor itself of the RD-41 itself)

(iv) Temperature and moisture are interpolated with values of theta and q, respectively. Thus, after interpolation, variables `temperature` and `relative_humidity` are recomputed from the interpolated values of `potential_temparature` and `specific_humidity`. The new values for T and RH will replace the previously interpolated T and RH variables from the raw sounding. Although, T and RH are the originally measured properties by the dropsonde sensors, for interpolation q and theta are preferred, as these variables are conserved.

(v) At this step of interpolation, the `time` variable is dropped from the dataset, since in the original files, time is an independent variable, and interpolating time here will not work in the same way as for all other variables. For the gridded product, time will be an artificial residue, if interpolated, and for almost all practical purposes, the sounding can be seen as a snapshot in time, with only the launch time being a relevant timestamp for every sonde.

<div id="flags"><h2><a href="#TOC">Flags</a></h2></div>

<h3>Low Flight Height Flag</h3>

- Since HALO and P3 had significantly different objectives and strategies for their respective flights in EUREC4A, they launched dropsondes from different altitudes. HALO typically had a flight altitude of ~10 km, when launching sondes, but for P3 this varied between 7.5 km and 2.5 km. The P3 dropped some sondes with its AXBT launches, which were at an altitude of ~2.5-3 km. This essentially means that these soundings sampled only the very low levels of the atmosphere, and had just half of P3's other sondes, and a third of HALO's typical sondes.

- The `low_height_flag` marks these sondes that have a launch altitude of less than 4 km, with a value of 1. Sondes with this flag's value as 0 have a launch altitude of at least 4 km. This flag is essential to put in to context estimates of integrated quantities such as total column moisture, as well as to act as an easy separator for users who want to look at profiles also above the typical inversion height. 

<h3>Cloud Flag</h3>

- There is considerable interest in classifying the soundings that passed through cloud/s and those that did not. For this purpose, a `cloud_flag` is part of the dataset. 

- This `cloud_flag` detects if there is a cloud present at each level, by passing through an algorithm (outlined in <a href="#zhang2010">Zhang (2010)</a>) that is mainly dependent upon RH thresholds as a function of altitude. There are some other checks also in place, to ensure that noisy signal don't lead to an inclusion (omission) of a non-cloud(potential cloud).

- Since the humidity sensor in the RD-41 dropsondes (deployed by HALO and P3 during EUREC<sup>4</sup>A) is much more reliable than previous sensors, on which the aforementioned algorithm has been tested, this provides some confidence in the `cloud_flag`. The flag can also be tested against (almost) collocated radar products from the airplane, which can show how well the flag is functioning. However, since the radar products will not be available till at least the end of this year (2020), it might be a good idea for the release of JOANNE. The satellite products might be faster in this sense, and can probably used to check to what extent the `cloud_flag` can be relied upon.

---
References:

<div id="vaisala_datasheet"><sup>1</sup> : <url>https://www.vaisala.com/sites/default/files/documents/RD41-Datasheet-B211706EN.pdf</url></div>
<br/>
<div id="hardy1998">Hardy (1998) : ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in the Range –100 to +100 C, Bob Hardy, Proceedings of the Third International Symposium on Humidity and Moisture, 1998</div>
<br/>
<div id="zhang2010">Zhang (2010) : <url>https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JD014030</url></div>
