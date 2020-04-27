# Level - 3A

<div style="text-align: justify">

Level-3A is a dataset where all Level-2 soundings are gridded to a uniform vertical resolution of 10 m, along the geopotential height dimension, up to an altitude of 10 km.

<div id="TOC">
<b>Contents :</b>
    <ul>
        <li>
            <a href="#file_str">File Structure</a>
        </li>
        <li>
            <a href="#variables">Added Variables</a>
        </li>
        <li>
            <a href="#interpolation">Interpolation</a>
        </li>
        <li>
            <a href="#gridding">Gridding</a>
        </li>
        <li>
            <a href="#flags">Flags</a>
        </li>
            <ul>
                <li>
                    <a href="#low_height_flag">Low Flight Height Flag</a>
                </li>
                <li>
                    <a href="#cloud_flag">Cloud Flag</a>
                </li>
            </ul>
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

<h4>Launch Time</h4>

Level-3A data are of the trajectory type with a single timestamp associated with each sounding, i.e. the launch time. The soundings in Level-3 do not have time at all levels, because for the original sounding, time is the independent variable, and interpolating that does not make sense. If it is essential for the user to obtain time for the relevant levels here, the data are still available in the Level-1 files, and thus, can still be retrieved. It would not be recommended though, to use "time of recording" in conjunction with interpolated variables.

<h4>Specific humidity</h4>

For the estimation of saturated vapour pressure, method by <a href="#hardy_1998">Hardy (1998)</a> is used.

<h4>Precipitable water</h4>

Precipitable water (PW) is computed before the interpolation of the soundings is carried out. So, values of PW is from measurements on the raw grid. Since, there will only be a single value per sounding, it makes sense to stay true to the raw sounding, in this case.

The value of PW is obtained by integrating from surface up to the top of the measured column and thus, will naturally depend on the height of the atmospheric column in consideration. This means that PW cannot be compared across soundings without considering flight altitude. Generally since the moisture in the upper layers of the atmosphere is very low, this does not cause a major difference in PW, even if the flight altitude varies by a couple of km, above ~6-7 km. However, the difference between flight altitudes during EUREC<sup>4</sup>A sonde launches was significant. Almost all of HALO's sondes were launched from ~10 km height, whereas for P3, this height varied, with values between 2 km (AXBT launches) and 8 km. See more in the discussion of the <a href="#low_height_flag">low_height_flag variable</a>.

<h4>Static stability</h4>

Static stability is considered here as the gradient of potential temperature with respect to pressure. This value is estimated after the interpolation of variables to the common grid.

<h4>Platform</h4>

Although all soundings are in a single file in Level-3, they can still be separated into HALO and P3 sondes, using this variable, which specifies the platform from which the dropsonde was launched. The values of the variable are strings, and have two possible values - "HALO" and "P3".


<div id="gridding"><h2><a href="#TOC">Gridding</a></h2></div>

<div id="interpolation"><h2><a href="#TOC">Interpolation</a></h2></div>

<div id="flags"><h2><a href="#TOC">Flags</a></h2></div>

<div id="low_height_flag"><h3><a href="#TOC">Low Flight Height Flag</a></h3></div>

<div id="cloud_flag"><h3><a href="#TOC">Cloud Flag</a></h3></div>





What should be the pressure level at upto which integration for PW is carried out? For now, no limit is specified. Thus, PW value is upto top of measured profile.

Static stability is estimated using method by Bluestein (1992).

For all estimations of saturated vapour pressure, method by Hardy (1998) is used.

Bluestein, H B. Synoptic-dynamic meteorology in midlatitudes: Volume 1, principles of kinematics and dynamics. United States: N. p., 1992. Web.

<div if="hardy1998">ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in the Range –100 to +100 C, Bob Hardy, Proceedings of the Third International Symposium on Humidity and Moisture, 1998</div>