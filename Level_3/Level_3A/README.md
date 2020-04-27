# Level - 3A

<div id="TOC">
Contents :
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


<div id="variables"><h2><a href="#TOC">Added Variables</a></h2></div>

<h4>Launch Time</h4>

Level-3A data are of the trajectory type with a single timestamp associated with each sounding, i.e. the launch time. The sounding does not have time at all levels, because for the original sounding, time is the independent variable, and interpolating that does not make sense. If it is essential for the user to obtain time for the relevant levels here, the data are still available in the Level-1 files, and thus, can still be retrieved. It would not be recommended though, to use "time of recording" in conjunction with interpolated variables.

<div id="interpolation"><h2><a href="#TOC">Interpolation</a></h2></div>

<div id="gridding"><h2><a href="#TOC">Gridding</a></h2></div>

<div id="flags"><h2><a href="#TOC">Flags</a></h2></div>

<div id="low_height_flag"><h3><a href="#TOC">Low Flight Height Flag</a></h3></div>

<div id="cloud_flag"><h3><a href="#TOC">Cloud Flag</a></h3></div>



Precipitable water is computed before the interpolation of the soundings is carried out. Values of precipitable water is from measurements on the raw grid.

What should be the pressure level at upto which integration for PW is carried out? For now, no limit is specified. Thus, PW value is upto top of measured profile.

Static stability is estimated using method by Bluestein (1992).

For all estimations of saturated vapour pressure, method by Hardy (1998) was used.

Bluestein, H B. Synoptic-dynamic meteorology in midlatitudes: Volume 1, principles of kinematics and dynamics. United States: N. p., 1992. Web.

ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in the Range –100 to +100 C, Bob Hardy, Proceedings of the Third International Symposium on Humidity and Moisture, 1998