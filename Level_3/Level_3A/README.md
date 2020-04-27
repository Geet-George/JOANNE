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