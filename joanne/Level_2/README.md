<div style="text-align: justify">

# Level 2

<div id="TOC">
Contents :
    <ul>
        <li>
            <a href="#file_str">File Structure</a>
        </li>
        <li>
            <a href="#quality_control">Quality Control</a>
        </li>
            <ul>
                <li>
                    <a href="#sonde_classification">Sonde Classification and QC Terminology</a>
                </li>
                <li>
                    <a href="#quality_checks">Quality Checks</a>
                </li>
                    <ul>
                        <li>
                            <a href="#ind_flags">ind_FLAG</a>
                        </li>
                        <li>
                            <a href="#srf_flags">srf_FLAG</a>
                        </li>
                        <li>
                            <a href="#ld_flag">ld_FLAG</a>
                        </li>
                        <li>
                            <a href="#flag">FLAG</a>
                        </li>
                    </ul>
            </ul>
            <li>
            <a href="#output">Output</a>
            </li>
                <ul>
                    <li>
                        <a href="#lv2_data_files">Level-2 Data Files</a>
                    </li>
                    <li>
                        <a href="#status_file">Status File</a>
                    </li>
                    <li>
                        <a href="#log_file">Log File</a>
                    </li>
                </ul>
    </ul>
</div>
<div id="file_str">
    <h1>
        <a href="#TOC">File Structure</a>
    </h1>
    <p>
        The Level-2 NC files are data from individual soundings, which passed the Level-2 QC checks applied after processing the raw file with ASPEN v3.4.3.

For Level-2, only variables that are measurements from the dropsonde sensors are included. There are no redundant state variables. 
    </p>

<table>
<tbody>
<tr>
<td><b>OBJECT</b></td>
<td><b>NAME</b></td>
<td><b>UNITS</b></td>
</tr>
<tr>
<td><b>Dimension</b></td>
<td>obs</td>
<td>&nbsp;</td>
</tr>
<tr>
<td><b>Coordinates</b></td>
<td>time (independent) </td>
<td>seconds since<br>1970-01-01 00:00:00 UTC</td>
</tr>
<tr>
<td>&nbsp;</td>
<td>height </td>
<td>m</td>
</tr>
<tr>
<td>&nbsp;</td>
<td>latitude </td>
<td>degree_north</td> 
</tr>
<tr>
<td>&nbsp;</td>
<td>longitude </td>
<td>degree_east</td> 
</tr>
<tr>
<td><b>Variables</b></td>
<td>pressure </td>
<td>hPa</td>
</tr>
<tr>
<td></td>
<td>temperature </td>
<td>&#x2103;</td>
</tr>
<tr>
<td></td>
<td>relative_humidity </td>
<td>%</td>
</tr>
<tr>
<td></td>
<td>wind_speed </td>
<td>m s-1</td>
</tr>
<tr>
<td></td>
<td>wind_direction </td>
<td>degree</td>
</tr>
<tr>
<td></td>
<td>GPS_altitude<sup>*</sup></td>
<td>m</td>
</tr>
<tr>
<td></td>
<td>vertical_velocity<sup>*</sup></td>
<td>m s-1</td>
</tr>
</tbody>
</table>

<sup>*</sup> - variables haven't been included in the data in the current version, as they need to go through further quality checks.

Note that file names of the soundings are not always indicative of launch times, although this will be the case for most soundings. The attribute `Launch_time_(UTC)` in every sounding file should be considered as the final authority on launch time.

</div>
<div id="quality_control">
    <h1>
        <a href="#TOC">Quality Control</a>
    </h1>
    <p style="text-align: justify;">
        The primary aim of the quality control process is to classify sondes based on their level of success in making measurements. We discard the sondes with zero or negligible usable data, and keep the remaining sondes as part of the Level-2 product.     </p>
    <p style="text-align: justify;">
        Most sondes have complete sounding profiles, while some have partially faulty profiles, thus allowing us to include some parts of the sounding that are still useful. An example is if the PTU sensor of a sonde failed, we will still include the sounding in Level-2 for the data that the GPS sensor provides.
    </p>
</div>

<div id="sonde_classification">
    <h2>
        <a href="#TOC">Sonde classification and QC terminology</a>
    </h2>

<b>`flag`</b> : record of result of local test  
<b>`good`</b> : passes local test  
<b>`bad`</b>  : fails local test  
<b>`ugly`</b> : fails preliminary local test, but some data may be salvaged  

<b>`FLAG`</b> : record of result of group of tests or final flag for sonde  
<b>`GOOD`</b> : sonde is good, passed all tests and can be used straightaway  
<b>`BAD`</b>  : sonde is bad, should not be considered for any further data processing or analysis  
<b>`UGLY`</b> : sonde has some data that may be salvaged later, but cannot be processed straightaway. Needs more QC.


<div id="quality_checks">
    <h2>
       <a href="#TOC"> Quality Checks</a>
    </h2>
</div>

There are 3 individual processes used to classify sondes, which are combined later to make a final grouping:

1. <div id="ind_flags">
    <h3>
       <a href="#TOC">ind_FLAGs:</a>
    </h3>
</div> 

- Ratio of individual parameter's count to the total count of records (time is used as proxy for measurement record. Multiple parameters can have values at a single record.)
    - There is one flag each for the following parameters:  
    `t_flag`   : *temperature (tdry)*  
    `p_flag`   : *pressure (pres)*  
    `rh_flag`  : *relative humidity (rh)*  
    `z_flag`   : *GPS altitude (gpsalt)*  
    `u_flag`   : *u wind (u_wind)*  
    `v_flag`   : *v wind (v_wind)*  

- Time values are recorded every 0.25 seconds. Although, the PTU and GPS sensors have a measurement frequency of 2 Hz and 4 Hz, respectively, the distribution of measurements vary slightly from the ideal case - which is for all parameters (except u,v) to have measurements at every other time record, and for u,v to have measurements at every time record. Since, the time records also include values during initialisation as well as during a little before and after the launch, when no signal can be sent back to the AVAPS PC, the actual ratio will always be lower than the ideal estimate of 1 (for u,v) and 0.5 (for the remaining parameters).

- The true distribution shows that peaks start to flatten around 0.8 and 0.4 for u,v and other parameters, respectively. Thus, sondes with ratios lower than these values are taken as not having a complete profile, and termed as `ugly` sondes. These ugly sondes still have data, but because they are expected to have more NaN fields than most sondes, they are kept for more QC and NaN-filling later, depending upon the extent of the dearth of measurements in that sonde.

- For each parameter, if this proportional ratio is higher than the set threshold, then the sonde is flagged as `good`. If the ratio is lower than set threshold, but non-zero, then the sonde is flagged `ugly` and if the ratio is zero, i.e. this particular parameter was not measured at all for the sonde, then, it is flagged as `bad`. (This value might later be changed from zero to 10% of the time counts).
  
  - If all individual ind_flags are `good`, the sonde is flagged as `GOOD` for `ind_FLAG`,  
  - if all individual ind_flags are `bad`, the sonde is flagged as `BAD` for `ind_FLAG`,  
  - if neither of these conditions is met, the sonde is flagged as `UGLY` for `ind_FLAG`.
    
<div style="text-align: justify">

----
2. <div id="srf_flags">
    <h3>
       <a href="#TOC">srf_FLAG:</a>
    </h3>
</div> 

- These flags act mostly to identify the sondes' behaviour in the lower layer of the atmosphere (mostly near surface, but also < 4 km for some tests) as sanity checks. Functions are defined to check for obvious errors in near-surface values of parameters, such as being out of bounds (limits of realistic values), poor agreement between estimates of altitude, etc.
  
- <div id="srf_p_flag">
    <h4>
       <a href="#TOC">srf_p_flag:</a>
    </h4>
  </div>
    
    - This flag checks if maximum pressure measured by sonde is within bounds: 1000 hPa - 1020 hPa.
  
    - If the value is higher than bound, it is unrealistic, and value lower than bound means sonde did not measure the bottommost levels of the atmosphere.
    
    - This flag does not check any GPS values. Even if there were no pressure values above 1000 hPa, there may still be GPS measurements in the lowest levels. Such sondes can still be useful for wind and wind-derived products.

- <div id="srf_t_flag">
    <h4>
       <a href="#TOC">srf_t_flag:</a>
    </h4>
  </div>

    - This flag checks if tdry (air temperature) is within bounds:
    
      1. Maximum air temperature recorded should not be greater than the upper limit, set to a default value of 30 deg C.
      2. Mean air temperature in the bottom 100 m (by gpsalt) should not be lesser than surface limit, set to a default of 20 deg C.
    
    - If any of the above limits is violated, the tdry for the sonde is considered out of bounds, and marked as `False`. The sonde is also marked `False`, if there are no measurements in the bottom 100 m (by gpsalt).

- <div id="srf_rh_flag">
    <h4>
       <a href="#TOC">srf_rh_flag:</a>
    </h4>
  </div>

    - Checking if rh (relative humidity) is within bounds:
    
      1. Mean RH in the bottom 100 m (by gpsalt) should not be lesser than srf_limit, set to a default of 50 %.
    
    - If the above limit is violated, the rh for the sonde is considered out of bounds, and marked as `False`. The sonde is also marked `False`, if there are no measurements in the bottom 100 m (by gpsalt). 
    
- <div id="srf_z_flag">
    <h4>
       <a href="#TOC">srf_z_flag:</a>
    </h4>
  </div>

  - This flag checks if maximum GPS altitude of sonde is within bounds: <= limit (default assigned as 30 m). Value higher than bound means there are no near-surface measurements
    
  - This flag does not include any pressure values. Even if there were no GPS values below 30 m,
    there may still be PTU measurements in the lowest levels. 

- <div id="palt_gpsalt_rms_flag">
    <h4>
       <a href="#TOC">palt_gpsalt_rms_flag:</a>
    </h4>
  </div>

    - This function estimates the root mean square (RMS) difference between geopotential altitude (palt) and the GPS altitude (gpsalt),for values below 4 km, and based on a limit (rms_limit; 
    set to a default value of 100 m), is flagged accordingly.
    
    - If the estimated RMS difference is below the limit, then the sonde is flagged as `True` for this test.
    
    - If the estimated RMS difference is greater than the limit, or if there are no values of either palt or gpsalt
    overlapping in the lower 4 km, then the sonde is flagged as `False` for this test. The lack of overlap could be 
    because either there are no palt values or no gpsalt values or both.
    
    - If the above limit is violated, the rh for the sonde is considered out of bounds,
    and marked as `False`. The sonde is also marked `False`, if there are no measurements in the 
    bottom 100 m (by gpsalt). 
---
- <div id="srf_flag">
    <h4>
       <a href="#TOC"><b>`srf_FLAG`:</b></a>
    </h4>
  </div>

  - For each of these srf_flags, if the sonde passes the test, it is marked as 1, else as 0, which stand for `good` and `bad` respectively. 
    - If all individual srf_flags are `good`, the sonde is flagged as `GOOD` for `srf_FLAG`,  
    - if all individual srf_flags are `bad`, the sonde is flagged as `BAD` for `srf_FLAG`,  
    - if neither of these conditions is met, the sonde is flagged as `UGLY` for `srf_FLAG`.  
 ---     
<div style="text-align: justify">

3. <div id="ld_flag">
    <h3>
       <a href="#TOC">ld_FLAG:</a>
    </h3>
</div> 

- Checking whether the sonde detected a launch automatically. If the sonde fails to do this, it does not switch to high-power signal transmission, and thus, stops sending data back to the AVAPS PC, after a short range. 
    
- The primary method to check launch detection is to parse through the log files of the sonde in raw format. These files have names starting with 'A' and are followed by the date and time of launch. The file extension is the number of the channel used to initialise the sonde and receive its signal, but for all practical purposes, it is a .txt file. (Note: For sondes that did not detect a launch, the file name has time when the sonde was initialised). The log file contains an internal record termed 'Launch Obs Done?'. If this value is 1, the launch was detected, else if it is 0, launch was not detected. The same values are used to mark the `ld_FLAG`.

- An alternative method to check launch detection is through an attribute of the PQC files, called 'reference_time'. If a launch was detected, this attribute shows the launch time. If a launch was not detected automatically, the attribute shows the correct date, but the time is 'T00:00:00.000'. This can also be used to mark the `ld_FLAG`. However, this should only be used as a quick-fix if the log files are not available, since this test is not fool-proof. Moreover, PQC files with no automatic launch detection do not have relevant sonde information such as SondeID added to them. 

<div style="text-align: justify">

---
<div id="flag">
    <h3>
       <a href="#TOC">FLAG:</a>
    </h3>
</div> 

- After the sondes pass through these three processes, the aforementioned `FLAG`s are used to determine the final `FLAG` value for the sonde.  
    - If `ld_FLAG` is 0 ==> the sonde FLAG is termed `BAD`  
    - If `ld_FLAG` is 1, and if both `srf_FLAG` and `ind_FLAG` are `GOOD`, then the sonde FLAG is termed `GOOD`.  
    - If `ld_FLAG` is 1, and if both `srf_FLAG` and `ind_FLAG` are `BAD`, then the sonde FLAG is termed `BAD`.  
    - If `ld_FLAG` is 1, and if `srf_FLAG` and `ind_FLAG` have different values, then the sonde FLAG is termed `UGLY`.   
      
- Although the process of classifying the sondes can be simplified by other combinations of the `ind_flag`s and `srf_flag`s, the current method ensures no good sondes are omitted, and no bad sondes are admitted. The rest of the sondes, the ugly sondes, still have data that can be salvaged, and after some QC and/or flagging, can be combined with the other good sondes.  
  
- The NC file generated as a product stores results for each individual test mentioned above, group of tests and the final classification. Thus, the user can still mould the classification based on their objectives or add/remove tests to the process and customise it for themselves.
---
<div id="output">
    <h1>
       <a href="#TOC">Output</a>
    </h1>
</div>

1. <div id="lv2_data_files">
    <h2>
       <a href="#TOC">Level-2 Data files</a>
    </h2>
</div>

- These files are the individual sounding data files, created in the format mentioned in <a href="file_str">file structure</a>. All files include flight variables such as position, height, speed, etc. in their metadata, among other information.
- The file names are in the format:
   
  ```
  [campaign]_[platform]_[instrument]_[launch date in the format YYYYMMDD]_[launch time in the format HHMMSS].nc

  e.g. EUREC4A_HALO_Dropsonde-RD41_20200128_195748.nc
    ```

2. <div id="status_file">
    <h2>
       <a href="#TOC">Status File</a>
    </h2>
</div>

- The status file is an NC file generated after the Level-2 QC process, and for every sonde launched during EUREC4A, it includes the following data:
  -  total measurement counts for all parameters (which goes into the estimate of `ind_flag`s)
  -  all individual ind_flag values and the combined `ind_FLAG` value
  -  all individual srf_flag values and the combined `srf_FLAG` value
  -  `ld_FLAG` value
  -  `FLAG` value  

- One status file each is generated  for HALO and P3.
  
3. <div id="log_file">
    <h2>
       <a href="#TOC">Log File</a>
    </h2>
</div>

- The log file is a .txt file which records the details of sondes that failed to detect an automatic launch. For the failed sondes, since these details are not available from the Level-1 .PQC files directly, this file is kept as a log. The file also mentions the total sondes that failed to detect a launch.
  
- The file also mentions how many sondes were classified as "good","bad" and "ugly" as per each test of `ind_FLAG`, `srf_FLAG` and `FLAG`.
  
- One log file each is generated for HALO and P3.
