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
            </ul>
    </ul>
</div>
<div id="file_str">
    <h2>
        <a href="#TOC">File Structure</a>
    </h2>
    <p>
        The Level-2 NC files are data from individual soundings, which passed the Level-2 QC checks applied after processing the raw file with ASPEN v3.4.3.

For Level-2, only variables that are measurements from the dropsonde sensors are included. There are no redundant state variables. 
    </p>

<table>
<tbody>
<tr>
<td><b>Dimension</b></td>
<td>obs</td>
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
<td>&#176N</td> 
</tr>
<tr>
<td>&nbsp;</td>
<td>longitude </td>
<td>&#176E</td> 
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
<td>relative humidity </td>
<td>%</td>
</tr>
<tr>
<td></td>
<td>wind speed </td>
<td>m/s</td>
</tr>
<tr>
<td></td>
<td>wind direction </td>
<td>&#176</td>
</tr>
<tr>
<td></td>
<td>GPS altitude<sup>*</sup></td>
<td>m</td>
</tr>
<tr>
<td></td>
<td>vertical velocity<sup>*</sup></td>
<td>m/s</td>
</tr>
</tbody>
</table>

<sup>*</sup> - variables haven't been included in the data in the current version, as they need to go through further quality checks.

</div>
<div id="quality_control">
    <h2>
        <a href="#TOC">Quality Control</a>
    </h2>
    <p style="text-align: justify;">
        The primary aim of the quality control process is to classify sondes based on their level of success in making measurements. We discard the sondes with zero or negligible usable data, and keep the remaining sondes as part of the Level-2 product.     </p>
    <p style="text-align: justify;">
        Most sondes have complete sounding profiles, while some have partially faulty profiles, thus allowing us to include some parts of the sounding that are still useful. An example is if the PTU sensor of a sonde failed, we will still include the sounding in Level-2 for the data that the GPS sensor provides.
    </p>
</div>

<div id="sonde_classification">
    <h3>
        <a href="#TOC">Sonde classification and QC terminology</a>
    </h3>
<p>
'flag' : record of result of local test  
'good' : passes local test  
'bad'  : fails local test  
'ugly' : fails preliminary local test, but some data may be salvaged  

'FLAG' : record of result of group of tests or final flag for sonde  
'GOOD' : sonde is good, passed all tests and can be used straightaway  
'BAD'  : sonde is bad, should not be considered for any further data processing or analysis  
'UGLY' : sonde has some data that may be salvaged later, but cannot be processed straightaway. Needs more QC.
</p>

<div id="quality_checks">
    <h3>
       <a href="#TOC"> Quality Checks</a>
    </h3>
</div>

There are 3 individual processes used to classify sondes, which are combined later to make a final grouping:

1. <div id="ind_flag">
    <h4>
       <a href="#TOC">ind_FLAG:</a>
    </h4>
</div> 

- Ratio of individual parameter's count to the total count of records (time is used as proxy for measurement record. Multiple parameters can have values at a single record.)
    - There is one flag each for the following parameters:  
    t   : *temperature (tdry)*  
    p   : *pressure (pres)*  
    rh  : *relative humidity (rh)*  
    z   : *GPS altitude (gpsalt)*  
    u   : *u wind (u_wind)*  
    v   : *v wind (v_wind)*  

- For each parameter, if this proportional ratio is higher than the set threshold (based on overall distribution, explained later), then the sonde is flagged as 'good'. If the ratio is lower than set threshold, but non-zero, then the sonde is flagged 'ugly' and if the ratio is zero, i.e. this particular parameter was not measured at all for the sonde, then, it is flagged as 'bad'. (This value might later be changed from zero to 10% of the time counts).
  
  - If all individual ind_flags are 'good', the sonde is flagged as 'GOOD' for ind_FLAG,  
  - if all individual ind_flags are 'bad', the sonde is flagged as 'BAD' for ind_FLAG,  
  - if neither of these conditions is met, the sonde is flagged as 'UGLY' for ind_FLAG.
    

2. <div id="srf_flag">
    <h4>
       <a href="#TOC">srf_FLAG:</a>
    </h4>
</div> 

- Values in the lower layer of the atmosphere (mostly near surface, but also < 4 km for some tests) as sanity checks.
  
- Functions are defined to check for obvious errors in near-surface values of parameters, such as being out of bounds (limits of realistic values), poor agreement between estimates of altitude, etc. These functions are described in greater detail where the function is defined.

- For each of these srf_flags, if the sonde passes the test, it is marked as 1, else as 0, which stand for 'good' and 'bad' respectively. 
  - If all individual srf_flags are 'good', the sonde is flagged as 'GOOD' for srf_FLAG,  
  - if all individual srf_flags are 'bad', the sonde is flagged as 'BAD' for srf_FLAG,  
  - if neither of these conditions is met, the sonde is flagged as 'UGLY' for srf_FLAG.  
      
    
3. <div id="ld_flag">
    <h4>
       <a href="#TOC">ld_FLAG:</a>
    </h4>
</div> 

- Checking whether the sonde detected a launch automatically. If the sonde fails to do this, it does not switch to high-power signal transmission, and thus, stops sending data back to the AVAPS PC, after a short range. 
    
- The primary method to check launch detection is to parse through the log files of the sonde in raw format. These files have names starting with 'A' and are followed by the date and time of launch. The file extension is the number of the channel used to initialise the sonde and receive its signal, but for all practical purposes, it is a .txt file. (Note: For sondes that did not detect a launch, the file name has time when the sonde was initialised). The log file contains an internal record termed 'Launch Obs Done?'. If this value is 1, the launch was detected, else if it is 0, launch was not detected. The same values are used to mark the ld_FLAG.
    
- An alternative method to check launch detection is through an attribute of the PQC files, called 'reference_time'. If a launch was detected, this attribute shows the launch time. If a launch was not detected automatically, the attribute shows the correct date, but the time is 'T00:00:00.000'. This can also be used to mark the ld_FLAG. However, this should only be used as a quick-fix if the log files are not available, since this test is not fool-proof. Moreover, PQC files with no automatic launch detection do not have relevant sonde information such as SondeID added to them. 

4. <div id="flag">
    <h4>
       <a href="#TOC">FLAG:</a>
    </h4>
</div> 

- After the sondes pass through these three processes, the FLAGs are used to determine the final FLAG value for the sonde.  
    - If ld_FLAG is 0, then the sonde FLAG is termed 'BAD'  
    - If ld_FLAG is 1, and if both srf_FLAG and ind_FLAG are 'GOOD', then the sonde FLAG is termed 'GOOD'.  
    - If ld_FLAG is 1, and if both srf_FLAG and ind_FLAG are 'BAD', then the sonde FLAG is termed 'BAD'.  
    - If ld_FLAG is 1, and if srf_FLAG and ind_FLAG have different values, then the sonde FLAG is termed 'UGLY'.   
      
- Although the process of classifying the sondes can be simplified by other combinations of the ind_flags and srf_flags, the current method ensures that no good sondes are omitted, and no bad sondes are admitted. The rest of the sondes, the ugly sondes, still have data that can be salvaged, and after some QC and/or flagging, can be combined with the other good sondes.  
  
- The NC file generated as a product stores results for each individual test mentioned above, group of tests and the final classification. Thus, the user can still mould the classification based on their objectives or add/remove tests to the process and customise it for themselves.