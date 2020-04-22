# Level-2 Quality Control

\textbf{How the sondes are classified :}

Terminology:  
'flag' : record of result of local test  
'good' : passes local test  
'bad'  : fails local test  
'ugly' : fails preliminary local test, but some data may be salvaged  

'FLAG' : record of result of group of tests or final flag for sonde  
'GOOD' : sonde is good, passed all tests and can be used straightaway  
'BAD'  : sonde is bad, should not be considered for any further data processing or analysis  
'UGLY' : sonde has some data that may be salvaged later, but cannot be processed straightaway. Needs more QC.


There are 3 individual processes used to classify sondes, which are combined later to make a final grouping:

1. \textbf{ind_FLAG:} Ratio of individual parameter's count to the total count of records (time is used as proxy for measurement record. Multiple parameters can have values at a single record.)
    - There is one flag each for the following parameters:  
    t : temperature (tdry)  
    p : pressure (pres)  
    rh: relative humidity (rh)  
    z : GPS altitude (gpsalt)  
    u : u wind (u_wind)  
    v : v wind (v_wind)  
    
    - For each parameter, if this proportional ratio is higher than the set threshold (based on overall distribution, explained later), then the sonde is flagged as 'good'. If the ratio is lower than set threshold, but non-zero, then the sonde is flagged 'ugly' and if the ratio is zero, i.e. this particular parameter was not measured at all for the sonde, then, it is flagged as 'bad'. (This value might later be changed from zero to 10% of the time counts).
    - If all individual ind_flags are 'good', the sonde is flagged as 'GOOD' for ind_FLAG,  
      if all individual ind_flags are 'bad', the sonde is flagged as 'BAD' for ind_FLAG,  
      if neither of these conditions is met, the sonde is flagged as 'UGLY' for ind_FLAG.
    

2. \textbf{srf_FLAG:} Values in the lower layer of the atmosphere (mostly near surface, but also < 4 km for some tests) as sanity checks
    - Functions are defined to check for obvious errors in near-surface values of parameters, such as being out of bounds (limits of realistic values), poor agreement between estimates of altitude, etc.
    - These functions are described in greater detail where the function is defined.
    - For each of these flags, if the sonde passes the test, it is marked as 1, else as 0, which stand for 'good' and 'bad' respectively. 
    - If all individual srf_flags are 'good', the sonde is flagged as 'GOOD' for srf_FLAG,  
      if all individual srf_flags are 'bad', the sonde is flagged as 'BAD' for srf_FLAG,  
      if neither of these conditions is met, the sonde is flagged as 'UGLY' for srf_FLAG.  
      
    
3. \textbf{ld_FLAG:} Checking whether the sonde detected a launch automatically. If the sonde fails to do this, it does not switch to high-power signal transmission, and thus, stops sending data back to the AVAPS PC, after a short range. 
    - The primary method to check launch detection is to parse through the log files of the sonde in raw format. These files have names starting with 'A' and are followed by the date and time of launch. The file extension is the number of the channel used to initialise the sonde and receive its signal, but for all practical purposes, it is a .txt file. (Note: For sondes that did not detect a launch, the file name has time when the sonde was initialised). The log file contains an internal record termed 'Launch Obs Done?'. If this value is 1, the launch was detected, else if it is 0, launch was not detected. The same values are used to mark the ld_FLAG.
    - An alternative method to check launch detection is through an attribute of the PQC files, called 'reference_time'. If a launch was detected, this attribute shows the launch time. If a launch was not detected automatically, the attribute shows the correct date, but the time is 'T00:00:00.000'. This can also be used to mark the ld_FLAG. However, this should only be used as a quick-fix if the log files are not available, since this test is not fool-proof. Moreover, PQC files with no automatic launch detection do not have relevant sonde information such as SondeID added to them. 
    
After the sondes pass through these three processes, the FLAGs are used to determine the final FLAG value for the sonde.  
    - If ld_FLAG is 0, then the sonde FLAG is termed 'BAD'  
    - If ld_FLAG is 1, and if both srf_FLAG and ind_FLAG are 'GOOD', then the sonde FLAG is termed 'GOOD'.  
    - If ld_FLAG is 1, and if both srf_FLAG and ind_FLAG are 'BAD', then the sonde FLAG is termed 'BAD'.  
    - If ld_FLAG is 1, and if srf_FLAG and ind_FLAG have different values, then the sonde FLAG is termed 'UGLY'.   
      
Although the process of classifying the sondes is simplified by summarising among the ind_flags and srf_flags, the method ensures that no good sondes are omitted, and no bad sondes are admitted. The rest of the sondes, the ugly sondes, still have data that can be salvaged, and after some QC and/or flagging, can be combined with the other good sondes.  
  
The NC file generated as a product stores results for each individual test, group of tests and the final classification. Thus, the user can still mould the classification based on their objectives or add/remove tests to the process and customise it for themselves.