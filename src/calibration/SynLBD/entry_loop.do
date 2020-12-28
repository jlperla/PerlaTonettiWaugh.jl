clear
********************************************************************************
set more off

*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
local pwd : pwd
global rootdir "`pwd'"
cd $rootdir

save temp_entry, emptyok replace
* creates an empty .dta file which we will fill in...

cd "/rdcprojects/tr/tr00612/data/synlbd/2.0.2/"
* where the files are located

clear

local files: dir "`c(pwd)'" files"*.dta"
*grab the lsit of the files in the directory

* Then the below commands will loop over each file, perform the computation, and 
* then append the results to the .dta file created above...
* This all works...

foreach file in `files'{
	di "`file'"
	use `file'
	
	*keep emp lastyear 
	* Take only employment and the last year
	
	keep if emp>19 & missing(emp)==0
	
	gen year = substr(`"`file'"',7,4)
	* This grabs the year from the file to record the year we are in. 
	
	destring(year), replace
	*Convert it too a float
	
	display year
	
	* Then the stuff below, compute all employment, then compute all the employment
	* based upon exit, and a size threshold.
	
	egen all_emp = sum(emp)
	
	egen entry = sum(emp) if firstyear == (year) 
	
	* Manipulate this stuff
	
	collapse all_emp entry

	keep all_emp entry

	* Generate the ratios...
	
	gen entry_rate = entry/all_emp
		
	display entry_rate 
	
	
	* Create the year again (some reason it gets lost)
	
	gen year = substr(`"`file'"',7,4)
	
	keep entry_rate year
	* This is the stuff we want
	
	*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
	cd $rootdir

	append using temp_entry.dta
	* Append it to the .dta file
	
	save temp_entry.dta, replace
	* Save the .dta file
	
	cd "/rdcprojects/tr/tr00612/data/synlbd/2.0.2/"
	* Go back to the directory where we started
	* and do this over again
	
	clear
	}
********************************************************************************
* Report the output...	
*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
cd $rootdir

use temp_entry.dta

sum entry_rate
