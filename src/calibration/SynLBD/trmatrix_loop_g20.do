clear
set more off
********************************************************************************
*cd "/home/spec676/"
*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
local pwd : pwd
global rootdir "`pwd'"
cd $rootdir


save temp_trmatrix_g20, emptyok replace
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
	* Shows the file name we are using
	use `file'
	
	gen year = substr(`"`file'"',7,4)
	* This grabs the year from the file to record the year we are in. 
	
	destring(year), replace
	*Convert it too a float
	
	display year
	
	* Only do this for years before 1996...
	
	if year < 1996 {
	
	keep if emp > 19 & missing(emp) == 0
	*keep if emp > 1 & missing(emp)==0
	* This is the Hurst Pugsley cutoff, take out employment.
	************************************************************************
	* This is if we want to look only at manufacturing...
	
	*gen nsic3 = real(sic3)

	*keep if nsic3 >199 & nsic3 < 400
	
	
	************************************************************************
	sum emp,detail 

	********************************************************************************
	* Now going to generate bins in the quantile, so the notation will be q1_y1 is 
	* bottom quintile in year 1 and so forth.
	
	************************************************************************
	* Figure out qurtiles...
	
	generate q1_y1=0 
	replace q1_y1=1 if emp<=r(p25)

	generate q2_y1=0
	replace q2_y1=1 if r(p25)<emp & emp<=r(p50)

	generate q3_y1=0 
	replace q3_y1=1 if r(p50)<emp & emp<=r(p75)

	generate q4_y1=0
	replace q4_y1=1 if r(p75)<emp & emp<=r(max)

	* Then keep only these values and save them...

********************************************************************************
	* Need to read in the new data set, this is a way to keep the name...
	* have to do this local stuff...
	local my_name = "synlbd"+ string((year+5))+"c.dta"
	
	*di my_name
	
	keep lbdnum q1_y1 q2_y1 q3_y1 q4_y1
	
	*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
	cd $rootdir

	save quantile_y1, replace
	
	************************************************************************
	
	cd "/rdcprojects/tr/tr00612/data/synlbd/2.0.2/"
	
	use `my_name', clear 
	* my_name has to be like this "local" variable, took half a day to figure
	* this out.
	
	keep if emp > 19 & missing(emp)==0
	*keep if emp > 1 & missing(emp)==0
	
	************************************************************************
	* This is if we want to look only at manufacturing...
	
	*gen nsic3 = real(sic3)

	*keep if nsic3 >199 & nsic3 < 400
	
	
	************************************************************************

	*sum emp,detail

	merge 1:1 lbdnum using $rootdir/quantile_y1.dta
	keep if _merge==3
	
	*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
	cd $rootdir

	do quartile_code.do
	* In quartile code, it works through all the different transitions...

	gen year = substr(`"`file'"',7,4)
	
	keep p11 p12 p13 p14 ///
	p21 p22 p23 p24 ///
	p31 p32 p33 p34 ///
	p41 p42 p43 p44 year
	
	append using temp_trmatrix_g20.dta
	* Append it to the .dta file
	
	save temp_trmatrix_g20.dta, replace
	* Save the .dta file
	clear
	}
	
	clear
	
	cd "/rdcprojects/tr/tr00612/data/synlbd/2.0.2/"
	* Go back to the directory where we started
	* and do this over again

}

********************************************************************************
*cd "/rdcprojects/tr/tr00612/programs/users/spec676/PTW_AER/"
cd $rootdir

use temp_trmatrix_g20.dta

********************************************************************************
* Report stuff...
sum p31

sum p32

sum p33

sum p34

sum p41

sum p42

sum p43

sum p44

