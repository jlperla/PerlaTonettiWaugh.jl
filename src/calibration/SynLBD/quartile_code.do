*sum emp,detail

* Compute the same kind of statistics as above. 
set more off

sum emp,detail

* Compute the same kind of statistics as above. 

generate q1_y2=0 
replace q1_y2=1 if emp<=r(p25)

generate q2_y2=0
replace q2_y2=1 if r(p25)<emp & emp<=r(p50)

generate q3_y2=0 
replace q3_y2=1 if r(p50)<emp & emp<=r(p75)

generate q4_y2=0
replace q4_y2=1 if r(p75)<emp & emp<=r(max)

********************************************************************************
* Now create some new variables, q11 which is the variable if a q1 to q1 transition
* occuerd.

generate q11=0 
replace q11=1 if q1_y1==1 & q1_y2==1

generate q12=0 
replace q12=1 if q1_y1==1 & q2_y2==1

generate q13=0 
replace q13=1 if q1_y1==1 & q3_y2==1

generate q14=0 
replace q14=1 if q1_y1==1 & q4_y2==1

generate q21=0 
replace q21=1 if q2_y1==1 & q1_y2==1

generate q22=0 
replace q22=1 if q2_y1==1 & q2_y2==1

generate q23=0 
replace q23=1 if q2_y1==1 & q3_y2==1

generate q24=0 
replace q24=1 if q2_y1==1 & q4_y2==1

generate q31=0 
replace q31=1 if q3_y1==1 & q1_y2==1

generate q32=0 
replace q32=1 if q3_y1==1 & q2_y2==1

generate q33=0 
replace q33=1 if q3_y1==1 & q3_y2==1

generate q34=0 
replace q34=1 if q3_y1==1 & q4_y2==1

generate q41=0 
replace q41=1 if q4_y1==1 & q1_y2==1

generate q42=0 
replace q42=1 if q4_y1==1 & q2_y2==1

generate q43=0 
replace q43=1 if q4_y1==1 & q3_y2==1

generate q44=0 
replace q44=1 if q4_y1==1 & q4_y2==1

egen sum_q1=sum(q1_y1)
egen sum_q2=sum(q2_y1)
egen sum_q3=sum(q3_y1)
egen sum_q4=sum(q4_y1)

egen sum_q11=sum(q11)
egen sum_q12=sum(q12)
egen sum_q13=sum(q13)
egen sum_q14=sum(q14)
egen sum_q21=sum(q21)
egen sum_q22=sum(q22)
egen sum_q23=sum(q23)
egen sum_q24=sum(q24)
egen sum_q31=sum(q31)
egen sum_q32=sum(q32)
egen sum_q33=sum(q33)
egen sum_q34=sum(q34)
egen sum_q41=sum(q41)
egen sum_q42=sum(q42)
egen sum_q43=sum(q43)
egen sum_q44=sum(q44)

collapse sum_q1 sum_q2 sum_q3 sum_q4 ///
sum_q11 sum_q12 sum_q13 sum_q14 ///
sum_q21 sum_q22 sum_q23 sum_q24 ///
sum_q31 sum_q32 sum_q33 sum_q34 ///
sum_q41 sum_q42 sum_q43 sum_q44

gen p11=sum_q11/sum_q1
gen p12=sum_q12/sum_q1
gen p13=sum_q13/sum_q1
gen p14=sum_q14/sum_q1
gen p21=sum_q21/sum_q2
gen p22=sum_q22/sum_q2
gen p23=sum_q23/sum_q2
gen p24=sum_q24/sum_q2
gen p31=sum_q31/sum_q3
gen p32=sum_q32/sum_q3
gen p33=sum_q33/sum_q3
gen p34=sum_q34/sum_q3
gen p41=sum_q41/sum_q4
gen p42=sum_q42/sum_q4
gen p43=sum_q43/sum_q4
gen p44=sum_q44/sum_q4

