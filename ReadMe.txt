“main.c”:
line 201 -> while (i_step <= TSTEPS):
before building the code, go to constants.c to adjust value of TSTEPS or TOTAL_TIME. Alternative statement should be “adv_time<TOTAL_TIME”,


“constants.h”
line 29 -> #define NEUTRON


“rates.c”

look at function:
void compute_sinks (double *s): to see the sink strength
double compute_diff_coeff (const int attr[LEVELS]): to change diffusivity of clusters
double *compute_bind_term(const int attr[LEVELS]): to change binding energy


“RateMatrix.c”
inline double computeDamages(struct RateMatrix* rateMatrix): to change damage rate, for example damage[0] is the neutron insertion rate.

IF Jaime want you to change any parameters in function compute_rates_objects in “rates.c”, just change them in  function computeDamage in “RateMatrix.c”. 


How to compile:
put in the right cpdf file:
>>make clean
>>make
>>job.q
……