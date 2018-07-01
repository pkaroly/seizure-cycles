# Seizure Cycles
Code for analysis of SeizureTracker data

All analsyis was undertaken in MATLAB R2017a (Mathworks). Note that all code requires access to the Circular Statistics Toolbox that can be downloaded online:
https://au.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-

## Rerefences

P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009 
http://www.jstatsoft.org/v31/i10

## File List

* **ST_cycles.m**: Main code

* **parse2matlab.m**: Convert excel data to matlab

* **ST_by_syndrome.m**: analyse the epilepsy types

* **ST_by_type.m**: analyse the seizure types

* **ST_simulation.m**: perform simulations for statistical significance

* **watson.m**: Compute watson's U2 test for circular data
