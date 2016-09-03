rem Calculate probability of pairwise correlation without genotype data, and with extensive log level
findr 12 0 0 pij_rank_a data/geuvadis/dt.dat data/geuvadis/dt2.dat 10 3000 360 data/geuvadis/dp.dat 0

rem Calculate probability of 5 tests on pairwise regulation with genotype data, and with default log level
findr 0 0 0 pijs_gassist_a data/geuvadis/dg.dat data/geuvadis/dt.dat data/geuvadis/dt2.dat 10 3000 360 data/geuvadis/dp1.dat data/geuvadis/dp2.dat data/geuvadis/dp3.dat data/geuvadis/dp4.dat data/geuvadis/dp5.dat 2 0

rem Calculate recommended combination of tests on pairwise regulation with genotype data, and with default log level
findr 0 0 0 pij_gassist_a data/geuvadis/dg.dat data/geuvadis/dt.dat data/geuvadis/dt2.dat 10 3000 360 data/geuvadis/dp.dat 2 0
