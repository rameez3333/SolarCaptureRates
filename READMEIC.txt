This code consists mainly of two scripts.

SolarCapRateSuite.py, which calculates the capture rate in the Sun by evaluating Eq. 3 of 1506.03386 for different velocity distribution functions, and 

PICOEventRateSuite.py, which calculates the event rates in the detector by evaluating Eq. 2 of  1506.03386 for different velocity distribution functions. This script relies on the file PICOEfficiencies.dat.

ResultsCombiner.py will combine the outputs from the other previous two scripts to produce the final figures.

batchsubmitter.py contains an example of how the first two scripts can be used on a slurm machine to compute all the rescaling factors.
