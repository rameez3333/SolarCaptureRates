from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplate.sh').readlines()


masses = [20., 35., 50., 100., 250., 500., 1000., 3000., 5000., 10000.]


jname = 'CapRate'

for mass in masses:
    jobname = jname+str(mass)+'job'
    fout = open(jobname+'.slurm', "w")
    #jobline = ' python SolarCapRateSuite.py -m ' +str(mass) + ' -s ' +str(1.0)
    jobline = ' python PICOEventRateSuite.py -m ' +str(mass)
    for line in flines:
        fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
    fout.close()
    os.system('chmod +x ' + jobname+'.slurm')
    os.system('sbatch -p icecube '+ jobname+'.slurm')
    raw_input('test')
    

