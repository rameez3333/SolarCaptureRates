from glob import glob
import os, sys
import numpy as np
import time
nside = str(32)

flines = open('submittemplate.sh').readlines()

#mode = sys.argv[1]

masses = [20., 35., 50., 100., 250., 500., 1000., 3000., 5000., 10000.]


jname = 'PICOEvRate'

if int(sys.argv[1])==1:
    for mass in masses:
      for O in [1, 2, 0]:
        jobname = jname+str(mass)+'job' + str(O)
        fout = open(jobname+'.slurm', "w")
    #jobline = ' python SolarCapRateSuite.py -m ' +str(mass) + ' -s ' +str(1.0)
        jobline = ' python PICOEventRateSuite.py -m ' +str(mass) + ' -o '+str(O)
        for line in flines:
            fout.write(line.replace('__NAME__', jobname).replace('__JOBLINE__', jobline))
        fout.close()
        os.system('chmod +x ' + jobname+'.slurm')
        os.system('sbatch -p icecube '+ jobname+'.slurm')
        #raw_input('test')
    
elif int(sys.argv[1])==2:
    for mass in masses:
        jobline = ' python PICOEventRateSuite.py -m ' +str(mass)
        os.system(jobline)
