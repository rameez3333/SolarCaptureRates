import numpy as np
from glob import glob
import matplotlib.pyplot as plt

nflist = sorted(glob('*_Summary.txt'))
flist = sorted(glob('*_DDSummary.txt'))

print flist

clist = ['red', 'green', 'blue', 'brown', 'black', 'cyan', 'magenta', 'violet', 'yellow', 'salmon']

for fname, fname2 in zip(flist, nflist):
    arr = np.genfromtxt(fname, skip_header=2, delimiter="|").transpose()
    arr2 = np.genfromtxt(fname2, skip_header=2, delimiter="|").transpose()
    mass = fname.split("_")[0]
    cl=clist.pop()
    print arr
    plt.plot(arr2[0], arr2[2], label=mass+' GeV', lw=2, color=cl)
    plt.plot(arr[0], arr[2], lw=2, color=cl, ls='--')
    #plt.title('DM Mass = ' +mass, fontsize=20.)
    
    
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Stream Peak Velocity', fontsize=15.)
plt.ylabel('Rescaling Factor', fontsize=15.)
plt.legend(loc='upper left', ncol=3, fontsize=10.)
plt.savefig('BothSummary.png')
plt.show()
del plt
