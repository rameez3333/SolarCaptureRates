import numpy as np
from glob import glob


flist = sorted(glob('*_Summary.txt'))
print flist

for fname in flist:
    arr = np.genfromtxt(fname, skip_header=2, delimiter="|").transpose()
    mass = fname.split("_")[0]
    import matplotlib.pyplot as plt
    print arr
    plt.plot(arr[0], arr[2])
    plt.title('DM Mass = ' +mass, fontsize=20.)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Stream Peak Velocity', fontsize=15.)
    plt.ylabel('Rescaling Factor', fontsize=15.)
    plt.show()
    plt.savefig(mass+'summary.png')
