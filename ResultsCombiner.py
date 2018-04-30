import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from scipy.optimize import fsolve, minimize

nflist = sorted(glob('*_Summary.txt'))
flist = sorted(glob('*_DDSummary.txt'))

print flist
colour={}
colour[11]='red'
colour[8]='green'
colour[5]='blue'
clist = ['red', 'green', 'blue', 'brown', 'black', 'cyan', 'magenta', 'violet', 'yellow', 'salmon']
dDD ={}
dDDx={}
dSC ={}
dSCx={}

for fname, fname2 in zip(flist, nflist):
    arr = np.genfromtxt(fname, skip_header=2, delimiter="|").transpose()
    arr2 = np.genfromtxt(fname2, skip_header=2, delimiter="|").transpose()
    mass = float(fname.split("_")[0])
    dDD[mass]=arr[2]
    dDDx[mass]=arr[0]
    dSC[mass]=arr2[2]
    dSCx[mass]=arr2[0]
    cl=clist.pop()
    #print arr, arr2
    #plt.plot(arr2[0], arr2[2], label=mass+' GeV', lw=2, color=cl)
    #plt.plot(arr[0], arr[2], lw=2, color=cl, ls='--')
    
    #plt.title('DM Mass = ' +mass, fontsize=20.)
    
Slines = open('LimitsToRescale/SolarWIMP_GvaSthcombnined_3yr.txt').readlines()[1:]
DDdat = np.genfromtxt('LimitsToRescale/PICO1702.txt', delimiter=',').transpose()

picollog = InterpolatedUnivariateSpline(np.log10(DDdat[0]), np.log10(DDdat[1]))
picol = InterpolatedUnivariateSpline(DDdat[0], DDdat[1])


channels = [11, 5, 8]

legendtxt={}
legendtxt[5] = r'$b \bar{b}$'
legendtxt[8] = r'$W^{+} W^{-}$'
legendtxt[11] = r'$\tau^{+} \tau^{-}$'
Smasses={}
Smassesall={}
SSDlimits={}
SSDlimitsall={}
Ssyst={}
PL={}



    
    



DDl={}
Sl={}

for ch in channels:
    Smasses[ch]=[]
    SSDlimits[ch]=[]
    Ssyst[ch]=[]
    DDl[ch]={}
    Sl[ch]={}

for line in Slines:
    Smasses[int(line.split()[1])].append(float(line.split()[0]))
    SSDlimits[int(line.split()[1])].append(float(line.split()[-3]))
    Ssyst[int(line.split()[1])].append(float(line.split()[-1]))
    
for ch in channels:
    Smassesall[ch] = np.asarray(Smasses[ch])
    Smasses[ch]=np.asarray(Smasses[ch])[0:-3]
    SSDlimitsall[ch]=np.asarray(SSDlimits[ch])*1.e-36
    SSDlimits[ch]=np.asarray(SSDlimits[ch])[0:-3]*1.e-36
    Ssyst[ch]=np.asarray(Ssyst[ch])[0:-3]
    PL[ch] = picol(Smasses[ch])
    for m, sl, dl in zip(Smasses[ch], SSDlimits[ch], PL[ch]):
        DDl[ch][m] = dl
        Sl[ch][m] = sl

def findIntersection(x0, fun1,fun2):
    return minimize(lambda x : np.absolute(fun1(x) - fun2(x)), x0, method='SLSQP')

def ResultsCombiner(m=100., ch=5, plots=False):
    DDL = DDl[ch][m]
    SL = Sl[ch][m]
    
    ddstreams, ddlims = dDDx[m], dDD[m]*DDL
    sstreams, slims = dSCx[m], dSC[m]*SL
    
    ddstreams, ddlims = ddstreams[np.isfinite(ddlims)], ddlims[np.isfinite(ddlims)]
    sstreams, slims = sstreams[np.isfinite(slims)], slims[np.isfinite(slims)]
    
    #print 'DDstreams', ddstreams, ddlims
    
    
    Df = InterpolatedUnivariateSpline(ddstreams, ddlims, ext=0, k=1)
    Sf = InterpolatedUnivariateSpline(sstreams, slims, ext=0, k=1)

    #pint = findIntersection(1.e-3, Df, Sf)
    
    #print 'This intersects at ', pint
    
    l, u = np.max([ddstreams[0], sstreams[0]]), np.min([ddstreams[-1], sstreams[-1]])
    
    intfound = True
    
    #try:

    print 'Edges', l, u
    xs = np.power(10., np.linspace(np.log10(l),np.log10(u), 50000))
    
    print np.min(xs), np.max(xs)
    print np.min(ddstreams), np.max(ddstreams)
    print np.min(sstreams), np.max(sstreams)
    
    
    diffs = np.absolute(Df(xs) - Sf(xs))
#xs = xs[diffs!=np.nan]
#diffs = diffs[diffs!=np.nan]
#print 'Try Now'

#for x, diff in zip(xs, diffs):
    #print x, Df(x), Sf(x), diff

    pint =  xs[np.argmin(diffs)]
    
    if ch==5:
        if m != 100:
            pint = np.min(xs)
    
    print 'Intersection at', pint
    #except:
        #intfound=False
    
    
    #plt.figure(0)
    #plt.plot(xs, diffs)
    #plt.plot([1.e-4])
    #plt.xscale('log')
    #plt.yscale('log')
    ##plt.show()
    
    
    
    
    plt.figure(1)
    plt.plot(ddstreams, ddlims, color='red', label='PICO', lw=3)
    plt.plot(sstreams, slims, color='green', label='IceCube', lw=3)
    if intfound:
        plt.plot([1.e-4, 1.e-1], [Sf(pint), Sf(pint)], color='black', label='Halo VDF Independent', lw=3, ls='--')
    plt.xscale('log')
    plt.xlabel('Stream Velocity[c]')
    plt.yscale('log')
    plt.ylabel('CS Limit')
    plt.legend(loc='best', fontsize=15.)
    plt.ylabel(r"$\sigma^{\mathrm{SD}}_{\chi\mathrm{-}p}$ [cm$^{2}$]", fontsize=20)
    plt.title('WIMP Mass = '+str(m) + ' GeV, ' + legendtxt[ch])
    plt.savefig('CombinedFigures/' + str(m) +'_'+ str(ch)+'.png')
    plt.show()
    return pint, Sf(pint)
    

CSsd={}
VP={}
for ch in channels:
    CSsd[ch]=[]
    VP[ch]=[]
    for m in Smasses[ch]:
        vp, cssd = ResultsCombiner(m, ch)
        CSsd[ch].append(cssd)
        VP[ch].append(vp)
        
plt.figure(0)
for ch in channels:
    plt.plot(Smasses[ch], CSsd[ch], label=legendtxt[ch], color=colour[ch], lw=4)
plt.xlabel(r"$m_{\chi}$ [GeV]", fontsize=15)
plt.ylabel(r"$\sigma^{\mathrm{SD}}_{HI}$ [cm$^{2}$]", fontsize=15)
plt.xscale('log')
plt.yscale('log')
plt.axis([5., 3.e3, 5.e-42, 1.e-36])
plt.legend(loc='best', fontsize=15.)
plt.savefig('OnlyHI.png')
#plt.show()

plt.figure(1)
lghand1={}
for ch in channels:
    plt.plot(Smasses[ch], CSsd[ch], label=legendtxt[ch], color=colour[ch], lw=4)
    plt.plot(Smassesall[ch], SSDlimitsall[ch], color=colour[ch], lw=2, ls='--')
plt.plot(DDdat[0], DDdat[1], lw=2, color='black', ls='--')
plt.xlabel(r"$m_{\chi}$ [GeV]", fontsize=15)
plt.ylabel(r"$\sigma^{\mathrm{SD}}_{HI}$ [cm$^{2}$]", fontsize=15)
plt.xscale('log')
plt.yscale('log')
plt.axis([5., 3.e3, 5.e-42, 1.e-36])
plt.legend(loc='best', fontsize=15.)
plt.savefig('OnlyHI.png')
plt.show()
    
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('Stream Peak Velocity', fontsize=15.)
#plt.ylabel('Rescaling Factor', fontsize=15.)
#plt.legend(loc='upper left', ncol=3, fontsize=10.)
#plt.show()
#plt.savefig('BothSummary.png')
#del plt
