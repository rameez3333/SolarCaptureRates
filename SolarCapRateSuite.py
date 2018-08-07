import numpy as np
import scipy.integrate as integrate
import sunpy.sun.models as sunmodel
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import InterpolatedUnivariateSpline
from optparse import OptionParser
import dmdd


usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-m", "--wimpmass", action="store", type="float", default=100., dest="WMASS", help="WIMP Mass")
parser.add_option("-l", "--velocitylow", action="store", type="float", default=5., dest="VLOW", help="Lowest Velocity Stream")
parser.add_option("-u", "--velocityhigh", action="store", type="float", default=2000., dest="VHIGH", help="Highest Velocity Stream")
parser.add_option("-s", "--velocitystep", action="store", type="float", default=30., dest="VSTEP", help="Step Size Velocity Stream")
parser.add_option("-f", "--formfactor", action = "store", default=0, type='int', dest="FF", help = "0: Use the traditional exponential form factor in energy 1: Use dmdd: https://github.com/veragluscevic/dmdd")
(options, args) = parser.parse_args()

wimpmass = options.WMASS
vlow = options.VLOW
vhigh = options.VHIGH
vstep = options.VSTEP

mpl.rcParams['text.usetex'] = True 

Avogadro = 6.0221409e+23
C = 299792.458

sundarr = np.genfromtxt('bp2000stdmodel.dat-org', skip_header=23, skip_footer=1).transpose()

solarradius = sunmodel.interior.units['radius'].value #solar radius in metres
solarmass = sunmodel.interior.units['mass'].value #solar mass in kg
Rarr = np.asarray(sundarr[1])*solarradius
Marr = np.asarray(sundarr[0])*solarmass
Darr = np.asarray(sundarr[3])
nucmfrac={}
nucmfrac[1] = sundarr[6]
nucmfrac[4] = sundarr[7]
nucmfrac[3] = sundarr[8]
nucmfrac[12] = sundarr[9]
nucmfrac[14] = sundarr[10]
nucmfrac[16] = sundarr[11]

MassFrac={}
for key in nucmfrac.keys():
    MassFrac[key] = InterpolatedUnivariateSpline(Rarr, nucmfrac[key], ext=1)


SolarMass = InterpolatedUnivariateSpline(Rarr, Marr, ext=1)
SolarDensity = InterpolatedUnivariateSpline(Rarr, Darr, ext=1) #in g/cm^3
UniGravConstant = 6.67e-11 #m^3kg^-1s^-2

class VDF(object):
    def __init__(self, vdftype='SMH', vesc=650., k=2.5, vpeak = 230., vsun=233., norm=1.):
        self.vdftype=vdftype
        self.norm = norm
        self.vesc = vesc
        self.k = k
        self.vpeak=vpeak
        self.vsun = vsun
        if self.vdftype=='SMH':
            def VelGalRestFrame(v, vesc=vesc, k=k, vpeak=vpeak):
                v = np.power((np.exp((vesc**2. - v**2.)/(k*vpeak**2.)) - 1.),k) * 0.5 * (np.sign(vesc - v) + 1.)
                if not np.isscalar(v):
                    v[np.isnan(v)]=0.
                elif np.isnan(v):
                    v=0.
                return v
            self.VelGalRestFrame = VelGalRestFrame
            def V2VelGalRestFrame(v, vesc=vesc, k=k, vpeak=vpeak):
                return np.power(v, 2.)*self.VelGalRestFrame(v, vesc=vesc, k=k, vpeak=vpeak)
            self.V2VelGalRestFrame=V2VelGalRestFrame
            self.integ =  integrate.quad(self.V2VelGalRestFrame, 0., vesc+200., (vesc, k, vpeak))[0]
            def NormalizedV2VelGalRestFrame(v, vesc=vesc, k=k, vpeak=vpeak):
                return self.V2VelGalRestFrame(v, vesc=vesc, k=k, vpeak=vpeak)/self.integ
            self.NormalizedV2VelGalRestFrame = NormalizedV2VelGalRestFrame
            if self.vdftype=='SMH':
                def VelSolarRestFrame(v, vsun = vsun, VF = self.V2VelGalRestFrame):
                    def ToIntegrate(c, v=300, vsun = self.vsun):
                        return VF(np.sqrt(v**2. +vsun**2. + 2.*v*vsun*c))
                    def Integrator(v, vsun=self.vsun):
                        return integrate.quad(ToIntegrate, -1., +1., (v, vsun))[0]*0.5
                    if not np.isscalar(v):
                        return np.asarray(map(Integrator, v))
                    else:
                        return Integrator(v)    
                self.VelSolarRestFrame = VelSolarRestFrame
        
        elif self.vdftype=='Dirac':
            print 'Constructing Solar Rest frame VDF directly'
            def VelSolarRestFrame(v, vpeak=self.vpeak):
                def Vel(v, vpeak=self.vpeak):
                    lmin = np.sqrt(v**2. + vsun**2. - 2.*v*vsun)
                    lmax = np.sqrt(v**2. + vsun**2. + 2.*v*vsun)
                    if (lmin < vpeak) and (lmax > vpeak):
                        return vpeak/v*vsun
                    else:
                        return 0.
                if not np.isscalar(v):
                    return np.asarray(map(Vel, v))
                else:
                    return Vel(v)
            self.VelSolarRestFrame = VelSolarRestFrame
        def V2VelSolarRestFrame(v):
            return v*v*VelSolarRestFrame(v)
        self.V2VelSolarRestFrame=V2VelSolarRestFrame
        self.integv2solarrestframe = integrate.quad(V2VelSolarRestFrame, 0, self.vpeak+2.*self.vsun)[0]
        def NormalizedV2VelSolarRestFrame(v):
            return V2VelSolarRestFrame(v)/self.integv2solarrestframe
        self.NormalizedV2VelSolarRestFrame = NormalizedV2VelSolarRestFrame
        def NormalizedVSolarRestFrame(v):
            return NormalizedV2VelSolarRestFrame(v)/v
        self.NormalizedVSolarRestFrame = NormalizedVSolarRestFrame
        
vels =  np.linspace(0.2, 1600, 7999)


def CapProb(v, vsesc, mDM, mi=0.938):
    E = 0.5*mDM*(v*v + vsesc*vsesc)/(C*C)
    delmax = 4.*mi*mDM/(mDM+mi)**2.
    delmin = v**2./(v**2. + vsesc**2.)
    
    #print E, delmax, delmin
    
    if delmax < delmin:
        return 0.
    
    #print E, delmin, delmax
    #Qmin = 0.5*mDM*v*v/(C*C)
    #betaplus = 4.*mDM*mi/(mDM+mi)**2.
    #Qmax = 0.5*betaplus*mDM*(v*v + vsesc*vsesc)/(C*C)
    r = 0.754*np.power(mi, 1./3.)*1.e3/197. # PPPC4DMnu paper paragraph before Eq 14, natural units GeV^-1
    E0 = 1.5/(mi*r**2.)
    def integrand(x):   #Form factor
        return np.exp(-1.*x/E0)
    integ = integrate.quad(integrand, E*delmin, E*delmax, epsabs = 1e-16, epsrel = 1.e-16, limit=5000, maxp1=5000, limlst=5000)[0]
    #integ = np.max([0, (delmax-delmin)/delmax])
    #return integ
    return 0.5*(np.sign(integ/(E*delmax))+1.)*integ/(E*delmax)
    
    

def SolarEscapeVelocity(R):
    if R < 4521000:
        return 40.
    elif R < 660064912.32:
        return np.sqrt(2.*UniGravConstant*SolarMass(R)/R)/1000.
    else:
        return np.sqrt(2.*UniGravConstant*solarmass/R)/1000.
    #Escape velocity at R(in metres) in km/s

def NumDensity(rad, i=1): # return in Num/m^3
    return SolarDensity(rad)*MassFrac[i](rad)/float(i)*Avogadro*1.e6
    

    
def VelIntegrand(v, vsesc, vdf, mdm, mi=0.938):
    return 4.*np.pi*v*v*vdf.NormalizedVSolarRestFrame(v)*(v*v + vsesc*vsesc)*CapProb(v, vsesc, mdm, mi)

def Radintegrand(rad, vdf, mdm, mi=0.938, i=1):
    if vdf.vdftype=='SMH':
        velinteg = integrate.quad(VelIntegrand, 0, 1000., (SolarEscapeVelocity(rad), vdf, mdm, mi))[0]
    elif vdf.vdftype=='Dirac':
        vsesc = SolarEscapeVelocity(rad)
        v = vdf.vpeak
        velinteg = 4.*np.pi*v*(v*v + vsesc*vsesc)*CapProb(v, vsesc, mdm, mi)
    return 4.*np.pi*rad*rad*NumDensity(rad, i)*velinteg


smh = VDF()
norm = smh.integv2solarrestframe
velocities = np.arange(vlow, vhigh, vstep)
streams={}

for vel in velocities:
    streams[vel] = VDF(vdftype='Dirac', vpeak=vel)

    
#d1 =VDF(vdftype='Dirac', vpeak=20.)
#d2 =VDF(vdftype='Dirac', vpeak=100.)
#d3 =VDF(vdftype='Dirac', vpeak=200.)
#d4 =VDF(vdftype='Dirac', vpeak=700.)

#d5 =VDF()


smhcaprate = 0.3/wimpmass*1.e-36*(integrate.quad(Radintegrand, 0.1, solarradius, (smh,wimpmass, 0.938, 1))[0] + integrate.quad(Radintegrand, 0.1, solarradius, (smh,wimpmass, 0.938*14., 14))[0])

print 'SMH Capture Rate', smhcaprate

fout = open(str(wimpmass)+'_Summary.txt', "w")
fout.write('SMH Capture Rate '+str(smhcaprate)+'\n')

print 'Stream Capture Rates'
fout.write('Stream Peak Velocity|Stream Capture Rates|Rescaling Factor\n')
streamcaprates={}
for vel in velocities:
    streamcaprates[vel] = 0.3/wimpmass*1.e-36*(integrate.quad(Radintegrand, 0.1, solarradius, (streams[vel],wimpmass))[0] + integrate.quad(Radintegrand, 0.1, solarradius, (streams[vel], wimpmass, 0.938*14., 14))[0])
    print vel, streamcaprates[vel]
    fout.write(str(vel/C)+'|'+str(streamcaprates[vel])+'|'+str(smhcaprate/streamcaprates[vel])+'\n')

fout.close()



#dist1 = np.asarray(d1.NormalizedVSolarRestFrame(vels))
#dist2 = np.asarray(d2.NormalizedVSolarRestFrame(vels))
#dist3 = np.asarray(d3.NormalizedVSolarRestFrame(vels))
#dist4 = np.asarray(d4.NormalizedVSolarRestFrame(vels))
#dist5 = np.asarray(d5.NormalizedVSolarRestFrame(vels))



#plt.plot(vels, dist1, color='red', linewidth=3)
#plt.plot(vels, dist2, color='green', linewidth=3)
#plt.plot(vels, dist3, color='blue', linewidth=3)
#plt.plot(vels, dist4, color='brown', linewidth=3)
#plt.plot(vels, dist5, color='black', linewidth=3)
#plt.xlabel('v(km/s)')
#plt.ylabel(r'$f_{\bigodot}(v)/v$')
##plt.yscale('log')
#plt.show()
