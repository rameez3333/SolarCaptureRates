import numpy as np
import scipy.integrate as integrate
import sunpy.sun.models as sunmodel
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d
from optparse import OptionParser



usage = 'usage: %prog [options]'
parser = OptionParser(usage)

parser.add_option("-m", "--wimpmass", action="store", type="float", default=100., dest="WMASS", help="WIMP Mass")
parser.add_option("-l", "--velocitylow", action="store", type="float", default=5., dest="VLOW", help="Lowest Velocity Stream")
parser.add_option("-u", "--velocityhigh", action="store", type="float", default=30000., dest="VHIGH", help="Highest Velocity Stream")
parser.add_option("-s", "--velocitysteps", action="store", type="float", default=3000., dest="VSTEP", help="number of velocity steps")

(options, args) = parser.parse_args()

wimpmass = options.WMASS
vlow = options.VLOW
vhigh = options.VHIGH
vstep = options.VSTEP

#Load PICO efficiency data, which corresponds to %of threshold as:
percthresh = np.asarray([0.0, 0.2, 0.5, 0.8, 1.0])
pdls = open('PICOEfficiencies.dat').readlines()
#C = Carbon, F=Flourine; b=best fit, o = optimistic, p = pessimistic (1 sigma)
C_b = np.asarray([float(k) for k in filter(None, pdls[0][pdls[0].find('[')+1:pdls[0].find(']')].split(' '))])
F_b = np.asarray([float(k) for k in filter(None, pdls[1][pdls[1].find('[')+1:pdls[1].find(']')].split(' '))])
C_p = np.asarray([float(k) for k in filter(None, pdls[2][pdls[2].find('[')+1:pdls[2].find(']')].split(' '))])
F_p = np.asarray([float(k) for k in filter(None, pdls[3][pdls[3].find('[')+1:pdls[3].find(']')].split(' '))])
C_o = np.asarray([float(k) for k in filter(None, pdls[4][pdls[4].find('[')+1:pdls[4].find(']')].split(' '))])
F_o = np.asarray([float(k) for k in filter(None, pdls[5][pdls[5].find('[')+1:pdls[5].find(']')].split(' '))])


Avogadro = 6.0221409e+23
C = 299792.458
masses = {}

masses['C'] = 12.0107
masses['F'] = 18.998

efficiencies={}
#efficiencies['C_b'] = interp1d(C_b, percthresh, fill_value=[0.0,1.0])
#efficiencies['F_b'] = interp1d(F_b, percthresh, fill_value=[0.0,1.0])
#efficiencies['C_p'] = interp1d(C_p, percthresh, fill_value=[0.0,1.0])
#efficiencies['F_p'] = interp1d(F_p, percthresh, fill_value=[0.0,1.0])
#efficiencies['C_o'] = interp1d(C_o, percthresh, fill_value=[0.0,1.0])
#efficiencies['F_o'] = interp1d(F_o, percthresh, fill_value=[0.0,1.0])

#scipy.interp1d should be faster than np.interp, but has a bug in many versions that doesnt allow it to take fill_value

efficiencies['C_b'] = lambda x: np.interp(x, C_b, percthresh, left = 0.0, right = 1.0)
efficiencies['F_b'] = lambda x: np.interp(x, F_b, percthresh, left = 0.0, right = 1.0)
efficiencies['C_p'] = lambda x: np.interp(x, C_p, percthresh, left = 0.0, right = 1.0)
efficiencies['F_p'] = lambda x: np.interp(x, F_p, percthresh, left = 0.0, right = 1.0)
efficiencies['C_o'] = lambda x: np.interp(x, C_o, percthresh, left = 0.0, right = 1.0)
efficiencies['F_o'] = lambda x: np.interp(x, F_o, percthresh, left = 0.0, right = 1.0)

abundances={}
abundances['C'] = masses['C']*3./(masses['C']*3. + masses['F']*8.)
abundances['F'] = masses['F']*8./(masses['C']*3. + masses['F']*8.)


class VDF(object):
    def __init__(self, vdftype='SMH', vesc=650., k=2.5, vpeak = 230., vsun=233., norm=1., earthvel=30.):
        self.vdftype=vdftype
        self.norm = norm
        self.vesc = vesc
        self.k = k
        self.vpeak=vpeak
        self.vsun = vsun
        self.earthvel = earthvel
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
        


def VelIntegrand(v, vdf, mdm, mi, cosangle=0):
    r = 0.754*np.power(mi, 1./3.)*1.e3/197.
    E0 = 1.5/(mi*r**2.)
    def integrand(x):  #Form factor goes in here
        return np.exp(-1.*x/E0)
    #E = 0.5*mdm*v*v/C*C
    return 4.*np.pi*vdf.NormalizedV2VelSolarRestFrame(v+cosangle*vdf.earthvel)/v#*integrand(0.5*mdm*(v*v)/(C*C))


def RecoilEnergyIntegrand(ER, vdf, mdm, nuc='C', err='b',  cosangle=0):
    r = 0.754*np.power(masses[nuc], 1./3.)*1.e3/197.
    E0 = 1.5/(masses[nuc]*r**2.)
    def integrand(x):  #Form factor goes in here
        return np.exp(-1.*x/E0)
    muM = mdm*masses[nuc]/(mdm + masses[nuc])
    vmin = np.sqrt(masses[nuc]*ER/1.e6/2.)/muM*C
    if vdf.vdftype=='SMH':
        integ = integrate.quad(VelIntegrand, vmin, 5000., (vdf, mdm, masses[nuc]))[0]
        retval = efficiencies[nuc+'_'+err](ER)*abundances[nuc]/masses[nuc]*integ*integrand(ER/1.e6)
        print nuc, ER, E0 ,retval, vmin
        return retval
    elif vdf.vdftype=='Dirac':
        v = vdf.vpeak
        if v>vmin:
            integ = 4.*np.pi*(v+ vdf.earthvel*cosangle)/v**2.
        else:
            integ = 0.
        retval = efficiencies[nuc+'_'+err](ER)*abundances[nuc]/masses[nuc]*integ*integrand(ER/1.e6)#*integrand(0.5*mdm*(v*v)/(C*C))
        print nuc, ER, E0 ,retval, vmin
        return retval




smh = VDF()
norm = smh.integv2solarrestframe
velocities = np.power(10., np.linspace(np.log10(vlow), np.log10(vhigh), vstep))

streams={}

for vel in velocities:
    streams[vel] = VDF(vdftype='Dirac', vpeak=vel)




smhevrate = 0.3/wimpmass*1.e-36*(((integrate.quad(RecoilEnergyIntegrand, 2.5, 500., (smh,wimpmass, 'F')))[0] + (integrate.quad(RecoilEnergyIntegrand, 2.5, 500., (smh,wimpmass, 'C')))[0]) )# + ((integrate.quad(RecoilEnergyIntegrand, 1000., 1000000., (smh,wimpmass, 'F')))[0] + (integrate.quad(RecoilEnergyIntegrand, 1000., 1000000., (smh,wimpmass, 'C')))[0]) +((integrate.quad(RecoilEnergyIntegrand, 1000000., np.inf, (smh,wimpmass, 'F')))[0] + (integrate.quad(RecoilEnergyIntegrand, 1000000., np.inf, (smh,wimpmass, 'C')))[0]))

print smhevrate


fout = open(str(wimpmass)+'_DDSummary.txt', "w")
fout.write('SMH Event Rate '+str(smhevrate)+'\n')

print 'Stream Event Rates'
fout.write('Stream Peak Velocity|Stream Event Rates|Rescaling Factor\n')
        
streamevrates={}
for vel in velocities:
    streamevrates[vel] = 0.3/wimpmass*1.e-36*( (integrate.quad(RecoilEnergyIntegrand, 2.5, 500., (streams[vel],wimpmass, 'F'))[0] + integrate.quad(RecoilEnergyIntegrand, 2.5, 500., (streams[vel], wimpmass, 'C'))[0])) #+ (integrate.quad(RecoilEnergyIntegrand, 1000., 1000000., (streams[vel],wimpmass, 'F'))[0] + integrate.quad(RecoilEnergyIntegrand, 1000., 1000000., (streams[vel], wimpmass, 'C'))[0]) + (integrate.quad(RecoilEnergyIntegrand, 1000000., np.inf, (streams[vel],wimpmass, 'F'))[0] + integrate.quad(RecoilEnergyIntegrand, 1000000., np.inf, (streams[vel], wimpmass, 'C'))[0]) )
    print vel, streamevrates[vel]
    if streamevrates[vel]:
        fout.write(str(vel/C)+'|'+str(streamevrates[vel])+'|'+str(smhevrate/streamevrates[vel])+'\n')
    else:
        fout.write(str(vel/C)+'|'+str(streamevrates[vel])+'|'+str(np.nan)+'\n')

fout.close()



















        
        
        
        
        
        
        
        
