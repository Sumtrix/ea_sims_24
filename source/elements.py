import numpy as np
from functools import reduce
import matplotlib.pyplot as plt

def zweitorQuer(Z):
    return [[1, Z], [0, 1]]

def zweitorLaengs(Z):
    return [[1, Z], [0, 1]]

# -------------------------

class phsyConst():
    def __init__(self):
        self.rho = 1.205
        self.c = 343
        self.mu = 1.8e-10
        self.p0 = 1013
        self.kappa = 1.4


class Zweitor():
    def __init__(self, Z, laengs=1, freqDep=None):
        self.val = Z
        self.freqDep = freqDep
        self.laengs = laengs
            
    def getM(self, w):
        if self.laengs:
            match self.freqDep:
                case "0100":
                    res = [[1, 1j*w*self.val],[0, 1]]
                case "0200":
                    res = [[1, self.val/(1j*w)],[0, 1]]
                case _:
                    res = np.matrix([[1, self.val], [0, 1]])
        else:
            match self.freqDep:
                case "0010":
                    res = [[1, 0],[(1/self.val)*(1/(1j*w)), 1]] # double inverse due to kehrwert
                case "0020":
                    res = [[1, 0],[(1/self.val)*(1j*w), 1]] # double inverse due to kehrwert
                case _:
                    res = np.matrix([[1, 0], [1/self.val, 1]])
        return res



class electrElem(Zweitor):
    def __init__(self, Z, laengs, freqDep=None):
        self.domein = "electrical"
        Zweitor.__init__(self, Z, laengs, freqDep)
        

class mechanElem(Zweitor):
    def __init__(self, Z, laengs, freqDep=None):
        self.domein = "mechanical"
        Zweitor.__init__(self, Z, laengs, freqDep)


class acoustElem():
    def __init__(self, Z, laengs, freqDep=None):
        super().__init__()
        self.domein = "acoustical"
        Zweitor.__init__(self, Z, laengs, freqDep)

class elecMechTransf():
    def __init__(self, Bl):
        self._val = Bl
        self._mat = np.matrix([[Bl, 0], [0, 1/Bl]])

    def getM(self):
        return self._mat
    
class mechAcoustTransf():
    def __init__(self, S):
        self._val = S
        self._mat = np.matrix([[0, 1/S], [S, 0]])

    def getM(self):
        return self._mat
        

class SoundField(phsyConst):
    def __init__(self, S):
        super().__init__()
        self.S = S
        self.a = np.sqrt(S/np.pi)
    
    def getLoad(self, w):
        k = w / self.c
        self.Zs = self.rho*self.c*(0.5*(k*self.a)**2 + 1j*(8/(3*np.pi))*k*self.a)
        self.Z = self.Zs / self.S
        return self.Z


class ElecDynSpeaker(phsyConst, electrElem, mechanElem, acoustElem, elecMechTransf):
    def __init__(self, TSParams, xmax):
        phsyConst.__init__(self)

        self.QES = TSParams["Qes"] 
        self.QMS = TSParams["Qms"]
        self.FS = TSParams["Fs"]
        self.VAS = TSParams["Vas"] * 1e-3
        self.MMS = TSParams["Mms"] * 1e-4 
        self.RE = TSParams["Re"]
        self.SD = TSParams["Sd"] * 1e-4 
        self.LE = TSParams["Le"] * 1e-3
        self.CMS = self.VAS / (self.rho * self.c**2 * self.SD**2)
        self.RMS = np.sqrt(self.MMS / self.CMS) / self.QMS
        self.BL = np.sqrt(self.RE * self.MMS * (2 * np.pi * self.FS)**2 * self.QES)
        self.CMES = 1 / ((2 * np.pi * self.FS)**2 * self.MMS)
        self.LCES = self.CMES / (self.BL**2)
        self.RES = self.BL**2 / self.RE
        self.QTS = (self.QMS * self.QES) / (self.QMS + self.QES)
        self.XMAX = xmax
        
        w = np.linspace(20, 20000, 19981)
        self.build_system()


    def build_system(self):
        self.R = electrElem(self.RE, 1)
        self.L = electrElem(1/self.LE, 1, "0200")
        self.BL = elecMechTransf(self.BL)
        self.m = mechanElem(1/self.MMS, 0, "0020")
        self.s = mechanElem(1/self.CMS, 0, "0010")
        self.r = mechanElem(1/self.RMS, 0)
        self.S = mechAcoustTransf(self.SD)
        #self.sf = SoundField(self.SD)

    def calculateSystem(self, w):
        matrices = [self.R.getM(w), self.L.getM(w), self.BL.getM(), self.m.getM(w), self.s.getM(w), self.r.getM(w), self.S.getM()]
        sys = reduce(np.matmul, matrices)
        return sys
    
    def calculateLoad(self, w):
        return self.sf.getLoad(w)

# ---------------------------------------------------------

TSParams = {"Qes":0.36, #
            "Qms":7.7, #
            "Fs":27, # Hz
            "Vas":37, # L
            "Mms":59, # g
            "Re":3.3, # ohm
            "Sd":210, # cm2
            "Le":0.73}  # mH

speaker = ElecDynSpeaker(TSParams, 0.006)

# frequency responce p_out/p_in with acostical load 
f = np.linspace(50, 10000)
omegas = f * 2 * np.pi
p_rel_u = np.empty(np.shape(omegas), dtype=complex)
for i, w in enumerate(omegas):
    sys = speaker.calculateSystem(w) 
    #load = speaker.calculateLoad(w)
    # schallharte wand 
    p_rel_u[i] = 1/sys[0,0] # assuming v = 0

#plt.plot(f, p_rel_u)
plt.plot(f, 20*np.log10(abs(p_rel_u)))
plt.show()

