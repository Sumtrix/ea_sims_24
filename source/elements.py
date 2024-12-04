import numpy as np
from functools import reduce
import matplotlib.pyplot as plt


class phsyConst():
    def __init__(self):
        self.rho = 1.205
        self.c = 343
        self.mu = 1.8e-10
        self.p0 = 1013
        self.kappa = 1.4


class Zweitor():
    def __init__(self, Z, type=None):
        self.Z = Z
        self.type = type
        
    
    def matrix(self, w):
        Z_w = self.Z(w) if callable(self.Z) else self.Z  # Z(w), falls Z frequenzabhängig ist
        if self.type == "laengs":
            return np.array([[1, Z_w],
                             [0, 1]])
        elif self.type == "quer":
            Y_w = 1 / Z_w
            return np.array([[1, 0],
                             [Y_w, 1]])
        else:
            raise ValueError("Ungültiger Typ. Bitte 'laengs' oder 'quer' wählen.")


    # def matrix(self, w):
    #     if self.laengs:
    #         match self.type:
    #             case "0100":
    #                 res = [[1, 1j*w*self.val],[0, 1]]
    #             case "0200":
    #                 res = [[1, self.val/(1j*w)],[0, 1]]
    #             case _:
    #                 res = np.matrix([[1, self.val], [0, 1]])
    #     else:
    #         match self.type:
    #             case "0010":
    #                 res = [[1, 0],[(1/self.val)*(1/(1j*w)), 1]] # double inverse due to kehrwert
    #             case "0020":
    #                 res = [[1, 0],[(1/self.val)*(1j*w), 1]] # double inverse due to kehrwert
    #             case _:
    #                 res = np.matrix([[1, 0], [1/self.val, 1]])
    #     return res


class elecMechTransf():
    def __init__(self, Bl):
        self._val = Bl
        self._mat = np.array([[Bl, 0], [0, 1/Bl]])

    def matrix(self):
        return self._mat
    

class mechAcoustTransf():
    def __init__(self, S):
        self._val = S
        self._mat = np.array([[0, 1/S], [S, 0]])

    def matrix(self):
        return self._mat
        

class SoundField(phsyConst):
    def __init__(self, S):
        super().__init__()
        self.S = S
        self.a = np.sqrt(S/np.pi)
    
    def getLoad(self, w):
        k = w / self.c
        #self.Zs = self.rho*self.c*((0.5*(k*self.a)**2) + (1j*(8/(3*np.pi))*k*self.a))      # false
        self.Zs = self.rho*self.c*((0.06*(k*self.a)**4) + (1j*(8/(3*np.pi))*k*self.a))    # full function
        #self.Zs = self.rho*self.c*(1j*(8/(3*np.pi))*k*self.a)      # ka<<1
        self.Z = self.Zs / self.S
        return self.Z


class ElecDynSpeaker(phsyConst, elecMechTransf):
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
        self.R = Zweitor(self.RE, "laengs")
        self.L = Zweitor(lambda w: 1j*w*self.LE, "laengs")
        self.Bl = elecMechTransf(self.BL)
        self.m = Zweitor(lambda w: 1j*w*self.MMS, "quer")
        self.s = Zweitor(lambda w: 1/(1j*w*self.CMS), "quer")
        self.r = Zweitor(self.RMS, "quer")
        self.S = mechAcoustTransf(self.SD)
        self.sf = SoundField(self.SD)


    def calculateSystem(self, w):
        #matrices = [self.R.getM(w), self.L.getM(w), self.BL.getM(), self.m.getM(w), self.s.getM(w), self.r.getM(w), self.S.getM()]
        
        matrices = [self.R.matrix(w), 
                    self.L.matrix(w), 
                    self.Bl.matrix(), 
                    self.m.matrix(w), 
                    self.s.matrix(w), 
                    self.r.matrix(w), 
                    self.S.matrix()]
        sys_matrix = np.eye(2)  # Start mit Identitätsmatrix
        for mat in matrices:
            sys_matrix = np.dot(sys_matrix, mat)  # Matrixmultiplikation
        return sys_matrix
    
    def calculateLoad(self, w):
        return self.sf.getLoad(w)
    
    def calculateImpedance(self, w):
        # transform acoustic imp to mechanical
        k = w / self.c
        a = np.sqrt(self.SD/np.pi)      # sqrt? script shows something else but i think it's an error
        rak = self.SD*self.rho*self.c*(0.006*(k*a)**4)
        mak = self.rho*(8/3)*a**3
        mstar = self.MMS + mak
        rstar = self.RMS + rak
        imp_mech_ac = 1 / ((1j*w*(mstar/(self.BL)**2))+((1/self.CMS)/(1j*w*(self.BL)**2))+(rstar/(self.BL)**2))
        self.imp = imp_mech_ac + self.R.Z + self.L.Z(w) 
        return self.imp


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
f = np.logspace(1, 4.4, 20000)
omegas = f * 2 * np.pi
pu_unterload = np.empty(np.shape(omegas), dtype=complex)
p_rel_u = np.empty(np.shape(omegas), dtype=complex)
imp = np.empty(np.shape(omegas), dtype=complex)

for i, w in enumerate(omegas):
    sys = speaker.calculateSystem(w) 
    load = speaker.calculateLoad(w)     #  imp = p / v
    imp[i] = speaker.calculateImpedance(w)
    pu_unterload[i] = 1 / (sys[0, 0] + sys[0, 1] * load)
    #p_rel_u[i] = 1/sys[0,0] # assuming v = 0
    

plt.semilogx(f, 20 * np.log10(abs(pu_unterload)))
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude (dB)")
plt.grid(True, which="both", linestyle="--")
plt.show()

plt.semilogx(f, abs(imp))
plt.xlabel("Frequency (Hz)")
plt.ylabel("Impedance")
plt.grid(True, which="both", linestyle="--")
plt.show()