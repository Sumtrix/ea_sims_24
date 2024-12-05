import numpy as np
from functools import reduce
import matplotlib.pyplot as plt


def multMatrices2x2(matrices):
    sys_matrix = np.eye(2)  # Start mit Identitätsmatrix
    for mat in matrices:
        sys_matrix = np.dot(sys_matrix, mat)  # Matrixmultiplikation
    return sys_matrix


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
        

class AcousticVolume(phsyConst, Zweitor):
    def __init__(self, V):
        phsyConst.__init__(self)
        self.V = V * 1e-3 # Liter to m3
        Na = V / (self.rho * self.c**2)
        Zweitor.__init__(self, lambda w : 1/(1j*w*Na), "quer")


class AcousticPort(phsyConst, Zweitor):
    def __init__(self, r, l, type="tube",  corr="wall"):
        '''
        r : floar
        Radious of port, if square port is used, then r is half the height/width
        '''
        phsyConst.__init__(self)
        match type:
            case "tube":
                a = r
                self.S = np.pi * a**2
                xi = (8*self.mu)/a**2
                match corr:
                    case "wall":
                        deltal = (8/(3*np.pi)) * a
                    case "freeair":
                        deltal = 0.6 * a
            case "square":
                h = r*2
                self.S = h*2
                xi = (12*self.mu)/self.S
                deltal = 0
        self.lstar = l + deltal
        self.Ma = (self.rho * self.lstar) / self.S
        self.Za = (xi * self.lstar) / self.S
        # coumpling to large volume basically no effect, thus not considered
        Zweitor.__init__(self, lambda w : 1j*w*self.Ma + self.Za, "laengs")


class PortedBox(phsyConst):
    def __init__(self, V, r, l, type="tube", corr="wall"):
        phsyConst.__init__(self)
        Vm3 = V * 1e-3
        self.encl = AcousticVolume(Vm3)
        self.port = AcousticPort(r, l, type, corr)
        self.sf = SoundField(self.port.S)
        self.f_tune = (self.c/(2*np.pi)) * np.sqrt(self.port.S / (self.port.lstar * Vm3))
        print(f"Tuned to {round(self.f_tune, 2)} Hz")
        
    def calculateSystem(self, w):
        matrices = [self.encl.matrix(w), 
                    self.port.matrix(w)]
        return multMatrices2x2(matrices)
        
    def calculateLoad(self, w):
        ZsS = self.sf.calculateLoad(w)
        load = ZsS
        #load = self.port.Z(w) + ZsS
        #load = 1 /((1/self.encl.Z(w))+(1/(self.port.Z(w)+ZsS)))
        return load
    
    def calculateImpedance(self, w):
        return self.encl.Z(w) + 1/self.port.Z(w)


class SoundField(phsyConst):
    def __init__(self, S):
        super().__init__()
        self.S = S
        self.a = np.sqrt(S/np.pi)
    
    def calculateLoad(self, w):
        k = w / self.c
        # elf.Zs = self.rho*self.c*((0.5*(k*self.a)**2) + (1j*(8/(3*np.pi))*k*self.a))      # Kolben in Schallwand
        self.Zs = self.rho*self.c*((0.06*(k*self.a)**4) + (1j*(8/(3*np.pi))*k*self.a))    # Kreisplatte in free air
        self.Z = self.Zs / self.S
        return self.Z


class ElecDynSpeaker(phsyConst, elecMechTransf):
    def __init__(self, TSParams, xmax):
        phsyConst.__init__(self)

        self.QES = TSParams["Qes"] 
        self.QMS = TSParams["Qms"]
        self.FS = TSParams["Fs"]
        self.VAS = TSParams["Vas"] * 1e-3
        self.MMS = TSParams["Mms"] * 1e-3
        self.RE = TSParams["Re"]
        self.SD = TSParams["Sd"] * 1e-4 
        self.LE = TSParams["Le"] * 1e-3
        self.RMS = TSParams["Rms"]
        self.CMS = TSParams["Cms"] * 1e-3
        self.BL = TSParams["Bl"]
        self.CMES = self.MMS / (self.BL)**2
        self.LCES = (self.BL**2) / self.SD
        self.RES = self.BL**2 / self.RE
        self.QTS = (self.QMS * self.QES) / (self.QMS + self.QES)
        self.XMAX = xmax
        self.build_system()


    def build_system(self):
        self.R = Zweitor(self.RE, "laengs")
        self.L = Zweitor(lambda w: 1/(1j*w*self.LE), "laengs")
        self.Bl = elecMechTransf(self.BL)
        self.m = Zweitor(lambda w: 1/(1j*w*self.MMS), "quer")
        self.s = Zweitor(lambda w: 1j*w*self.CMS, "quer")
        self.r = Zweitor(1/self.RMS, "quer")
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
        return multMatrices2x2(matrices)
    
    
    def calculateLoad(self, w):
        return self.sf.calculateLoad(w)
    
    
    def calculateImpedance(self, w):
        # transform acoustic imp to mechanical
        # k = w / self.c
        # a = np.sqrt(self.SD/np.pi)      # sqrt? script shows something else but i think it's an error
        # rak = self.SD*self.rho*self.c*(0.006*(k*a)**4)
        # mak = self.rho*(8/3)*a**3
        # mstar = self.MMS + mak
        # rstar = self.RMS + rak 
        # s = 1 / self.CMS
        # imp_mech_ac = 1 / ((1j*w*(mstar/self.BL**2))+\
        #                    (s/(1j*w*self.BL**2))+\
        #                    (rstar/(self.BL**2)))
        # m* and r* from TSParameters are already the transformed acoustical + mechanical masses and frictions
        imp_mech_ac = 1 / ((1j*w*(self.MMS/(self.BL)**2))+((1/self.CMS)/(1j*w*(self.BL)**2))+(self.RMS/(self.BL)**2))
        self.imp = imp_mech_ac + self.R.Z + 1/self.L.Z(w) 
        return self.imp


# ---------------------------------------------------------

TSParams = {"Qes":0.36, #
            "Qms":7.7, #
            "Fs":27, # Hz
            "Vas":37, # L
            "Mms":59, # g
            "Re":3.3, # ohm
            "Sd":210, # cm2
            "Le":0.73,  # mH
            "Rms":1.3,
            "Cms":0.59,  # mm/n
            "Bl":9.6}    # Tm
            

speaker = ElecDynSpeaker(TSParams, 0.006)
# enclosure = AcousticVolume(30) # L
# port = AcousticPort(r=0.07, l=37.5)
#port_sf = SoundField(BR_Box.port.S)
V = 30 # L
r = 0.035 # m
l = 0.375 # m
BR_Box = PortedBox(V, r, l, type="tube", corr="wall")


# frequency responce p_out/p_in with acostical load 
f = np.logspace(1, 4.4, 20000)
omegas = f * 2 * np.pi

speaker_imp = np.empty(np.shape(omegas), dtype=complex)
speaker_tf = np.empty(np.shape(omegas), dtype=complex)
br_tf = np.empty(np.shape(omegas), dtype=complex)
brbox_imp = np.empty(np.shape(omegas), dtype=complex)
pu_speaker_front = np.empty(np.shape(omegas), dtype=complex)
pu_port_out = np.empty(np.shape(omegas), dtype=complex)
pu = np.empty(np.shape(omegas), dtype=complex)
sfload = np.empty(np.shape(omegas), dtype=complex)

for i, w in enumerate(omegas):
    speaker_imp[i] = speaker.calculateImpedance(w)
    speak_sys = speaker.calculateSystem(w) 
    speaker_load = speaker.calculateLoad(2)
    speaker_tf[i] = 1/speak_sys[0,0]        # assuming v = 0       # ! check whether i changed that this does not include S mat - for only speaker the excursion of speaker is more interesting. 
    
    brbox_imp[i] = BR_Box.calculateImpedance(w)
    brbox_sys = BR_Box.calculateSystem(w)
    brbox_load = BR_Box.calculateLoad(w) 
    br_tf[i] = 1/ (brbox_sys[0,0] + brbox_sys[0, 1] * (1/brbox_load))
    
    comb_sys = multMatrices2x2([speak_sys, brbox_sys])

    pu_speaker_front[i] = 1/(speak_sys[0,0] + (speak_sys[0, 1] * (1/speaker_load)))
    pu_port_out[i] =  1/(comb_sys[0,0] + (comb_sys[0, 1] * (1/brbox_load)))
    
plt.semilogx(f, 20 * np.log10(abs(pu_speaker_front)), linestyle="--", color="k", alpha=0.3)
plt.semilogx(f, 20 * np.log10(abs(pu_port_out)), linestyle=":", color="k", alpha=0.3)
plt.semilogx(f, 20 * np.log10(abs(pu_speaker_front+pu_port_out)), linestyle="-", color="k")
plt.legend(["speaker front p", "port out p", "total"])

# plt.semilogx(f, 20 * np.log10(abs(speaker_tf)), linestyle="--", color="b", alpha=0.3)
plt.semilogx(f, 20 * np.log10(abs(br_tf)), linestyle=":", color="b", alpha=0.3)

plt.xlabel("Frequency (Hz)")
plt.ylabel(r"$\vert\dfrac{\hat{p}}{\hat{u}}\vert$ (dB)")
plt.grid(True, which="both", linestyle="--")
plt.show()

plt.semilogx(f, abs(speaker_imp), linestyle="-", color="k",)
plt.semilogx(f, abs(brbox_imp), linestyle="-", color="k",)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Impedance")
plt.grid(True, which="both", linestyle="--")
plt.show()