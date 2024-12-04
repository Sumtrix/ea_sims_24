from numpy import array, dot, pi, log10, matrix
import numpy as np
import matplotlib.pyplot as plt
from matrizen import *


roh = 1.205
c = 343
mu = 1.8e-10
p0 = 1013
kappa = 1.4


class Speaker():
    def __init__(self, Re, Le, Bl, SD, Cms, Mms, Fs, Vas, Qms, Qes):
        self.tsp = {"Re" : Re,
               "Le" : Le * 0.001,            # mH
               "Bl" : Bl,
               "SD" : SD * 0.0001,           # cm2
               "Cms" : Cms * 0.001, 
               "Mms" : Mms * 0.001,          # gramm
               "Rms" : (1/Cms) / ((Fs*2*pi) * Qms),
            #    "Cmes" : Cmes,
            #    "Lces" : Lces,
            #    "Res" : Res,
               "Fs" : Fs,
               "Vas" : Vas * 0.001,          # Liter
               "Qms" : Qms,
               "Qes" : Qes,
               "Qts" : 1/((1/Qms)+(1/Qes))}
                

    def update(self, w):
        '''
        Build the entire transfer matrix with current values.
        '''
        self.A = np.empty((2, 2, len(w)), dtype=complex)
        self.v_over_u = np.empty(len(w), dtype=complex)
        for ww_idx, ww in enumerate(w):
            # electrical
            R = laengs_mat(self.tsp["Re"])
            Le = quer_mat(1/self.tsp["Le"])
            Bl = array([[self.tsp["Bl"], 0], [0, 1/self.tsp["Bl"]]])
            self.tot_el = np.matmul(R, Le)
            self.bl_out = np.matmul(self.tot_el, Bl)
            
            # mechanical
            m = quer_mat(1/(1j*ww*self.tsp["Mms"]))
            r = quer_mat(1/(self.tsp["Rms"]))
            s = quer_mat(1j*ww*self.tsp["Cms"])
            S = array([[0, 1/self.tsp["SD"]], [self.tsp["SD"], 0]])
            self.tot_mech = np.matmul(np.matmul(m, r), s) 
            self.Mech_out = np.matmul(self.bl_out, self.tot_mech)
            self.S_out = np.matmul(self.Mech_out, S)
                        
            # acoustic impedance
            k = ww / c
            a = np.sqrt(self.tsp["SD"] / pi)
            # ea lecture script p. 53 eqation (5.23) - impedance in inf. wall
            self.ac_imp = 1 / (self.tsp["SD"] * roh * c * (0.5 * (k * a)**2 + 1j * (8 / 3 * pi) * k * a))
            # ea lecture script p. 52 eqation (5.22)
            #self.ac_load = 
            
            
            # parallel restiances 
            m_star = self.tsp["Mms"] + roh*(8/3)*a**3
            m = array([[1, 1j*ww*m_star], [0, 1]])

            s = array([[1, 0], [1/(1j*ww*self.tsp["Cms"]), 1]])
            
            r_star = self.tsp["Rms"] + (0.5*roh*c*self.tsp["SD"]*(k*a)**2)
            r_el = self.tsp["Bl"]**2 / (self.tsp["Re"] + self.tsp["Le"])
            r = array([[1, 0], [r_star+r_el, 1]])

            temp = np.matmul(m, s)
            total = np.matmul(temp, r)
            
            # model on p. 54
            #source = self.tsp["Bl"] / self.tsp["Re"]
            self.v_over_u[ww_idx] = 1 / (total[0, 0] * (self.tsp["Bl"]/self.tsp["Re"]))
            
            # frequency responce
            #Zs = roh * c * (0.06 * (k * a)**4 + 1j * (8 / 3 * pi) * k * a)
            #freq_resp = V**2 / u**2 * self.tsp["SB"] * np.real(Zs)
            
    def plot__p_over_u(self, w):
        self.update(w)
        plt.semilogx(w/(2*pi),log10(abs(self.v_over_u)))
        plt.grid(which="both")
        
        
        

class PassiveRadiatorEnclosure():
    def __init__(self):
        tsp = {"Vas" : None,
               "Vas" : None,
               "Vas" : None,
               "Vas" : None,
               "Vas" : None,}
    
    def update(self):
        pass
    
    
class BassReflexEnclosure():
    def __init__(self):
        tsp = {"Vas" : None,
               "Vas" : None,
               "Vas" : None,
               "Vas" : None,
               "Vas" : None,}
    
    def update(self):
        pass
    
        
class Volume():
    def __init__(self, V):
        self.V = V
        self.Na = self.V / (roh * (c**2))
        
    def update(self):
        pass
        
        
class Port():
    def __init__(self, l, d):
        l_corr = l + l * (8/(3*pi) * (d/2))
        S = pi * (d/2)**2
        self.Ma = (roh * l_corr) / S
    
    def update(self):
        pass



if __name__ == "__main__":

    import sim_elements as se
    import numpy as np

    # Dayton ND654 21/2

    sp = se.Speaker(Re=3.5,
                Le=0.44, 
                Bl=2.68, 
                SD=15.6,
                Cms=1.26,
                Mms=2.5,
                Fs=89.5,
                Vas=0.53,
                Qms=5.56,
                Qes=0.68)
    print("yes")
    f = np.linspace(10, 20000)
    w = 2 * pi * f
    sp.plot__p_over_u(w)
    plt.ylim([-10, 0])
    plt.show()
