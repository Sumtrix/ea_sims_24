from numpy import array, dot, pi, log10, matrix
import numpy as np
import matplotlib.pyplot as plt


# -------------------------------------------
# speaker and sim parameters

Re=3.5
Le=0.44 
Bl=2.68 
SD=15.6
Cms=1.26
Mms=2.5
Fs=89.5
Vas=0.53
Qms=5.56
Qes=0.68


f = np.linspace(10, 20000)
w = 2 * pi * f


# -------------------------------------------
# constants

roh = 1.205
c = 343
mu = 1.8e-10
p0 = 1013
kappa = 1.4

tsp = {"Re" : Re,
        "Le" : Le ,            # mH
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
        

# -------------------------------------------

A = np.empty((2, 2, len(w)), dtype=complex)
v_over_u = np.empty(len(w), dtype=complex)
for ww_idx, ww in enumerate(w):
    # electrical
    R = array([[1, tsp["Re"]], [0, 1]])
    Le = array([[1, tsp["Le"]], [0, 1]])
    Bl = array([[tsp["Bl"], 0], [0, 1/tsp["Bl"]]])
    tot_el = np.matmul(R, Le)
    bl_out = np.matmul(tot_el, Bl)
    
    # mechanical
    m = array([[1, 0], [1j*ww*tsp["Mms"], 1]])
    r = array([[1, 0],[tsp["Rms"], 1]])
    s = array([[1, 0], [1j*ww*tsp["Cms"], 1]])
    S = array([[0, 1/tsp["SD"]], [tsp["SD"], 0]])
    tot_mech = np.matmul(np.matmul(m, r), s) 
    Mech_out = np.matmul(bl_out, tot_mech)
    S_out = np.matmul(Mech_out, S)
                
    # acoustic impedance
    k = ww / c
    a = np.sqrt(tsp["SD"] / pi)
    
    box_volume = 0.5 * 0.001

    # dayton passive radiator : DMA70-PR 2-1/2" 
    fs_pm = 34.5
    Cms = 1
    
    nv = box_volume / (roh * c**2) 
    mpr = 1
    0.00177



    NV = np.array([1,0],[1j*ww*nv,1])
    MPM = np.array([1,1j*ww*mpr],[0,1])
    NPM = np.array([1,0],[1j*ww*npm,1])
    ZSS = np.array([1,0],[1j*ww*zss,1])
    acoustic_load = np.matmul(np.matmul(np.matmul(NV, MPM), NPM), ZSS) 




    

# -------------------------------------------
# plot results

plt.semilogx(w/(2*pi), log10(abs(v_over_u)))
plt.grid(which="both")
