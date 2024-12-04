# sim_test.py

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
