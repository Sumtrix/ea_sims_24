import numpy as np
import matplotlib.pyplot as plt


class Zweitor:
    def __init__(self, Z, typ="längs"):
        """
        Initialisiert ein Zweitor.
        
        Parameters:
        Z: Funktion oder konstante Impedanz (z.B. lambda w: 1/(1j*w*m) für Masse)
        typ: "längs" oder "quer", um die Anordnung des Zweitors anzugeben
        """
        self.Z = Z
        self.typ = typ.lower()

    def matrix(self, w):
        """
        Gibt die Zweitormatrix in Abhängigkeit von der Kreisfrequenz w zurück.
        
        Parameters:
        w: Kreisfrequenz (rad/s)
        
        Returns:
        2x2-Matrix der Zweitor-Übertragungsfunktion
        """
        Z_w = self.Z(w) if callable(self.Z) else self.Z  # Z(w), falls Z frequenzabhängig ist
        if self.typ == "längs":
            # Längselement: Z als Serienimpedanz
            return np.array([[1, Z_w],
                             [0, 1]])
        elif self.typ == "quer":
            # Querelement: Z als Parallelimpedanz
            Y_w = 1 / Z_w
            return np.array([[1, 0],
                             [Y_w, 1]])
        else:
            raise ValueError("Ungültiger Typ. Bitte 'längs' oder 'quer' wählen.")

class Speaker:
    def __init__(self, thiele_small_parameters):
        """
        Initialisiert den Lautsprecher mit Thiele-Small-Parametern.
        """
        self.ts = thiele_small_parameters  # Dictionary der TS-Parameter
        self.build_system()

    def build_system(self):
        """
        Erstellt die Zweitore des Lautsprechersystems basierend auf den TS-Parametern.
        """
        R = self.ts['R']
        L = self.ts['L']
        Bl = self.ts['Bl']
        m = self.ts['m']
        s = self.ts['s']
        r = self.ts['r']
        S = self.ts['S']
        
        # Elemente als Zweitore initialisieren
        self.R = Zweitor(lambda w: R, typ="längs")
        self.L = Zweitor(lambda w: 1j * w * L, typ="längs")
        self.Bl = Zweitor(lambda w: 1 / (Bl), typ="quer")
        self.m = Zweitor(lambda w: 1j * w * m, typ="längs")
        self.s = Zweitor(lambda w: 1 / (1j * w * s), typ="längs")
        self.r = Zweitor(lambda w: r, typ="längs")
        self.S = Zweitor(lambda w: 1 / S, typ="quer")  # Akustische Impedanz

    def calculate_sys(self, w):
        """
        Berechnet die Gesamtsystemmatrix bei einer bestimmten Kreisfrequenz.
        
        Parameters:
        w: Kreisfrequenz (rad/s)
        
        Returns:
        Gesamtsystemmatrix
        """
        # Matrizen berechnen
        matrices = [self.R.matrix(w), 
                    self.L.matrix(w), 
                    self.Bl.matrix(w), 
                    self.m.matrix(w), 
                    self.s.matrix(w), 
                    self.r.matrix(w), 
                    self.S.matrix(w)]
        
        # Matrizen multiplizieren
        sys_matrix = np.eye(2)  # Start mit Identitätsmatrix
        for mat in matrices:
            sys_matrix = np.dot(sys_matrix, mat)  # Matrixmultiplikation
        
        return sys_matrix

# Beispiel: Thiele-Small-Parameter definieren
thiele_small_parameters = {
    'R': 6.0,       # Spulenwiderstand in Ohm
    'L': 0.001,     # Spuleninduktivität in H
    'Bl': 5.0,      # Elektro-mechanischer Wandler in T·m
    'm': 0.02,      # Bewegte Masse in kg
    's': 0.0001,    # Mechanische Aufhängung in N/m
    'r': 1.0,       # Mechanische Reibung in Ns/m
    'S': 0.01       # Mechanisch-akustischer Wandler in m²/Pa
}

# Lautsprecher instanziieren
speaker = Speaker(thiele_small_parameters)

# Beispiel: Gesamtsystemmatrix bei w = 100 rad/s berechnen
# w = 100
# sys_matrix = speaker.calculate_sys(w)
# print("Systemmatrix bei w = 100 rad/s:")
# print(sys_matrix)


f = np.linspace(50, 10000)
omegas = f * 2 * np.pi
p_rel_u = np.empty(np.shape(omegas), dtype=complex)
for i, w in enumerate(omegas):
    sys = speaker.calculate_sys(w) 
    #load = speaker.calculateLoad(w)
    # schallharte wand 
    p_rel_u[i] = 1/sys[0,0] # assuming v = 0

#plt.plot(f, p_rel_u)
plt.semilogx(f, 20*np.log10(abs(p_rel_u)))
plt.show()
