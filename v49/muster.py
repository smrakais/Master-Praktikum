import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from scipy.signal import peak_widths
from texutils.table import TexTable

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18


#########################################################################################
# objektorientierter Plot

z, z_intensity = np.genfromtxt('data/z_scan_test.UXD', unpack = True)

fig, ax = plt.subplots()
ax.set_xlabel(r'$z$ / mm')
ax.set_ylabel(r'Intensität')
ax.grid(ls =  '--')
ax.minorticks_on()        # x,    y,  Breite, Höhe
rect = patches.Rectangle((-0.152, 0), 0.30, 0.99,   linewidth=1, 
                                                    edgecolor='g', 
                                                    facecolor='none',
                                                    label = 'Strahlbreite')
ax.add_patch(rect)
ax.plot(z,z_intensity*1e-6,'k-', label = 'Messwerte' )
ax.legend(loc ='best')
plt.tight_layout()
ax.axvline(x = 0.0263184, linestyle="--", color="r")
ax.axhline(y = 0.497296, linestyle="--", color="r")
plt.show()

#########################################################################################
#Fehlerrechnung mit curve fit und plot
########################################################################################################
# read data
detector_phi, detector_intensity = np.genfromtxt('data/Detektorscan.UXD', unpack=True)

# fitting function 
def gauss(x, amp, mu, sigma):
    return amp/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x- mu)**2/(2*sigma**2))

# fit and errors
detector_params, detector_pcov = curve_fit(gauss, detector_phi, detector_intensity)
detector_err = np.sqrt(np.diag(detector_pcov))
print('amp =', detector_params[0], '±', detector_err[0])
print('alpha_0 =', detector_params[1], '±', detector_err[1])
print('sigma =', detector_params[2], '±', detector_err[2],'\n')

# Der Plot zum curve fit
detector_phi_new = np.linspace(detector_phi[0]-0.05, detector_phi[-1]+0.05, 10000)
plt.figure()                        #create new fig
plt.xlabel(r"$\alpha_{i}$ / °")
plt.ylabel(r"Intensität")
plt.plot(detector_phi, detector_intensity, "x", label="Messwerte")
plt.plot(detector_phi_new, gauss(detector_phi_new, *detector_params), label="Ausgleichskurve")
plt.legend(loc="best") 
plt.grid(ls = '--')
plt.tight_layout()

plt.show()
########################################################################################################

##########
# Tabelle
##########
t = TexTable([dt, N], [r"$t$ / ns",r"Zählrate(10s) "], 
            label='tab:Koinzidenz',
            caption='Messwerte der Zählraten in Abhängigkeit der Verzögerungsleitungen.')
t.set_row_rounding(0, 1) #reihe und rundung
t.set_row_rounding(1, 0)
t.write_file('build_new/tabKoinzidenz_new.tex')
print('Die Tabelle der Koinzidenz wurde erzeugt!\n')