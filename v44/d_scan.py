import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
#from scipy import stats
from uncertainties import ufloat 
from scipy.optimize import curve_fit
import matplotlib.patches as patches

# z_scan
 
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

# Geometriewinkel auch aus Strahldicke errechenbar.
d = 0.3
D = 20
print('Die Strahlbreite beträgt ' + str(d) +'mm.')
print('Aus der Strahlbreite ' + str(d) +'mm'+' und der Probendicke ' + str(D)+ 'mm' )
print('ergibt sich mit arcsin(d/D) ein Geometriewinkel von ' + str(np.degrees(np.arcsin(d/D))) + '°.')
fig.savefig('build/z_scan.pdf')


# Detektor Scan

# read data
detector_phi, detector_intensity = np.genfromtxt('data/Detektorscan.UXD', unpack=True)

# fitting function x = theta in degree
def gauss(x, amp, mu, sigma):
    return amp/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x- mu)**2/(2*sigma**2))

# fit and errors
detector_params, detector_pcov = curve_fit(gauss, detector_phi, detector_intensity)
detector_err = np.sqrt(np.diag(detector_pcov))

# Plot
detector_phi_new = np.linspace(detector_phi[0]-0.05, detector_phi[-1]+0.05, 10000)
plt.figure()
plt.xlabel(r"$\theta$ / °")
plt.ylabel(r"Intensität")
plt.plot(detector_phi, detector_intensity*1e-6, ".", label="Messwerte")
plt.plot(detector_phi_new, gauss(detector_phi_new, *detector_params)*1e-6, label="Ausgleichskurve")
plt.legend(loc="best") 
plt.grid(ls = '--')
plt.tight_layout()
plt.axhline(y = 0.483409,xmin=0.40, xmax=0.66 ,linestyle="--", color="r")
plt.axvline(x = -0.04,ymin=0, ymax=0.5 ,linestyle="--", color="r")
plt.axvline(x = 0.062,ymin=0, ymax=0.5 ,linestyle="--", color="r")
plt.xlim(-0.2,0.2) 

# max_value
print('Die maximale Detektorintensität beträgt: I = ' + str( np.max(detector_intensity)*1e-6))
print('Die halbe maximale Detektorintensität beträgt: I = ' + str( np.max((detector_intensity)*1e-6)/2))

# Halbwertsbreite
print("Die Halbwertsbreite beträgt: " + str(round(0.062-0.004,3)) +  "°")

plt.savefig("build/detector_scan.pdf")
#plt.show()



# Reflektivitätsscan und ideale Kurve
# read data
messung_phi, messung_intensity = np.genfromtxt("data/reflektivitatscan.UXD", unpack=True)

diffus_phi, diffus_intensity = np.genfromtxt("data/reflektivitatscan_2.UXD", unpack=True)

# relative data
rel_intensity = messung_intensity - diffus_intensity
rel_phi = diffus_phi

# ideale Kurve
alpha_si= 0.223

def ideale_kurve(x):
    R_f =(alpha_si/messung_phi)**4 
    return R_f

# Plot
plt.figure()
plt.xlabel(r"$\theta$ / °")
plt.yscale("log")
plt.ylabel(r"Reflektivität")
plt.plot(messung_phi, messung_intensity/np.max(detector_intensity), label="Messwerte")
plt.plot(diffus_phi, diffus_intensity/np.max(detector_intensity), label="Diffuser Scan")
plt.plot(rel_phi, rel_intensity/np.max(detector_intensity), label="Korrigierte Messwerte")
plt.plot(messung_phi,ideale_kurve(messung_phi), label = 'Theoriekurve glattes Si')
plt.grid(ls= '--')
plt.minorticks_on()
plt.tight_layout()
plt.ylim(top = 5)
plt.xlim(right = 1.5)
plt.axvline( x = 0.195,linewidth = 0.9,linestyle= '--',color ='k', label = r'$\alpha_c$ für PS')
plt.legend(loc="best") 
plt.savefig("build/messwerte_relativ.pdf")

alpha_crit= 0.195

# Korrektur und Geometriewinkel
# Geometriewinkel aus Rockingscan ablesen

theta , intensity = np.genfromtxt('data/dreieck_1.UXD', unpack = 'True')

plt.figure()
plt.xlabel(r"$\theta$ / °")
#plt.yscale("log")
plt.ylabel(r"Intensität ")
plt.plot(theta, intensity, "x", label="Messwerte")
plt.scatter(0.68,418,color= 'r') # g_winkel_1
plt.scatter(-0.76,778,color= 'r') # g_winkel_2
plt.legend(loc="best") 
plt.grid(ls= '--')
plt.tight_layout()
print('Der Geometriewinkel beträgt: ' + str((0.68 + abs(-0.76))/2 )+ '°.')
plt.savefig("build/dreieck.pdf")
