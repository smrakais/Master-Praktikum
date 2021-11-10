import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
#from scipy import stats
from uncertainties import ufloat 
from scipy.optimize import curve_fit

# z_scan
 
z, z_intensity = np.genfromtxt('data/z_scan_test.UXD', unpack = True)

plt.figure()
plt.xlabel(r'$z$ / mm')
plt.ylabel(r'Intensität')
plt.grid()
plt.minorticks_on
plt.plot(z,z_intensity*1e-6,'k-', label = 'Messwerte' )
plt.legend(loc ='best')
plt.tight_layout()
plt.axvline(x = 0.0263184, linestyle="--", color="r")
plt.axhline(y = 0.497296, linestyle="--", color="r")
plt.savefig('build/z_scan.pdf')



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
plt.grid()
plt.tight_layout()
plt.axhline(y = 0.483409,xmin=0.45, xmax=0.58 ,linestyle="--", color="r")
plt.axvline(x = -0.04,ymin=0, ymax=0.5 ,linestyle="--", color="r")
plt.axvline(x = 0.062,ymin=0, ymax=0.5 ,linestyle="--", color="r")

# max_value
print(np.max(detector_intensity)*1e-6)
print(np.max((detector_intensity)*1e-6)/2)

# Halbwertsbreite
print("Die Halbwertsbreite beträgt: " + str(round(0.062-0.004,3)) +  "°")

plt.savefig("build/detector_scan.pdf")
#plt.show()



# Messung und Diffusionskorrektur
#read data
messung_phi, messung_intensity = np.genfromtxt("data/reflektivitatscan.UXD", unpack=True)

diffus_phi, diffus_intensity = np.genfromtxt("data/reflektivitatscan_2.UXD", unpack=True)

#relative data
rel_intensity = messung_intensity - diffus_intensity
rel_phi = diffus_phi

#Plot
plt.figure()
plt.xlabel(r"$\theta$ / °")
plt.yscale("log")
plt.ylabel(r"Intensität")
plt.plot(messung_phi, messung_intensity, label="Messwerte")
plt.plot(diffus_phi, diffus_intensity, label="Diffuser Scan")
plt.plot(rel_phi, rel_intensity, label="Korrigierte Messwerte")
plt.legend(loc="best") 
plt.grid()
plt.tight_layout()
plt.savefig("build/messwerte_relativ.pdf")


#korrektur und geometriewinkel
