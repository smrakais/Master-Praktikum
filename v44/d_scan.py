from abc import abstractstaticmethod
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat 
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from scipy.signal import peak_widths

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
print('Die Strahlbreite beträgt d = ' + str(d) +'mm.')
print('Aus der Strahlbreite d = ' + str(d) +'mm'+' und der Probendicke D = ' + str(D)+ 'mm,' )
print('ergibt sich mit arcsin(d/D) ein Geometriewinkel von ' + str(np.degrees(np.arcsin(d/D))) + '°.''\n')
fig.savefig('build/z_scan.pdf')

# Geometriewinkel
def geo(x):
    return D*np.degrees(np.sin(x)/d)

# Detektor Scan

# read data
detector_phi, detector_intensity = np.genfromtxt('data/Detektorscan.UXD', unpack=True)

# fitting function 
def gauss(x, amp, mu, sigma):
    return amp/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x- mu)**2/(2*sigma**2))

# fit and errors
detector_params, detector_pcov = curve_fit(gauss, detector_phi, detector_intensity)
detector_err = np.sqrt(np.diag(detector_pcov))

# Plot
detector_phi_new = np.linspace(detector_phi[0]-0.05, detector_phi[-1]+0.05, 10000)
plt.figure()
plt.xlabel(r"$\alpha_{i}$ / °")
plt.ylabel(r"Intensität")
plt.plot(detector_phi, detector_intensity, "x", label="Messwerte")
plt.plot(detector_phi_new, gauss(detector_phi_new, *detector_params), label="Ausgleichskurve")
plt.legend(loc="best") 
plt.grid(ls = '--')
plt.tight_layout()
plt.axhline(y = 483409,xmin=0.40, xmax=0.66 ,linestyle="--", color="r")
plt.axvline(x = -0.04,ymin=0, ymax=0.5 ,linestyle="--", color="r")
plt.axvline(x = 0.062,ymin=0, ymax=0.5 ,linestyle="--", color="r")
plt.xlim(-0.2,0.2) 

# max_value
print('Die maximale Detektorintensität beträgt: I = ' + str( np.max(detector_intensity)))
#print('Die halbe maximale Detektorintensität beträgt: I = ' + str( np.max((detector_intensity))/2))

# Halbwertsbreite
print("Die Halbwertsbreite(FWHM) beträgt: " + str(round(0.062-0.004,3)) +  "°" '\n')

plt.savefig("build/detector_scan.pdf")
#plt.show()



# Reflektivitätsscan und ideale Kurve
# read data
messung_phi, messung_intensity = np.genfromtxt("data/reflektivitatscan.UXD", unpack=True)

diffus_phi, diffus_intensity = np.genfromtxt("data/reflektivitatscan_2.UXD", unpack=True)

# relative data

i = 4 #damit ist der Graph etwas schöner

rel_phi = diffus_phi #[i:]  # der Winkel ändert sich ja nicht.
rel_phi = rel_phi[i:] 

messung_phi = messung_phi[i:]
messung_intensity = messung_intensity[i:]

diffus_phi = diffus_phi[i:]
diffus_intensity = diffus_intensity[i:]

rel_intensity = messung_intensity - diffus_intensity


# ideale Kurve der Reflektivität von Silizium

# kritischer Winkel
alpha_si_krit= 0.223 

def ideale_kurve(messung_phi):
    R_f =(alpha_si_krit/(2*messung_phi))**4 
    return R_f

####TEST####  JAN für ideale Kurve
#r_lambda = 1.54e-10
#k = 2*np.pi / r_lambda
#n = 1 - 7.6e-6 + 1.54e-8j*141/(4*np.pi)
## Ideale Kurve
#def ideal(alpha):
#    return (np.abs((k * np.sin(alpha)- k*np.sqrt(n**2-np.cos(alpha)**2))/(k * np.sin(alpha)+ k*np.sqrt(n**2-np.cos(alpha)**2))))**2
#############

# Plot
plt.figure()
plt.xlabel(r"$\alpha_{i}$ / °")
plt.yscale("log")
plt.ylabel(r"Reflektivität")
plt.plot(messung_phi, messung_intensity/np.max(detector_intensity),color = 'k', label="Messwerte")
plt.plot(diffus_phi, diffus_intensity/np.max(detector_intensity),color = 'r', label="Diffuser Scan")
plt.plot(rel_phi, rel_intensity/np.max(detector_intensity),color='orange', label="Korrigierte Messwerte")
plt.plot(messung_phi[35:],ideale_kurve(messung_phi)[35:],color='green', label = 'Theoriekurve glattes Si') # bei 35 abgeschnitten da dort der krit winkel
####TEST#### JAN für ideale Kurve
#plt.plot(rel_phi, ideal(np.deg2rad(rel_phi)), label="Ideale Siliziumoberfläche")
############
plt.plot(rel_phi, (rel_intensity/np.max(detector_intensity)*geo(rel_phi)),linewidth = 1.3,color ='b',label = 'Geometriefaktor Korrektur')


# Schichtdicke von PS abschätzen
lambda_ = 1.53*10**(-10)

def schichtdicke(delta):
    d = lambda_/(2*delta)
    return print('Die Schichtdicke beträgt ' + str(d) +'m.')

# Peaks finden
x_peak_factor = 2.5/484 # 0.0056338 für 484 ist es gut. eig müsste das 497 sein das 
y = rel_intensity/np.max(detector_intensity)*geo(rel_phi)

# OLD
#peaks, _  = find_peaks(y, height = 2 )
#plt.plot((peaks*x_peak_factor)[2:],y[peaks][2:],'x')

# NEW 
low  = argrelextrema(y, np.less)
x_new = x_peak_factor*low[0][2:9]
y_new = y[low[0][2:9]]
plt.plot(x_new,y_new,'x')

# Breite der Peaks
print(abs(x_new[0]-x_new[1]))
abstaende = np.array([abs(x_new[0]-x_new[1]),
                      abs(x_new[1]-x_new[2]),
                      abs(x_new[2]-x_new[3]),
                      abs(x_new[3]-x_new[4]),
                      abs(x_new[4]-x_new[5]),
                      abs(x_new[5]-x_new[6]),
                      ])
delta = np.mean(abstaende)                       
print('Der Mittelwert ist '+ str(delta) +'.')

# Ausgabe Schichtdicke
schichtdicke(delta)


plt.grid(ls= '--')
plt.minorticks_on()
plt.tight_layout()
plt.ylim(top = 1400)
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