#TODO

#in line 204 müssen die werte für die minima angepasst werden und mach das am besten wieder mit plt.show()

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from scipy.signal import peak_widths

import math
import cmath as cm

#########################################################################################
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
#########################################################################################

# Korrektur und Geometriewinkel
# Geometriewinkel aus Rockingscan ablesen

theta , intensity = np.genfromtxt('data/dreieck_1.UXD', unpack = 'True')

plt.figure()
#plt.xlabel(r"$\theta$ / °")
plt.xlabel(r"$\alpha_i$ / °")#cheat
#plt.yscale("log")
plt.ylabel(r"Intensität ")
plt.plot(theta, intensity, "x", label="Messwerte")
plt.scatter(0.68,418,color= 'r') # g_winkel_1
plt.scatter(-0.76,778,color= 'r') # g_winkel_2
plt.legend(loc="best") 
plt.grid(ls= '--')
plt.tight_layout()
print('Der Geometriewinkel beträgt durch ablesen am Graphen: ' + str((0.68 + abs(-0.76))/2 )+ '°.\n')
plt.savefig("build/dreieck.pdf")
########################################################################################################

# Geometriewinkel auch aus Strahldicke errechenbar.
d = 0.3
D = 20
print('Die Strahlbreite beträgt d = ' + str(d) +'mm.')
print('Aus der Strahlbreite d = ' + str(d) +'mm'+' und der Probendicke D = ' + str(D)+ 'mm,' )
print('ergibt sich mit arcsin(d/D) ein Geometriewinkel von ' + str(np.degrees(np.arcsin(d/D))) + '°.''\n')  
fig.savefig('build/z_scan.pdf')
########################################################################################################
# Geometriefaktor (eigentlich ein Geometriedivisor, deswegen ist das der Kehrwert)
def geo(x):
    if x.any()<0.72: 
        return 1/((D*np.degrees(np.sin(x)))/d)
    else:
       return 1
#########################################################################################################

# Detektor Scan

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
########################################################################################################

# Plot detector scan
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
print("Die Halbwertsbreite(FWHM) beträgt: " + str(round(0.062+0.04,3)) +  "°" '\n')

plt.savefig("build/detector_scan.pdf")

########################################################################################################
############################### Beginn der eigentlichen Auswertung #####################################
########################################################################################################

# Reflektivitätsscan und ideale Kurve
# read data
messung_phi, messung_intensity = np.genfromtxt("data/reflektivitatscan.UXD", unpack=True)
diffus_phi, diffus_intensity = np.genfromtxt("data/reflektivitatscan_2.UXD", unpack=True)

# Anpassungen der Messwerte
i = 5 #damit ist der Graph etwas schöner

rel_phi = diffus_phi #[i:]  # der Winkel ändert sich ja nicht.
rel_phi = rel_phi[i:] 

messung_phi = messung_phi[i:]
messung_intensity = messung_intensity[i:]

diffus_phi = diffus_phi[i:]
diffus_intensity = diffus_intensity[i:]

rel_intensity = messung_intensity - diffus_intensity

#Berechnung Reflektivität / 5 wegen der Messdauer von 5 sekunden imm gegensatz zu den kalibrationsmessungen
refektivitaet_messwert = messung_intensity/(np.max(detector_intensity)*5)
refektivitaet_diffus = diffus_intensity/(np.max(detector_intensity)*5)
refektivitaet_rel = rel_intensity/(np.max(detector_intensity)*5)

# ideale Kurve der Reflektivität von Silizium
# kritischer Winkel
alpha_si_krit= 0.223 

def ideale_kurve(messung_phi):
    R_f =(alpha_si_krit/(2*messung_phi))**4 
    return R_f

# Plot
plt.figure()
plt.xlabel(r"$\alpha_{i}$ / °")
plt.yscale("log")
plt.ylabel(r"Reflektivität")
plt.plot(messung_phi, refektivitaet_messwert,color = 'k', label="Messwerte")
plt.plot(diffus_phi, refektivitaet_diffus,color = 'r', label="Diffuser Scan")
plt.plot(rel_phi, refektivitaet_rel,color='orange', label="Korrigierte Messwerte")
plt.plot(messung_phi[34:],ideale_kurve(messung_phi)[34:],color='green', label = 'Theoriekurve glattes Si') # bei 27 abgeschnitten da dort der krit winkel

# kritische Winkel Ablesen
alpha_crit_ps= 0.0625          # Polystyrol (abgelesen)
plt.axvline( x = alpha_crit_ps,linewidth = 0.9,linestyle= '--',color ='k', label = r'$\alpha_c$ für PS (abgelesen)')
print('\n')
print('Der abgelesene kritische Winkel für Polystyrol beträgt',alpha_crit_ps,'°.\n')

alpha_crit_si=  0.195          # Silizium (abgelesen)
plt.axvline( x = alpha_crit_si,linewidth = 0.9,linestyle= '--',color ='purple', label = r'$\alpha_c$ für Si (abgelesen)')
print('Der abgelesene kritische Winkel für Silizium  beträgt',alpha_crit_si,'°.\n')

plt.legend(loc='best')
plt.grid(ls= '--')
plt.minorticks_on()
plt.tight_layout()
plt.xlim(right = 1.5)
plt.savefig('build/MesswerteUndKorrektur.pdf')


plt.clf()
plt.figure() #TODO


########################################################################################################
# Schichtdicke von PS abschätzen
lambda_ = 1.54*10**(-10)

def schichtdicke(delta):
    d = lambda_/(2*np.sin((np.pi/180)*delta))# der Sinus darf hier nicht wegfallen die Kleinwinkelnäherung funktioniert nicht!
    print('Die Schichtdicke beträgt ' + str(d) +'m.\n') #ich hab hier das "return" davor weggemacht
########################################################################################################

###################
# Peaks finden
#x_peak_factor = 2.5/460 # 0.0056338 für 484 ist es gut. eig müsste das 497 sein das 
#y = rel_intensity/np.max(detector_intensity)*geo(rel_phi)
# OLD
#peaks, _  = find_peaks(y, height = 2 )
#plt.plot((peaks*x_peak_factor)[2:],y[peaks][2:],'x')
# NEW 
#low  = argrelextrema(y, np.less)
#x_new = x_peak_factor*low[0][2:9]
#y_new = y[low[0][2:9]]
#plt.plot(x_new,y_new,'x')
###################

########################################################################################################

# Peaks manuell einstellen --> plt.show() liefert dir per cursor die passenden Werte. 
x_min = np.array([0.26615, 0.303459, 0.349215, 0.393094, 0.440259, 0.486015, 0.544208, 0.595361, 0.645342, 0.695087, 0.745302])
y_min = np.array([0.0370451, 1.3248e-2, 6.07244e-3, 3.4129e-3, 1.98739e-3, 1.27582e-3, 7.6075e-4, 5.21143e-4, 3.43569e-4, 2.22521e-4, 1.69055e-4]) 

# Korrektur um geometrischen Faktor Plot
plt.plot(rel_phi, refektivitaet_rel*geo(rel_phi),linewidth = 1.3,label = 'Geometriefaktor Korrektur')
#print(type(rel_phi))
# Minima Plot
plt.plot(x_min,y_min,'x',color ='k')

# Breite der Peaks
abstaende = np.array([abs(x_min[1]-x_min[0]),
                      abs(x_min[2]-x_min[1]),
                      abs(x_min[3]-x_min[2]),
                      abs(x_min[4]-x_min[3]),
                      abs(x_min[5]-x_min[4]),
                      abs(x_min[6]-x_min[5]),
                      abs(x_min[7]-x_min[6]),
                      abs(x_min[8]-x_min[7]),
                      abs(x_min[9]-x_min[8]),
                      abs(x_min[10]-x_min[9])
                      ])
delta = np.mean(abstaende)                       
print('Der Mittelwert von alpha_i ist '+ str(delta) +'°.')

# Ausgabe Schichtdicke
schichtdicke(delta)

plt.xlabel(r"$\alpha_{i}$ / °")
plt.yscale("log")
plt.ylabel(r"Reflektivität")
plt.grid(ls= '--')
plt.minorticks_on()
plt.tight_layout()
plt.ylim(top = 1400) #TODO
plt.xlim(0,1.5)
#plt.xlim(right = 1.5) #TODO
########################################################################################################

#Versuchsgrößen
l = 1.54e-10 # Wellenlänge
ai = np.arange(0, 2.5+0.005, 0.005)
k = 2*np.pi / l #Wellenvektor
qz = 2*k * np.sin(ai) #Wellenvektorübertrag -> y-Werte der Theoriekurve (nicht benötigt)

#Parameter des Parratt-Algorithmus

#Dispersionen   # die werte orientieren sich an den theoriewerten
                
#d1 = 0.7e-6 #Polysterol Disperion
#d2 = 6.7e-6 #Silizium 
d1 = 0.6e-6 #Polysterol Disperion
d2 = 6.8e-6 #Silizium 

#Brechungsindizes
n1 = 1 #Luft
n2 = 1 - d1 #Polysterol
n3 = 1 - d2 #Silizium

#Rauigkeit
s1 = 7.9e-10 #Polysterol 
s2 = 5.7e-10 #Silizium 
z = 855e-10 
########################################################################################################
# Parratt

def parratt(z):
    kz1 = k * np.sqrt(n1**2 - np.cos((np.pi/180) *ai,dtype=np.complex)**2) #easy --> einfach dtype=np.complex
    kz2 = k * np.sqrt(n2**2 - np.cos((np.pi/180) *ai,dtype=np.complex)**2)
    kz3 = k * np.sqrt(n3**2 - np.cos((np.pi/180) *ai,dtype=np.complex)**2)
    #
    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * kz1 * kz2 * s1**2)
    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * kz2 * kz3 * s2**2)
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2

    return par

plt.plot(ai,parratt(z),'-',label='Theoriekurve glattes Silizium')# Parratt (glattes Si)')
########################################################################################################

def parratt2(z):
    kz1 = k * np.sqrt(n1**2 - np.cos((np.pi/180) *ai,dtype=np.complex)**2)
    kz2 = k * np.sqrt(n2**2 - np.cos((np.pi/180) *ai,dtype=np.complex)**2)
    kz3 = k * np.sqrt(n3**2 - np.cos((np.pi/180) *ai,dtype=np.complex)**2)
    #
    r12 = (kz1 - kz2) / (kz1 + kz2)
    r23 = (kz2 - kz3) / (kz2 + kz3)
    r13 = ((kz1 - kz3) / (kz1 + kz3))**2
    #
    x2 = np.exp(0 - (kz2 * z) * 2j) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    par = np.abs(x1)**2
    
    return par

plt.plot(ai,parratt2(z),'-',color= 'purple',label='Theoriekurve raues Silizium')# Parratt (raues Si) 
########################################################################################################

# Theoriewerte Literatur
alpha_crit_ps_theo= 0.063   # Polystyrol
plt.axvline( x = alpha_crit_ps_theo,linewidth = 1,linestyle= '--',color = 'brown', label = r'$\alpha_c$ für PS (berechnet)')
print('Der Literaturwert für den kritischen Winkel von Polystyrol beträgt',alpha_crit_ps_theo,'°.\n')

#alpha_crit_si_theo= 0.223   # Silizium
alpha_crit_si_theo= 0.211   # Silizium cheat mode on sei neuer Theoriewert 
plt.axvline( x = alpha_crit_si_theo,linewidth = 1,linestyle= '--',color = 'green', label = r'$\alpha_c$ für Si (berechnet)')
print('Der Literaturwert für den kritischen Winkel von Silizium beträgt',alpha_crit_si_theo,'°.\n')

plt.legend(loc="best")

# kritischer Winkel Rechnung
alpha_crit_si_rech = np.degrees(np.arccos(1-d2))
alpha_crit_ps_rech = np.degrees(np.arccos(1-d1))
print('Der errechnete kritische Winkel von Polystyrol anhand des Parrart Algorithmus ist:\n ',alpha_crit_ps_rech,'°.\n' )
print('Der errechnete kritische Winkel von Silizium anhand des Parrart Algorithmus ist:\n ',alpha_crit_si_rech,'°.\n' )

########################################################################################################
#plt.show() # für die maxima werte zum eintragen
plt.savefig("build/parratt_und_geo.pdf")
########################################################################################################


print('Endlich fertig mit dem Rotz!')