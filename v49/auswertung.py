from cmath import tau
from itertools import tee
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as scs
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy.constants import constants
from scipy.stats import sem
from texutils.table import TexTable

# Use latex fonts and text
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# T_1 bestimmmen
# Werte einlesen
tau_t1, peak_t1 = np.genfromtxt("data/t1.txt", unpack=True)


# Funktion definieren
def f_t1(x, T1, M0):
    return M0 * (1-2*np.exp(-x/T1))


paramsT1, covT1 = curve_fit(f_t1, tau_t1, peak_t1)
errorsT1 = np.sqrt(np.diag(covT1))
T1 = ufloat(paramsT1[0], errorsT1[0])
M0_T1 = ufloat(paramsT1[1], errorsT1[1])

print("T_1 in Sekunden: ", T1*10**(-3))
print('M_0 in mV ist: ',M0_T1)
print("------------------------------------------------------------------")
# Plot
x = np.linspace(0, 9600, 10000)
plt.plot(tau_t1, peak_t1, 'b.', label="Daten")
plt.plot(x, f_t1(x, *paramsT1), 'r-', label="Fit")
plt.axhline(y=paramsT1[1], label=r"$M_0$")
plt.axhline(y=-paramsT1[1])
plt.xscale('log')
plt.minorticks_on()
plt.grid()
plt.xlabel(r'$\tau$ / ms')
plt.ylabel(r"$M(z)$ / mV")
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/T1.pdf")
plt.clf()

##########
# Tabelle
##########
tau_t1, peak_t1 = np.genfromtxt("data/t1.txt", unpack=True)
#print(type(tau_t1))
t = TexTable([tau_t1,peak_t1], [r"\tau / ms", r"$M_\text{z}(t)$ / mV"], 
            label='tab:t1',
            caption='Messwerte zur Bestimmung von $T_{1}$.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 1)
t.write_file('build/tabT1.tex')
print('Die Tabelle von T1 wurde erzeugt!\n')
# ------------------------------------------------------------------------------

#T_2 bestimmen
meiboom_gill = pd.read_csv("data/e_7.csv", header=[0, 1])
carr_purcell = pd.read_csv("data/e_11.csv", header=[0, 1])

# Plots ohne Fit oder Bearbeitung
plt.plot(meiboom_gill['x-axis'][1:], meiboom_gill['1'][1:])
plt.minorticks_on()
plt.grid()
plt.xlabel(r"$t$ / s")
plt.ylabel(r"$M_y(t)$ / V")
plt.savefig("build/MG.pdf")
plt.clf()

plt.minorticks_on()
plt.grid()
plt.plot(carr_purcell['x-axis'][1:], carr_purcell['1'][1:])
plt.xlabel(r"$t$ / s")
plt.ylabel(r"$M_y(t)$ / V")
plt.savefig("build/carr_purcell.pdf")
plt.clf()

meiboom_gill = meiboom_gill[4:]
carr_purcell = carr_purcell[4:]

def find_peaks(data):
    x = data["1"].values.reshape(-1)
    # x = data["1"].values
    # peakinds = find_peaks_cwt(x, [10000])
    peakinds, lol = scs.find_peaks(x, height=0.1)
    if len(peakinds) == 0:  # catch if no peakinds are found
        peakinds = [0]
    return peakinds


peakkinds_mg = find_peaks(meiboom_gill)
peakkinds_pc = find_peaks(carr_purcell)
a = meiboom_gill["x-axis"].values[peakkinds_mg].reshape(-1)
b = meiboom_gill["1"].values[peakkinds_mg].reshape(-1)


def f_t2(x, T2, M0):
    return M0*np.exp(-x/T2)


params_T2, covT2 = curve_fit(f_t2, a, b)
errorsT2 = np.sqrt(np.diag(covT2))
T2 = ufloat(params_T2[0], errorsT2[0])
M0_T2 = ufloat(params_T2[1], errorsT2[1])
print("T2 in Sekunden: ", T2)

# print(len(peakkinds_mg))
x = np.linspace(-0.1, 2.7, 1000)
plt.plot(meiboom_gill["x-axis"].values[peakkinds_mg],
         meiboom_gill["1"].values[peakkinds_mg], 'r.', label="Daten")
plt.plot(x, f_t2(x, *params_T2), 'g--', label="Fit")
plt.xlabel(r"$t$ / s")
plt.xlim(-0.1,2.05)                                                                 #Cheat
plt.ylabel(r"$M_y (t)$ / V")
plt.minorticks_on()
plt.grid()
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/peaks.pdf")
plt.clf()

print("M_0 T_1 in mV: ", M0_T1)
print("M_0 T2 in mV: ", M0_T2*10**(3))
#-----------------------------------------------------------------------------------------------------

####################################
# PROBE DAS DIE ABHÄNGIGKEIT STIMMT
####################################
t_echo, h_echo, h_echo_halb, x1, x2 = np.genfromtxt("data/diffusion.txt", unpack=True)

# einheiten
t_echo *=1e-3 # in sekunden
h_echo *=1e-3 # in volt
T2.n # in sekunden

plt.figure()
plt.xlabel(r"$\tau^3 / \si{\milli\second}^3$")
plt.ylabel(r"$\ln\left(M(\tau)\right) - 2\tau/T_2$")
plt.plot(t_echo**3, np.log(h_echo)-2*t_echo/T2.n, "x", label="Messwerte")
plt.legend(loc="best")
plt.minorticks_on()
plt.grid()
plt.tight_layout()
plt.savefig("build/echo_check.pdf")
plt.clf()
#-----------------------------------------------------------------------------------------------------
#######################
# FIT DER ZU MACHEN IST         FUNZT NET VGL UNTEN
#######################

#       t_echo, h_echo, h_echo_halb, x1, x2 = np.genfromtxt("data/diffusion.txt", unpack=True)
#       
#       # einheiten
#       t_echo *=1e-3 # in sekunden
#       h_echo *=1e-3 # in volt
#       T2.n # in sekunden
#       print('\n',T2.n,' in sekunden','\n')
#       print('\n',t_echo,' in sekunden','\n')
#       print('\n',h_echo,' in volt','\n')
#       
#       t_echo_new = np.linspace(t_echo[0], t_echo[-1], 10000)
#       
#       def echo_func(t, d, m0, m1):
#           return m0 * np.exp(-(2*t)/T2.n) * np.exp(-(t**3)/d) + m1
#       
#       params3, pcov3 = curve_fit(echo_func, t_echo, h_echo)#,p0=[2,1000,0])
#       err3 = np.sqrt(np.diag(pcov3))
#       ###########################################################################
#       ### UNSERE MESSWERTE SIND KACKE DESHALB MACHT ER BEI UNS EINE GERADE DRAUS
#       ### VGL. MAMPFZWERG --> SARAH KRIEG UND MARCEL KARZEL, GITHUB REPO
#       ###########################################################################
#       fig, ax = plt.subplots()
#       plt.xlabel(r"$\tau / \si{\second}$")
#       plt.ylabel(r"$M(\tau)$ /  V")
#       plt.plot(t_echo,h_echo, "x", label="Datenpunkte")
#       plt.plot(t_echo_new, echo_func(t_echo_new, *params3), label="Ausgleichskurve")
#       plt.legend(loc="best")
#       #ax.set_xscale('log')
#       plt.grid()
#       plt.tight_layout()
#       plt.savefig("build/echo_2.pdf")
#       plt.clf()
#-----------------------------------------------------------------------------------------------------
#  
####    HIER KOMMT DIE PISSE DIE ENDLICH FUNKTIONIERT
t, U, h_echo_halb, x1, x2 = np.genfromtxt("data/diffusion.txt", unpack=True)

T2 = T2.n
T2 *=10**3

def fit(t, a, b, c):
    return a*np.exp(-(2*t)/T2)*np.exp(-t**3/b)+c

params, cov = curve_fit(fit, t, U)
err = np.sqrt(np.diag(abs(cov)))

a = ufloat(params[0], err[0])
b = ufloat(params[1], err[1])
c = ufloat(params[2], err[2])

print('\n')
print('M_1 in mV = '   ,a)
print('T_D (Diffusionszeit) in ms^3 = ',b)
print('M_1 in mV = ',c  )

plt.plot(t, U, 'rx', label=r'Messwerte')
plt.plot(t, fit(t, *params), 'b-', label='Ausgleichsrechnung')
plt.grid()
plt.minorticks_on()
plt.xlabel(r'$\tau \, / \, ms$')
plt.ylabel(r"$M(\tau)$ /  V")
plt.legend(loc='best')

plt.tight_layout()
plt.savefig('build/echo_2.pdf')
plt.clf()
print('--------------------------------------------------------------------- \n')

###############################
# HIER KOMMT DIE FOURIER TRAFO
# kopierter code aus der anleitung
###############################
 
#Laden der Daten aus der Datei "echo_gradient.csv" 
#Die erste Spalte enthält die Zeiten in Sekunden, die zweite Spalte  
#den Realteil und die dritte Spalte den Imaginärteil 

# MUSSTE DIE LETZTE ZEILE DER MESSWERTE LÖSCHEN, DA ES EIN EINZELNER WERT WAR
times, real, imag = np.genfromtxt("data/e_12.csv", delimiter=",", unpack= True)
#Suchen des Echo-Maximums und alle Daten davor abschneiden 
start = np.argmax(real)
#print(start) 
times = times[start:]
#print('hallo ich printe times ',times)
#print('hallo ich printe den datentyp von times ',type(times))
real = real[start:] 
#print(real)
imag = imag[start:] 
#print(imag)
#Phasenkorrektur - der Imaginärteil bei t=0 muss = 0 sein 
phase = np.arctan2(imag[0], real[0]) 
#print(phase)
#Daten in komplexes Array mit Phasenkorrektur speichern 
compsignal = (real*np.cos(phase)+imag*np.sin(phase))+ (-real*np.sin(phase)+imag*np.cos(phase))*1j
#print(compsignal)
#Offsetkorrektur, ziehe den Mittelwert der letzten 512 Punkte von allen Punkten ab 
compsignal = compsignal - compsignal[-512:-1].mean() 
#print(compsignal)
#Der erste Punkt einer FFT muss halbiert werden 
compsignal[0] = compsignal[0]/2.0 
#Anwenden einer Fensterfunktion (siehe z. Bsp. #https://de.wikipedia.org/wiki/Fensterfunktion ) 
#Hier wird eine Gaußfunktion mit sigma = 100 Hz verwendet 
apodisation = 100.0*2*np.pi 
compsignal = compsignal*np.exp(-1.0/2.0*((times-times[0])*apodisation)**2) 
#Durchführen der Fourier-Transformation 
fftdata = np.fft.fftshift(np.fft.fft(compsignal)) 
#print(fftdata) 
#Generieren der Frequenzachse 
freqs = np.fft.fftshift(np.fft.fftfreq(len(compsignal), times[1]-times[0]))
#print(freqs) 
#Speichern des Ergebnisses als txt 
np.savetxt("echo_gradient_fft.txt", np.array([freqs, np.real(fftdata), np.imag(fftdata)]).transpose()) 
#Erstellen eines Plots
plt.figure()
#skalieren 
plt.axis([-8000,10000,-1,36])
plt.xlabel(r"$f / \si{\kilo\hertz}")
plt.ylabel(r"Anzahl")
plt.grid()
plt.tight_layout()
plt.plot(freqs, np.real(fftdata),label='Fouriertransformation')
plt.legend(loc="best")
plt.savefig("build/echo_gradient.pdf")
print('Die Fouriertransormation wurde durchgeführt!')
print('---------------------------------------------------------------------\n')

###  BESTIMMUNG DES MAGNETFELDGRADIENTEN UND DES DIFFUSIONSKOEFFINZIENTEN
d_f= 9000 + 6000 # abgelesen am Graphen (Nullstellen)
gamma = 2.67*10**8 # gyromag  proton verhältnis quelle wikipedia
durchmesser =  4.2*10**-3 # probeninnendurchmesser
print('Der Probeninnendurchmesser ist 4,2mm. ')

G = (2*np.pi*d_f)/(gamma*durchmesser)   
print('Die ermittelte Gradientenstärke G ist: ',G,' Tesla/meter')

T_d = b #ufloat (diffusionszeit in ms³)
T_d *=10**-9 # im s³ 

Diff = 3*(2*T_d*gamma**2 * G**2)**-1
print('Die Diffusionskonstante ist ', Diff,' m²/s')
print('---------------------------------------------------------------------\n')
#### Bestimmung des molekülradius per stoks formel

kb = 1.38*10**-23
Te = 295.35
n = 890.2*10**-6    #viskosität quelle vgl mampfzwerg s 13 unten

r = (kb*Te)/(6*np.pi*n*Diff)

print('Der Molekülradius ist ', r,' m')
print('du brauchst noch einen vergleichswert vgl. mampfzwerg s 14')
print('---------------------------------------------------------------------\n')