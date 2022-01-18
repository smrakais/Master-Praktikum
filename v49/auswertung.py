import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.signal as scs
# import uncertainties.unumpy as unp
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
plt.xlabel(r"$t$ / s")
plt.ylabel(r"$M_y(t)$ / V")
plt.savefig("build/MG.pdf")
plt.clf()

plt.plot(carr_purcell['x-axis'][1:], carr_purcell['1'][1:])
plt.xlabel(r"$t$ / s")
plt.ylabel(r"$M_y(t)$ / V")
plt.savefig("build/carr_purcell.pdf")
plt.clf()

meiboom_gill = meiboom_gill[4:]
carr_purcell = carr_purcell[4:]

# print(carr_purcell.head())
# print(meiboom_gill.head())


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
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("build/peaks.pdf")
plt.clf()

print("M_0 T_1 in mV: ", M0_T1)
print("M_0 T2 in mV: ", M0_T2*10**(3))
#   ----------------------------------------------------------------------------
## Viskosität und Diffusion bestimmen
#tau_d, peak_d, x_1, x_2 = np.genfromtxt("Auswertung/Daten/d.txt", unpack=True)
#tau_d = tau_d * 10**(-3)    # in Sekunden
#peak_d = peak_d * 10**(-3)  # Volt
#x_1 = x_1 * 10**(-3)    # Sekunden
#x_2 = x_2 * 10**(-3)    # Sekunden
## Halbwertsbreite
#halbwertsbreite = np.mean(np.abs(x_1-x_2))
#print("Halbwertsbreite: ", round(halbwertsbreite*10**6, 2),
#      round(sem(np.abs(x_1-x_2))*10**6, 2))
#
## Tabellen
#np.savetxt('Auswertung/Tabellen/t12.txt',
#           np.column_stack([x_1*10**3, x_2*10**3, np.abs(x_1-x_2)*10**3]),
#           delimiter=' & ', newline=r' \\'+'\n', fmt="%.3f")
#np.savetxt('Auswertung/Tabellen/D.txt',
#           np.column_stack([tau_d * 10**3, peak_d*10**3]),
#           delimiter=' & ', newline=r' \\'+'\n', fmt="%.1f")
#
## Größen
#gyro = 42.576*10**6*2*np.pi     # Hz pro Tesla
#d = 4.4*10**(-3)    # Probenduchmesser in meter
#G = 2.2*4/(d * gyro * halbwertsbreite)
#print("G: ", G)
#
#
#def f_d(x, M0, D, A):
#    return M0 * np.exp(-x/T2.n) * np.exp(-D*gyro**2*G**2*x**3/12) + A
#
#
#def log(x, M0, D):
#    return np.log(M0) - x/T2.n - (D*gyro**2*G**2/12)*x**3
#
#
## params_d, cov_d = curve_fit(log, 2*tau_d, np.log(-peak_d),
##                             p0=[0.5, 1.6*10**(-9)])
#params_d, cov_d = curve_fit(f_d, 2*tau_d, -peak_d,
#                            p0=[0.5, 1.6*10**(-9), 1])
#error_d = np.sqrt(np.diag(cov_d))
#M0_d = ufloat(params_d[0], error_d[0])
#D = ufloat(params_d[1], error_d[1])
#M_1 = ufloat(params_d[2], error_d[2])
#print("-------------------------------------------------------------------")
#print("M0 aus Viskosität: ", M0_d*10**(3))
#print("Diffusionskoeffizient: ", D)
#print("y-Achsenabschnitt: ", M_1)
#
#x = np.linspace(np.min(2*tau_d)-0.001, np.max(2*tau_d)+0.001, 1000)
## plt.plot(2*tau_d, np.log(-peak_d), 'r.', label="Daten")
## plt.plot(x, log(x, *params_d), 'b--', label="Fit")
#plt.plot(2*tau_d, -peak_d, 'r.', label="Daten")
#plt.plot(x, f_d(x, *params_d), 'b--', label="Fit")
#plt.xlabel(r"$2\tau$ / s")
#plt.ylabel(r"$\ln(M_y (t))$ / ln(V)")
#plt.legend(loc='best')
#plt.tight_layout()
#plt.savefig("Auswertung/Plots/Diffusion.pdf")
#plt.clf()
#
#
## Viskosität einlesen
#t_visko, delta = np.genfromtxt('Auswertung/Daten/viskosität.txt', unpack=True)
#delta_t = 15*60+35+42/60
#
#np.savetxt('Auswertung/Tabellen/viskositaet.txt',
#           np.column_stack([t_visko, delta]),
#           delimiter=' & ', newline=r' \\'+'\n', fmt="%.1f")
#
#
#def linear(x, a, b):
#    return a*x + b
#
#
#params_visko, cov_visko = curve_fit(linear, t_visko[5:9], delta[5:9])
#errors_visko = np.sqrt(np.diag(cov_visko))
#
#x = np.linspace(700, 1040, 1000)
#plt.plot(t_visko[5:9], delta[5:9], 'r.')
#plt.plot(x, linear(x, *params_visko), 'b--')
#
#delta_1 = linear(delta_t, *params_visko)
#print("Delta: ", delta_1)
#eta = 998.2*1.024*10**(-9)*(delta_t - delta_1)
#print("Eta: ", eta)
#
#r = (constants.k * 293.15)/(6*np.pi*D*eta)
#print("---------------------------------------------------------------------")
#print("Molekülradius :", r)
#
## -----------------------------------------------------------------------------
## Theoretische Vergleichwerte
#print("Hexagonaler Wert: ",
#      (28.89*10**(-27)*0.74/(4/3 * np.pi * 998.2))**(1/3))
#print("Hexagonal :",
#      (1-r/(28.89*10**(-27)*0.74/(4/3 * np.pi * 998.2))**(1/3))*100)
#print("Van-der-Waal-Wert: ",
#      (3*constants.k*647.05/(128*np.pi*22.04*10**6))**(1/3))
#print("Van-der-Waal: ",
#      (1-r/(3*constants.k*647.05/(128*np.pi*22.04*10**6))**(1/3))*100)
#

# meiboom_gill = pd.read_csv("Auswertung/Daten/pc_t2_csv.csv")
#
# print(meiboom_gill["x-axis"][1:])
# print(meiboom_gill["1"][1:])
# plt.plot(meiboom_gill["x-axis"], meiboom_gill["1"])
# plt.savefig("Auswertung/Plots/Test_2.pdf")
# plt.clf()