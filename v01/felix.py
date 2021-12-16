# Header
import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from texutils.table import TexTable

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18


def linear(x, m, b):
    return m*x + b


# Daten einlesen
dt, N = np.genfromtxt('data_new/plateau_daten.txt', unpack=True)
daten = np.genfromtxt('data_new/ereignisse.Spe', unpack=True)
print("Summe Stopp: ", ufloat(np.sum(daten), np.sum(np.sqrt(daten))))
dt = dt - 20    # Zentrieren
N *= 10
hoehe = np.mean(N[4:16])
links = N[0:4]
rechts = N[16:]
dt_rechts = dt[16:]

# linearer Fit links und rechts
params_links, cov_links = curve_fit(linear, dt[0:4], links)
errors_links = np.sqrt(np.diag(cov_links))
m = ufloat(params_links[0], errors_links[0])
b = ufloat(params_links[1], errors_links[1])
print("Steigung links: ", m)
print("y-Achsenabschnitt rechts: ", b)
print("Halbe Höhe: ", hoehe/2)

params_rechts, cov_rechts = curve_fit(linear, dt_rechts, rechts)
errors_rechts = np.sqrt(np.diag(cov_rechts))
m = ufloat(params_rechts[0], errors_rechts[0])
b = ufloat(params_rechts[1], errors_rechts[1])
print("Steigung rechts: ", m)
print("y-Achsenabschnitt rechts: ", b)

# Berechnen des Schnittpunktes
x_links = np.linspace(-25, -17)
x_rechts = np.linspace(6, 20)

links_w = linear(x_links, *params_links)
rechts_w = linear(x_rechts, *params_rechts)

plt.figure(1)
plt.ylabel(r"$N(t) \, / \, (10\mathrm{s})^{-1}$")
plt.xlabel(r"$\mathrm{d}t \, / \, \mathrm{ns}$")
plt.errorbar(dt, N, yerr=np.sqrt(N), fmt='kx', label="Messwerte")
plt.plot(x_links, linear(x_links, *params_links), 'r',
         label="Regression links")
plt.plot(x_rechts, linear(x_rechts, *params_rechts), 'r',
         label="Regression rechts")
# plt.plot(dt, N, 'r.', label="Messwerte")
plt.axhline(y=hoehe, xmin=0.15, xmax=0.75, label="Plateau")
plt.axhline(y=hoehe/2, xmin=0.11, xmax=0.87, color="green",
            label="Halbwertsbreite")
# plt.ylim(0, 130)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.minorticks_on()
plt.savefig("build_new/Plateau.pdf")
plt.clf()

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
################################################################################################

###################
#  Kanalauswertung
###################
kanal,kal_t, hits = np.genfromtxt("data_new/kalibration.txt", unpack=True)

params_kal, cov_kal = curve_fit(linear, kanal, kal_t)
errors_kal = np.sqrt(np.diag(cov_kal))
m = ufloat(params_kal[0], errors_kal[0])
b = ufloat(params_kal[1], errors_kal[1])
print("Steigung: ", m)
print("y-Achsenabschnitt: ", b)


# Plot dazu
x = np.linspace(0, 210)
plt.plot(kanal, kal_t, 'r+', label="Daten", markersize=25)
plt.plot(x, linear(x, *params_kal), 'b--', label="Regression")
plt.xlabel("Kanal")
plt.ylabel(r"$T_{VZ} \, / \, \mathrm{\mu s}$")
plt.xlim(0, 210)
plt.tight_layout()
plt.minorticks_on()
plt.grid()
plt.legend(loc="best")
plt.savefig("build_new/kal.pdf")
plt.clf()

##########
# Tabelle
##########
t = TexTable([kal_t, kanal], [r"Pulsabstand / \si{\micro\second}",r"Kanal"], 
            label='tab:Kalibration',
            caption='Messwerte der Kalibration des MCA.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 0)
t.write_file('build_new/tabKalibration_new.tex')
print('Die Tabelle der Kalibration wurde erzeugt!\n')
################################################################################################

# Bestimmung Untergrundrate
messdauer = 147182  # Sekunden
Startimpulse = 3061879
Startimpulse = ufloat(Startimpulse, np.sqrt(Startimpulse))
n = Startimpulse/messdauer
Ts = 20*10**(-6)    # Sekunden
Nf = Startimpulse*n*Ts*unp.exp(-n*Ts)
Nf_kanal = Nf/450
print("-------------------")
print(Startimpulse)
print("Fehlmessungen: ", Nf)
print("Untergrundrate: ", Nf_kanal)

# Umrechnung Kanäle in Zeit
kanaele = np.arange(0, 512, 1)
zeiten1 = linear(kanaele, *params_kal)
# Rausnehmen von komischen Werten
zeiten = zeiten1[3:6]
zeiten = np.append(zeiten, zeiten1[7:14])
zeiten = np.append(zeiten, zeiten1[15:455])

daten_ang = daten[3:6]
daten_ang = np.append(daten_ang, daten[7:14])
daten_ang = np.append(daten_ang, daten[15:455])


# Nullenverarbeitung
for i in range(len(daten_ang)):
    if daten_ang[i] == 0:
        a = np.array([daten_ang[i-1], daten_ang[i+1]])
        daten_ang[i] = np.round(np.mean(a))


# Entfernte Daten und Zeiten
zeiten_ent = zeiten1[:3]
zeiten_ent = np.append(zeiten_ent, zeiten1[6])
zeiten_ent = np.append(zeiten_ent, zeiten1[14])
daten_ent = daten[:3]
daten_ent = np.append(daten_ent, daten[6])
daten_ent = np.append(daten_ent, daten[14])

# Definition der exp-Funktion
def e(x, N_0, l, U):
    return N_0*np.exp(-l*x) + U

params_fit, cov_fit = curve_fit(e, zeiten, daten_ang,
                                sigma=1/np.sqrt(daten_ang))
errors_fit = np.sqrt(np.diag(cov_fit))
N_0 = ufloat(params_fit[0], errors_fit[0])
l = ufloat(params_fit[1], errors_fit[1])
U = ufloat(params_fit[2], errors_fit[2])

print("-------------------")
print("N_0: ", N_0)
print("Lambda: ", l)
print("Lebensdauer: ", 1/l)
print("Untergrundrate: ", U)
tau = 2.2
ta_fit = 1/l
print("Verhältnis: ", (1-ta_fit/tau)*100)
print("Verhältnis Untergrund: ", (1-U/Nf_kanal)*100)

plt.plot(zeiten, daten_ang, 'r.', label="Daten")
#plt.plot(zeiten_ent, daten_ent, 'gx', label="Nicht-betrachtete Daten")
plt.plot(zeiten, e(zeiten, *params_fit), 'b--', label="Fit")
plt.xlabel(r"$\tau \,  / \, \mathrm{\mu s}$")
plt.ylabel(r"$N(t)$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.minorticks_on()
plt.savefig('build_new/fit.pdf')
plt.clf()


print("Alles ohne Probleme ausgeführt!")