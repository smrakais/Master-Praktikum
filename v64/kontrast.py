import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from texutils.table import TexTable


#theta_P ist der polarisationswinkel
theta_P, U_max, U_min = np.genfromtxt('data/kontrast.txt', unpack=True)

#damit in radiant
theta_P = theta_P / 360 * (2*np.pi) 

Z = U_max-U_min
N = U_max+U_min


def Fitf(theta, a, b, c, d):
    return np.abs(a*np.sin(b*theta + c)) + d


Kontr = Z/N

params, cov = curve_fit(Fitf, theta_P, Kontr, )
errors = np.sqrt(np.diag(cov))
a = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
c = ufloat(params[2], errors[2])
d = ufloat(params[3], errors[3])


x = np.linspace(-0.4, 3.5, 1000)

theta = (np.pi/2 - c)/b

print('jetzt kommen die fitparameter')
print(a)
print(b)
print(c)
print(d,'\n')

print('---------------')
print(theta * 360 / (2 * np.pi))
print('das ist der beste winkel mit dem größsten kontrast')
print('den wir im experiment haben. wir haben 130° genommen.')
print('da wir das am anfang nicht ausgerechnet haben')

plt.plot(theta_P, Kontr, 'r+', label="Daten")
plt.plot(x, Fitf(x, *params), 'b', label="Regression")
plt.xlabel(r"$\theta_P \, / \, \mathrm{rad}$")
plt.ylabel('K')
plt.xticks([0, 0.5*np.pi, np.pi], ['0', r'$\frac{\pi}{2}$', r'$\pi$'])
plt.xlim(-0.4, 3.5)
plt.ylim(0, 1)
plt.tight_layout()
plt.legend(loc="best")
plt.savefig("build/Kontrast.pdf")
plt.clf()


##########
# Tabelle
##########
winkel, umax ,umin, kontrast = np.genfromtxt('data/gesicherte_messwerte/kontrast_(Kopie).txt',unpack = True)

print('-------------')
t = TexTable([winkel, umax,umin,kontrast], [r"Winkel / $°$", r"$U_\text{max}$ / mV",r"$U_\text{min}$/ mv",r"Kontrast"], 
            label='tab:Kontrast',
            caption='Messwerte der unterschiedlichen Spannungen.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 0)
t.set_row_rounding(2, 0)
t.set_row_rounding(3, 4)

t.write_file('build/tabKontrast.tex')
print('Die Tabelle des Kontrastes wurde erzeugt!\n')