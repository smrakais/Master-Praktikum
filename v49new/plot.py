import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from scipy import stats
from uncertainties import ufloat
from scipy.optimize import curve_fit
from something import some 
#Generate data 
p1, pulse1, max1= some.neueWerte(file_name="data/dataa.txt", finished_file="build/taba.tex",  vars_name=[r"$p / \si{\milli\bar}$", r"$\frac{N}{\SI{120}{\second}}$", r"$\text{Maximum Position}$"], label_text="taba", caption_text=r"Die Werte für den Druck in dem Glaszylinder, die Anzahl der Pulse und die Position des Maximums bei einem Abstand von $d_1 = \SI{2.7}{\centi\meter}$." , precision=1)
p2, pulse2, max2= some.neueWerte(file_name="data/datab.txt", finished_file="build/tabb.tex",  vars_name=[r"$p / \si{\milli\bar}$", r"$\frac{N}{\SI{120}{\second}}$", r"$\text{Maximum Position}$"], label_text="tabb", caption_text=r"Die Werte für den Druck in dem Glaszylinder, die Anzahl der Pulse und die Position des Maximums bei einem Abstand $d = \SI{2}{\centi\meter}$." , precision=1)
twelve, pulse3 = some.neueWerte(file_name="data/datac.txt", finished_file="build/tab100.tex",  vars_name=[r"$\text{Anzahl}$"], label_text="tab100", caption_text=r"Die Werte für den Druck in dem Glaszylinder, die Anzahl der Pulse und die Position des Maximums." , precision=1)
#Generate table with calculated data
some.tabelle([pulse3[0:19],pulse3[20:39],pulse3[40:59],pulse3[60:79],pulse3[80:99]], finished_file="build/tabc.tex", vars_name=[r"Pulse", r"Pulse", r"Pulse", r"Pulse", r"Pulse"], label_text="tabc", caption_text=r"Die Pulse wurden zur Analyse der Statistik des Radioaktiven Zerfalls bestimmt.", precision=1) 

#extra values
p0 = 1013 
x01 = 0.027
x02 = 0.02
#functions 

x1 =  x01* p1/p0
x2 =  x02* p2/p0

max1 = max1/1120 * 4e6
max2 = max2/1099 * 4e6

Steigung1, yAbschnitt1, err1 = some.linReg(x=x1*1e3, y=pulse1, p=x1[16:19]*1e3, q=pulse1[16:19], x_name=r"$x_1 / \si{\milli\meter}$", y_name=r"$\frac{N}{\SI{120}{\second}}$", num=1,  x_add=1, file_name="build/plota.pdf")
Steigung1b, yAbschnitt1b, err1b = some.linReg(x=x1*1e3, y=max1*1e-6, p=x1[:17]*1e3, q=max1[:17]*1e-6, x_name=r"$x_1 / \si{\milli\meter}$", y_name=r"$E / \si{\mega\electronvolt}$", num=2,  x_add=1, file_name="build/plotb.pdf")

steigungerr1 = ufloat(Steigung1, err1)*1e3
mitReichw1 = (1/2 * pulse1[0] - yAbschnitt1)/steigungerr1
EnergiemitReichw1 = (mitReichw1*1e3/3.1)**(2/3)

steig1 = ufloat(Steigung1, err1)
E_R1 = mitReichw1*1000*Steigung1+yAbschnitt1

Steigung2, yAbschnitt2, err2 = some.linReg(x=x2*1e3, y=pulse2, p=x2[5:18]*1e3, q=pulse2[5:18], x_name=r"$x_2 / \si{\milli\meter}$", y_name=r"$\frac{N}{\SI{120}{\second}}$", num=3,  x_add=1, file_name="build/plotc.pdf")
Steigung2b, yAbschnitt2b, err2b = some.linReg(x=x2*1e3, y=max2*1e-6, p=x2[:]*1e3, q=max2[:]*1e-6, x_name=r"$x_2 / \si{\milli\meter}$", y_name=r"$E / \si{\mega\electronvolt}$", num=4,  x_add=1, file_name="build/plotd.pdf")

steigungerr2 = ufloat(Steigung2, err2)*1e3
mitReichw2 = (1/2 * pulse2[0] - yAbschnitt2)/steigungerr2
EnergiemitReichw2 = (mitReichw2*1e3/3.1)**(2/3)

steig2 = ufloat(Steigung2, err2)
E_R2 = mitReichw2*1000*Steigung2+yAbschnitt2

some.tabelle([x1*1e3, pulse1, max1*1e-6], finished_file="build/tab1.tex", vars_name=[r"$x_1 / \si{\milli\meter}$", r"$\frac{N}{\SI{120}{\second}}$", r"$E /\si{\mega\electronvolt} $"], label_text="tab1", caption_text=r"Die Reichweite $x_1$, die Anzahl der Impulse und die Position des Maximums.", precision=2) 
some.tabelle([x2*1e3, pulse2, max2*1e-6], finished_file="build/tab2.tex", vars_name=[r"$x_2 / \si{\milli\meter}$", r"$\frac{N}{\SI{120}{\second}}$", r"$E / \si{\mega\electronvolt}$"], label_text="tab2", caption_text=r"Die Reichweite $x_2$, die Anzahl der Impulse und die Position des Maximums.", precision=2) 
#some.tabelle(, finished_file="tab<++>.tex", vars_name=[r"<++>", r"<++>"], label_text="tab<++>", caption_text=r"<++>", precision=2) 
#Generate linReg-Plot
#some.linReg(x=<++>, y=<++>, x_name=r"<++>", y_name=r"<++>", num=<++>,  x_add=<++>, file_name="build/plot<++>.pdf")
#Generate curve-fit-Plot 
#some.curvefit(x=<++>, y=<++>, num=<++>, x_add=<++>, function=<++>, x_name=r"<++>", y_name=r"<++>", file_name="build/plot<++>.pdf")

#def gauß(x, mean, std):
#    return 1/np.sqrt(4*np.pi*std**2) *np.exp(-(x-mean)**2/(2*std**2)) 
pulse3 = (pulse3- np.min(pulse3))/100
mean = pulse3.mean()
std = pulse3.std()
np.random.seed(42)
gauß = np.random.normal(mean, std, 10000)
poisson = np.random.poisson(mean, 10000)

plt.figure(5) 
plt.hist(pulse3, bins=20, label="Daten", density=True)
plt.hist(gauß, bins=20, label="Gauß", color="green", density= True, histtype='step')
plt.ylabel("Wahrscheinlichkeit")
plt.xlabel("Pulse")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("build/plotf.pdf") 
plt.figure(6) 
plt.hist(pulse3, bins=20, label="Daten", density=True)
plt.hist(poisson, bins=20, label="Poisson", density = True, histtype='step')
plt.ylabel("Wahrscheinlichkeit")
plt.xlabel("Pulse")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("build/plote.pdf") 


#save solution
file = open("build/solution.txt", "w")
file.write(f"V701\n\nWerte aus Fit 1\nSteigung 1 = {Steigung1}+-{err1}\nyAbschnitt1 = {yAbschnitt1}\nmittlere Reichweite 1= {mitReichw1}\nEnergie bei mittlerer Reichweite = {EnergiemitReichw1}\ndE/dx = {Steigung1b}+-{err1b}\nyAbschnitt1b = {yAbschnitt1b}\n\nWerte aus Fit 2\nSteigung 2 = {Steigung2}+-{err2}\nyAbschnitt2 = {yAbschnitt2}\nmittlere Reichweite 2= {mitReichw2}\nEnergie bei mittlerer Reichweite = {EnergiemitReichw2}\ndE/dx = {Steigung2b}+-{err2b}\nyAbschnitt2b = {yAbschnitt2b}\n\nZählrate:\nMittelwert = {mean}\nStandardabweichung={std}\nVarianz=Std^2={std**2}\nSeeed = 42\n\n\nEnergie1: {E_R1}MeV\nEnergie2: {E_R2} MeV")
file.close()

