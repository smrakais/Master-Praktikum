import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from texutils.table import TexTable
import uncertainties.unumpy as unp
from uncertainties import ufloat
####################################################
####Aufgabe bestimme die Lebensdauer eines Myons####
####################################################


#########################################################################################################
# Bestimmung der Halbwertsbreite der Koinzidenz
#########################################################################################################

# lin reg funktion
def geradengleichung(x,a,b):
    return a*x+b

# Daten einlesen
time, intensity = np.genfromtxt('data/plateau_daten.txt', unpack = True)
time = time[:-5] # die letzen 5 messwerte werden verworfen
intensity = intensity[:-5]

#create new fig
plt.figure()    

plt.xlabel(r"$t$ / ns")
plt.ylabel(r"Ereignisrate $(10s)$")
plt.grid(ls = '--')
plt.tight_layout()
plt.minorticks_on()

###################################
# plot für die ersten 11 Messwerte
##################################
time3 = time[:12]
intensity3 = intensity[:12]

params, cov = np.polyfit(time3,intensity3,deg=1,cov=True)
errors = np.sqrt(np.diag(cov))
#print('a = {:.2f} ± {:.2f}'.format(params[0], errors[0]))
#print('b = {:.2f} ± {:.2f}'.format(params[1], errors[1]))

time3_new=np.linspace(time3[0], time3[-1])
plt.plot(time3_new, geradengleichung(time3_new,*params),'b-',label = 'Regressionsgerade')

###################################
# plot für die ersten 11 Messwerte
###################################
time1 = time[15:]
intensity1 = intensity[15:]

params, cov = np.polyfit(time1,intensity1,deg=1,cov=True)
errors = np.sqrt(np.diag(cov))
#print('a = {:.2f} ± {:.2f}'.format(params[0], errors[0]))
#print('b = {:.2f} ± {:.2f}'.format(params[1], errors[1]))

plt.plot(time,intensity ,"x", label="Messwerte")
time1_new=np.linspace(time1[0], time1[-1])
plt.plot(time1_new, geradengleichung(time1_new,*params),'b-',label = 'Regressionsgerade')

###################################
# plot für die Messerte in der Mitte
###################################
plt.hlines(y=116, xmin=3, xmax=9, linewidth=1, color='r',label='Plateau')

########
# FWHM 
########
y = 58
plt.hlines(y, xmin=0.05, xmax=10.6, linewidth=1, color='g',label='Halbwertsbreite')
print('Die FWHM ist ',y,'ns groß.')
plt.legend(loc="best") 
plt.savefig('build/koinzidenz.pdf')
print('Der Plot der Koinzidenz wurde erzeugt!')

##########
# Tabelle
##########
t = TexTable([time, intensity], [r"$t$ / ns",r"Zählrate"], 
            label='tab:Koinzidenz',
            caption='Messwerte der Zählraten in Abhängigkeit der Verzögerungsleitungen.')
t.set_row_rounding(0, 1) #reihe und rundung
t.set_row_rounding(1, 0)
t.write_file('build/tabKoinzidenz.tex')
print('Die Tabelle der Koinzidenz wurde erzeugt!\n')
#########################################################################################################

###########################################################################################
###########                 KEINE AHNUNG OB DAS SINNVOLLER IST ?
# Gauss Fit
###########
#
#plt.figure()
## fitting function 
#def gauss(x, amp, mu, sigma):
#    return amp/(sigma*np.sqrt(2*np.pi)) * np.exp(-(x- mu)**2/(2*sigma**2))
#
## fit and errors
#time_params, time_pcov = curve_fit(gauss, time, intensity)
#time_err = np.sqrt(np.diag(time_pcov))
#print('Intensität_peak =',      time_params[0], '±', time_err[0])
#print('time =',  time_params[1], '±', time_err[1])
#print('sigma =',    time_params[2], '±', time_err[2],'\n')
#time_new = np.linspace(time[0], time[-1])
#
#
#plt.hlines(y=116/2, xmin=0.3, xmax=10.6, linewidth=1, color='g',label='Halbwertsbreite')
#plt.plot(time, intensity, "x", label="Messwerte")
#plt.plot(time_new, gauss(time_new, *time_params), label="Ausgleichskurve")
#plt.legend(loc="best") 
#plt.grid(ls = '--')
#plt.xlabel(r"$t$ / ns")
#plt.ylabel(r"Impulsrate $(10s)$")
#plt.tight_layout()
#plt.minorticks_on()
#plt.savefig('build/gauss.pdf')
###########################################################################################


#############################################
# Kalibration des MCA (Multichannel Analyzer)
#############################################

kanal, time  = np.genfromtxt('data/kalibration.txt',unpack = True)

#create new fig
plt.figure()    

plt.xlabel(r"Kanal")
plt.ylabel(r"Pulsabstand / \si{\micro\second} ")
plt.grid(ls = '--')
plt.tight_layout()
plt.minorticks_on()

params, cov = np.polyfit(kanal,time,deg=1,cov=True)
errors = np.sqrt(np.diag(cov))
print('Die Regression liefert die Daten:')
print('m = {:.6f} ± {:.6f}'.format(params[0], errors[0]))
print('b = {:.6f} ± {:.6f}'.format(params[1], errors[1]),'\n')

plt.plot(kanal,time,'x', label='Messwerte',color ='r')
kanal_new=np.linspace(kanal[0], kanal[-1])

plt.plot(kanal_new, geradengleichung(kanal_new,*params),'b-',label = 'Regressionsgerade')
plt.legend(loc="best") 

plt.savefig('build/kalibration.pdf')
print('Der Plot der Kalibration wurde erzeugt!')

##########
# Tabelle
##########
t = TexTable([time, kanal], [r"Pulsabstand / \si{\micro\second}",r"Kanal"], 
            label='tab:Kalibration',
            caption='Messwerte der Kalibration des MCA.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 0)
t.write_file('build/tabKalibration.tex')
print('Die Tabelle der Kalibration wurde erzeugt!\n')
#########################################################################################################

###########################
#Bestimmung des Untergrunds
###########################

plt.figure()

anz_kanal= 8191-6453+1+50 #1789 nur die gefüllten channel, weil bei den leeren haben wir auch keinen Untergrund
messdauer = 256355
Startimpulse = 2912226
Startimpulse = ufloat(Startimpulse, np.sqrt(Startimpulse))
N_quer = Startimpulse/messdauer # N_quer
Ts = 10*10**(-6)    # Suchzeit (Sekunden)
Nf = N_quer*Ts*unp.exp(-N_quer*Ts)*Startimpulse
Nf_kanal = Nf/anz_kanal  # Untergrund pro Kanal
print("-------------------")
print(Startimpulse)
print("Fehlmessungen: ", Nf)
print("Untergrundrate: ", Nf_kanal)
print('Die Untergrundrate wurde bestimmt!')

y  = np.genfromtxt('data/ereignisse.Spe',unpack = True)
x=np.arange(0,8192)
plt.plot(x*params[0],y,'x')
plt.savefig('build/plot.pdf')
plt.clf()
#########################################################################################################

############################
#Bestimmung der Lebensdauer
############################

# Löscht alle Nullen aus den Daten raus
data = np.genfromtxt('data/ereignisse.Spe',unpack = True)
valueToBeRemoved = 0.0 # die Nachkommastelle ist wichtig
data = list(filter((valueToBeRemoved).__ne__, data))
#print(np.shape(data))
#print(np.shape(data)[0])

plt.figure()
#t = np.arange(0,np.shape(data)[0]) # Anzahl der befüllten Kanäle
t=np.arange(0,8192)                 # Anzahl aller Kanäle
t = t * params[0]
print(t)



data = np.genfromtxt('data/test.Spe',unpack = True)
t=np.arange(0,512)

plt.plot(t* params[0],data,'x')
plt.savefig('build/ereignisse.pdf')
