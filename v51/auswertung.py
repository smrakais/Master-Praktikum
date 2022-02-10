import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from texutils.table import TexTable
from texutils.table import Combined

def linear(x, m, b):
    return m*x + b

def exponent(x,c,b ):
    return c * np.exp(b*-x)

###########################################################
################# Plot 1 ## Verstärkung = 100 #############
###########################################################
f, ua, t1, t2 = np.genfromtxt('data/A1_linear100.txt', unpack = True) 

ue = 50e-3 # Eingangsspannung Volt

# Phase berechnen
dt= (t2-t1)*1e-6 # in sekunden
phase_rad = 2 * np.pi *f*dt # radiant
phase = np.degrees(phase_rad) # degrees

# verstärkung
gain = ua/ ue

# plot
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'$log(f)$ / log(Hz)')
ax1.set_ylabel(r'$log(V)$/ Verstärkung')
ax1.grid(ls =  '--')
ax1.minorticks_on()
#ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.plot(np.log(f),np.log(gain),'x',label= 'Messwerte')

#fit
params, covariance_matrix = curve_fit(linear, np.log(f[0:6]), np.log(gain[0:6]))
errors = np.sqrt(np.diag(covariance_matrix))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
#print(type(b))
#print(b)
Verstärkungsfaktor_100 = unp.exp(b)
#print(type(unp.exp(b)))
#print(np.shape(Verstärkungsfaktor_100))
#print(type(float(unp.nominal_values(Verstärkungsfaktor_100)))) #convert to float
print('\n')
print("Steigung Plateau: ", m)
print("y-Achsenabschnitt Plateau: ", b)
print('Die Verstäkung ist ',Verstärkungsfaktor_100)
f_new =f[0:6]
f_new = np.linspace(f_new[0], f_new[-1], 100)
ax1.plot(np.log(f_new), linear(np.log(f_new), *params), 'r', label="lineare Regression")

#rausgenommen
plt.plot(np.log(f[6:9]),np.log(gain[6:9]),'x',color= 'black',label= 'rausgenommen')

##fit 2
params2, covariance_matrix2 = curve_fit(linear, np.log(f[9:]), np.log(gain[9:]))
errors2 = np.sqrt(np.diag(covariance_matrix2))
m = ufloat(params2[0], errors2[0])
b = ufloat(params2[1], errors2[1])
print('\n')
print("Steigung am Abfall: ", m)
print("y-Achsenabschnitt am Abfall: ", b,'\n')
f_new2 =f[9:]
f_new2= np.linspace(f_new2[0],f_new2[-1],100)
plt.plot(np.log(f_new2),linear(np.log(f_new2),*params2),'b',label= 'lineare Regression')
plt.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/a.pdf')
plt.clf()

#save 
array_log_verstärkung  = linear(np.log(f_new2),*params2)
array_log_frequenz = np.log(f_new2)

# Brechnungen
v_gr= unp.log((unp.exp(ufloat(params[1], errors[1]))/np.sqrt(2))) # logarithmischer wert der verstärkung bei der grenzfrequenz
print('log(V/sqrt(2)), Grenzverstärkung = ',v_gr) #np.log(x) = ln(x)

#f_gr Berechnung
f_gr_100 =unp.exp((v_gr-ufloat(params2[1], errors2[1]))/ufloat(params2[0], errors2[0]))
print('Die logarithmierte Grenzfrequenz ist ',unp.log(f_gr_100), 'log(Hz)\n')
print('Die Grenzfrequenz ist ',f_gr_100,'Hz')

# Berechnung Bandbreitenprodukt (BWP)
BWP_100 = Verstärkungsfaktor_100 * f_gr_100
print('Das Bandbreitenprodukt ist ',BWP_100,'Hz')

##########
# Tabelle
##########
t = TexTable([f, ua, phase], [r"$f$ / Hz",r"$U_a$ / V",r"$\phi$ / °"], 
            label='tab:Linearverstärker_100 ',
            caption='Messwerte des Linearverstärkers bei einem Verstärkungsfaktor von 100.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 2)
t.set_row_rounding(2, 0)

t.write_file('build/Linearverstärker_100.tex')
print('Die Tabelle des Linearverstärkers bei einem Verstärkungsfaktor von 100 wurde erzeugt!\n')
print('--------------------done-------------------------')
#############################################################################################################

###########################################################
################# Plot 2 ## Verstärkung = 10 ##############
###########################################################
f, ua, t1, t2 = np.genfromtxt('data/A2_linear10.txt', unpack = True) 

ue = 50e-3 # Eingangsspannung Volt

# Phase berechnen
dt= (t2-t1)*1e-6 # in sekunden
phase_rad = 2 * np.pi *f*dt # radiant
phase = np.degrees(phase_rad) # degrees

# verstärkung
gain = ua/ ue

# plot
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'$log(f)$ / log(Hz)')
ax1.set_ylabel(r'$log(V)$/ Verstärkung')
ax1.grid(ls =  '--')
ax1.minorticks_on()
#ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.plot(np.log(f),np.log(gain),'x',label= 'Messwerte')

#fit
params, covariance_matrix = curve_fit(linear, np.log(f[0:10]), np.log(gain[0:10]))
errors = np.sqrt(np.diag(covariance_matrix))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
Verstärkungsfaktor_10 = unp.exp(b)
print('\n')
print("Steigung Plateau: ", m)
print("y-Achsenabschnitt Plateau: ", b)
print('Die Verstäkung ist ',Verstärkungsfaktor_10)
f_new =f[0:10]
f_new = np.linspace(f_new[0], f_new[-1], 100)
ax1.plot(np.log(f_new), linear(np.log(f_new), *params), 'r', label="lineare Regression")

#rausgenommen
#plt.plot(np.log(f[6:9]),np.log(gain[6:9]),'x',color= 'black',label= 'rausgenommen')

##fit 2
params2, covariance_matrix2 = curve_fit(linear, np.log(f[10:]), np.log(gain[10:]))
errors2 = np.sqrt(np.diag(covariance_matrix2))
m = ufloat(params2[0], errors2[0])
b = ufloat(params2[1], errors2[1])
print('\n')
print("Steigung am Abfall: ", m)
print("y-Achsenabschnitt am Abfall: ", b,'\n')
f_new2 =f[10:]
f_new2= np.linspace(f_new2[0],f_new2[-1],100)
plt.plot(np.log(f_new2),linear(np.log(f_new2),*params2),'b',label= 'lineare Regression')
plt.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/b.pdf')
plt.clf()

#save 
array_log_verstärkung  = linear(np.log(f_new2),*params2)
array_log_frequenz = np.log(f_new2)

# Brechnungen
v_gr= unp.log((unp.exp(ufloat(params[1], errors[1]))/np.sqrt(2))) # logarithmischer wert der verstärkung bei der grenzfrequenz
print('log(V/sqrt(2)), Grenzverstärkung = ',v_gr) #np.log(x) = ln(x)

#f_gr Berechnung
f_gr_10 =unp.exp((v_gr-ufloat(params2[1], errors2[1]))/ufloat(params2[0], errors2[0]))
print('Die logarithmierte Grenzfrequenz ist ',unp.log(f_gr_10), 'log(Hz)\n')
print('Die Grenzfrequenz ist ',f_gr_10,'Hz')

# Berechnung Bandbreitenprodukt (BWP)
BWP_10 = Verstärkungsfaktor_10 * f_gr_10
print('Das Bandbreitenprodukt ist ',BWP_10,'Hz')

##########
# Tabelle
##########
t = TexTable([f, ua, phase], [r"$f$ / Hz",r"$U_a$ / V",r"$\phi$ / °"], 
            label='tab:Linearverstärker_10 ',
            caption='Messwerte des Linearverstärkers bei einem Verstärkungsfaktor von 10.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 2)
t.set_row_rounding(2, 0)

t.write_file('build/Linearverstärker_10.tex')
print('Die Tabelle des Linearverstärkers bei einem Verstärkungsfaktor von 10 wurde erzeugt!\n')
print('--------------------done-------------------------')
#########################################################################################################
###########################################################
################# Plot 3 ## Verstärkung = 1000 ############
###########################################################
f, ua, t1, t2 = np.genfromtxt('data/A3_linear1000.txt', unpack = True) 

ue = 50e-3 # Eingangsspannung Volt

# Phase berechnen
dt= (t2-t1)*1e-6 # in sekunden
phase_rad = 2 * np.pi *f*dt # radiant
phase = np.degrees(phase_rad) # degrees

# verstärkung
gain = ua/ ue

# plot
fig, ax1 = plt.subplots()
ax1.set_xlabel(r'$log(f)$ / log(Hz)')
ax1.set_ylabel(r'$log(V)$/ Verstärkung')
ax1.grid(ls =  '--')
ax1.minorticks_on()
#ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.plot(np.log(f),np.log(gain),'x',label= 'Messwerte')

#fit
params, covariance_matrix = curve_fit(linear, np.log(f[0:6]), np.log(gain[0:6]))
errors = np.sqrt(np.diag(covariance_matrix))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
Verstärkungsfaktor_1000 = unp.exp(b)
print('\n')
print("Steigung Plateau: ", m)
print("y-Achsenabschnitt Plateau: ", b)
print('Die Verstäkung ist ',Verstärkungsfaktor_1000)
f_new =f[0:6]
f_new = np.linspace(f_new[0], f_new[-1], 100)
ax1.plot(np.log(f_new), linear(np.log(f_new), *params), 'r', label="lineare Regression")

#rausgenommen
plt.plot(np.log(f[6:9]),np.log(gain[6:9]),'x',color= 'black',label= 'rausgenommen')

##fit 2
params2, covariance_matrix2 = curve_fit(linear, np.log(f[9:]), np.log(gain[9:]))
errors2 = np.sqrt(np.diag(covariance_matrix2))
m = ufloat(params2[0], errors2[0])
b = ufloat(params2[1], errors2[1])
print('\n')
print("Steigung am Abfall: ", m)
print("y-Achsenabschnitt am Abfall: ", b,'\n')
f_new2 =f[9:]
f_new2= np.linspace(f_new2[0],f_new2[-1],100)
plt.plot(np.log(f_new2),linear(np.log(f_new2),*params2),'b',label= 'lineare Regression')
plt.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/c.pdf')
plt.clf()

#save 
array_log_verstärkung  = linear(np.log(f_new2),*params2)
array_log_frequenz = np.log(f_new2)

# Brechnungen
v_gr= unp.log((unp.exp(ufloat(params[1], errors[1]))/np.sqrt(2))) # logarithmischer wert der verstärkung bei der grenzfrequenz
print('log(V/sqrt(2)), Grenzverstärkung = ',v_gr) #np.log(x) = ln(x)

#f_gr Berechnung
f_gr_1000 =unp.exp((v_gr-ufloat(params2[1], errors2[1]))/ufloat(params2[0], errors2[0]))
print('Die logarithmierte Grenzfrequenz ist ',unp.log(f_gr_1000), 'log(Hz)\n')
print('Die Grenzfrequenz ist ',f_gr_1000,'Hz')

# Berechnung Bandbreitenprodukt (BWP)
BWP_1000 = Verstärkungsfaktor_1000 * f_gr_1000
print('Das Bandbreitenprodukt ist ',BWP_1000,'Hz')

##########
# Tabelle
##########
t = TexTable([f, ua, phase], [r"$f$ / Hz",r"$U_a$ / V",r"$\phi$ / °"], 
            label='tab:Linearverstärker_1000 ',
            caption='Messwerte des Linearverstärkers bei einem Verstärkungsfaktor von 1000.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 1)
t.set_row_rounding(2, 0)

t.write_file('build/Linearverstärker_1000.tex')
print('Die Tabelle des Linearverstärkers bei einem Verstärkungsfaktor von 1000 wurde erzeugt!\n')
print('--------------------done-------------------------')
#########################################################################################################

##############################
# Tabelle GESAMMELTE PARAMETER
##############################

V_theo = np.array([100,10,1000])
V_berech = np.array([float(unp.nominal_values(Verstärkungsfaktor_100)),float(unp.nominal_values(Verstärkungsfaktor_10)),float(unp.nominal_values(Verstärkungsfaktor_1000))]) 
V_log_berech =np.log(V_berech) 
#check
print(V_theo,'\n', V_berech,'\n',V_log_berech,'\n')

f_gr = np.array([float(unp.nominal_values(f_gr_100)),float(unp.nominal_values(f_gr_10)),float(unp.nominal_values(f_gr_1000))])
f_gr_log = np.log(f_gr)
GBP = np.array([BWP_100.n,BWP_10.n,BWP_1000.n])
#check
print(f_gr,'\n',f_gr_log,'\n',GBP,'\n')

# Tabelle 
t = TexTable([V_theo, V_berech, V_log_berech,f_gr, f_gr_log, GBP], [r"$V_\text{theo}$ ",r"$V_\text{b}$ ",r"ln($V_\text{b}$) ",r"$f_\text{gr}$ / Hz ",r"ln($f_\text{gr}$) / ln(Hz)",r"$GBP (V_\text{b} \cdot f_\text{gr})$ / Hz"], 
            label='tab:Verstärkungen ',
            caption='Tabelle der berechneten Verstärkungen, Grenzfrequenzen und des Bandbreitenprodukts (GBP).'
            'bla')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 2)
t.set_row_rounding(2, 2)
t.set_row_rounding(3, 0)
t.set_row_rounding(4, 2)
t.set_row_rounding(5, 0)
t.write_file('build/tabParameter_Linearverstärker.tex')
print('Die Tabelle des Linearverstärkers bei einem Verstärkungsfaktor von 1000 wurde erzeugt!\n')

###########################
# Für kombinierte Tabellen
###########################
# Kombinierte Tabelle
# t=Combined([t1, t2], label= 'tab_1_2',caption='Kombinierte Tabelle der berechneten Verstärkungen, Grenzfrequenzen und des Bandbreitenprodukts (GBP)')
# t.write_file('build/tabKombiniert_1_2.tex')
###########################
#########################################################################################################