import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from texutils.table import TexTable
from texutils.table import Combined

def linear(x, m, b):
    return m*x + b

def log_linear(x,m,b ):
    return m * np.log(x) + b # --> log(y) = m * log(x) + b | du musst die y-werte im log plotten!

def fx(x,m,b):
    return np.exp(np.log(x)*m+b)

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

##############################################
# Phasenbeziehung in Abhängigkeit der Frequenz
##############################################
fig, ax = plt.subplots()
plt.ylabel(r"$\phi \, / \, \mathrm{rad}$")
plt.xlabel(r"$f \, / \, \mathrm{Hz}$")
ax.grid(ls =  '--')
ax.minorticks_on()
ax.set_xscale('log')
ax.plot(f,phase,'x', label = 'Messwerte' )
ax.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/Phasenbeziehung_100.pdf')

print('Der Plot der Phasenbeziehung <-> Frequenz wurde erstellt!')
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
##############################################
# Phasenbeziehung in Abhängigkeit der Frequenz
##############################################
fig, ax = plt.subplots()
plt.ylabel(r"$\phi \, / \, \mathrm{rad}$")
plt.xlabel(r"$f \, / \, \mathrm{Hz}$")
ax.grid(ls =  '--')
ax.minorticks_on()
ax.set_xscale('log')
ax.plot(f,phase,'x', label = 'Messwerte' )
ax.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/Phasenbeziehung_10.pdf')

print('Der Plot der Phasenbeziehung <-> Frequenz wurde erstellt!')
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
##############################################
# Phasenbeziehung in Abhängigkeit der Frequenz
##############################################
fig, ax = plt.subplots()
plt.ylabel(r"$\phi \, / \, \mathrm{rad}$")
plt.xlabel(r"$f \, / \, \mathrm{Hz}$")
ax.grid(ls =  '--')
ax.minorticks_on()
ax.set_xscale('log')
ax.plot(f,phase,'x', label = 'Messwerte' )
ax.plot(f[4],phase[4],'x',color='black', label = 'rausgenommen' )
ax.plot(f[8],phase[8],'x',color='black' )
ax.plot(f[13],phase[13],'x',color='black' )
ax.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/Phasenbeziehung_1000.pdf')

print('Der Plot der Phasenbeziehung <-> Frequenz wurde erstellt!')
print('--------------------done-------------------------')
#########################################################################################################

##############################
# Tabelle GESAMMELTE PARAMETER
##############################

V_theo = np.array([100,10,1000])
V_berech = np.array([float(unp.nominal_values(Verstärkungsfaktor_100)),float(unp.nominal_values(Verstärkungsfaktor_10)),float(unp.nominal_values(Verstärkungsfaktor_1000))]) 
V_log_berech =np.log(V_berech) 
#check
#print(V_theo,'\n', V_berech,'\n',V_log_berech,'\n')

f_gr = np.array([float(unp.nominal_values(f_gr_100)),float(unp.nominal_values(f_gr_10)),float(unp.nominal_values(f_gr_1000))])
f_gr_log = np.log(f_gr)
GBP = np.array([BWP_100.n,BWP_10.n,BWP_1000.n])
#check
#print(f_gr,'\n',f_gr_log,'\n',GBP,'\n')

# Tabelle 
t = TexTable([V_theo, V_berech, V_log_berech,f_gr, f_gr_log, GBP], [r"$V_\text{theo}$ ",r"$V_\text{b}$ ",r"ln($V_\text{b}$) ",r"$f_\text{gr}$ / Hz ",r"ln($f_\text{gr}$) / ln(Hz)",r"$GBP (V_\text{b} \cdot f_\text{gr})$ / Hz"], 
            label='tab:Verstärkungen ',
            caption='Tabelle der berechneten Verstärkungen, Grenzfrequenzen und des Bandbreitenprodukts (GBP).'
            'bla')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 0)
t.set_row_rounding(2, 2)
t.set_row_rounding(3, 0)
t.set_row_rounding(4, 2)
t.set_row_rounding(5, 0)
t.write_file('build/tabParameter_Linearverstärker.tex')
print('Die Tabelle der berechneten Paramter wurde erzeugt!\n')
print('--------------------done-------------------------')
###########################
# Für kombinierte Tabellen
###########################
# Kombinierte Tabelle
# t=Combined([t1, t2], label= 'tab_1_2',caption='Kombinierte Tabelle der berechneten Verstärkungen, Grenzfrequenzen und des Bandbreitenprodukts (GBP)')
# t.write_file('build/tabKombiniert_1_2.tex')
###########################
#########################################################################################################
#############################################
# Der invertierende Integrator mit Fremddaten
#############################################
ue, ua , f = np.genfromtxt('data_of_others_cuz_ours_suck/int/data_int.txt',unpack = True)
f *= 1e3 # in Hz

#plot
fig, ax = plt.subplots()
plt.ylabel(r"$U_a \, / \, \mathrm{V}$")
plt.xlabel(r"$f \, / \, \mathrm{Hz}$")
ax.grid(ls =  '--')
ax.minorticks_on()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(f,ua,'x', label = 'Messwerte' )

# fit 
params1, covariance_matrix1 = curve_fit(fx, f[:5], ua[:5])
errors1 = np.sqrt(np.diag(covariance_matrix1))
m = ufloat(params1[0], errors1[0])
b = ufloat(params1[1], errors1[1])
print("Steigung am Integrator 1: ", m)
print("y-Achsenabschnitt am Integrator 1 ", b,'\n')
plt.plot(f[:5],fx(f,*params1)[:5],'b',label= 'lineare Regression')

# fit teil 2
params2, covariance_matrix2 = curve_fit(fx, f[10:], ua[10:])
errors2 = np.sqrt(np.diag(covariance_matrix1))
m = ufloat(params2[0], errors2[0])
b = ufloat(params2[1], errors2[1])
print("Steigung am Integrator 2: ", m)
print("y-Achsenabschnitt am Integrator 2: ", b,'\n')
plt.plot(f[10:],fx(f,*params2)[10:],'r',label= 'lineare Regression 2' )

plt.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/integrator.pdf')
plt.clf()

print('Der Plot des Integrators wurde erstellt!','\n')
##########
# Tabelle
##########
t = TexTable([f, ua, ue], [r"$f$ / Hz",r"$U_a$ / V",r"$U_e$ / V"], 
            label='tab:integrator ',
            caption='Messwerte des invertierenden Integrators.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 1)
t.set_row_rounding(2, 1)

t.write_file('build/tab_integrator.tex')
print('Die Tabelle des invertierenden Integrators wurde erzeugt!\n')
print('--------------------done-------------------------')
#########################################################################################################
##################################################
# Der invertierende Differenzierer mit Fremddaten   # verlauf oke.. aber ein paar werte seltsam über 250 Volt ...?
##################################################
ue, ua , f = np.genfromtxt('data_of_others_cuz_ours_suck/diff/data_diff.txt',unpack = True)
f *= 1e3 # in Hz

#plot
fig, ax = plt.subplots()
plt.ylabel(r"$U_a \, / \, \mathrm{V}$")
plt.xlabel(r"$f \, / \, \mathrm{Hz}$")
ax.grid(ls =  '--')
ax.minorticks_on()
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(f,ua,'x', label = 'Messwerte' )
ax.plot(f[2],ua[2],'x',color='black', label = 'rausgenommen' )
ax.plot(f[4],ua[4],'x',color='black' )

# fit 
int_pos = np.array([0,1,3,5,6]) # damit wähle ich in ua und f nur bestimmte array werte aus, da ich nicht alle will. 
params1, covariance_matrix1 = curve_fit(fx, f[int_pos], ua[int_pos])
errors1 = np.sqrt(np.diag(covariance_matrix1))
m = ufloat(params1[0], errors1[0])
b = ufloat(params1[1], errors1[1])
print("Steigung am Differenzierer: ", m)
print("y-Achsenabschnitt am Differenzierer ", b,'\n')
plt.plot(f[int_pos],fx(f,*params1)[int_pos],'b',label= 'lineare Regression')

ax.legend(loc ='best')
plt.tight_layout()
plt.savefig('build/differenzierer.pdf')
plt.clf()

print('Der Plot des Differenzierers wurde erstellt!\n')
##########
# Tabelle
##########
t = TexTable([f, ua, ue], [r"$f$ / Hz",r"$U_a$ / V",r"$U_e$ / V"], 
            label='tab:differenzierer ',
            caption='Messwerte des invertierenden Differenzierers.')
t.set_row_rounding(0, 0) #reihe und rundung
t.set_row_rounding(1, 1)
t.set_row_rounding(2, 1)

t.write_file('build/tab_differenzierer.tex')
print('Die Tabelle des invertierenden Differenzierers wurde erzeugt!\n')
print('--------------------done-------------------------')
#########################################################################################################
#############################################
# Der invertierende Schmitt Trigger
#############################################
# Eingangsignal ist Dreieckspannnung
R1 = 1000 # Ohm
R2 = 100e3 
Ub = 15
U_thres_min = -(R1/R2) * 15
U_thres_max = (R1/R2) * 15

print('Der theoretische Wert der maximalen Schwellspannung ist ' ,U_thres_max,'V')
print('Der theoretische Wert der minimalen Schwellspannung ist ',U_thres_min,'V','\n')

print('Abgelesen anhand des Graphen wird U_thres_max = 151mV')
print('Abgelesen anhand des Graphen wird U_thres_min = −168mV ','\n')

Abweichung_max = (1-(150/151))*100
Abweichung_min = (1-(150/168))*100

print('Abweichung ist von U_thres_max_theo ist',Abweichung_max,'%')
print('Abweichung ist von U_thres_min_theo ist',Abweichung_min,'%')
print('--------------------done-------------------------')
#########################################################################################################
######################
# Der Signalgenerator
######################
# Schmitt Trigger und dann eine Ivertierender Integrator --> Rechteckspannung

R1 = 1000
R2 = 100000  
R3 = 100000
C = 1e-6

f = R2/(4*C*R1*R3) 
f_gem = 216
print('Die zu erwartende Freuquenz der Rechteckspannung beträgt ', f, 'Hz\n')
print('Die gemessene Freuquenz der Rechteckspannung beträgt ', f_gem, 'Hz\n')
print('Damit ergibt sich eine Abweichung von ',round((1-(f_gem/f))*100,2),'% \n')

print('Arschloch, Wixxer, Hurensohn deine Mutter hatte schon!')
print('Endlich ist dieses kack Praktikum vorbei! Meine Fresse hat mich die Scheiße angepisst.')
print('Was für eine elendige, reudige Kacke das doch war!')
print('Hurensohn!')