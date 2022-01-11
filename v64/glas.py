import numpy as np
from uncertainties import ufloat
from scipy.stats import stats
from texutils.table import TexTable

M = np.genfromtxt('data/durchgaenge_glas.txt', unpack=True) 

T = 1 * 10**-3
lam = 632.990 * 10**(-9)           
theta_0 = 10*2 * np.pi / 180        # der faktor muss noch bearbeitet werden es sind nicht ganz 2 TODO
theta_kipp = 10 * np.pi / 180
#theta = D_theta * (np.pi)/180 #in radiant

def n(M):
    return (1-( (lam*M) / (T*theta_0*theta_kipp) ))**(-1)

print('\nArray der Brechungsindexe')
print(n(M),'\n')

M_mean= np.mean(M)
print('Die mittlere Anzahl der Maxima ist ')
print(M_mean)

b_index = ufloat(np.mean(n(M)), stats.sem(n(M)))
print('----------')
print('jetzt mit fehler des mittelwerts')
print('--> ',b_index,' Brechungsindex Glas.')

print('-------------')
##########
# Tabelle
##########
t = TexTable([M], [r"Anzahl der Maxima"], 
            label='tab:Maxima',
            caption='Gemessene Anzahl der Maxima bei Rotation des Glases.')
t.set_row_rounding(0,0) #reihe und rundung
t.write_file('build/tabMaxima_Glas.tex')
print('Die Tabelle der Maxima bei Glas wurde erzeugt!\n')