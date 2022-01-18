import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from table import TexTable
from scipy import stats
from scipy.optimize import curve_fit


def neueWerte(file_name="data/dataa.txt", finished_file="build/taba.tex",  vars_name=[r"t/\si{\second}", r"s/\si{\meter}"], label_text="taba", caption_text=r"Neue Tabelle." , precision=2):
    vars = np.genfromtxt(file_name, unpack=True)
    tab_name = TexTable([*vars], vars_name, label= label_text , caption= caption_text, roundPrecision= precision)
    tab_name.writeFile(finished_file)
    return vars

def tabelle(vars, finished_file="taba.tex", vars_name=[r"t/\si{\second}", r"s/\si{\meter}"], label_text="taba", caption_text=r"Eine neue Tablle", precision=2):
    tab_name = TexTable([*vars], vars_name, label= label_text , caption= caption_text, roundPrecision= precision)
    tab_name.writeFile(finished_file)

def gerade(x, m, n):
    return m*x + n

def linReg(x, y, p, q, x_name=r"t/\si{\second}", y_name=r"x/\si{\meter}", num=1,  x_add=5, file_name="build/plota.pdf"):
    Steigung1, yAbschnitt1, r_value1, p_value1, std_err1= stats.linregress(p,q)
    plt.figure(num) 
    newx= np.linspace(p[0]-x_add,p[-1]+x_add, num=1000)
    plt.plot(x, y, "xr", label="Daten")
    plt.plot(newx, gerade(newx, Steigung1, yAbschnitt1), "r", label="Fit", linewidth=1.0)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(file_name) 
    return Steigung1, yAbschnitt1, std_err1

def plot(x, y, x_name=r"t/\si{\second}", y_name=r"x/\si{\meter}", num=1, file_name="build/plota.pdf"):
    plt.figure(num) 
    plt.plot(x, y, "xr", label="Daten")
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(file_name) 

def curvefit(x, y, num=1, x_add=5, function=gerade, x_name=r"t/\si{\second}", y_name=r"s/\si{\meter}", file_name="build/plota.pdf"):
    params, pcov = curve_fit(function, x, y)
    plt.figure(num) 
    newx= np.linspace(x[0]-x_add,x[-1]+x_add, num=1000)
    plt.plot(x, y, "xr", label="Daten")
    plt.plot(newx, function(newx, *params), "r", label="Fit", linewidth=1.0)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(file_name)
    err = np.sqrt(np.diag(pcov))
    return params, err


