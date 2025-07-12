import seaborn as sns
import random, itertools, math, pickle, joblib, scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import t
from scipy.linalg import sqrtm
from scipy.stats import norm

from operator import itemgetter
from itertools import groupby
from scipy.interpolate import pade
from IPython.display import display, Latex
import statsmodels.api as sm
from scipy.stats import iqr, genexpon
from matplotlib.ticker import FuncFormatter

plt.rcParams['figure.figsize'] = (10.0, 10.0) 
plt.rcParams['font.family'] = ['Times New Roman']
e = math.e

# some auxiliary functions
def I(x, v):
    return 1 if x <= v else 0

def Gaussian_Kernel(u):
    return math.exp(-1/2*u**2)/(math.sqrt(2*math.pi))

def Gaussian_Kernel_h(x, X, h):
    return Gaussian_Kernel((x-X)/h)/h

def scale(c, D):
    M = max(1, np.linalg.norm(D)/len(D)**.5)
    return c/M

def generate_delta():  
    u = random.uniform(0,1)
    if u<= 0.5:
        rst = 1
    else: 
        rst = -1
    return rst

def project(x_old, x_new, feasible_region):
    lb, ub = feasible_region
    # if lb <= x_new <= ub: rst = x_new
    # else: rst = x_old
    rst = min(ub, max(lb, x_new))
    return rst

def project_vec(x_new, ub_norm):
    rst = x_new/np.linalg.norm(x_new) *min(np.linalg.norm(x_new), ub_norm)
    return rst

# the main procedure
def sensitivity(N, psi, phi, d, theta, generate_rv,
                init_v,
                cst, cst_h1, cst_h2, cst_c,
                set_v, set_co, bound_D,
                crn=False,T=None):
    co = [init_v]; v = [init_v]; 
    Dv = [np.ones(d)]; Dco = [np.ones(d)]     
    for i in range(int(N)):
        h1 = cst_h1/(i+1)**(1/5)
        h2 = cst_h2/(i+1)**(1/7)    
        
        g = cst/(i+1)**.998
        a = cst/(i+1)**.997
        e = cst/(i+1)**.998
        b = cst/(i+1)**.999
        
        c = cst_c/(i+1)**(1/7)
        c = scale(c, Dv[-1])
        c = scale(c, Dco[-1])        
        
        delta = np.array([generate_delta() for j in range(d)])
        if crn==True:
            Z, L = generate_rv(theta, i, T)
            Zp, Lp = generate_rv(theta+c*delta, i, T)
            Zn, Ln = generate_rv(theta-c*delta, i, T)
        else:
            Z, L = generate_rv(theta)
            Zp, Lp = generate_rv(theta+c*delta)
            Zn, Ln = generate_rv(theta-c*delta)

        H = -(I(Lp, v[-1]+c*np.inner(delta, Dv[-1]))-I(Ln, v[-1]-c*np.inner(delta, Dv[-1])))/(2*c)*delta
        Dv_ = project_vec(Dv[-1]+g*H, bound_D)
        
        Kp = Gaussian_Kernel_h(v[-1]+c*np.inner(delta, Dv[-1]), Lp, h2)
        Kn = Gaussian_Kernel_h(v[-1]-c*np.inner(delta, Dv[-1]), Ln, h2)        
        H = (Kp*(psi-I(Zp, co[-1]+c*np.inner(delta, Dco[-1])))-Kn*(psi-I(Zn, co[-1]-c*np.inner(delta, Dco[-1]))))/(2*c)*delta
        Dco_ = project_vec(Dco[-1]+a*H, bound_D)
        
        co_ = project(co[-1],co[-1]+e*Gaussian_Kernel_h(v[-1], L, h1)*(psi-I(Z, co[-1])), set_co)
        v_ = project(v[-1],v[-1]+b*(phi-I(L, v[-1])), set_v)
        
        v.append(v_)
        Dv.append(Dv_)
        Dco.append(Dco_)      
        co.append(co_)
    return co, v, Dv, Dco

# independent replications
def replication(N, psi, phi, d, theta, generate_rv,
                init_v,
                cst, cst_h1, cst_h2, cst_c,
                set_v, set_co, bound_D,crn=False,df=None):
    co_lst = []; v_lst = []; Dv_lst = []; Dco_lst = []
    r = 0
    while True:
        with np.errstate(invalid='raise'):
            try:
                if crn==True:
                    T=t.rvs(df=df,size=2*int(N))
                else:
                    T=None
                co, v, Dv, Dco = sensitivity(N, psi, phi, d, theta, generate_rv,init_v,
                                             cst, cst_h1, cst_h2, cst_c,
                                             set_v, set_co, bound_D,crn=crn,T=T)
                co_lst.append(co); v_lst.append(v); 
                Dv_lst.append(Dv); Dco_lst.append(Dco); 
                r+=1
            except (FloatingPointError, ZeroDivisionError):
                print('Error: Division by Zero')
        if r>=100:
            break
    return co_lst, v_lst, Dv_lst, Dco_lst