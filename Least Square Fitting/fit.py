import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import os
import get_data


countryname = "India"

df = get_data.getJHU(countryname,thresh=100,end="04/18/2020")
N = get_data.getPopulationData(countryname)

a_i = 0.25     # relative infectivity of infected who are asymptomatics
p_a = 0.6      # fraction of infected who are asymptomatics
f_a  = 0.1     # fraction of asymtomatics that get detected

if not os.path.exists(countryname):
    os.makedirs(countryname)

folder = '/a_i_{}_p_a_{}_f_a_{}/'.format(a_i,p_a,f_a)
os.makedirs(countryname+folder)
os.makedirs(countryname+folder+"bootstrap/")

predictionsize = 5                                 # Leave Last 'predictionsize' days out

time = np.arange(0,df.shape[0]-predictionsize,0.1) # Time Series
w_alpha= 0.3                                       # Weightage Parameter
g = 1/5.1                                          # Mean Gamma
r = p_a/(1-p_a)                                    # Ratio of Aymptomatics to Symptomatics
d_a = 1/8                                         # recovery rate for Asymptomatics
a_test = f_a/(1-f_a)                              # asyptomatic detection rate

# params are a,b,logE0,logIS0
guess = np.array([0.5,0.5,2,2,0.06,0.01,0.01])
param_bounds= ((0,1/8,0,0,0,0,0),(10,1/2,np.log10(N),np.log10(N),1,1,1))
data_weightage = w_alpha*np.power((1-w_alpha),np.tile(np.arange(df.shape[0]-predictionsize),3))

def curvefit_model(time,*p):

    E0 = 10**p[2]
    IS0 = 10**p[3]
    A0 = r*IS0
    QS0 = df.active[0]
    RS0 = df.recovered[0]
    D0 = df.dead[0]
    RU0 = r*df.recovered[0]
    SI0 = 0
    QA0 = 0
    RA0 = 0
    b = p[0]
    d = p[1]
    s = 0
    a = p[4]
    l = p[5]
    k = p[6]

    def dX(X,t):
        a_t=a_test
        M = np.array([[ -a,       0,     0,              0,     0,     0,     0,     0,     s,     0,     0],
                      [  0,-(1+r)*g,     0,              0,     0,     0,     0,     0,     0,     0,     0],
                      [  0,       g,    -d,              0,     0,     0,     0,     0,     0,     0,     0],
                      [  0,     r*g,     0,   -d_a-d_a*a_t,     0,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     d,              0,  -l-k,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,              0,     l,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,              0,     k,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,            d_a,     0,     0,     0,     0,     0,     0,     0],
                      [  a,       0,     0,              0,     0,     0,     0,     0,    -s,     0,     0],
                      [  0,       0,     0,        d_a*a_t,     0,     0,     0,     0,     0,    -l,     0],
                      [  0,       0,     0,              0,     0,     0,     0,     0,     0,     l,     0]])
        C = X[0]*X[2]*np.array([[  -b/N],
                                [   b/N],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0]])
        C = C + X[0]*X[3]*np.array([[-a_i*b/N],
                                    [ a_i*b/N],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0]])
        return (np.matmul(M,X.reshape(1,-1).T)+C).flatten()

    output = odeint(dX,[N-E0-IS0-QS0-RS0-D0-RU0-SI0-QA0-RA0,E0,IS0,A0,QS0,RS0,D0,RU0,SI0,QA0,RA0],time)
    return np.array([output[::10,4]+output[::10,9],output[::10,5],output[::10,6]]).T.flatten()


#print("Starting Initial Fit.....",end="")
try:
    popt,pcov = curve_fit(curvefit_model,time,df[["active","recovered","dead"]][:-predictionsize].values.flatten(),SI0=guess,method='trf',bounds=param_bounds,sigma=data_weightage)
except:
    print("failed")
#print("Fitted. \nCalculating Goodness of Fit for Uncertainly Analysis.....",end="")

def model(time,*p,changeATE=None,changeS=None,changeDA=None,lag=None):

    E0 = 10**p[2]
    IS0 = 10**p[3]
    A0 = r*IS0
    QS0 = df.active[0]
    RS0 = df.recovered[0]
    D0 = df.dead[0]
    RU0 = r*df.recovered[0]
    RA0 = 0
    SI0 = 0
    QA0 = 0
    b = p[0]
    d = p[1]
    ss = 0
    aa = p[4]
    l = p[5]
    k = p[6]

    def dX(X,t):
        if (lag is None and t>df.shape[0]) or (lag is not None and t>df.shape[0]+lag):
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
            if changeATE is not None:
                a_t=changeATE/(1-changeATE)
                s = 0
                a = aa
                d_a_ = d_a
            elif changeS is not None:
                a_t=a_test
                s=changeS
                a = aa
                d_a_ = d_a
            elif changeDA is not None:
                a_t=a_test
                s= ss
                a = aa
                d_a_ = changeDA
        else:
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
        M = np.array([[ -a,       0,     0,              0,     0,     0,     0,     0,     s,     0,     0],
                      [  0,-(1+r)*g,     0,              0,     0,     0,     0,     0,     0,     0,     0],
                      [  0,       g,    -d,              0,     0,     0,     0,     0,     0,     0,     0],
                      [  0,     r*g,     0, -d_a_-d_a_*a_t,     0,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     d,              0,  -l-k,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,              0,     l,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,              0,     k,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,           d_a_,     0,     0,     0,     0,     0,     0,     0],
                      [  a,       0,     0,              0,     0,     0,     0,     0,    -s,     0,     0],
                      [  0,       0,     0,       d_a_*a_t,     0,     0,     0,     0,     0,    -l,     0],
                      [  0,       0,     0,              0,     0,     0,     0,     0,     0,     l,     0]])
        C = X[0]*X[2]*np.array([[  -b/N],
                                [   b/N],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0]])
        C = C + X[0]*X[3]*np.array([[-a_i*b/N],
                                    [ a_i*b/N],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0]])
        return (np.matmul(M,X.reshape(1,-1).T)+C).flatten()

    return odeint(dX,[N-E0-IS0-QS0-RS0-D0-RU0-SI0-QA0-RA0,E0,IS0,A0,QS0,RS0,D0,RU0,SI0,QA0,RA0],time)

n = (df.shape[0]-predictionsize)*3                      # Number of Data Points
p = len(popt)                          # Number of Parameters
dof = max(0, n - p)                    # Degrees of Freedom

prediction = model(np.arange(0,df.shape[0],0.1),*popt)[::10,:]

RSS = np.sum((np.array([prediction[:-predictionsize,4]+prediction[:-predictionsize,9],prediction[:-predictionsize,5],prediction[:-predictionsize,6]]).T-df[["active","recovered","dead"]][:-predictionsize].values)**2)
AIC = 2*p+n*np.log(RSS/n)
print("RMS:",np.sqrt(RSS/n),"AIC:",AIC)
prederror = np.sum((np.array([prediction[-predictionsize:,4]+prediction[-predictionsize:,9],prediction[-predictionsize:,5],prediction[-predictionsize:,6]]).T-df[["active","recovered","dead"]][-predictionsize:].values)**2)

fig= plt.figure(figsize=(12,5))
ax1= fig.add_subplot(1,1,1)
ax1.plot(df.date[:-predictionsize],prediction[:-predictionsize,2]/1000,'-',color='gray',label="Symptomatic Infected")
ax1.plot(df.date[:-predictionsize],prediction[:-predictionsize,3]/1000,'-',color='salmon',label="Asymptomatic Infected")
ax1.plot(df.date[:-predictionsize],(prediction[:-predictionsize,4]+prediction[:-predictionsize,9])/1000,'r-',label="Quarantined")
ax1.plot(df.date[:-predictionsize],prediction[:-predictionsize,5]/1000,'g-',label="Recovered")
ax1.plot(df.date[:-predictionsize],prediction[:-predictionsize,6]/1000,'k-',label="Dead")
ax1.plot(df.date,df[["active","recovered","dead"]]/1000,'.')
ax1.plot(df.date[-predictionsize-1:],prediction[-predictionsize-1:,2]/1000,'--',color='gray',label="(P)Symptomatic Infected")
ax1.plot(df.date[-predictionsize-1:],prediction[-predictionsize-1:,3]/1000,'--',color='salmon',label="(P)Asymptomatic Infected")
ax1.plot(df.date[-predictionsize-1:],prediction[-predictionsize-1:,4]/1000+prediction[-predictionsize-1:,9]/1000,'r--',label="(P)Quarantined")
ax1.plot(df.date[-predictionsize-1:],prediction[-predictionsize-1:,5]/1000,'g--',label="(P)Recovered")
ax1.plot(df.date[-predictionsize-1:],prediction[-predictionsize-1:,6]/1000,'k--',label="(P)Dead")
plt.title("Leave Last {:0.1f} days: F(Asymp) = {:0.1f} Pop {:0.1f} Qua A = {:0.0f}% of S".format(predictionsize,p_a,f_a,a_i*100))
plt.legend()
plt.savefig(countryname+folder+"Leave{:0.1f}days.png".format(predictionsize))
plt.close(fig)
RS0 = (d_a*popt[0]*(a_test+1)+popt[1]*a_i*popt[0]*r)/(d_a*popt[1]*(r+1)*(a_test+1))
print(a_i,predictionsize,p_a,f_a,popt[0],1/g,1/popt[1],RS0,np.sqrt(RSS/n),np.sqrt(prederror/n),popt[5],popt[6],10**popt[2],10**popt[3],popt[4],1/d_a,sep=',')
fig= plt.figure(figsize=(12,2))
ax1= fig.add_subplot(1,1,1)
rt = RS0*prediction[:,0]/N
plt.plot(df.date,rt)
plt.savefig(countryname+folder+"RT_Leave{:0.1f}days.png".format(predictionsize))
plt.close(fig)

np.save(countryname+folder+"optimal_prediction.npy",model(np.arange(0,df.shape[0]+365,0.1),*popt)[::10,:])

F = np.cumsum(np.array([prediction[:-predictionsize,4]+prediction[:-predictionsize,9],prediction[:-predictionsize,5],prediction[:-predictionsize,6]]).T,axis=0)
Po_mean = np.concatenate([[[df.active[0],df.recovered[0],df.dead[0]]],np.diff(F,axis=0)],axis=0)

n_bootstrap = 1000
popts= []
fit_fail=0
for i in range(n_bootstrap):
    resample = np.random.poisson(Po_mean)
    try:
        popt_b,pcov_b = curve_fit(curvefit_model,time,resample.flatten(),SI0=guess,method='trf',bounds=param_bounds,sigma=data_weightage)
        popts.append(popt_b)
        np.save(countryname+folder+"bootstrap/optimal_prediction_{}.npy".format(i),model(np.arange(0,df.shape[0],0.1)+365,*popt)[::10,:])
    except:
        fit_fail+=1
    print("MC iteration: {} of {}".format(i+1,n_bootstrap-fit_fail),end="\r")
popts = np.array(popts)

np.save(countryname+folder+"parameters.npy",popts)
np.save(countryname+folder+"data.npy",df[["active","recovered","dead"]].values)
