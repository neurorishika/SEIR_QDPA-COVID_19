import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from subprocess import call
from scipy.stats.distributions import gamma
import json
import wget
import os

def pooled_SD(sds,means):
    return np.sqrt((np.sum(sds**2,axis=0)+np.sum(means-np.mean(means,axis=0)))/sds.shape[0])

data = pd.read_csv("Data/India.csv")
data.date = pd.to_datetime(data.date)

wget.download('https://api.covid19india.org/data.json', os.getcwd()+"\\Data\\covid19india.json")

dates = np.array([pd.to_datetime(i['date']+"2020") for i in json.load(open('Data/covid19india.json',))['cases_time_series']])
confirmed = np.array([int(i['dailyconfirmed'])for i in json.load(open('Data/covid19india.json',))['cases_time_series']])

confirmed = confirmed[data.date[0]<=dates]
dates = dates[data.date[0]<=dates]

real_data = confirmed[data.date[0]<=dates]
df = pd.DataFrame()
df['active'] = real_data
df['imported'] = np.zeros(real_data.shape[0])
lol = [15,2,1,2,3,6,2,7,10,5,11,5,12,14,21,27,34,35,45,55,38,34,24,40,16,20,11,20,7,8,6,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
df['imported'][:len(lol)] = lol

df.to_csv('temp_files/dataset.csv',index=False)

rt = []
for n in range(1000):
    print("Iteration: ",n+1,end='\r')
    dataset = np.copy(real_data)
    G = gamma(3.325+0.616*np.random.normal(),0.979+0.195*np.random.normal())
    for i in range(len(dataset)):
        send_back = np.clip(np.round(G.rvs(dataset[i])),0,10)
        send_back = send_back[i-send_back>=0]
        dataset[i] = 0
        for j in np.unique(np.int32(send_back)):
            dataset[i-j] += np.sum(send_back==j)
    df = pd.DataFrame()
    df['active'] = dataset[:-10]
    df.to_csv('temp_files/dataset.csv',index=False)
    call(['RScript.exe','Rscripts/Rt-bootstrap (Uncorrected).R'])
    rt.append(pd.read_csv('temp_files/rtoutput.csv'))

means = np.array([x["Mean(R)"].values for x in rt])
sds = np.array([x["Std(R)"].values for x in rt])
np.save('Output/Rt/India/India_rt_means.npy',means)
np.save('Output/Rt/India/India_rt_sd.npy',sds)
np.save('Output/Rt/India/India_rt_dates.npy',pd.Series(dates))
