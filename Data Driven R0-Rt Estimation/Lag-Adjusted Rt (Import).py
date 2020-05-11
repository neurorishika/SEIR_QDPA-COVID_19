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

real_data = confirmed
df = pd.DataFrame()
df['active'] = real_data
df.to_csv('temp_files/dataset.csv',index=False)

rt1 = []
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
    rt1.append(dataset[:-10])

real_data = data["imported"].values
df = pd.DataFrame()
df['active'] = real_data
df.to_csv('temp_files/dataset.csv',index=False)

rt2 = []
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
    rt2.append(dataset[:-10])

np.save("Output/Rt/India/India_Onset_means_total.npy",np.array(rt1).mean(axis=0))
np.save("Output/Rt/India/India_Onset_sds_total.npy",np.array(rt1).std(axis=0))
np.save("Output/Rt/India/India_Onset_means_import.npy",np.array(rt2).mean(axis=0))
np.save("Output/Rt/India/India_Onset_sds_import.npy",np.array(rt2).std(axis=0))

local = np.max([np.array(rt1)[:,:-1],np.array(rt2)],axis=0)-np.array(rt2)
imported = np.array(rt2)

rt = []
for i in range(1000):
    print("{}/{}".format(i+1,local.shape[0]),end='\r')
    df_import = pd.DataFrame()
    df_import["local"] = local[i,:]
    df_import["imported"] = imported[i,:]
    df_import.to_csv("temp_files/imported.csv",index=False)
    call(['RScript.exe','Rscripts/Rt-bootstrap (Import Corrected).R'])
    rt.append(pd.read_csv('temp_files/rtoutput.csv'))

means = np.array([x["Mean(R)"].values for x in rt])
sds = np.array([x["Std(R)"].values for x in rt])
np.save('Output/Rt/India/India_rt_import_means.npy',means)
np.save('Output/Rt/India/India_rt_import_sd.npy',sds)
