import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from subprocess import call
from scipy.stats.distributions import gamma
import os

state_id = {"KL":"Kerala",
"MH":"Maharashtra",
"GJ":"Gujarat",
"DL":"Delhi",
"RJ":"Rajasthan",
"TN":"Tamil Nadu",
"MP":"Madhya Pradesh",
"UP":"Uttar Pradesh",
"TG":"Telangana",
"AP":"Andhra Pradesh",
"KA":"Karnataka",
"WB":"West Bengal",
"JK":"Jammu and Kashmir",
"HR":"Haryana",
"PB":"Punjab",
"BR":"Bihar" }

data = pd.read_csv("Data/States.csv")
for state in state_id.keys():
    boots = 1000
    real_data = data[state].values
    df = pd.DataFrame()
    df['active'] = real_data
    df.to_csv('temp_files/dataset.csv',index=False)
    try:
        means = np.load("Output/Rt/States/Rt_"+state+"_means.npy")
        sds = np.load("Output/Rt/States/Rt_"+state+"_sds.npy")
        dat_means = np.load("Output/Rt/States/"+state+"_means.npy")
        dat_sds = np.load("Output/Rt/States/"+state+"_sds.npy")
    except:
        rt = []
        dats = []
        for n in range(boots):
            print("Iteration: ",n+1,end='\r')
            G = gamma(3.325+0.616*np.random.normal(),0.979+0.195*np.random.normal())
            dataset = np.copy(real_data)
            for i in range(len(dataset)):
                send_back = np.clip(np.round(G.rvs(int(dataset[i]))),0,10)
                send_back = send_back[i-send_back>=0]
                dataset[i] = 0
                for j in np.unique(np.int32(send_back)):
                    dataset[i-j] += np.sum(send_back==j)
            df = pd.DataFrame()
            df['active'] = dataset[:-10]
            dats.append(dataset[:-10])
            df.to_csv('temp_files/dataset.csv',index=False)
            call(['RScript.exe','Rscripts/Rt-bootstrap (Uncorrected).R'])
            rt.append(pd.read_csv('temp_files/rtoutput.csv'))

        means = np.array([x["Mean(R)"].values for x in rt])
        sds = np.array([x["Std(R)"].values for x in rt])
        np.save("Output/Rt/States/Rt_"+state+"_means.npy",means)
        np.save("Output/Rt/States/Rt_"+state+"_sds.npy",sds)

        dat_means = np.mean(dats,axis=0)
        dat_sds = np.std(dats,axis=0)
        np.save("Output/Rt/States/"+state+"_means.npy",dat_means)
        np.save("Output/Rt/States/"+state+"_sds.npy",dat_sds)
