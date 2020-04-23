import numpy as np
import pandas as pd

def getJHU(country,start=None,end=None,thresh= None):
    """
    return Pandas Dataframe with columns 'date','confirmed','recovered','dead','new_cases','active'
    Function to Get Data from JHU Database on GitHub
    url://https://github.com/CSSEGISandData/COVID-19

    country: str; Country/Region in JHU Database (see above)
    start: str; Start Date (MM/DD/YYYY)
    end: str; End Date (MM/DD/YYYY)

    """
    dfr = pd.read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
    dfd = pd.read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
    dfc = pd.read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

    dfc.head()

    dfc = dfc[dfc['Country/Region'] == country]
    dfr = dfr[dfr['Country/Region'] == country]
    dfd = dfd[dfd['Country/Region'] == country]
    dfr = dfr.drop(['Province/State','Lat','Long'],axis=1)
    dfd = dfd.drop(['Province/State','Lat','Long'],axis=1)
    dfc = dfc.drop(['Province/State','Lat','Long'],axis=1)
    dfc=dfc.groupby(by=dfc['Country/Region'],axis=0).sum().reset_index(drop=True).transpose()
    dfr=dfr.groupby(by=dfr['Country/Region'],axis=0).sum().reset_index(drop=True).transpose()
    dfd=dfd.groupby(by=dfd['Country/Region'],axis=0).sum().reset_index(drop=True).transpose()
    df= pd.concat([dfc,dfr,dfd],axis=1)

    df.reset_index(inplace=True)
    df.columns=['date','confirmed','recovered','dead']
    df['new_cases'] = np.concatenate([[0],np.diff(df.confirmed)])
    df['active'] = df['confirmed']-df['dead']-df['recovered']
    df = df.iloc[1:]
    df.date = pd.to_datetime(df.date)
    if start is not None:
        df = df[df.date>=pd.to_datetime(start)]
    if end is not None:
        df = df[df.date<=pd.to_datetime(end)]
    if thresh is not None:
        df = df[df.active>=thresh]
    df.to_csv("{}.csv".format(country),index=False)
    return df.reset_index(drop=True)

def getPopulationData(country):
    """
    return dtype(int)
    Function to Get Data from WorldBank Population Dataset hosted on GitHub
    url://https://github.com/technosap/technosap.github.io/blob/master/population.csv

    country: str; Country Name in World Bank Database (see above)

    """
    pop = pd.read_csv("https://raw.githubusercontent.com/technosap/technosap.github.io/master/population.csv")
    pop = pop[pop["Country Name"]==country]
    return int(pop[pop["Year"]==pop["Year"].max()].Value)
