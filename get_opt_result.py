# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 00:41:44 2019

@author: meiru
"""
import pandas as pd
import numpy as np
import datetime
import matplotlib.pyplot as plt

portHY = pd.read_csv('Cov_Input_ret_mat_HY.csv')
portIG = pd.read_csv('Cov_Input_ret_mat_IG.csv')
baseHY = pd.read_csv("Cov_Input_benchmark_HY.csv")
baseIG = pd.read_csv("Cov_Input_benchmark_IG.csv")

wtHY = pd.read_csv("rolling_opt_result_HY.csv")
wtIG = pd.read_csv("rolling_opt_result_IG.csv")

portHY['Date'] = pd.to_datetime(portHY['Date']).dt.date
dateHY = portHY['Date']
portHY = portHY.set_index(['Date'])
portIG['Date'] = pd.to_datetime(portIG['Date']).dt.date
dateIG = portIG['Date']
portIG = portIG.set_index(['Date'])
wtHY['Date'] = pd.to_datetime(wtHY['Date']).dt.date
wtHY = wtHY.set_index(['Date'])
wtIG['Date'] = pd.to_datetime(wtIG['Date']).dt.date
wtIG = wtIG.set_index(['Date'])
baseHY['Date'] = pd.to_datetime(baseHY['Date']).dt.date
baseHY = baseHY.set_index(['Date'], drop = False)
baseIG['Date'] = pd.to_datetime(baseIG['Date']).dt.date
baseIG = baseIG.set_index(['Date'], drop = False)

def getReturn(port, wt_df, base):
    ret = pd.DataFrame(index=np.arange(len(port)), columns = ['Date', 'portfolio', 'benchmark'])
    for row in range(len(port)):
        returns = np.array(port.iloc[row])
        wt = np.array(wt_df.iloc[row])
        wt_returns = (np.dot(returns, wt)).sum()
        ret.iloc[row]['portfolio'] = wt_returns
        ret.iloc[row]['benchmark'] = base.iloc[row]['benchmark']
        ret.iloc[row]['Date'] = base.iloc[row]['Date']
    return ret

def getError(df):
    res = pd.DataFrame(columns = ['Date', 'Error'])
    df = df[df['portfolio'] != float(0)] 
    res['Error'] = df['benchmark'] - df['portfolio']
    res['Date'] = df['Date']

    errorMean = res['Error'].mean()
    errorSqW = ((res['Error'] - errorMean)**2).sum()
    errorW = (errorSqW / len(res['Error'])) ** (1/2)

    errorSq = (res['Error']**2).sum()
    errorR = (errorSq / (len(res['Error']) - 1)) **(1/2)
    return res, errorR, errorW

def getTurnover(wt, date_df):
    turnoverL = pd.DataFrame(columns = ['Date', 'Turnover'])
    wt = wt[(wt.T != float(0)).any()]
    date_df = date_df[-len(wt):]
    for month in range(1, len(wt)):
        turnover = abs(wt.iloc[month] - wt.iloc[month - 1]).sum()
        date = date_df.iloc[month]
        turnoverL.loc[month] = [date, turnover]
    
    turn_avg = []
    turnoverL['Year'] = pd.to_datetime(turnoverL['Date']).dt.year
    for y, g in turnoverL.groupby(turnoverL['Year']):
        turn_avg.append(g['Turnover'].sum())
    turn_avg_annual = sum(turn_avg) / len(turn_avg)
    return turnoverL, turn_avg_annual



return_HY = getReturn(portHY, wtHY, baseHY)
Error_HY, errHYR, errHYW = getError(return_HY)
return_IG = getReturn(portIG, wtIG, baseIG)
Error_IG, errIGR, errIGW = getError(return_IG)

turnover_HY, turnover_Annual_HY = getTurnover(wtHY, dateHY)
turnover_IG, turnover_Annual_IG = getTurnover(wtIG, dateIG)

return_HY.plot()
return_IG.plot()

return_HY.to_csv('Opt_return_HY.csv')
return_IG.to_csv('Opt_return_IG.csv')
