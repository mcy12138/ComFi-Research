#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 13:26:34 2017

@author: burtonh
"""

##########################################################################################
#### Calling Modules #####################################################################

import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas_datareader.data as web
import pandas as pd

from cvxopt import matrix
from cvxopt.solvers import qp
from cvxopt.blas import dot
from cvxopt.solvers import options


##########################################################################################
##########################################################################################
## Functions #############################################################################

def gmv_wts(return_mat):
    #Returns and Variances
    avg_ret  = return_mat.mean(axis=0);
    sigma    = np.cov(return_mat.T);           n        = np.size(avg_ret);
    
    #Creating CVXOPT Variables
    avg_ret  = matrix(avg_ret);                sigma    = matrix(sigma);
    P        = sigma;                          q        = matrix(np.zeros((n,1)));

    #Inequality Constraint: Gpx <= h captures the constraints ((x >= 0)
    Gp    = matrix(-np.identity(n));            hp    = matrix(np.zeros((n,1)));

    #Global Minimum Variance Portfolio
    A        = matrix(1.0, (1,n));
    b        = matrix(1.0)
    options['show_progress'] = False
    sol      = qp(P, q, Gp, hp, A, b);
    wgmv     = np.array(sol['x'])
    return wgmv

##########################################################################################

def gmv_rtn(return_mat):

    #Returns and Variances
    avg_ret  = return_mat.mean(axis=0);
    sigma    = np.cov(return_mat.T);           n        = np.size(avg_ret);

    #Creating CVXOPT Variables
    avg_ret  = matrix(avg_ret);                sigma    = matrix(sigma);
    P        = sigma;                          q        = matrix(np.zeros((n,1)));

    #Inequality Constraint: Gpx <= h captures the constraints ((x >= 0)
    Gp    = matrix(-np.identity(n));            hp    = matrix(np.zeros((n,1)));

    #Global Minimum Variance Portfolio
    A        = matrix(1.0, (1,n));
    b        = matrix(1.0)
    sol      = qp(P, q, Gp, hp, A, b);
    wgmv     = sol['x']
    rtn      = dot(avg_ret,wgmv);
    r_min    = np.array(rtn)
    App      = matrix(np.concatenate((np.array(A),np.array(avg_ret).T)))
    bp       = matrix(np.concatenate((b,r_min*b)))
    options['show_progress'] = False
    sol      = qp(P, q, Gp, hp, App, bp)
    return r_min

#######################################################################################
    
# Rolling Samples
def rolling(return_mat,assets,num_obs):

    #Creating place holder for weights in return_mat
    n_return_mat = return_mat
    cov_mat = []
    n = num_obs

    for i in assets:
        n_return_mat[i+"_w_GMV"]  = 0.0
        l      = len(n_return_mat)
        s      = len(assets)
        w_GMV  = [i+"_w_GMV" for i in assets]

    for i in range(n,l):   
        #rolling window
#        w_gmv       = gmv_wts(n_return_mat[assets].iloc[i-n:i])
#        cov_mat.append(n_return_mat[assets].iloc[i-n:i].cov())
        
        #expanding window
        w_gmv = gmv_wts(n_return_mat[assets].iloc[0:i])
        cov_mat.append(n_return_mat[assets].iloc[0:i].cov())
        
        for j in range(s):
            n_return_mat[w_GMV[j]].iloc[i] = w_gmv[j]

    ## Final Weights
#    w_gmv       = gmv_wts(n_return_mat[assets].iloc[l-n:l])
#    for j in range(s):
#        n_return_mat[w_GMV[j]].iloc[l-1]  = w_gmv[j]
#    n_return_mat  = n_return_mat.dropna()
    return n_return_mat, cov_mat


##########################################################################################
##########################################################################################
##########################################################################################
##### 1. Data ############################################################################

##### HY #####
sample_returns_HY = pd.read_csv("Cov_Input_ret_mat_HY.csv")
sample_returns_HY['Date'] = pd.to_datetime(sample_returns_HY.Date)
sample_returns_HY.set_index('Date', inplace=True)

bench_returns_HY = pd.read_csv("Cov_Input_benchmark_HY.csv")
bench_returns_HY['Date'] = pd.to_datetime(bench_returns_HY.Date)
bench_returns_HY.set_index('Date', inplace=True)
samples_HY = sample_returns_HY.columns.values

for i in range(len(samples_HY)):
    name = samples_HY[i] + "_dif"
    bench_returns_HY[name] = sample_returns_HY.iloc[:,i] - bench_returns_HY["benchmark"]
return_mat_HY = bench_returns_HY
return_mat_HY = return_mat_HY.drop('benchmark',axis=1)
assets_HY = return_mat_HY.columns.values    


##### IG #####
sample_returns_IG= pd.read_csv("Cov_Input_ret_mat_IG.csv")
sample_returns_IG['Date'] = pd.to_datetime(sample_returns_IG.Date)
sample_returns_IG.set_index('Date', inplace=True)

bench_returns_IG = pd.read_csv("Cov_Input_benchmark_IG.csv")
bench_returns_IG['Date'] = pd.to_datetime(bench_returns_IG.Date)
bench_returns_IG.set_index('Date', inplace=True)
samples_IG = sample_returns_IG.columns.values

for i in range(len(samples_IG)):
    name = samples_IG[i] + "_dif"
    bench_returns_IG[name] = sample_returns_IG.iloc[:,i] - bench_returns_IG["benchmark"]
return_mat_IG = bench_returns_IG
return_mat_IG = return_mat_IG.drop('benchmark',axis=1)
assets_IG = return_mat_IG.columns.values    


## Rolling Windows on Actual Data
#l       = len(return_mat)
num_obs = 36
return_mat_roll_HY, cov_mat_roll_HY = rolling(return_mat_HY,assets_HY,num_obs)
return_mat_roll_HY.loc[:, '0_dif_w_GMV':].to_csv('rolling_opt_result_HY.csv')

return_mat_roll_IG, cov_mat_roll_IG = rolling(return_mat_IG,assets_IG,num_obs)
return_mat_roll_IG.loc[:, '0_dif_w_GMV':].to_csv('rolling_opt_result_IG.csv')
#bench_returns['benchmark'].plot()
#return_mat_roll['GMV'].plot()



##########################################################################################
##########################################################################################
















