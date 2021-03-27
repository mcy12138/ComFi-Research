import pandas as pd
import numpy as np
import datetime

portHY = pd.read_csv('C:\\Users\\meiru\\Documents\\mary\\university\\CUSR\\Alternative Index Replication in Fixed Income\\Data\\Portfolio_wt_HY.csv')
portIG = pd.read_csv('C:\\Users\\meiru\\Documents\\mary\\university\\CUSR\\Alternative Index Replication in Fixed Income\\Data\\Portfolio_wt_IG.csv')
baseHY = pd.read_csv("C:\\Users\\meiru\\Documents\\mary\\university\\CUSR\\Alternative Index Replication in Fixed Income\\Data\\BaselineHY.csv")
baseIG = pd.read_csv("C:\\Users\\meiru\\Documents\\mary\\university\\CUSR\\Alternative Index Replication in Fixed Income\\Data\\BaselineIG.csv")

baseHY['Date'] = pd.to_datetime(baseHY['prd_formatted']).dt.to_period('D')
baseHY = baseHY.set_index(['Date'])
baseHY = baseHY.drop(['prd_formatted', 'universe'], axis = 1)

baseIG['Date'] = pd.to_datetime(baseIG['prd_formatted']).dt.to_period('D')
baseIG = baseIG.set_index(['Date'])
baseIG = baseIG.drop(['prd_formatted', 'universe'], axis = 1)

portHY['Date'] = pd.to_datetime(portHY['Date']).dt.to_period('D')
portHY_by_month = []
for _, g in portHY.groupby(portHY['Date']):
    portHY_by_month.append(g)
    
portIG['Date'] = pd.to_datetime(portIG['Date']).dt.to_period('D')
portIG_by_month = []
for _, g in portIG.groupby(portIG['Date']):
    portIG_by_month.append(g)
    
def getReturnbyBucket(Port, numBucket):
    res = pd.DataFrame(index=np.arange(len(Port)), columns = ['Date'])
    for month in range(len(Port)):
        res.iloc[month]['Date'] = Port[month].iloc[0]['Date']
    for b in range(numBucket):
        ret_bucket = []
        for month in range(len(Port)):
            portfolio = Port[month]
            portfolio_bucket = portfolio.loc[(portfolio['bucket'] == b)]
#            tt_mv = portfolio_bucket['market_value'].sum()
#            wt = portfolio_bucket['market_value'] / tt_mv
#            tt_ret = (portfolio_bucket['total_return_mtd'] * wt).sum()
            
            tt_ret = (portfolio_bucket['total_return_mtd']).sum()
            ret_bucket_month = tt_ret / len(portfolio_bucket)
            ret_bucket.append(ret_bucket_month)
        res[str(b)] = pd.Series(ret_bucket).values
    return res

retHY = getReturnbyBucket(portHY_by_month, 27)
retHY = retHY.set_index(['Date'])
retHY.to_csv('Cov_Input_ret_mat_HY.csv')

retIG = getReturnbyBucket(portIG_by_month, 27)
retIG = retIG.set_index(['Date'])
retIG.to_csv('Cov_Input_ret_mat_IG.csv')

baseHY.to_csv('Cov_Input_benchmark_HY.csv')
baseIG.to_csv('Cov_Input_benchmark_IG.csv')
