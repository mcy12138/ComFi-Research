# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:09:24 2019

@author: 12833
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
#from sklearn import preprocessing

HY = pd.read_csv('data/HY_TEST.csv')
HY['date'] = pd.to_datetime(HY['prd_formatted'])
HY = HY.set_index(['date'])


HY_2007 = HY.loc['2007', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2008 = HY.loc['2008', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2009 = HY.loc['2009', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2010 = HY.loc['2010', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2011 = HY.loc['2011', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2012 = HY.loc['2012', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2013 = HY.loc['2013', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2014 = HY.loc['2014', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2015 = HY.loc['2015', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2016 = HY.loc['2016', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2017 = HY.loc['2017', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2018 = HY.loc['2018', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]
HY_2019 = HY.loc['2019', ['class_2', 'cusip','total_return_volatility', 'market_value', 
                        'prev_6m_total_return', 'coupon', 'oas', 'duration']]

UTIL_2007 = HY_2007.loc[HY_2007['class_2']=='UTILITY']
del UTIL_2007['class_2']
UTIL_2007 = UTIL_2007.reset_index()
UTIL_2007.set_index(['date', 'cusip'], inplace = True)

UTIL_2008 = HY_2008.loc[HY_2008['class_2']=='UTILITY']
del UTIL_2008['class_2']
UTIL_2008 = UTIL_2008.reset_index()
UTIL_2008.set_index(['date', 'cusip'], inplace = True)

UTIL_2009 = HY_2009.loc[HY_2009['class_2']=='UTILITY']
del UTIL_2009['class_2']
UTIL_2009 = UTIL_2009.reset_index()
UTIL_2009.set_index(['date', 'cusip'], inplace = True)

UTIL_2010 = HY_2010.loc[HY_2010['class_2']=='UTILITY']
del UTIL_2010['class_2']
UTIL_2010 = UTIL_2010.reset_index()
UTIL_2010.set_index(['date', 'cusip'], inplace = True)

UTIL_2011 = HY_2011.loc[HY_2011['class_2']=='UTILITY']
del UTIL_2011['class_2']
UTIL_2011 = UTIL_2011.reset_index()
UTIL_2011.set_index(['date', 'cusip'], inplace = True)

UTIL_2012 = HY_2012.loc[HY_2012['class_2']=='UTILITY']
del UTIL_2012['class_2']
UTIL_2012 = UTIL_2012.reset_index()
UTIL_2012.set_index(['date', 'cusip'], inplace = True)

UTIL_2013 = HY_2013.loc[HY_2013['class_2']=='UTILITY']
del UTIL_2013['class_2']
UTIL_2013 = UTIL_2013.reset_index()
UTIL_2013.set_index(['date', 'cusip'], inplace = True)

UTIL_2014 = HY_2014.loc[HY_2014['class_2']=='UTILITY']
del UTIL_2014['class_2']
UTIL_2014 = UTIL_2014.reset_index()
UTIL_2014.set_index(['date', 'cusip'], inplace = True)

UTIL_2015 = HY_2015.loc[HY_2015['class_2']=='UTILITY']
del UTIL_2015['class_2']
UTIL_2015 = UTIL_2015.reset_index()
UTIL_2015.set_index(['date', 'cusip'], inplace = True)

UTIL_2016 = HY_2016.loc[HY_2016['class_2']=='UTILITY']
del UTIL_2016['class_2']
UTIL_2016 = UTIL_2016.reset_index()
UTIL_2016.set_index(['date', 'cusip'], inplace = True)

UTIL_2017 = HY_2017.loc[HY_2017['class_2']=='UTILITY']
del UTIL_2017['class_2']
UTIL_2017 = UTIL_2017.reset_index()
UTIL_2017.set_index(['date', 'cusip'], inplace = True)

UTIL_2018 = HY_2018.loc[HY_2018['class_2']=='UTILITY']
del UTIL_2018['class_2']
UTIL_2018 = UTIL_2018.reset_index()
UTIL_2018.set_index(['date', 'cusip'], inplace = True)

UTIL_2019 = HY_2019.loc[HY_2019['class_2']=='UTILITY']
del UTIL_2019['class_2']
UTIL_2019 = UTIL_2019.reset_index()
UTIL_2019.set_index(['date', 'cusip'], inplace = True)

UTIL_mon = []
#2007
UTIL_mon.append(UTIL_2007.loc['2007-01', :])
UTIL_mon.append(UTIL_2007.loc['2007-02', :])
UTIL_mon.append(UTIL_2007.loc['2007-03', :])
UTIL_mon.append(UTIL_2007.loc['2007-04', :])
UTIL_mon.append(UTIL_2007.loc['2007-05', :])
UTIL_mon.append(UTIL_2007.loc['2007-06', :])
UTIL_mon.append(UTIL_2007.loc['2007-07', :])
UTIL_mon.append(UTIL_2007.loc['2007-08', :])
UTIL_mon.append(UTIL_2007.loc['2007-09', :])
UTIL_mon.append(UTIL_2007.loc['2007-10', :])
UTIL_mon.append(UTIL_2007.loc['2007-11', :])
UTIL_mon.append(UTIL_2007.loc['2007-12', :])
#2008
UTIL_mon.append(UTIL_2008.loc['2008-01', :])
UTIL_mon.append(UTIL_2008.loc['2008-02', :])
UTIL_mon.append(UTIL_2008.loc['2008-03', :])
UTIL_mon.append(UTIL_2008.loc['2008-04', :])
UTIL_mon.append(UTIL_2008.loc['2008-05', :])
UTIL_mon.append(UTIL_2008.loc['2008-06', :])
UTIL_mon.append(UTIL_2008.loc['2008-07', :])
UTIL_mon.append(UTIL_2008.loc['2008-08', :])
UTIL_mon.append(UTIL_2008.loc['2008-09', :])
UTIL_mon.append(UTIL_2008.loc['2008-10', :])
UTIL_mon.append(UTIL_2008.loc['2008-11', :])
UTIL_mon.append(UTIL_2008.loc['2008-12', :])
#2009
UTIL_mon.append(UTIL_2009.loc['2009-01', :])
UTIL_mon.append(UTIL_2009.loc['2009-02', :])
UTIL_mon.append(UTIL_2009.loc['2009-03', :])
UTIL_mon.append(UTIL_2009.loc['2009-04', :])
UTIL_mon.append(UTIL_2009.loc['2009-05', :])
UTIL_mon.append(UTIL_2009.loc['2009-06', :])
UTIL_mon.append(UTIL_2009.loc['2009-07', :])
UTIL_mon.append(UTIL_2009.loc['2009-08', :])
UTIL_mon.append(UTIL_2009.loc['2009-09', :])
UTIL_mon.append(UTIL_2009.loc['2009-10', :])
UTIL_mon.append(UTIL_2009.loc['2009-11', :])
UTIL_mon.append(UTIL_2009.loc['2009-12', :])
#2010
UTIL_mon.append(UTIL_2010.loc['2010-01', :])
UTIL_mon.append(UTIL_2010.loc['2010-02', :])
UTIL_mon.append(UTIL_2010.loc['2010-03', :])
UTIL_mon.append(UTIL_2010.loc['2010-04', :])
UTIL_mon.append(UTIL_2010.loc['2010-05', :])
UTIL_mon.append(UTIL_2010.loc['2010-06', :])
UTIL_mon.append(UTIL_2010.loc['2010-07', :])
UTIL_mon.append(UTIL_2010.loc['2010-08', :])
UTIL_mon.append(UTIL_2010.loc['2010-09', :])
UTIL_mon.append(UTIL_2010.loc['2010-10', :])
UTIL_mon.append(UTIL_2010.loc['2010-11', :])
UTIL_mon.append(UTIL_2010.loc['2010-12', :])
#2011
UTIL_mon.append(UTIL_2011.loc['2011-01', :])
UTIL_mon.append(UTIL_2011.loc['2011-02', :])
UTIL_mon.append(UTIL_2011.loc['2011-03', :])
UTIL_mon.append(UTIL_2011.loc['2011-04', :])
UTIL_mon.append(UTIL_2011.loc['2011-05', :])
UTIL_mon.append(UTIL_2011.loc['2011-06', :])
UTIL_mon.append(UTIL_2011.loc['2011-07', :])
UTIL_mon.append(UTIL_2011.loc['2011-08', :])
UTIL_mon.append(UTIL_2011.loc['2011-09', :])
UTIL_mon.append(UTIL_2011.loc['2011-10', :])
UTIL_mon.append(UTIL_2011.loc['2011-11', :])
UTIL_mon.append(UTIL_2011.loc['2011-12', :])
#2012
UTIL_mon.append(UTIL_2012.loc['2012-01', :])
UTIL_mon.append(UTIL_2012.loc['2012-02', :])
UTIL_mon.append(UTIL_2012.loc['2012-03', :])
UTIL_mon.append(UTIL_2012.loc['2012-04', :])
UTIL_mon.append(UTIL_2012.loc['2012-05', :])
UTIL_mon.append(UTIL_2012.loc['2012-06', :])
UTIL_mon.append(UTIL_2012.loc['2012-07', :])
UTIL_mon.append(UTIL_2012.loc['2012-08', :])
UTIL_mon.append(UTIL_2012.loc['2012-09', :])
UTIL_mon.append(UTIL_2012.loc['2012-10', :])
UTIL_mon.append(UTIL_2012.loc['2012-11', :])
UTIL_mon.append(UTIL_2012.loc['2012-12', :])
#2013
UTIL_mon.append(UTIL_2013.loc['2013-01', :])
UTIL_mon.append(UTIL_2013.loc['2013-02', :])
UTIL_mon.append(UTIL_2013.loc['2013-03', :])
UTIL_mon.append(UTIL_2013.loc['2013-04', :])
UTIL_mon.append(UTIL_2013.loc['2013-05', :])
UTIL_mon.append(UTIL_2013.loc['2013-06', :])
UTIL_mon.append(UTIL_2013.loc['2013-07', :])
UTIL_mon.append(UTIL_2013.loc['2013-08', :])
UTIL_mon.append(UTIL_2013.loc['2013-09', :])
UTIL_mon.append(UTIL_2013.loc['2013-10', :])
UTIL_mon.append(UTIL_2013.loc['2013-11', :])
UTIL_mon.append(UTIL_2013.loc['2013-12', :])
#2014
UTIL_mon.append(UTIL_2014.loc['2014-01', :])
UTIL_mon.append(UTIL_2014.loc['2014-02', :])
UTIL_mon.append(UTIL_2014.loc['2014-03', :])
UTIL_mon.append(UTIL_2014.loc['2014-04', :])
UTIL_mon.append(UTIL_2014.loc['2014-05', :])
UTIL_mon.append(UTIL_2014.loc['2014-06', :])
UTIL_mon.append(UTIL_2014.loc['2014-07', :])
UTIL_mon.append(UTIL_2014.loc['2014-08', :])
UTIL_mon.append(UTIL_2014.loc['2014-09', :])
UTIL_mon.append(UTIL_2014.loc['2014-10', :])
UTIL_mon.append(UTIL_2014.loc['2014-11', :])
UTIL_mon.append(UTIL_2014.loc['2014-12', :])
#2015
UTIL_mon.append(UTIL_2015.loc['2015-01', :])
UTIL_mon.append(UTIL_2015.loc['2015-02', :])
UTIL_mon.append(UTIL_2015.loc['2015-03', :])
UTIL_mon.append(UTIL_2015.loc['2015-04', :])
UTIL_mon.append(UTIL_2015.loc['2015-05', :])
UTIL_mon.append(UTIL_2015.loc['2015-06', :])
UTIL_mon.append(UTIL_2015.loc['2015-07', :])
UTIL_mon.append(UTIL_2015.loc['2015-08', :])
UTIL_mon.append(UTIL_2015.loc['2015-09', :])
UTIL_mon.append(UTIL_2015.loc['2015-10', :])
UTIL_mon.append(UTIL_2015.loc['2015-11', :])
UTIL_mon.append(UTIL_2015.loc['2015-12', :])
#2016
UTIL_mon.append(UTIL_2016.loc['2016-01', :])
UTIL_mon.append(UTIL_2016.loc['2016-02', :])
UTIL_mon.append(UTIL_2016.loc['2016-03', :])
UTIL_mon.append(UTIL_2016.loc['2016-04', :])
UTIL_mon.append(UTIL_2016.loc['2016-05', :])
UTIL_mon.append(UTIL_2016.loc['2016-06', :])
UTIL_mon.append(UTIL_2016.loc['2016-07', :])
UTIL_mon.append(UTIL_2016.loc['2016-08', :])
UTIL_mon.append(UTIL_2016.loc['2016-09', :])
UTIL_mon.append(UTIL_2016.loc['2016-10', :])
UTIL_mon.append(UTIL_2016.loc['2016-11', :])
UTIL_mon.append(UTIL_2016.loc['2016-12', :])
#2017
UTIL_mon.append(UTIL_2017.loc['2017-01', :])
UTIL_mon.append(UTIL_2017.loc['2017-02', :])
UTIL_mon.append(UTIL_2017.loc['2017-03', :])
UTIL_mon.append(UTIL_2017.loc['2017-04', :])
UTIL_mon.append(UTIL_2017.loc['2017-05', :])
UTIL_mon.append(UTIL_2017.loc['2017-06', :])
UTIL_mon.append(UTIL_2017.loc['2017-07', :])
UTIL_mon.append(UTIL_2017.loc['2017-08', :])
UTIL_mon.append(UTIL_2017.loc['2017-09', :])
UTIL_mon.append(UTIL_2017.loc['2017-10', :])
UTIL_mon.append(UTIL_2017.loc['2017-11', :])
UTIL_mon.append(UTIL_2017.loc['2017-12', :])
#2018
UTIL_mon.append(UTIL_2018.loc['2018-01', :])
UTIL_mon.append(UTIL_2018.loc['2018-02', :])
UTIL_mon.append(UTIL_2018.loc['2018-03', :])
UTIL_mon.append(UTIL_2018.loc['2018-04', :])
UTIL_mon.append(UTIL_2018.loc['2018-05', :])
UTIL_mon.append(UTIL_2018.loc['2018-06', :])
UTIL_mon.append(UTIL_2018.loc['2018-07', :])
UTIL_mon.append(UTIL_2018.loc['2018-08', :])
UTIL_mon.append(UTIL_2018.loc['2018-09', :])
UTIL_mon.append(UTIL_2018.loc['2018-10', :])
UTIL_mon.append(UTIL_2018.loc['2018-11', :])
UTIL_mon.append(UTIL_2018.loc['2018-12', :])
#2019
UTIL_mon.append(UTIL_2019.loc['2019-01', :])
UTIL_mon.append(UTIL_2019.loc['2019-02', :])
#UTIL_mon.append(UTIL_2019.loc['2019-03', :])
UTIL_mon.append(UTIL_2019.loc['2019-04', :])


#####turnover
def turnover_rate(x,y): #x, y are this month and next month data
    n_x = len(x.index)
    n_y = len(y.index)
    mv_x = x['market_value'].sum()
    mv_y = y['market_value'].sum()
    #num = min(len(x.index), len(y.index))
    
    a = [] #previous month bonds cusip
    for i in range(n_x):
        a.append(x.index[i][1])    
    b = [] #current month bonds cusip
    for i in range(n_y):
        b.append(y.index[i][1])
    
    outflow = list(set(a)-set(b))
    num_outflow = len(outflow)
    mv_outflow = 0
    for i in range(num_outflow):
        ip = outflow[i]
        mv = x.xs(ip, level = 'cusip')['market_value'][0]
        mv_outflow += mv
    
    inflow = list(set(b)-set(a)) 
    num_inflow  = len(inflow)
    mv_inflow = 0
    for i in range(num_inflow):
        ip = inflow[i]
        mv = y.xs(ip, level = 'cusip')['market_value'][0]
        mv_inflow += mv
        
    #total = num_inflow + num_outflow     
    return [num_outflow, num_inflow, round(mv_outflow/mv_x, 4), round(mv_inflow/mv_y, 4)]

diff = []    
for i in range(len(UTIL_mon)-1):
    rate = turnover_rate(UTIL_mon[i], UTIL_mon[i+1])
    text = str(UTIL_mon[i].index[0][0]) + ' ' + str(UTIL_mon[i+1].index[0][0])
    rate.insert(0,text)
    diff.append(rate)

'''
df = pd.DataFrame(diff)
df.to_csv('Turnover_Utility.csv', sep='\t')
'''

