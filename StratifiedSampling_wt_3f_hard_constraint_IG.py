# -*- coding: utf-8 -*-
import pandas as pd
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import matplotlib.pyplot as plt
from dateutil import parser
import numpy as np
import datetime

IG = pd.read_csv('IG_UPDATED_7-1.csv')
base = pd.read_csv('BaselineIG.csv')
baseAvg = pd.read_csv("IG_Index_Averages_7-8.csv")
'''
Clearn up data: add Date column
'''
IG['Date'] = pd.to_datetime(IG['prd_formatted']).dt.to_period('D')
base['Date'] = pd.to_datetime(base['prd_formatted']).dt.to_period('D')
base.drop('prd_formatted', axis = 1)
base = base.set_index(['Date'], drop = False)

baseAvg['Date'] = pd.to_datetime(baseAvg['prd_formatted']).dt.to_period('D')
baseAvg.drop('prd_formatted', axis = 1)
baseAvg = baseAvg.set_index(['Date'], drop = False)
 
'''
IG_by_month = a list of dataframes of each month
'''
IG_by_month = []
for _, g in IG.groupby(IG['Date']):
    IG_by_month.append(g)
    
'''
getQuants:
    Input: df = dataframe, factor = a factor by which the dataframe will be splitted
    Output: A list of dataframes split into tertiles if the factor is quantitative,
            or the number of categories if the factor is qualitative.
'''    
def getCutoffs(df, factor):
    quantiles = pd.qcut(df[factor], 3)    
    valCount = list((quantiles.value_counts(sort = False)).index)
    C1 = valCount[0].right
    C2 = valCount[1].right
    return (C1, C2)

def getQuants(df, factor):
    if is_numeric_dtype(df[factor]):
        (C1, C2) = getCutoffs(df, factor)
        Q1 = df.loc[(df[factor] <= C1)]
        Q2 = df.loc[(df[factor] <= C2) & (df[factor] > C1)]
        Q3 = df.loc[(df[factor] > C2)]
        res = [Q1, Q2, Q3]
    else:
        res = []
        for _, g in df.groupby(df[factor]):
            res.append(g)
    return res

'''
getBuckets:
    Inputs: df = dataframe, f1, f2, f3 =  3 factors of the stratified sampling in order
    Outputs: A 2-D list with each sub-list contains the buckets for each month as dataframes.
             len(Buckets_by_month) = # of months in the original dataframe
             len(Buckets_by_month[i]) = # of buckets for the ith month
'''
def getBuckets(df, f1, f2, f3):
  Buckets_by_month = []
    
  for month in range(len(df)): 
    Buckets_each_month = []
    d0 = df[month].dropna(subset=[f1, f2, f3])#drop rows that miss values for the factors
    #split by factor 1
    if len(d0) < 3: #if the df has fewer than 3 rows, append the bucket, and stop further splitting
        Buckets_each_month.append(d0)
        continue
    
    dataByF1 = getQuants(d0, f1)
    for i in range(len(dataByF1)):
      d1 = dataByF1[i]
      if len(d1) < 3:
        Buckets_each_month.append(d1)
        continue
            
      #split by factor 2
      dataByF2 = getQuants(d1, f2)    
      for j in range(len(dataByF2)):
        d2 = dataByF2[j]
        if len(d2) < 3:
          Buckets_each_month.append(d2)
          continue
                
        #split by factor 3
        dataByF3 = getQuants(d2, f3)
        for k in range(len(dataByF3)):
            d3 = dataByF3[k]
            Buckets_each_month.append(d3)
            
    Buckets_by_month.append(Buckets_each_month)
  return Buckets_by_month                              

'''
labelBuckets:
    This function takes the output from getBuckets and add a column that 
    records the bucket index for each bond.
    (This function is slow, taking about 3 minutes)
'''
def labelBuckets(Buckets_by_month):
    for month in range(len(Buckets_by_month)):
        for b in range(len(Buckets_by_month[month])):
            Buckets_by_month[month][b]['bucket'] = b
    return Buckets_by_month

'''
getNumBonds:
    Input: data = a list of buckets for one month, numBond = the avg. # of bonds per bucket
    Output: A list of the number of bonds for each bucket.
'''
def getNumBonds(data, numBond):
    numBuckets = len(data)
    tt_numNeeded = numBond * numBuckets
    tt_mv = 0
    numBondsL = pd.DataFrame(columns = ['market_value'])
    for b in range(numBuckets):
        bucket = data[b]
        mv = bucket['market_value'].sum()
        numBondsL.loc[b] = mv
        tt_mv += mv

    numBondsL['market_value_perc'] = numBondsL['market_value'].div(tt_mv)
    numBondsL['numBonds'] = numBondsL['market_value_perc'] * tt_numNeeded
    numBondsL['numBonds'] = numBondsL['numBonds'].astype(int)
    return numBondsL

'''
getPorfolio:
    Input: df = A 2-D list of buckets for each month (Output from getBuckets)
           numBonds = avg # of bonds from each bucket
    Output: A 2-D list where each sublist has the bonds chosen for that month
    This function returns the desired portfolio of each month.
'''
def getPortfolio(df, numBonds):
    tt_portfolio = [] #A list of portfolio for each month
    for month in range(len(df)):
        data = df[month]
        numBuckets = len(data)
        tt_numWanted = numBonds * numBuckets
        if numBuckets == 0: #Edge case: has no bonds/buckets in that month
            continue
        
        if len(tt_portfolio) == 0: #address the starting case
            portfolio_prev = pd.DataFrame()
        else:         
            portfolio_prev = tt_portfolio[month - 1]
        
        numBondsL = getNumBonds(data, numBonds) #a list that contains # of bonds to choose from each bucket
        portfolio = pd.DataFrame() #a dataframe of bonds for current month
        mv_bucket_max = 0
        for b in range(numBuckets): 
            bucket = data[b]
            bucket_curr = bucket #bond candidates in this bucket
            mv_bucket = bucket_curr['market_value'].sum()
            numNeeded = int(numBondsL.loc[b]['numBonds'])

            #sort bonds from previous buckets by market value (decreasing)
            #keep bonds from previous buckets that remain in current bucket
            #buy additional bonds if too few bonds remain
            #sell additional bonds if too many bonds remain
            if len(portfolio_prev) != 0: 
                
                # STEP 1: keep bonds from previous buckets that remain in current bucket
                portfolio_per_bucket_prev = portfolio_prev.loc[(portfolio_prev['bucket'] == b)]
                staying = list(set(portfolio_per_bucket_prev['cusip']).intersection(set(bucket_curr['cusip'])))
                
                portfolio_per_bucket_prev = portfolio_per_bucket_prev.set_index(['cusip'], drop = False)
                bucket_curr = bucket_curr.set_index(['cusip'], drop = False)
                staying_in_portfolio = bucket_curr.loc[staying]
                staying_in_portfolio = staying_in_portfolio.sort_values(by = ['market_value'], ascending = False)
                
                if len(staying_in_portfolio) > numNeeded:
                    #if number remained > number needed, buy number of needed bonds
                    portfolio = portfolio.append(staying_in_portfolio[:numNeeded])
                    continue
                #if not enough bonds remained from last month, append all staying bonds
                portfolio = portfolio.append(staying_in_portfolio)
                bucket_curr = bucket_curr.drop(staying, axis = 0)
                portfolio_per_bucket_prev = portfolio_per_bucket_prev.drop(staying, axis = 0)

                numNeeded -= len(staying)   
                
            # STEP 2: buy additional bonds if too few bonds remain
            if len(bucket_curr) < numNeeded: #not enough bonds to pick
                portfolio = portfolio.append(bucket_curr)
            else:
                bucket_curr = bucket_curr.sort_values(by = ['market_value'], ascending = False)
                portfolio = portfolio.append(bucket_curr.iloc[:numNeeded])
                bucket_curr = bucket_curr[numNeeded:]
            
            if mv_bucket > mv_bucket_max:
                mv_bucket_max = mv_bucket
                maxBucket = bucket_curr 
                
        #buy additional bonds from largest MV bucket to keep numBonds fixed
        if len(portfolio) < tt_numWanted:
            numNeed = tt_numWanted - len(portfolio)
            portfolio = portfolio.append(maxBucket.iloc[:numNeed]) 
        
        tt_portfolio.append(portfolio)
    return tt_portfolio

'''
getPorfolioActual:
    Input: data = IG_by_month, tt_portfolio = the desired portfolio (output from
           getPortfolio), turnover = the turnover constraint
    Output: A 2-D list where each sublist has the bonds chosen for that month,
           and a list of # of changing bonds for each month
    This function takes in the desired portfolio and returns the actual portfolio
    under the turnover constraint.
'''
def getPortfolioActual(data, tt_portfolio, turnover):
    tt_port_final = []
    tt_turnover = pd.DataFrame(columns = ['Date', 'Inflow', 'Outflow', 'Same'])
    
    for month in range(len(tt_portfolio)):
        portfolio = tt_portfolio[month]
        tt_numWanted = len(portfolio)
        if month == 0:
            portfolio_actual = portfolio
            Inflow = 0
            Outflow = 0
            Same = 0
        else:
            #Initialize
            Inflow = 0
            Outflow = 0
            Same = 0
            df_month = data[month]
            portfolio_prev = tt_portfolio[month - 1]
            portfolio_actual = pd.DataFrame()
            numTransactions = int(turnover * tt_numWanted / 2)
            
            # STEP 1: keep remaining bonds
            staying = list(set(portfolio_prev['cusip']).intersection(set(portfolio['cusip'])))
            Same += len(staying)

            portfolio = portfolio.set_index(['cusip'], drop = False)
            portfolio_prev = portfolio_prev.set_index(['cusip'], drop = False)
            staying_in_portfolio = portfolio.loc[staying]

            portfolio_actual = portfolio_actual.append(staying_in_portfolio)
            portfolio = portfolio.drop(staying, axis = 0) 
            portfolio_prev = portfolio_prev.drop(staying, axis = 0)
            
            # STEP 2: Buy bonds & Sell bonds by # of numTransactions
            if len(portfolio) <= numTransactions: #the # to buy < # transactions
                print("happen", len(portfolio))
                portfolio_actual.append(portfolio)
                Inflow += len(portfolio)
            else:
                portfolio = portfolio.sort_values(by = ['market_value'], ascending = False)
                portfolio_prev = portfolio_prev.sort_values(by = ['market_value'], ascending = True)
                
                portfolio = portfolio.reset_index(drop = True)
                portfolio_prev = portfolio_prev.reset_index(drop = True)
                
                #buy maxValue bonds
                portfolio_actual = portfolio_actual.append(portfolio[:numTransactions]) 
                Inflow += numTransactions
                portfolio = portfolio[numTransactions:]
                
                #sell minValue bonds from prev
                sold = list(set(portfolio_prev['cusip'][:numTransactions]).intersection(set(df_month['cusip'])))
                Outflow += len(sold)
               
                portfolio_prev = portfolio_prev[numTransactions:]
                
            # STEP 3: Get additional bonds from the remains of previous month's portfolio
                i = 0
                while len(portfolio_actual) < tt_numWanted and i < len(portfolio_prev):
                    portfolio_prev = portfolio_prev.reset_index(drop = True)
                    if portfolio_prev['cusip'][i] in df_month['cusip'].values: #append bonds that did not mature
                        row = np.where(df_month['cusip'] == portfolio_prev['cusip'][i])[0][0]
                        portfolio_actual = portfolio_actual.append(df_month.iloc[row,:])
                        Same += 1
                    i += 1
                   
                if len(portfolio_actual) < tt_numWanted: #too many bonds mature
                    numNeed = tt_numWanted - len(portfolio_actual)
                    portfolio_actual = portfolio_actual.append(portfolio[:numNeed])
                    Inflow += numNeed
        
        date = portfolio_actual.iloc[0]['prd_formatted']
        tt_turnover.loc[month] = [date, Inflow, Outflow, Same]                
        tt_port_final.append(portfolio_actual)
    return tt_port_final, tt_turnover

'''
getReturn:
    Input: A 2-D list where each sublist has the bonds chosen for that month
           (Output from getPortfolioActual)
    Output: A dataframe with col1 is date and col2 is the return for that month
'''
def getReturn(df): #Unweighted return
    avg_returnL = pd.DataFrame(columns = ['Date', 'Return'])
    
    for month in range(len(df)):
        portfolio = df[month]
        return_wt = portfolio['total_return_mtd'] * 0.9 + (portfolio['total_return_mtd'] - 0.5) * 0.1
        tt_return = return_wt.sum()
        avg_return = tt_return / len(portfolio)
        date = portfolio.iloc[0]['prd_formatted']
        avg_returnL.loc[month] = [date, avg_return]    

    avg_returnL['Date'] = pd.to_datetime(avg_returnL['Date']).dt.to_period('D')
    avg_returnL= avg_returnL.set_index(['Date'])
    return avg_returnL
'''
getError:
    This function takes in the portfolio dataframe and benchmark dataframe and
    returns a tuple of dataframe with the error and return for each month, the 
    tracking error calculatred as std., and the tracking error calculated in the
    right way.
'''
def getError(portfolio, base):
    myReturn = getReturn(portfolio)
    dfinal = myReturn.join(base)
    dfinal['Error'] = abs(dfinal['Return'] - dfinal['benchmark'])
    errorMean = dfinal['Error'].mean()
    errorSqW = ((dfinal['Error'] - errorMean)**2).sum()
    errorW = (errorSqW / len(dfinal['Error'])) ** (1/2)
    
    errorSq = (dfinal['Error']**2).sum()
    errorR = (errorSq / (len(dfinal['Error']) - 1)) **(1/2)
    return (dfinal, errorW, errorR)

'''
getTurnoverActual:
    Input: tt_turnover = a list of # of changing bonds of each month
           tt_turnover and portfolio are outputs from getPortfolioActual
    Output: A list of turnover for each month
'''
def getTurnoverActual(tt_turnover, portfolio):
    turnoverL = pd.DataFrame(columns = ['Date', 'Turnover'])
    for month in range(1, len(tt_turnover)):
        tt_numWanted = len(portfolio[month])
        incoming = tt_turnover['Inflow'][month]
        outgoing = tt_turnover['Outflow'][month]
        turnover = (incoming + outgoing) / tt_numWanted
        date = tt_turnover['Date'][month]
        turnoverL.loc[month] = [date, turnover]
    return turnoverL

'''
getAnnTurnover:
    Input: a list of turnover for each month
    Output: the average of annual turnover
'''
def getAnnTurnover(turnoverL):
    turnoverL['Year'] = pd.to_datetime(turnoverL['Date']).dt.to_period('A')
    turn_avg = []
    for y, g in turnoverL.groupby(turnoverL['Year']):
        turn_avg.append(g['Turnover'].sum())
    turn_avg_annual = sum(turn_avg) / len(turn_avg)
    return turn_avg_annual

def writeFile(path, contents):
    with open(path, "wt") as f:
        f.write(contents)    

def ss(df, base, f1, f2, f3, numBond, numBucket, turnover): 
    Buckets_by_month = getBuckets(df, f1, f2, f3)
    print("Buckets done")
    Buckets_by_month = labelBuckets(Buckets_by_month)
    print('labeling done')
    
    Port1 = getPortfolio(Buckets_by_month, numBond)
    numBonds = numBond * numBucket
    Port1A, turnover1 = getPortfolioActual(df, Port1, turnover)
    print("portfolio done")
    
    Error1, error1W, error1R = getError(Port1A, base)
    Error2, error2W, error2R = getError(Port1A[36:], base[36:])
    Turnover1 = getTurnoverActual(turnover1, Port1A) 
    Turnover1A = getAnnTurnover(Turnover1) 
    
    contents = str(numBucket) + " Buckets * " + str(numBond) + " Bond; " + "Turnover Constraint: " + str(turnover) + "\n"
    contents += f1 + " " + f2 + " " + f3 + "\n"
    contents += "2007 - 2019: " + "Error: " + str(error1R) + " Error(sd): " + str(error1W) + "\n"
    contents += "2010 - 2019: " + "Error: " + str(error2R) + " Error(sd): " + str(error2W) + "\n"
    contents += "Turnover: Min:" + str(Turnover1['Turnover'].min()) + " Median: " + str(Turnover1['Turnover'].median()) + " Max: " + str(Turnover1['Turnover'].max()) + "\n"
    contents += "Annualized Turnover: " + str(Turnover1A) + "\n"
    
    writeFile('Results_wt_IG_%.2f.txt' %(turnover), contents)
    return Port1A, Error1, Turnover1

#calling main function
Port, Error, Turnover = ss(IG_by_month, base, 'coupon', 'oad', 'oas', 11, 27, 0.04)

#Writing output files: 
Error.to_csv('TrackingError_wt_IG.csv')
Turnover.to_csv('Turnover_wt_IG.csv')
pd.concat(Port).to_csv('Portfolio_wt_IG.csv')

'''
getAvg: This function returns a dataframe with monthly average of a variable of
        both porfolio and benchark, a list of annual average for portfolio, and
        a list of annual average for benchmark.
'''
def getVar(df, var): #Unweighted return
    avg_L = pd.DataFrame(columns = ['Date', var])
    
    for month in range(len(df)):
        portfolio = df[month]
        tt = portfolio[var].sum()
        avg = tt / len(portfolio)
        date = portfolio.iloc[0]['prd_formatted']
        avg_L.loc[month] = [date, avg]    
    #avg = avg_returnL['return'].mean()
    avg_L['Date'] = pd.to_datetime(avg_L['Date']).dt.to_period('D')
    avg_L= avg_L.set_index(['Date'], drop = False)
    return avg_L

def getAvg(portfolio, base, var, baseVar):
    myVar = getVar(portfolio, var)
    dfinal = myVar.join(base[[baseVar]])
    dfinal['Year'] = dfinal['Date'].dt.year
    avgL = dict()
    baseAvgL = dict()
    for y, g in dfinal.groupby(dfinal['Year']):
        avgL[y]=g[var].mean()
        baseAvgL[y]=g[baseVar].mean()
    return dfinal, avgL, baseAvgL

OAD, oad_avg, oad_avg_b = getAvg(Port, baseAvg, 'oad', 'duration_average')
OAS, oas_avg, oas_avg_b = getAvg(Port, baseAvg, 'oas', 'spread_average')
Rating, rat_avg, rat_avg_b = getAvg(Port, baseAvg, 'sp_rating_number', 'credit_average')
MV, mv_avg, mv_avg_b = getAvg(Port, baseAvg, 'market_value', 'market_value_average')

Error['Year'] = Error['Date'].dt.year
ret_avg = dict()
ret_avg_b = dict()
for y, g in Error.groupby(Error['Year']):
        ret_avg[y]=12*(g['Return'].mean())
        ret_avg_b[y]=12*(g['benchmark'].mean())

OAD.to_csv('OAD_wt_IG.csv')
OAS.to_csv('OAS_wt_IG.csv')
Rating.to_csv('Rating_wt_IG.csv')
MV.to_csv('MV_wt_IG.csv')   

#annualRes is a dataframe of annual results of portfolio and benchmark oad, oas, mv, rating, return.
annualRes = pd.DataFrame([oad_avg, oad_avg_b, oas_avg, oas_avg_b,  rat_avg, rat_avg_b, mv_avg, mv_avg_b, ret_avg, ret_avg_b]).transpose()
annualRes.columns = ['oad_p', 'oad', 'oas_p', 'oas', 'rating_p', 'rating', 'mv_p', 'mv', 'return_p', 'return']
annualRes.to_csv('Annual_Results_wt_IG.csv')   
