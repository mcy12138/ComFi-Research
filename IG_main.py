

import pandas as pd
import IG_fun_library as lib
import math
import numpy as np

IG = pd.read_csv('data/IG_UPDATED_7-1.csv')
IG['date'] = pd.to_datetime(IG['prd_formatted'])
IG = IG.set_index(['date'])
IG = IG.sort_index()
start_yr = IG.index.year[0]
end_yr = IG.index.year[-1]
num_of_yrs = end_yr - start_yr + 1

ind_yr, ind_mon, ind_tot_mv = lib.buckets('INDUSTRIAL', IG, 
                                          start_yr, num_of_yrs)
util_yr, util_mon, util_tot_mv = lib.buckets('UTILITY', IG,
                                             start_yr, num_of_yrs)
fin_yr, fin_mon, fin_tot_mv = lib.buckets('FINANCIAL_INSTITUTIONS', IG,
                                          start_yr, num_of_yrs)


# compute the total market value of all the bonds in the index for each month
index_mv = []
for i in range(len(ind_tot_mv)):
    cur = 0
    cur += ind_tot_mv[i]
    cur += util_tot_mv[i]
    cur += fin_tot_mv[i]
    index_mv.append(cur)
    
n = 12
total_num = 320
to_rate = 0.10

ind_map = lib.cluster(ind_yr, ind_mon, n, index_mv, total_num, to_rate)
util_map = lib.cluster(util_yr, util_mon, n, index_mv, total_num, to_rate)
fin_map = lib.cluster(fin_yr, fin_mon, n, index_mv, total_num, to_rate)

ind_pf, ind_in, ind_out, ind_mat = lib.get_portfolio(ind_map, n, index_mv, total_num, to_rate)
util_pf, util_in, util_out, util_mat = lib.get_portfolio(util_map, n, index_mv, total_num, to_rate)
fin_pf, fin_in, fin_out, fin_mat = lib.get_portfolio(fin_map, n, index_mv, total_num, to_rate)

# compute the total number of bonds in each month's portfolio
# compute the total return for each month's portfolio
amt = [] # actual number of bonds we choose each month
for i in range(len(ind_pf)):
    amt.append(len(util_pf[i].index)+len(ind_pf[i].index)+len(fin_pf[i].index))
 

def cluster_mv(indMap):# clusters mv of by 149 months
    cluster_mv = []
    for i in range(len(indMap)):
        df = indMap[i]
        sum_mv = []
        for j in range(n):
            tmp = df.loc[df['cluster'] == j]
            sum_mv.append(tmp['market_value'].sum())
        cluster_mv.append(sum_mv)
    return cluster_mv

ind_cluster_mv = cluster_mv(ind_map)   
fin_cluster_mv = cluster_mv(fin_map)
util_cluster_mv = cluster_mv(util_map)

def portfolio_cluster_mv(indPf):
    cluster_mv = []
    for i in range(len(indPf)):
        df = indPf[i]
        sum_mv = [0]*n
        for index, row in df.iterrows():
            cluster_id = int(row['cluster'])
            sum_mv[cluster_id] += row['market_value']
        cluster_mv.append(sum_mv)
    return cluster_mv

ind_portfolio_cluster_mv = portfolio_cluster_mv(ind_pf)
fin_portfolio_cluster_mv = portfolio_cluster_mv(fin_pf)
util_portfolio_cluster_mv = portfolio_cluster_mv(util_pf)

    
def eq_weight(col):
    eq_weight_ret = []
    for i in range(len(ind_pf)):
        num = len(ind_pf[i]) + len(util_pf[i]) + len(fin_pf[i])
        if len(ind_pf[i]) > 0:
            ind_ret = ind_pf[i][col].sum()
        if len(util_pf[i]) > 0:
            util_ret = util_pf[i][col].sum()
        if len(fin_pf[i]) > 0:
            fin_ret = fin_pf[i][col].sum()
        ret_sum = ind_ret + util_ret + fin_ret
        eq_ret = ret_sum / num
        eq_weight_ret.append(eq_ret)
    return eq_weight_ret

eq_weight_ret = eq_weight('total_return_mtd')
eq_weight_oad = eq_weight('oad')
eq_weight_oas = eq_weight('oas')
eq_weight_credit = eq_weight('sp_rating_number')
eq_weight_mv = eq_weight('market_value') 
    

def turnover_sector(sec_in, sec_out):
    sec_in = pd.Series(sec_in)
    sec_out = pd.Series(sec_out)
    sec_to = sec_in.add(sec_out)
    return sec_to

def turnover_all(ind_to, util_to, fin_to):
    all_to = (ind_to.add(fin_to)).add(util_to)
    return all_to

ind_to_1 = lib.turnover_sector(ind_in, ind_out)
util_to_1 = lib.turnover_sector(util_in, util_out)
fin_to_1 = lib.turnover_sector(fin_in, fin_out)

turnover_1 = lib.turnover_all(ind_to_1, util_to_1, fin_to_1)
turnover_rate = turnover_1 / pd.Series(amt)
avg_turnover = np.mean(turnover_rate)

#annual turnover; sum all
t = turnover_rate[1:]
annual_turnover= [sum(t[i:i+12]) for i in range(0, len(t), 12)]
annual_turnover.insert(0, turnover_rate[0])


def annualize(T): #take the AVERAGE of annual data
    t = T[1:]
    annual_t = [np.mean(t[i:i+12]) for i in range(0, len(t), 12)]
    annual_t.insert(0, T[0])
    return annual_t

annual_dur = annualize(eq_weight_oad)
annual_spr = annualize(eq_weight_oas)
annual_credit = annualize(eq_weight_credit)
annual_mv = annualize(eq_weight_mv)


# return after accounting for the assumed 50 bps transaction costs
# !!!
for i in range(len(eq_weight_ret)):
    eq_weight_ret[i] = eq_weight_ret[i] * (1-turnover_rate[i])
    + (eq_weight_ret[i] - 0.5) * turnover_rate[i]
    
annual_ret = annualize(eq_weight_ret)

baseline = pd.read_csv('data/Baseline.csv')
baseline_IG = baseline.loc[baseline['universe'] == 'Investment Grade']
baseline_IG['date'] = pd.to_datetime(baseline['prd_formatted'])
baseline_IG = baseline_IG.set_index('date')
baseline_IG = baseline_IG.sort_index()
benchmark_ret = baseline_IG['benchmark'].tolist()

def tracking_error_1(P, B): # portfolio and benchmark
    error_square_sum = 0
    N = 0
    for i in range(len(B)):
        if not pd.isna(P[i]):
            error_square_sum += (P[i] - B[i]) ** 2
            N += 1
    TE = math.sqrt(error_square_sum / (N-1))
    return TE

TE_1 = tracking_error_1(eq_weight_ret, benchmark_ret)
    
wr1 = eq_weight_ret[25:]
br1 = benchmark_ret[25:]
TE_1_2010 = tracking_error_1(wr1, br1) # 2010-present
annual_TE = TE_1_2010*np.sqrt(12)

def tracking_error_2(P, B): # portfolio and benchmark
    error = []
    N = 0
    for i in range(len(B)):
        if not pd.isna(P[i]):
            error.append(P[i] - B[i])
            N += 1
    err = pd.Series(error)
    avg = err.mean()
    err_sum_sq = ((err-avg)**2).sum()
    TE = math.sqrt(err_sum_sq / (N))
    return TE

TE_2 = tracking_error_2(eq_weight_ret, benchmark_ret) #2007-2019

wr2 = eq_weight_ret[25:]
br2 = benchmark_ret[25:]
TE_2_2010 = tracking_error_2(wr2, br2) # 2010-present

ann = pd.DataFrame() #annualized data for 2007-2019
ann['credit'] = annual_credit # sp rating number
ann['dur'] = annual_dur
ann['mv'] = annual_mv
ann['ret'] = annual_ret #return
ann['oas'] = annual_spr
ann['to'] = annual_turnover
ann.to_excel('IG_Annualized.xlsx')

result = pd.DataFrame()
result['oad_pf'] = eq_weight_oad
result['oas_pf'] = eq_weight_oas
result['credit_pf'] = eq_weight_credit
result['market_value'] = eq_weight_mv
result['ret_pf'] = eq_weight_ret
result['turnover'] = turnover_rate
result.to_excel('IG_Result_pf.xlsx')



