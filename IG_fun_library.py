
import pandas as pd
from sklearn.cluster import KMeans

def buckets(sector, IG, start_yr, num_of_yrs):
    IG_yr = [] # a list of dataframes for each year
    sector_yr = [] # a list of dataframes for each year's INDUSity sector
    sector_mon = [] # a 2d list, row is year, col is month    
    sector_mv = [] # market value of the sector each month   
    for i in range(num_of_yrs):
        year = start_yr+i
        IG_cur_yr = IG.loc['%d' %year, ['class_2', 'cusip', 'total_return_volatility', 
                                 'market_value', 'prev_6m_total_return', 'coupon', 
                                 'oas', 'oad', 'total_return_mtd', 'sp_rating_number']]
        IG_yr.append(IG_cur_yr)
        sector_cur_yr = IG_cur_yr.loc[IG_cur_yr['class_2'] == sector]

        #change the index
        sector_cur_yr = sector_cur_yr.reset_index()
        # del sector_cur_yr['class_2']
        sector_cur_yr = sector_cur_yr.set_index(['date', 'cusip'], drop=False)
        sector_yr.append(sector_cur_yr)
        
        for _, g in sector_cur_yr.groupby(sector_cur_yr['date']):
            del g['class_2']
            del g['date']
            del g['cusip']
            sector_mon.append(g)
            cur_mv = g.loc[:, 'market_value'].sum()
            sector_mv.append(cur_mv)
            
    return (sector_yr, sector_mon, sector_mv)


def cluster(sector_yr, sector_mon, n, index_mv, total_num, to_rate):
# n is the number of clusters
# index_mv is the list of total market value of all the bonds in the index for each moth
# total_num is the total number of bonds we want to be in monthly portfolio
    for j in range(len(sector_mon)):
        D = sector_mon[j]

        for factor in list(D):
            col = D[factor]
            norm = (col - col.mean()) / col.std()
            sector_mon[j][factor] = norm
        #INDUS_mon[j] = INDUS_mon[j].dropna()
        sector_mon[j] = sector_mon[j].fillna(0) # fill all NAs with 0
    
    # n = 12
    # !!! Performing kmeans clustering
    clusters = []
    for D in sector_mon:
        km = KMeans(n_clusters=n, random_state = 170).fit(D)
        clusters.append(km)
    
    cluster_map = []    
    
    for i in range(len(sector_mon)):
        df = pd.DataFrame()
        df['date'] = sector_mon[i].index.get_level_values(0)
        df['cusip'] = sector_mon[i].index.get_level_values(1)
        df['cluster'] = clusters[i].labels_
        cluster_map.append(df)
    
    
    
    for i in range(13): # line 63 in backet_clustering
        year = 2007+i
        cur_yr = sector_yr[i]
        
        if year == 2019:
            for j in range(6):
                idx = (i-1) * 12 + j + 1
                D = cur_yr.loc['%d-0%d' %(year, j+1), ['market_value', 'total_return_mtd', 'oad', 'oas', 'sp_rating_number']].reset_index()
                cluster_map[idx]['market_value'] = D['market_value']
                cluster_map[idx]['total_return_mtd'] = D['total_return_mtd']
                cluster_map[idx]['oad'] = D['oad']
                cluster_map[idx]['oas'] = D['oas']
                cluster_map[idx]['sp_rating_number'] = D['sp_rating_number']
    
        else:
            for j in range(12):
                idx = (i-1) * 12 + j + 1
                if j < 9:
                    D = cur_yr.loc['%d-0%d' %(year, j+1), ['market_value', 'total_return_mtd', 'oad', 'oas', 'sp_rating_number']].reset_index()
                else:
                    D = cur_yr.loc['%d-%d' %(year, j+1), ['market_value', 'total_return_mtd', 'oad', 'oas', 'sp_rating_number']].reset_index()
                cluster_map[idx]['market_value'] = D['market_value']
                cluster_map[idx]['total_return_mtd'] = D['total_return_mtd']
                cluster_map[idx]['oad'] = D['oad']
                cluster_map[idx]['oas'] = D['oas']
                cluster_map[idx]['sp_rating_number'] = D['sp_rating_number']
                
    return cluster_map

def get_exp_portfolio(curmap, prev_pf, n, cur_index_mv, total_num, to_rate):
    if prev_pf.empty == False:
        prev_pf_cusip = prev_pf.loc[:, 'cusip'].tolist()
    mv_wt = (curmap['market_value'].sum()) / cur_index_mv
    new_pf = pd.DataFrame()
    for j in range(n):
        C = curmap.loc[curmap['cluster'] == j]
        C_cusip = C.loc[:, 'cusip'].tolist()
        C = C.set_index('cusip', drop=False)
        C = C.sort_values(by='market_value', ascending=False)
        cluster_mv = C['market_value'].sum()
        cluster_wt = cluster_mv / cur_index_mv
        # compute the number of bonds we want from this cluster
        N = int(cluster_wt * total_num)
        if N == 0: continue
        # simply market value
        new_pf_mv = pd.DataFrame()
        length = len(C.index)
        if length < N:
            new_pf_mv = new_pf_mv.append(C)
        else:
            new_pf_mv = new_pf_mv.append(C.nlargest(N, 'market_value'))
        # try to reduce turnover by comparing with prev_pf
        intersection = list(set(C_cusip) & set(prev_pf_cusip))
        I = len(intersection)
        if I < N:
            new_pf = new_pf.append(C.reindex(intersection))
            
            num_to_buy = N - I
            bonds_diff = C.drop(intersection)
            new_num = min(num_to_buy, len(bonds_diff.index))
            new_bonds = bonds_diff.nlargest(new_num, 'market_value')
            new_pf = new_pf.append(new_bonds.reset_index(drop=True), sort=False)
        elif I > N:
            num_to_sell = I - N
            sell = C.loc[intersection].nsmallest(num_to_sell, 'market_value')
            sell_cusip = sell.loc[:,'cusip'].tolist()
            prev_pf_idx = prev_pf.set_index('cusip', drop=False)
            remaining = prev_pf_idx.drop(sell_cusip)
            remaining_cusip = remaining.loc[:,'cusip'].tolist()
            new_pf = new_pf.append(C.reindex(remaining_cusip).dropna(), sort=False)
        else: 
            new_pf = new_pf.append(C.reindex(intersection))
                
    expected = int(mv_wt * total_num)
    actual = len(new_pf.index)
    if actual < expected:
        more = expected - actual
        curmap_idx = curmap.set_index('cusip', drop=False)
        cur_cusip = curmap.loc[:, 'cusip'].tolist()
        # get the bonds in the previous pf that remains in current month
        remaining = set(cur_cusip) & set(prev_pf_cusip)
        # remove the bonds that are already selected to be in current pf
        new_pf_cusip = set(new_pf.loc[:, 'cusip'].tolist())
        pool_cusip = list(remaining.difference(new_pf_cusip))
        pool = curmap_idx.loc[pool_cusip, :]
        # pool = pool.sort_values(by='market_value', ascending=False)   
        if pool.empty:
            pass
        elif len(pool_cusip) < more:
            new_pf = new_pf.append(pool)
        else:
            new_pf = new_pf.append(pool.nlargest(more, 'market_value'))
    
    return (mv_wt, new_pf_mv, new_pf)

def get_portfolio(cluster_map, n, index_mv, total_num, to_rate):   
    bonds = [pd.DataFrame()]*len(cluster_map) # portfolio list
    inflow = [0]*len(cluster_map)
    outflow = [0]*len(cluster_map)
    matured = [0]*len(cluster_map)
    mv_wt = [0]*len(cluster_map)
    
    # get the portfolio for the first month
    collection = pd.DataFrame()
    all_mv_0 = 0
    for j in range(n): #access cluster
        #print(str(i) + ',' + str(j))
        C = cluster_map[0][cluster_map[0].cluster == j] # cluster 
        cluster_mv = C['market_value'].sum()
        all_mv_0 += cluster_mv
        cluster_wt = cluster_mv / index_mv[0]
        # compute the number of bonds we want from this cluster
        num = int(cluster_wt * total_num)
        # get the portfolio (6 bonds with largest market values)
        length = len(C.index) # total number of bonds in the cluster
        if length < num:
            pf = C
        else:
            pf = C.nlargest(num, 'market_value')
        collection = collection.append(pf)
        
    all_wt_0 = all_mv_0 / index_mv[0]
    mv_wt[0] = (cluster_map[0]['market_value'].sum()) / index_mv[0]
    map_0 = cluster_map[0]
    expected = int(all_wt_0 * total_num)
    actual = len(collection.index)
    if actual < expected:
        more = expected - actual
        map_0_idx = map_0.set_index('cusip', drop=False)
        # get all the bonds in the first month
        cusip_0 = map_0.loc[:, 'cusip'].tolist()
        # get the cusips of the bonds already in pf
        cusip_col = collection.loc[:, 'cusip'].tolist()
        # delete the bonds that are already in pf
        delete = list(set(cusip_0) & set(cusip_col))
        pool = map_0_idx.drop(delete)
        add = pool.nlargest(more, 'market_value')
        collection = collection.append(add)

    bonds[0] = collection
    #print(bonds[0])
    
    # input: portfolio of the previous month as a list of tuples
    #        cluster map of the current month for a specific sector
    # group the bonds into 12 dataframes based on cluster
    # compute the number of bonds we want from each cluster: N
    # compute the intersection of each cluster with previous month's portfolio: I
    # if I < N, buy N-I bonds with largest market value
    # if I = N, remain the same
    # if I > N, sell I-N bonds with smallest market value
    for i in range(1, len(cluster_map)): #access month
        #print('i is: ', i)
        prev_pf = bonds[i-1]
        if prev_pf.empty == False:
            prev_pf_cusip = prev_pf.loc[:, 'cusip'].tolist()
        curmap = cluster_map[i]
        mv_wt[i] = (cluster_map[i]['market_value'].sum()) / index_mv[i]
        new_pf = pd.DataFrame()
        for j in range(n):
            C = curmap.loc[curmap['cluster'] == j]
            C_cusip = C.loc[:, 'cusip'].tolist()
            C = C.set_index('cusip', drop=False)
            C = C.sort_values(by='market_value', ascending=False)
            cluster_mv = C['market_value'].sum()
            cluster_wt = cluster_mv / index_mv[i]
            # compute the number of bonds we want from this cluster
            N = int(cluster_wt * total_num)
            if N == 0: continue
            # simply market value
            new_pf_mv = pd.DataFrame()
            length = len(C.index)
            if length < N:
                new_pf_mv = new_pf_mv.append(C)
            else:
                new_pf_mv = new_pf_mv.append(C.nlargest(N, 'market_value'))
            # try to reduce turnover by comparing with prev_pf
            intersection = list(set(C_cusip) & set(prev_pf_cusip))
            I = len(intersection)
            if I < N:
                new_pf = new_pf.append(C.reindex(intersection))
                
                num_to_buy = N - I
                bonds_diff = C.drop(intersection)
                new_num = min(num_to_buy, len(bonds_diff.index))
                new_bonds = bonds_diff.nlargest(new_num, 'market_value')
                new_pf = new_pf.append(new_bonds.reset_index(drop=True), sort=False)
            elif I > N:
                num_to_sell = I - N
                sell = C.loc[intersection].nsmallest(num_to_sell, 'market_value')
                sell_cusip = sell.loc[:,'cusip'].tolist()
                prev_pf_idx = prev_pf.set_index('cusip', drop=False)
                remaining = prev_pf_idx.drop(sell_cusip)
                remaining_cusip = remaining.loc[:,'cusip'].tolist()
                new_pf = new_pf.append(C.reindex(remaining_cusip).dropna(), sort=False)
            else: 
                new_pf = new_pf.append(C.reindex(intersection))
                    
        expected = int(mv_wt[i] * total_num)
        actual = len(new_pf.index)
        if actual < expected:
            more = expected - actual
            curmap_idx = curmap.set_index('cusip', drop=False)
            cur_cusip = curmap.loc[:, 'cusip'].tolist()
            # get the bonds in the previous pf that remains in current month
            remaining = set(cur_cusip) & set(prev_pf_cusip)
            # remove the bonds that are already selected to be in current pf
            new_pf_cusip = set(new_pf.loc[:, 'cusip'].tolist())
            pool_cusip = list(remaining.difference(new_pf_cusip))
            pool = curmap_idx.loc[pool_cusip, :]
            # pool = pool.sort_values(by='market_value', ascending=False)   
            if pool.empty:
                pass
            elif len(pool_cusip) < more:
                new_pf = new_pf.append(pool)
            else:
                new_pf = new_pf.append(pool.nlargest(more, 'market_value'))
        
        #print('i is:', i)
        curmap_idx = curmap.set_index('cusip', drop=False)
        num_bonds = mv_wt[i] * total_num
        num_of_tran = int(to_rate * num_bonds / 2)
        #transactions[i] = num_of_tran
        new_pf_cusip = set(new_pf.loc[:, 'cusip'].tolist())
        out_cusip = list((set(prev_pf_cusip)).difference(new_pf_cusip))
        in_cusip = list(new_pf_cusip.difference(set(prev_pf_cusip)))
        prev_pf_idx = prev_pf.set_index('cusip', drop=False)
        out_pool = prev_pf_idx.reindex(out_cusip)
        out_pool_cp = out_pool.copy(deep=True)
        #  print('outpool len', len(out_pool))
        new_pf_idx = new_pf.set_index('cusip', drop=False)
        in_pool = new_pf_idx.reindex(in_cusip)
        # print('inpool len', len(in_pool))
        mature = pd.DataFrame()
        in_num = 0
        out_num = 0
        quota = 0
        prev_pf_idx_drop = prev_pf_idx.copy(deep=True)
        # sell the bonds in out_pool with smallest market value
        if len(out_cusip) == 0:
            pass
        else:
            for index, row in out_pool.iterrows():
                if index in curmap_idx.index: 
                    new_mv = curmap_idx.loc[index, 'market_value']
                    out_pool_cp.loc[index, 'market_value'] = new_mv
                else:
                    mature = mature.append(row)
                    # matured bonds are not counted in transactions
                    out_pool_cp = out_pool_cp.drop(index, axis=0)
                    prev_pf_idx_drop = prev_pf_idx_drop.drop(index)
            # count mature bonds as outflow, so the actual number of bonds to sell decreases
            prev_pf_idx_1 = prev_pf_idx_drop.copy(deep=True)
            if num_of_tran - len(mature) <= 0:
                #print('num_of_tran - len(mature) <= 0')
                out_num = 0
            else:
                out_num = num_of_tran - len(mature)
            if out_num == 0: pass
            elif len(out_pool_cp) < out_num:
                #print('short of bonds to SELL')
                quota += out_num - len(out_pool_cp)
                out_num = len(out_pool_cp)
                out_ser = out_pool_cp.loc[:, 'cusip']
                prev_pf_idx_1 = prev_pf_idx_1.drop(out_ser)
            else:
                # out_num = num_of_tran
                out_ser = out_pool_cp.nsmallest(out_num, 'market_value').loc[:, 'cusip']
                prev_pf_idx_1 = prev_pf_idx_1.drop(out_ser)
            prev_pf_idx_2 = prev_pf_idx_1.copy(deep=True)
            outflow[i] = out_num
            matured[i] = len(mature)
        # buy the bonds in in_pool with largest market value
        if len(in_cusip) == 0:
            pass
        else:
            if out_num == 0:
                in_num = len(mature) + quota
            else: 
                in_num = num_of_tran + quota
            if quota != 0:
                #print('left over quota')
                quota = 0 # reset the quota
            if len(in_pool) < in_num:
                #print('not enough bonds to BUY')
                prev_pf_idx_2 = prev_pf_idx_2.append(in_pool)
                more = in_num - len(in_pool)
                prev_pf_idx_2_cusip = prev_pf_idx_2['cusip']
                diff = set(new_pf_mv['cusip'].tolist()) - set(prev_pf_idx_2_cusip.tolist())
                pool = new_pf_mv.reindex(pd.Series(list(diff)))
                if len(pool) < more:
                    rest = curmap_idx.drop(prev_pf_idx_2_cusip)
                    prev_pf_idx_2 = prev_pf_idx_2.append(rest.nlargest(more, 'market_value'))
                else:
                    prev_pf_idx_2 = prev_pf_idx_2.append(pool.nlargest(more, 'market_value'))
                
            else:
                # in_num = num_of_tran
                in_bonds = in_pool.nlargest(in_num, 'market_value')
                prev_pf_idx_2 = prev_pf_idx_2.append(in_bonds)
        
        inflow[i] = in_num
        idx = list(prev_pf_idx_2.index.values)
        actual_pf = curmap_idx.loc[idx, :]

        
        if 'cusip' in actual_pf.columns:
            actual_pf.reset_index(drop=True)
        else: 
            actual_pf.reset_index(drop=False)
        
        bonds[i] = actual_pf


    return (bonds, inflow, outflow, matured)

def compute_mv(pflist):
    # takes in the list of portfolio for each moth and cluster map
        
    pf_mv = [] 
    for i in range(len(pflist)):
    # compute the total market value of bonds in monthly portfolio
        cur_mv = 0
        cur_pf = pflist[i]
#        portfolio_mv = cluster_map[i]['market_value'].sum()
        for index, row in cur_pf.iterrows():
            # loop over bonds in the current portfolio
            cur_mv += row['market_value']
            
        pf_mv.append(cur_mv)
    
    return pf_mv

def compute_ret(pflist):
    
    pf_return = []
    for i in range(len(pflist)): 
        total_return = 0
        cur_pf = pflist[i]
        for index, row in cur_pf.iterrows():
            ret = row['total_return_mtd']
            if pd.isna(ret):
                total_return += 0
            else:
                total_return += ret
        # save the total return and number of bonds in the current portfolio
        pf_return.append((total_return, len(pflist[i].index)))
        
    return pf_return
 


def turnover_sector(sec_in, sec_out):
    sec_in = pd.Series(sec_in)
    sec_out = pd.Series(sec_out)
    sec_to = sec_in.add(sec_out)
    return sec_to

def turnover_all(ind_to, util_to, fin_to):
    all_to = (ind_to.add(fin_to)).add(util_to)
    return all_to
    
    # Important variables
    # INDUS_mon: list of dataframes by month (normalized)
    # cluster_map:  list of all bonds with corresponding index, cluster, and mv
    # bonds: 2d list: protoflio with max mv bonds for each month
    # sample_mv_percent: percentage of mv of our bond portfolio with respect to all bonds
    # bonds_diff: number of bonds different for two monthly portfolios
        # you can divide by 2 and get inflow = outflow
            
    
        
        