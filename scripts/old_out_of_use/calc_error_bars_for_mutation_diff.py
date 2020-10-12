# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 18:05:42 2020

@author: shosh
"""
#import scipy.stats as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats


def calc_strict_model_distribuation(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n):
    
    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
    p=float(strong_syn_sites)/syn_a
    z = 1-(1-confidence_level)/2
    
    expected_non_syn_sites_dist = []
    for i in range(n):
        p_rand = np.random.normal(p,z*np.sqrt(p*(1-p)/syn_a))
        lam = non_syn_a*p_rand
        expected_non_syn_sites_dist.append(np.random.poisson(lam))
    
    expected_non_syn_eg_mut_dist = []
    for i in range(n):
        expected_eg_mut = (strong_non_syn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
        expected_non_syn_eg_mut_dist.append(np.random.poisson(expected_eg_mut))    
    
    return expected_non_syn_eg_mut_dist


def calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,confidence_level,n):
    
    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
    non_syn_ag_mut_rate = float(non_syn_ag_mut)/non_syn_a
    z = 1-(1-confidence_level)/2
    p_weak_nss=float(weak_non_syn_sites)/non_syn_a
    p_weak_ss=float(weak_syn_sites)/syn_a
    p_strong_ss=float(strong_syn_sites)/syn_a
    
    
    expected_non_syn_sites_dist = []
    for i in range(n):
        p_weak_nss_rand = np.random.normal(p_weak_nss,z*np.sqrt(p_weak_nss*(1-p_weak_nss)/non_syn_a))
        p_weak_ss_rand = np.random.normal(p_weak_ss,z*np.sqrt(p_weak_ss*(1-p_weak_ss)/syn_a))
        p_strong_ss_rand = np.random.normal(p_strong_ss,z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a))
        lam = (p_weak_nss_rand/p_weak_ss_rand)*p_strong_ss_rand*non_syn_a
        expected_non_syn_sites_dist.append(np.random.poisson(lam))
        
    expected_non_syn_eg_mut_dist = []
    for i in range(n):
        expected_eg_mut_from_excess = (strong_non_syn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
        p_non_syn_ag_mut_rate_rand = np.random.normal(non_syn_ag_mut_rate,z*np.sqrt(non_syn_ag_mut_rate*(1-non_syn_ag_mut_rate)/non_syn_a))
        expected_eg_mut_from_expected_nss = expected_non_syn_sites_dist[i]*p_non_syn_ag_mut_rate_rand
        expected_non_syn_eg_mut_dist.append(np.random.poisson(expected_eg_mut_from_excess)+np.random.poisson(expected_eg_mut_from_expected_nss))
        
    return expected_non_syn_eg_mut_dist


def calc_strict_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n):
    
    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
    z = 1-(1-confidence_level)/2
    p=float(strong_syn_sites)/syn_a
    p_rand = stats.truncnorm(-2,2,loc=p,scale=z*np.sqrt(p*(1-p)/syn_a)).rvs(n)
    
    expected_non_syn_sites_dist = []
    for i in range(n):
        lam = non_syn_a*p_rand[i]
        expected_non_syn_sites_dist.append(np.random.poisson(lam))
    
    expected_non_syn_eg_mut_dist = []
    for i in range(n):
        expected_eg_mut_lam = (strong_non_syn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
        expected_non_syn_eg_mut_dist.append(np.random.poisson(expected_eg_mut_lam))
    
    return expected_non_syn_eg_mut_dist


def calc_liberal_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n):
    
    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
    non_syn_ag_mut_rate = float(non_syn_ag_mut)/non_syn_a
    z = 1-(1-confidence_level)/2
    p_weak_nss=float(weak_non_syn_sites)/non_syn_a
    p_weak_ss=float(weak_syn_sites)/syn_a
    p_strong_ss=float(strong_syn_sites)/syn_a
    p_weak_nss_rand = stats.truncnorm(-2,2,loc=p_weak_nss,scale=z*np.sqrt(p_weak_nss*(1-p_weak_nss)/non_syn_a)).rvs(n)
    p_weak_ss_rand = stats.truncnorm(-2,2,loc=p_weak_ss,scale=z*np.sqrt(p_weak_ss*(1-p_weak_ss)/syn_a)).rvs(n)
    p_strong_ss_rand = stats.truncnorm(-2,2,loc=p_strong_ss,scale=z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a)).rvs(n)
    p_non_syn_ag_mut_rate_rand = stats.truncnorm(-2,2,loc=non_syn_ag_mut_rate,scale=z*np.sqrt(non_syn_ag_mut_rate*(1-non_syn_ag_mut_rate)/non_syn_a)).rvs(n)
    
    expected_non_syn_sites_dist = []
    for i in range(n):
        lam = (p_weak_nss_rand[i]/p_weak_ss_rand[i])*p_strong_ss_rand[i]*non_syn_a
        expected_non_syn_sites_dist.append(np.random.poisson(lam))
        
    expected_non_syn_eg_mut_dist = []
    for i in range(n):
        expected_eg_mut_lam = (strong_non_syn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
        expected_eg_mut_from_expected_nss = expected_non_syn_sites_dist[i]*p_non_syn_ag_mut_rate_rand[i]
        expected_non_syn_eg_mut_dist.append(np.random.poisson(expected_eg_mut_lam)+np.random.poisson(expected_eg_mut_from_expected_nss))
        
    return expected_non_syn_eg_mut_dist
    
    


syn_a = 269168
non_syn_a = 859775
syn_ag_mut = 23058
non_syn_ag_mut = 5446
weak_syn_sites = 1170
weak_non_syn_sites = 2118
strong_syn_sites = 347
strong_non_syn_sites = 1971
eg_non_syn_mutations = 30
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)
    

syn_a = 267167
non_syn_a = 862246
syn_ag_mut = 19937
non_syn_ag_mut = 6325
weak_syn_sites = 1696
weak_non_syn_sites = 3287
strong_syn_sites = 341
strong_non_syn_sites = 1850
eg_non_syn_mutations = 60
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)



syn_a = 261427
non_syn_a = 854796
syn_ag_mut = 22756
non_syn_ag_mut = 6509
weak_syn_sites = 694
weak_non_syn_sites = 1374
strong_syn_sites = 245
strong_non_syn_sites = 1552
eg_non_syn_mutations = 49
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)


syn_a = 276896
non_syn_a = 862909
syn_ag_mut = 33014
non_syn_ag_mut = 8956
weak_syn_sites = 1577
weak_non_syn_sites = 2801
strong_syn_sites = 375
strong_non_syn_sites = 1909
eg_non_syn_mutations = 120
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)






syn_a = 269168
non_syn_a = 859775
syn_ag_mut = 23058
non_syn_ag_mut = 5446
weak_syn_sites = 1170
weak_non_syn_sites = 2118
strong_syn_sites = 347
strong_non_syn_sites = 1971
eg_non_syn_mutations = 30
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)
    

syn_a = 267167
non_syn_a = 862246
syn_ag_mut = 19937
non_syn_ag_mut = 6325
weak_syn_sites = 1696
weak_non_syn_sites = 3287
strong_syn_sites = 341
strong_non_syn_sites = 1850
eg_non_syn_mutations = 60
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)



syn_a = 261427
non_syn_a = 854796
syn_ag_mut = 22756
non_syn_ag_mut = 6509
weak_syn_sites = 694
weak_non_syn_sites = 1374
strong_syn_sites = 245
strong_non_syn_sites = 1552
eg_non_syn_mutations = 49
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)


syn_a = 276896
non_syn_a = 862909
syn_ag_mut = 33014
non_syn_ag_mut = 8956
weak_syn_sites = 1577
weak_non_syn_sites = 2801
strong_syn_sites = 375
strong_non_syn_sites = 1909
eg_non_syn_mutations = 120
n=1000000
confidence_level = 0.95

strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
strict_p = len([j for j in strict_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(strict_std)
print(strict_p)

liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation_trunc(syn_a,non_syn_a,syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,non_syn_ag_mut,confidence_level,n)
liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
liberal_p = len([j for j in liberal_expected_non_syn_eg_mut_dist if j<eg_non_syn_mutations])/n
print(liberal_std)
print(liberal_p)

# =============================================================================
# path = 'C:/Users/shosh/OneDrive/Desktop/test/'
# edited_msa_columns = pd.read_csv(path+'edited_msa_columns',sep='\t', index_col = None)
# 
# col_names = ['component', 'protein', 'location', 'mm_type', 'DNA_A', 'DNA_T', 'DNA_G', 'DNA_C',
#              'RNA_A', 'RNA_T', 'RNA_G', 'RNA_C', 'Trinity', 'RNA_coverage', 'DNA_coverage',
#              'p_val', 'AA_before', 'AA_after', 'type', 'protein_length', 'editing_level', 'strand']
#    
# def get_number_of_conserved_animsls(row,a):
#     return edited_msa_columns.loc[edited_msa_columns[a+'_site_key']==row['site_key'],'edited_animals'].tolist()[0]
# 
# animals = ['oct','bim','squ','sep','bob','lin']
# for a in animals:
#     sites_df = pd.read_csv(path+a+'_comps_strongest_is_syn', sep = '\t', names = col_names, index_col = None)
#     sites_df['site_key'] = sites_df.apply(lambda row: a+'|'+row['component']+';'+str(row['location']), axis=1)
#     sites_df['edited_animals'] = sites_df.apply(lambda row: get_number_of_conserved_animsls(row,a), axis=1)
#     sites_df.to_csv(path+a, sep='\t', index=False)
# =============================================================================









