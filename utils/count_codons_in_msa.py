# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 15:24:21 2021

@author: shosh
"""

import os
import sys
import warnings
import copy
import pandas as pd
import numpy as np
import itertools as it
import statsmodels.stats as smstats
import glob
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
    
    
try:
    import matplotlib.pyplot as plt
except ImportError:
    print('Could not Import data-viz libraries')
    
# from plots_for_hypotheses_test_results import pval_str, animals_names_dict
from build_full_tree_nucl_matrix import create_super_orthologs_df_by_orthomcl_hits, natural_sort


def trinity_comps_stops(trinity_file):
    """
    read relevat information from transcriptome file
    """
    def stop_codon(row):
        if row['strand']=="+":
            if len(row['sequence'])-3<row['end']:
                stop = '---'
            else:
                stop = str(row['sequence'][row['end']:row['end']+3])
        elif row['strand']=="-":
            if row['start']<=3:
                stop = '---'
            else:
                stop = str(row['sequence'][row['start']-4:row['start']-1].reverse_complement())
        else:
            stop = "unknown_strand"
    
        if stop not in ["---","TAG","TGA","TAA"]:
            print(row["id"]+': unrecognized stop')
        
        return stop
    
    data = []
    col_names = ['id','start','end','strand','sequence']
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split('\t')
        rec_data[0] = rec_data[0].rstrip()  #some fastafiles have spaces after each id, so fixing it here.
        strand = rec_data[6]
        start = int(rec_data[2])
        end = int(rec_data[4])
        sequence = record.seq        
        data.append((record.id,start,end,strand,sequence))
    df = pd.DataFrame(data = data, columns = col_names)
    df['stop']=df.apply(lambda row: stop_codon(row), axis=1)
    return df



def count_stops_for_super_orthologs(super_orthologs, trinity_stops_dfs_dict, animals):
    
    stops_dict = {}
    for a in animals:
        stops_dict.update({a:{'TAA':0,'TAG':0,'TGA':0}})
    
    for i, row in super_orthologs.iterrows():
        temp_stops_dict_per_row={}
        count_codons=True
        for a in animals:
            stop = trinity_stops_dfs_dict[a].loc[row[a],"stop"]
            if stop=="---":
                count_codons=False
                break
            else:
                temp_stops_dict_per_row.update({a:stop})
        if count_codons:
            for a in animals:
                stops_dict[a][temp_stops_dict_per_row[a]]+=1

    return stops_dict
     

def create_concatenated_msa_with_stops(msas_path,trinity_stops_dfs_dict,outdir,animals):

    def add_stop_to_align_object(file,trinity_stops_dfs_dict,outpath):
        
        out_file_name = file.split('/')[-1].split('_')[-1]
        with open(file,"r") as f: data = f.read()      
        alignment=AlignIO.read(StringIO(data), 'fasta')
        
        with open(outpath+out_file_name,"w") as fout:
            for row in alignment:
                stop = trinity_stops_dfs_dict[row.id.split('|')[0]].loc[row.id,'stop']    
                seq = str(row.seq)
                last_ungapped_col = max([str(seq).rindex(n) for n in ['A','T','G','C']])
                new_seq = seq[:last_ungapped_col+1]+stop+seq[last_ungapped_col+1:]
                fout.write(">"+row.id+'\n')
                fout.write(new_seq+'\n')
        fout.close()
        
    def concatMSAfiles(animals,msas_path,files_prefix=''):
        joined_msas_dict={}
        for a in animals:
            joined_msas_dict.update({a:''})
        for f in glob.glob(os.path.join(msas_path,files_prefix+'*.fasta')):
            alignment = AlignIO.read(open(f), "fasta")
            for i in range(len(animals)):
                current_id = alignment[i].id
                current_animal = current_id.split('|')[0]
                current_seq = str(alignment[i].seq)
                joined_msas_dict[current_animal]+=current_seq 
        
        with open(msas_path+'concatinated_msas.fasta',"w") as handle:
            for k,v in joined_msas_dict.items():
                handle.write('>'+k+'\n')
                handle.write(v+'\n')

    
    for f in natural_sort(glob.glob(os.path.join(msas_path,'*/codons_msa_for_super_orthologs_*.fasta'))):        
        add_stop_to_align_object(f,trinity_stops_dfs_dict,outdir)    
    concatMSAfiles(animals,outdir)


def count_codons_in_alignment(align, remove_gaped_columns=True):
    """
    Codons distribution in MSA
    Accept
    file - MSA file to count codods per each species.
    remove_gaped_columns - A flag for filtering MSA for columns that not containing gaps
    use_alignment_obj - load and iterate alignment AlignIO alignment object instead of SeqIO to parse fasta content
    
    Returns
    Dataframe with codons as indeces (including codons containing gap in not remove_gaped_columns)
    and columns for codons count (per each species).
    """
    nucl = ['A','T','G','C','-']
    codons=[''.join(x) for x in it.product(nucl,nucl,nucl)] #all posible codons (triplets) from nucl
    codons_dict=dict.fromkeys(codons,0)
    results = []
    
    #read alignment if align is alignment file path
    if type(align) is str:
        with open(align,"r") as f: data = f.read()
        alignment=AlignIO.read(StringIO(data), 'fasta')
    else: #else, continue with alignment object
        alignment=align
    
    length = alignment.get_alignment_length()
    if length%3:
        warnings.warn('Alignment length is not a multiple of 3 !!!') 
    
    # initialize dict for alignment ids with initialized codons count dicts as values  
    alignment_animals=[alignment[i].id for i in range(len(alignment))]
    all_animals_codons_dicts={}
    for a in alignment_animals:
        all_animals_codons_dicts.update({a:copy.deepcopy(codons_dict)})            
    
    for i in np.arange(0,length,3):
        tmp_dict = {}
        codon_align=alignment[:,i:i+3]
        count_column=True
        for j in range(len(alignment_animals)):
            codon_id= codon_align[j].id
            codon_seq=str(codon_align[j].seq)
            if remove_gaped_columns:
                if '-' in codon_seq:
                    count_column=False
                    break
            tmp_dict.update({codon_id:codon_seq})
        
        if count_column:
            for k,v in tmp_dict.items():
                all_animals_codons_dicts[k][v]+=1   
    
    for a, codons_cnt_dict in all_animals_codons_dicts.items():
        results.append(pd.Series(codons_cnt_dict, name=a))
            
    return pd.concat(results,axis=1)
    

def calc_t_for_codons_rates(row,animals,controls,df):
        
    group = [a for a in animals if a not in controls]
    control_rates = [row[a]/sum(df[a].values) for a in controls]
    group_rates = [row[a]/sum(df[a].values) for a in group]        
    control_rates_per_aa = [row[a]/sum(df[df['aa']==row['aa']][a].values) for a in controls]   
    group_rates_per_aa = [row[a]/sum(df[df['aa']==row['aa']][a].values) for a in group]   
    t,p = stats.ttest_ind(group_rates,control_rates)
    x,q = stats.ttest_ind(group_rates_per_aa,control_rates_per_aa)
    row['t_all']=t
    row['p_all']=p
    row['t_aa']=x
    row['p_aa']=q
    
    return row

def calc_t_for_aa_codons(row,animals,controls,df):
        
    group = [a for a in animals if a not in controls]
    aa_codons = list(df[df['aa']==row['aa']]['codon'].values)
    group_rates = []
    control_rates = []
    print(aa_codons)
    if len(aa_codons)==1:
        row['t_aa']=np.nan
        row['p_aa']=np.nan
        row['codon_group_average_rate']=1
        row['codon_control_average_rate']=1

    else:
        try:
            for c in aa_codons:
                codon_row = df[df['codon']==c].squeeze()
                codon_group_average_rate = np.mean([codon_row[a]/sum(df[df['aa']==row['aa']][a].values) for a in group])
                codon_control_average_rate = np.mean([codon_row[a]/sum(df[df['aa']==row['aa']][a].values) for a in controls])
                
                group_rates.append(np.mean([codon_row[a]/sum(df[df['aa']==row['aa']][a].values) for a in group]))
                control_rates.append(np.mean([codon_row[a]/sum(df[df['aa']==row['aa']][a].values) for a in controls]))
        
                row['codon_group_average_rate']=codon_group_average_rate
                row['codon_control_average_rate']=codon_control_average_rate
            
                t,p = stats.ttest_ind(group_rates,control_rates)
            row['t_aa']=t
            row['p_aa']=p
        
        except ValueError:
            row['t_aa']=np.nan
            row['p_aa']=np.nan
            row['codon_group_average_rate']=0
            row['codon_control_average_rate']=0
        
    return row

def calc_chi2_for_aa_codons(row,animals,controls,df):
        
    group = [a for a in animals if a not in controls]
    aa_codons = list(df[df['aa']==row['aa']]['codon'].values)
    group_rates = []
    control_rates = []
    print(aa_codons)
    if len(aa_codons)==1:
        row['chi2_aa']=np.nan
        row['p_chi2_aa']=np.nan
    else:
        try:
            for c in aa_codons:
                codon_row = df[df['codon']==c].squeeze()
                group_rates.append(np.mean([codon_row[a]/sum(df[df['aa']==row['aa']][a].values) for a in group]))
                control_rates.append(np.mean([codon_row[a]/sum(df[df['aa']==row['aa']][a].values) for a in controls]))
            
            chi,p,dof,expected = stats.chi2_contingency([group_rates,control_rates])
            row['chi2_aa']=chi
            row['p_chi2_aa']=p
        
        except ValueError:
            row['chi2_aa']=np.nan
            row['p_chi2_aa']=np.nan
        
    return row


def plot_codons_dist_for_multiple_animals(df,animals,outpath,title):
      
    sorted_df=df.sort_values(['aa','codon'],ascending=[True,True])
    fig, ax = plt.subplots(figsize=(20, 10))
    colors = ['red','blue','yellow','grey','green','pink','black','orange']    
    bar_width=0.4
    codons=list(sorted_df['codon'].values)
    ticks_steps=6
    locs = np.arange(start=0,stop=len(codons)*ticks_steps,step=ticks_steps)   
    all_values=[]
    for i,a in enumerate(animals):
        bars_locs=[x+i*bar_width for x in locs]
        events = list(sorted_df[a].values)
        values = [x/sum(events) for x in events]
        plt.bar(bars_locs,values,width=bar_width,color=colors[i],label=a)
        all_values+=values 
    max_val = max(all_values)
    
    plt.xticks([x+bar_width*(len(animals)-1)/2 for x in locs],codons,rotation=90)   
    plt.ylabel('Codon fraction')
    plt.ylim((0,max_val*1.9))
    aas=sorted(list(set(list(sorted_df['aa'].values))))
    ticks=ax.get_xticks()

    j=0
    k=0
    ticks_vals = list(enumerate(ticks))
    trans = ax.get_xaxis_transform()
    for i,t in ticks_vals:
        if i==j:
            aa=aas[k]
            ticks_n=len(sorted_df[sorted_df['aa']==aa])
            last_tick_for_aa=ticks_vals[i+ticks_n-1]
            plt.text((t+float(ticks_n-1)*ticks_steps/1.8-2.5),-0.1,aa,size=10,transform=trans)    
            plt.plot([t-2.5,last_tick_for_aa[1]+1.2],[-0.07,-0.07],color = 'black',transform=trans,clip_on=False)        
            j+=ticks_n
            k+=1
        else:
            pass
    
    plt.text(-10,max_val*1.07,"T per aa",size=10,rotation=90)
    plt.text(-10,max_val*1.24,"T all 64",size=10,rotation=90)    
        
    p_val_aa_list = list(sorted_df['q_aa'].values)
    p_val_all_list = list(sorted_df['q_all'].values)
    for i,t in enumerate(ticks):
        p_val_aa=pval_str(p_val_aa_list[i])
        plt.text((t-2.5),max_val*1.07,p_val_aa,size=10,rotation=90)    
        plt.plot([t-2.5,t+1.2],[max_val*1.05,max_val*1.05],color = 'black',clip_on=False)        
        p_val_all=pval_str(p_val_all_list[i])
        plt.text((t-2.5),max_val*1.24,p_val_all,size=10,rotation=90)    
        plt.plot([t-2.5,t+1.2],[max_val*1.22,max_val*1.22],color = 'black',clip_on=False)        
        
    
    plt.legend(loc=(0.9,0.75))
    plt.savefig(outpath+title+'.jpg')
    plt.close()



def plot_codons_dist_for_multiple_animals_horizontal(df,animals,outpath,title):
      
    sorted_df=df.sort_values(['aa','codon'],ascending=[True,True])
    fig, ax = plt.subplots(figsize=(20, 10))
    colors = ['red','blue','yellow','grey','green','pink','black','orange']    
    bar_width=0.4
    codons=list(sorted_df['codon'].values)
    ticks_steps=6
    locs = np.arange(start=0,stop=len(codons)*ticks_steps,step=ticks_steps)   
    all_values=[]
    sorted_df.plot.barh()
    for i,a in enumerate(animals):
        bars_locs=[x+i*bar_width for x in locs]
        events = list(sorted_df[a].values)
        values = [x/sum(events) for x in events]
        plt.barh(bars_locs,values,width=bar_width,color=colors[i],label=a)
        all_values+=values 
    max_val = max(all_values)
    
    plt.xticks([x+bar_width*(len(animals)-1)/2 for x in locs],codons,rotation=90)   
    plt.ylabel('Codon fraction')
    plt.ylim((0,max_val*1.9))
    aas=sorted(list(set(list(sorted_df['aa'].values))))
    ticks=ax.get_xticks()

    j=0
    k=0
    ticks_vals = list(enumerate(ticks))
    trans = ax.get_xaxis_transform()
    for i,t in ticks_vals:
        if i==j:
            aa=aas[k]
            ticks_n=len(sorted_df[sorted_df['aa']==aa])
            last_tick_for_aa=ticks_vals[i+ticks_n-1]
            plt.text((t+float(ticks_n-1)*ticks_steps/1.8-2.5),-0.1,aa,size=10,transform=trans)    
            plt.plot([t-2.5,last_tick_for_aa[1]+1.2],[-0.07,-0.07],color = 'black',transform=trans,clip_on=False)        
            j+=ticks_n
            k+=1
        else:
            pass
    
    plt.text(-10,max_val*1.07,"T per aa",size=10,rotation=90)
    plt.text(-10,max_val*1.24,"T all 64",size=10,rotation=90)    
        
    p_val_aa_list = list(sorted_df['q_aa'].values)
    p_val_all_list = list(sorted_df['q_all'].values)
    for i,t in enumerate(ticks):
        p_val_aa=pval_str(p_val_aa_list[i])
        plt.text((t-2.5),max_val*1.07,p_val_aa,size=10,rotation=90)    
        plt.plot([t-2.5,t+1.2],[max_val*1.05,max_val*1.05],color = 'black',clip_on=False)        
        p_val_all=pval_str(p_val_all_list[i])
        plt.text((t-2.5),max_val*1.24,p_val_all,size=10,rotation=90)    
        plt.plot([t-2.5,t+1.2],[max_val*1.22,max_val*1.22],color = 'black',clip_on=False)        
        
    
    plt.legend(loc=(0.9,0.75))
    plt.savefig(outpath+title+'.jpg')
    plt.close()


if __name__=='__main__':
    
    count_stops = True
    count_codons=False
    calc_results=False
    
    if count_stops:
        
        # all_best_hits_file=sys.argv[1]
        # trinity_path=sys.argv[2]
        # msas_path=sys.argv[3]
        # out_dir = sys.argv[4]
        all_best_hits_file = "C:/Users/shosh/Google_Drive/adaptive recoding in cephalopods/new_native_transcriptomes/all_best_hits.txt"
        trinity_path = "C:/Users/shosh/Google_Drive/adaptive recoding in cephalopods/new_native_transcriptomes/"
        msas_path = "D:/RNA_Editing_large_files_Backup_20201205/test/"
        out_dir = 'count_codons'
        animals = ['apl','nau','oct','bim','sep','squ','bob','lin']
        
        
        out_dir = '/'+'/'.join(all_best_hits_file.split('/')[:-1])+'/'+out_dir+'/'
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score']) 
        all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
        all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
        super_orthologs = create_super_orthologs_df_by_orthomcl_hits(all_best_hits,animals)
        
        trinity_stops_dfs_dict = {}
        for a in animals:
            print("Reading stops for "+a)
            trinity_file = trinity_path + "orfs_"+a+".fa"
            trinity_df = trinity_comps_stops(trinity_file)
            trinity_df["id"]= trinity_df.apply(lambda row: a+"|"+row["id"], axis=1)
            trinity_df.set_index("id",inplace=True)
            trinity_stops_dfs_dict.update({a:trinity_df})


        stops_dict = count_stops_for_super_orthologs(super_orthologs, trinity_stops_dfs_dict, animals)                
        data_dict = {k:(v['TAA'],v['TAG'],v['TGA']) for k,v in stops_dict.items()}
        data_dict.update({'codon':('TAA','TAG','TGA')})
        stops_cnt_df = pd.DataFrame.from_dict(data_dict,orient='columns')
        stops_cnt_df.set_index('codon',inplace=True)
        stops_cnt_df.to_excel(out_dir+'stops_cnt_per_orthologs.xlsx')
    
        create_concatenated_msa_with_stops(msas_path,trinity_stops_dfs_dict,out_dir,animals)
        remove_gaps = True
        concat_msa_file = out_dir+'concatinated_msas.fasta'
        codons_df = count_codons_in_alignment(concat_msa_file,remove_gaped_columns=remove_gaps)
        codons_df.to_csv(out_dir+'/'+"all_codons_cnt_from_msa.txt",sep='\t',index=True)
    
    
    if count_codons:
        file = sys.argv[1]
        remove_gaps = sys.argv[2]
        out_file = sys.argv[3]
        
        # file = 'D:/RNA_Editing_large_files_Backup_20201205/codons_msa_for_super_orthologs_2613.fasta'
        # remove_gaps=True
        # out_file='codons_no_gaps_test'
        concat_msa_file = msas_path
        codons_df = count_codons_in_alignment(file,remove_gaped_columns=remove_gaps)
        codons_df.to_csv(os.getcwd()+'/'+out_file,sep='\t',index=True)

    if calc_results:
        # file='D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Sanchez/all8/codons_cnt_tbl'
        file='D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Raxml/all8/codons_count_tbl_no_gaps'
        
        outpath='/'.join(file.split('/')[:-1])+'/'
        title=file.split('/')[-1]+'_dist'
        raw_animals_order = ['apl','nau','oct','bim','sep','squ','bob','lin']
        codons_df = pd.read_csv(file,sep='\t',index_col=0)
        codons_df=codons_df[raw_animals_order]
        codons_df=codons_df.rename(columns={key:val for key,val in [(a,animals_names_dict[a]) for a in codons_df.columns.values]})
        codons_df=codons_df[[False if '-' in c else True for c in codons_df.index]]
        animals=codons_df.columns.values
        codons_df.reset_index(inplace=True)
        codons_df=codons_df.rename(columns={'index':'codon'})
        codons_df['aa']=codons_df.apply(lambda row: Seq(row.codon).translate(),axis=1)
        # codons_df = codons_df.apply(lambda row: calc_t_for_codons_rates(row,animals,['Apl','Nau'],codons_df.copy()),axis=1)
        # codons_df.sort_values('p_all',ascending=False,inplace=True)
        # codons_df['q_all'] = np.hstack((fdrcorrection([p for p in codons_df['p_all'].values if p>=0], alpha=0.05, method='indep', is_sorted=False)[1], 
        #                                 [np.nan for p in codons_df['p_all'].values if not p>=0]))         
        # codons_df.sort_values('p_aa',ascending=False,inplace=True)
        # codons_df['q_aa'] = np.hstack((fdrcorrection([p for p in codons_df['p_aa'].values if p>=0], alpha=0.05, method='indep', is_sorted=False)[1], 
        #                                 [np.nan for p in codons_df['p_aa'].values if not p>=0])) 
        
        # codons_df = codons_df.apply(lambda row: calc_fisher_for_aa_codons(row,animals,['Apl','Nau'],codons_df.copy()),axis=1)
        # codons_df.sort_values('p_chi2_aa',ascending=False,inplace=True)
        # codons_df['q_chi2_aa'] = np.hstack((fdrcorrection([p for p in codons_df['p_chi2_aa'].values if p>=0], alpha=0.05, method='indep', is_sorted=False)[1], 
        #                                [np.nan for p in codons_df['p_chi2_aa'].values if not p>=0])) 
        
        codons_df = codons_df.apply(lambda row: calc_t_for_aa_codons(row,animals,['Apl','Nau'],codons_df.copy()),axis=1)
        codons_df.sort_values('p_aa',ascending=False,inplace=True)
        codons_df['q_aa'] = np.hstack((fdrcorrection([p for p in codons_df['p_aa'].values if p>=0], alpha=0.05, method='indep', is_sorted=False)[1], 
                                       [np.nan for p in codons_df['p_aa'].values if not p>=0])) 
        



        codons_df.to_excel(outpath+file.split('/')[-1].split('.')[0]+'_Tstats'+'.xlsx')
        # plot_codons_dist_for_multiple_animals_horizontal(codons_df,animals,outpath,title)