# -*- coding: utf-8 -*-
"""
Created on Wed May 17 18:09:01 2020

@author: shosh

Thie file finds all super orthologs - genes that are the two-way best hits in each pair of animals among a group of several animals
It does so based on the all_best_hits file containing pair of orthologs that are the 2-way best hits in the orthomcl pipeline that uses all vs all blast
"""
import os
import pandas as pd
import numpy as np
import itertools as it
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, IUPAC
import argparse
import time


def read_trinity_mrna_files(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = ['component','protein','protein_name','orfs_start','orfs_end','strand','sequence']
    
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split('\t')
        if rec_data[0][-1] == ' ':  #some fastafiles have spaces after each id, so fixing it here.
            rec_data[0] = rec_data[0].replace(' ','')
#        protein = rec_data[-1].split('|')[1]   #reading proteing from description assuming it was added to header using the transcriptome built pipeline we have for trinity
        protein = rec_data[-1].split('|')[2].split(' ')[0]
        rec_data = (rec_data[0],protein,protein.split('_')[0],int(rec_data[2]),int(rec_data[4]),rec_data[6],record.seq)
        data.append(rec_data)
    
    df = pd.DataFrame(data = data, columns = col_names)
    
    return df

def retrive_native_coding_sequence_in_orfs(trinity_row, editing_sites_tbl):
    """
    this function create the native coding sequence sequence
    given a sequence that is edited (in all strong sites, but anyway the function iterate over all sites) 
    and a dataframe of edited sites (with mm_type column)
    """
    sites_in_sequence = editing_sites_tbl[editing_sites_tbl['component']==trinity_row.name]
    complementaries = {'A':'T','T':'A','G':'C','C':'G'}
    sequence = str(trinity_row['sequence'])
    
    for i, es_row in sites_in_sequence.iterrows():
        
        #determine native and edited nucs in input seq
        native_coding_nuc = es_row['mm_type'][0]
        if es_row['strand'] == '-':
            native_nuc_in_seq = complementaries[native_coding_nuc]
        elif es_row['strand'] == '+':
            native_nuc_in_seq = native_coding_nuc
                
        sequence = sequence[:es_row['location']-1] + native_nuc_in_seq + sequence[es_row['location']:]
        
    if trinity_row['strand'] == '-':
        coding_sequence = str(Seq(sequence[trinity_row['orfs_start']-1:trinity_row['orfs_end']], generic_dna).reverse_complement())
    elif trinity_row['strand'] == '+':
        coding_sequence = sequence[trinity_row['orfs_start']-1:trinity_row['orfs_end']]

    return coding_sequence


def read_editing_sites_tbl(editing_sites_file, mm_type = 'AG'):
    """
    read the editing sites tabel
    """
    col_names = ['component', 'protein', 'location', 'mm_type', 'DNA_A', 'DNA_T', 'DNA_G', 'DNA_C',
                 'RNA_A', 'RNA_T', 'RNA_G', 'RNA_C', 'Trinity', 'RNA_coverage', 'DNA_coverage', 
                 'p_val', 'AA_before', 'AA_after', 'type', 'protein_length', 'editing_level', 'strand']
    sites_df = pd.read_csv(editing_sites_file, sep = '\t', names = col_names, index_col = None)
    sites_df = sites_df[sites_df['mm_type'] == mm_type]
    return sites_df


def create_super_orthologs_df(transcripts_dict,animals):
    """
    group super orthologs by protein name
    """
    proteins = []
    for k,v in transcripts_dict.items():
        proteins+=list(v['protein_name'])
        
    super_orthologs_list=[]
    previous_prots_names=[]
    for prot in set(proteins):    
        if prot not in previous_prots_names:
            group_genes=True
            super_orthologs=(prot,)        
                
            for a in animals:
                df=transcripts_dict[a]
                prot_in_a = df[df['protein_name']==prot].copy()
                if len(prot_in_a)==1:
                    super_orthologs+=(prot_in_a.squeeze()['component'],)
                elif len(prot_in_a)>1:
                    prot_in_a['length'] = df.apply(lambda row: len(row['sequence']), axis=1)
#                    df['length'] = df.apply(lambda row: row['orfs_end']-row['orfs_start'], axis=1)
                    longest_prot_row=prot_in_a.loc[prot_in_a['length'].idxmax(),:]
                    super_orthologs+=(longest_prot_row.squeeze()['component'],)
                else:
                    group_genes=False
                    
            if group_genes:
                super_orthologs_list.append(super_orthologs)
            previous_prots_names.append(prot)
        
    columns=['protein_name']+animals
    super_orthologs_df=pd.DataFrame(data=super_orthologs_list, columns=columns)
        
    return super_orthologs_df



def create_super_ortholog_fasta(super_orthologs_df, super_ortholog_transcripts_dict, animals, output_dir,proteins_msa=False):

    n_files=0
    n_ambiguos_mappings = 0
    
    for index, row in super_orthologs_df.iterrows():
         
        n_files+=1
        related_proteins = {}
        for a in animals:
            related_proteins.update({a:super_ortholog_transcripts_dict[a].loc[row[a],'protein']})
        if len(set(list(related_proteins.values()))) > 1:
            n_ambiguos_mappings+=1
#            print('ambigouos protein:')
#            print(related_proteins)
#
        file_name = 'super_orthologs_'+str(n_files)+'.fasta'
        writer =  FastaWriter(open(output_dir + file_name, 'w'), wrap=None)
        writer.write_header()
        
        for a in animals:
            trinity_data = super_ortholog_transcripts_dict[a].loc[row[a],:]
            sequence = Seq(trinity_data['native_coding_sequence'], generic_dna)
            if proteins_msa:
                sequence = str(sequence.translate())
                if '*' in sequence:
                    print(file_name + ' has a stop')
                sequence = Seq('W'.join(sequence.split('*')), IUPAC) #assuming all stop codonds in sequences are lost due to editing and translated as W
            writer.write_record(SeqRecord(sequence, id = row[a], description=''))
        writer.write_footer()
        
    print('\n'+str(n_ambiguos_mappings) + ' ambigouos uniprot mappings within super ortholog genes')

def filter_editing_sites_beyond_orfs(es_df,trinity_df):
    
    es_dfs_list = []
    for comp, row in trinity_df.iterrows():
        es = es_df[es_df['component']==comp]
        if len(es):
            es = es[np.logical_and(es['location']>=row['orfs_start'], es['location']<=row['orfs_end'])]
            es_dfs_list.append(es)
    if len(es_dfs_list):
        return pd.concat(es_dfs_list)
    else:
        return None

def count_stop_codonds_in_super_orthologs_orfs(row, transcripts_dict, animals):
    
    stop_codons=0
    for a in animals:
        native_coding_sequence = transcripts_dict[a].loc[row[a],'native_coding_sequence']
        protein = str(Seq(native_coding_sequence, generic_dna).translate())
        stop_codons+=protein.count('*')
    return stop_codons
        
    

if __name__ == '__main__':
    
# =============================================================================
#     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Finding ortholog sites for 2 species based on orthomcl and blastp results')
#     run_parser = parser.add_argument_group('Run muscle for a list of super orthologs genes')
#     run_parser.add_argument('-parent_dir', dest='parent_dir', action='store', required = True, help='path to parent directory in which all_best_hits file and trinity transcriptomes directory resides')
#     run_parser.add_argument('-o', dest='super_orthologs_output_dir', action='store', required = True, help='name of super ortholog fastas directory in which to create all fasta files and run msa processes')
#     run_parser.add_argument('-animals', dest='animals', action='store', nargs = '+', default = ['oct','squ','bim','nau','sep','apl'], help='subset of animals for finding super orthologs and perform MSA')
#     run_parser.add_argument('-program', dest='program', action='store', default = 'clu', help='select either clu or mus (clustal omega or muscle)')
#     run_parser.add_argument('-outfmt', dest='outfmt', action='store', default = 'clu', help='select out format for specified program')
#     run_parser.add_argument('-msa_processes', dest='msa_proc', action='store', default = '30', help='number of msa processes to run in parallel')
#     run_parser.add_argument('-gapopen', dest='gapopen', action='store', default = None, help='gapopen penalty for muscle msa program')
#     run_parser.add_argument('-proteins_msa', dest='proteins_msa', action='store', default = 'True', help='run msa for proteins sequences')
#     run_parser.add_argument('-only_edited_super_orthologs', dest='only_edited_super_orthologs', action='store', default = 'True', help='only keep super orthologs with editing sites')
#     
#     arguments = parser.parse_args()
#     
#     parent_dir = arguments.parent_dir
#     super_orthologs_output_dir = arguments.super_orthologs_output_dir
#     animals = arguments.animals
#     program = arguments.program
#     outfmt = arguments.outfmt
#     msa_proc = int(arguments.msa_proc)
#     gapopen = arguments.gapopen
#     proteins_msa = eval(arguments.proteins_msa)
#     only_edited_super_orthologs=arguments.only_edited_super_orthologs
#     
#     assert (program in ['clu','mus']), "program must be either clu or mus"
#     if gapopen is not None:
#         assert (program=='clu'), "cant specify gapopn penalty for clustalo"
#         
#     editing_sites_path = parent_dir + 'editing_sites/'
#     trinity_files_path = parent_dir + 'trinity_comps/'
#     ortholog_fastas_output_dir = parent_dir+super_orthologs_output_dir+'/'
# =============================================================================
    
    parent_dir = 'E:/RNA_editing_Large_files/Phylogeny/hpm_reconstruction/' 
    animals = ['oct','squ','bim','nau','sep','apl']
    proteins_msa = True
    super_orthologs_output_dir = 'super_orthologs_by_names'
    only_edited_super_orthologs = True
    editing_sites_path = parent_dir + 'editing_sites_chinese_fix/'
    trinity_files_path = parent_dir + 'new_native_transcriptomes_chinese_fix/'
    ortholog_fastas_output_dir = parent_dir+super_orthologs_output_dir+'/'
    
    print('\nMSA for super orthologs in ' + str(animals) + '\n')
    if not os.path.exists(ortholog_fastas_output_dir):
        os.makedirs(ortholog_fastas_output_dir)
    else:
        import warnings
        warnings.warn('Output dir already exists and existing fasta files will be overwritten')
       
    #collencting fasta records for super ortholog sequences only
    print('Collecting super ortholog records from fasta transcriptomes')
    transcripts_dict = {}    
    for a in animals:
        trinity_df = read_trinity_mrna_files(trinity_files_path+'orfs_'+a+'.fa')
        trinity_df['component'] = trinity_df.apply(lambda row: a+'|'+row['component'], axis=1)
        transcripts_dict.update({a:trinity_df})
    
    super_ortholog_genes = create_super_orthologs_df(transcripts_dict,animals)
    print('Super ortholog genes: ' + str(len(super_ortholog_genes)))
    
    #collencting editing sites for super ortholog sequences only
    print('Collecting editing sites for super ortholog sequences')
    super_ortholog_editing_sites_dict = {}
    for a in animals:
        es_df = read_editing_sites_tbl(editing_sites_path+a+'_editing_site_info.txt')
        es_df['component'] = es_df.apply(lambda row: a+'|'+row['component'], axis=1)
        super_orthologs_es = es_df[es_df['component'].isin(super_ortholog_genes[a])]
        super_ortholog_editing_sites_dict.update({a:super_orthologs_es})
        print(a + ' editing sites: ' + str(len(super_orthologs_es)) + ' out of ' + str(len(es_df)))
         
    super_ortholog_transcripts_dict={}
    for k,trinity_df in transcripts_dict.items():
        trinity_df.set_index('component', inplace = True)
        super_orthologs_trinity_df = trinity_df[trinity_df.index.isin(super_ortholog_genes[k])].copy()
        super_orthologs_trinity_df['native_coding_sequence'] = super_orthologs_trinity_df.apply(lambda row: retrive_native_coding_sequence_in_orfs(row,super_ortholog_editing_sites_dict[k]), axis = 1)
        super_ortholog_editing_sites_dict[k] = filter_editing_sites_beyond_orfs(super_ortholog_editing_sites_dict[k],super_orthologs_trinity_df)
        super_ortholog_transcripts_dict.update({k:super_orthologs_trinity_df})
    
    super_ortholog_genes['editing_sites'] = super_ortholog_genes.apply(lambda row: sum([len(super_ortholog_editing_sites_dict[a][super_ortholog_editing_sites_dict[a]['component']==row[a]]) for a in animals if super_ortholog_editing_sites_dict[a] is not None]), axis=1)
    super_ortholog_genes['stop_codons'] = super_ortholog_genes.apply(lambda row : count_stop_codonds_in_super_orthologs_orfs(row, super_ortholog_transcripts_dict, animals), axis=1)
    if only_edited_super_orthologs:
        super_ortholog_genes = super_ortholog_genes[super_ortholog_genes['editing_sites']>0]
    super_ortholog_genes.to_excel('E:/RNA_editing_Large_files/Phylogeny/hpm_reconstruction/orthologs_by_name.xlsx', index=False)
    
# =============================================================================
#     print('Writing super ortholog sequences to seperate fasta files')
#     create_super_ortholog_fasta(super_ortholog_genes,super_ortholog_transcripts_dict,animals, ortholog_fastas_output_dir, proteins_msa = proteins_msa)  
#     
#     msa_results_dir = ortholog_fastas_output_dir+'msa_results/'
#     if not os.path.exists(msa_results_dir):
#         os.makedirs(msa_results_dir)
#         msa_logs_dir = msa_results_dir+'msa_logs/'
#         os.makedirs(msa_logs_dir)
#     else:
#         import warnings
#         warnings.warn('MSA results dir already exists and existing files will be overwritten')
#         msa_logs_dir = msa_results_dir+'msa_logs/'
#         if not os.path.exists(msa_logs_dir):
#             os.makedirs(msa_logs_dir)
#         else:
#             warnings.warn('MSA log dir already exists and existing files will be overwritten')
#     
#     print('Running MSA for super orthologs fasta files')
#     processes = []
#     files_list = glob.glob(os.path.join(ortholog_fastas_output_dir,'super_orthologs*.fasta'))
#     pennding_processes = True
#     current_proc = 0
#     while pennding_processes:
#         if sum([p.poll() is None for p in processes]) < msa_proc:
#             file_path = files_list[current_proc]
#             file_name = file_path.split('/')[-1].split('.')[0]
#             print('starting MSA for ' + file_name)
#             if program == 'mus':
#                 if outfmt=='clu':
#                     outfmt_pass = ' --clwout '
#                     suffix='aln'
#                 elif outfmt=='fasta':
#                     outfmt_pass = ' --fastaout '
#                     suffix='fasta'
#                 if gapopen is None:
#                     msa_cmd = 'nohup muscle -in '+file_path + outfmt_pass+msa_results_dir+'msa_results_for_'+file_name+'.'+suffix + ' > ' + msa_logs_dir+'msa_log_for_'+file_name+'.txt'
#                 else:
#                     msa_cmd = 'nohup muscle -in '+file_path + outfmt_pass+msa_results_dir+'msa_results_for_'+file_name+'.'+suffix + ' -gapopen '+ gapopen + ' > ' + msa_logs_dir+'msa_log_for_'+file_name+'.txt'
#             elif program == 'clu':
#                 if proteins_msa:
#                     seq_type = 'Protein'
#                 else:
#                     seq_type = 'DNA'
#                 msa_cmd = 'nohup clustalo -i '+file_path + ' --outfile='+msa_results_dir+'msa_results_for_'+file_name+'.aln' + ' --seqtype='+seq_type + ' --outfmt '+outfmt + ' > ' + msa_logs_dir+'msa_log_for_'+file_name+'.txt'
# 
#             processes.append(subprocess.Popen(msa_cmd, shell = True, universal_newlines = True))
#             current_proc+=1
#         if len(files_list)==current_proc:
#             pennding_processes = False
#             
#     #wait until all processes are finished        
#     while sum([p.poll() is None for p in processes]):
#         time.sleep(1)
#     
#     print('Finished running msa for all super orthologs genes')
# =============================================================================


# =============================================================================
# aln = 'E:/RNA_editing_Large_files/Phylogeny/hpm_reconstruction/cephalopod_alignment.txt'
# 
# proteins_in_aln = []
# for i in open(aln, "r").readlines():
#     if ">" in i:
#         proteins_in_aln.append(i.split('\t')[1].rstrip())
# proteins_in_aln = list(set(proteins_in_aln))
# 
# proteins_super_orthologs = list(set(list(super_ortholog_genes['protein_name'])))
# 
# not_in_super = []
# not_in_aln = []
# 
# for i in proteins_super_orthologs:
#     if i not in proteins_in_aln:
#         not_in_aln.append(i)
#         
# for i in proteins_in_aln:
#     if i not in proteins_super_orthologs:
#         not_in_super.append(i)
# 
# =============================================================================
