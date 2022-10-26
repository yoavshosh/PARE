# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:43:25 2019

@author: shosh

This screpts is creating a list of ortholog editing sites based on the folowing data of the 2 species: 
the transctiptome fasta files (from trinity pipline)
the editing sites lists - created by the perl scripts as a proceeding piplen after trinity transcriptomes are built
a list of the 2-way best hits of the potential orthologs created using the orthomcl pipline
the blastp results (based on which orthomcl pipline was done)
"""
import re
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
import argparse


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



def read_trinity_mrna_files(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = ['comp_id','orfs_start','orfs_end','strand','sequence']
    
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split('\t')
        if rec_data[0][-1] == ' ':  #some fastafiles have spaces after each id, so fixing it here.
            rec_data[0] = rec_data[0].replace(' ','')
        rec_data = (rec_data[0],int(rec_data[2]),int(rec_data[4]),rec_data[6],record.seq)
        data.append(rec_data)
    
    df = pd.DataFrame(data = data, columns = col_names)
    df.set_index('comp_id', inplace = True)
    
    return df


def read_blast_results(blast_file, a1, a2):
    """
    read blast resutls file (all-vs_all) and keep only records for a1 as query animal and a2 as subject animal
    """    
    data = []
    col_names = ['comp_query', 'comp_subject','identity','alignment_length','mismatch','gapopen',
               'query_start','query_end','subject_start','subject_end','e_value','bitscore','btop']
    
    with open(blast_file, "r") as bf:
        for i, line in enumerate(bf.readlines()):
            fields = line.split('\t')
            if (fields[0][:len(a1)]==a1 and fields[1][:len(a2)]==a2) or (fields[0][:len(a2)]==a2 and fields[1][:len(a1)]==a1):
                fields[-1] = fields[-1].replace('\n','')
                data.append(fields)
    blast_df = pd.DataFrame(data = data, columns = col_names)
    
    for c in col_names:
        blast_df[c] = pd.to_numeric(blast_df[c],errors = 'ignore')
    
    return blast_df

    
def calculate_site_mrna_motif_data(row, trinity_df, nucs_around_site):
    """
    calculate site location with respect to coding sequence
    the codon (nucleotides), and location within the codon
    """
    trinity_data = trinity_df.loc[row['component'],:].squeeze()
    coding_sequence = trinity_data['native_coding_sequence']
    
    if row['strand'] == '+':
        site_coding_loc_base_0 = int(row['location']) - int(trinity_data['orfs_start'])
    elif row['strand'] == '-':
        site_coding_loc_base_0 = int(int(trinity_data['orfs_end'] - row['location']))
    
    site_pos_in_codon = site_coding_loc_base_0%3  
    codon = coding_sequence[site_coding_loc_base_0-site_pos_in_codon:site_coding_loc_base_0+3-site_pos_in_codon]
        
    row['nuc'] = coding_sequence[site_coding_loc_base_0]
    if site_coding_loc_base_0-nucs_around_site < 0:
        row['location_in_area'] = site_coding_loc_base_0
    else:
        row['location_in_area'] = nucs_around_site
    row['site_area'] = coding_sequence[max(0,site_coding_loc_base_0-nucs_around_site):min(len(coding_sequence),site_coding_loc_base_0+nucs_around_site)]
    row['site_coding_loc_base0'] = site_coding_loc_base_0
    row['site_loc_in_codon_base0'] = site_pos_in_codon
    row['original_codon'] = codon
    row['site_key'] = row['component'] + ';' + str(row['location'])
    
    return row


def find_corresponding_location_in_subject_coding_mrna(blastp_data, q_trinity_comp, q_es, reversed_query = False):
    """
    This function pars the btop results for query and subject sequences 
    and for each editing site within q_es it findes the corresponding location within the mrna of the subject sequence
    for sites that lies beyond the alignment or that
    """
    btop = re.findall(r'\d+|[A-Z-+]{2}', str(blastp_data['btop']))  #a list of numbers representing matching regions in blast and pairs of letters representing aa mismatches (or deletions) 

    #determining the query and subject sequences (blastp resutls are specific way but i use them to find all orthologs for the animal defined as subject in blastp and this script)
    query_start = 'query_start'
    query_end = 'query_end'
    subject_start = 'subject_start'
    if reversed_query:
        query_start = 'subject_start'
        query_end = 'subject_end'
        subject_start = 'query_start'    
        
    current_q_loc = blastp_data[query_start]*3 - 3
    current_s_loc = blastp_data[subject_start]*3 - 3
    
    q_site_keys_to_subject_coding_location_dict = {}
    q_es_no_orthologs_before_alignment = q_es[q_es['site_coding_loc_base0']<current_q_loc]
    q_es_no_orthologs_after_alignment = q_es[q_es['site_coding_loc_base0']>(blastp_data[query_end])*3-1] #this was checked varified to give the last base0 nucleotide of the alignment, if this cause a problem, may be easier to asign those after loop iterate over all sites in alignment range
    
    #for sites that lie beyond the alignment region - there are no ortholgs
    for index, row in pd.concat([q_es_no_orthologs_before_alignment,q_es_no_orthologs_after_alignment]).iterrows():
        q_site_keys_to_subject_coding_location_dict.update({row['site_key']:'site_is_beyond_alignment'})
        
    for i,b in enumerate(btop):
        #mrna window start is previous mrna window end
        previous_q_loc = current_q_loc
        previous_s_loc = current_s_loc
        
        #move number of nucleotides in queiry and subject coding mrna sequences according to btop results to determine mrna window end as the current position acording to the btop results
        try:
            nucs_in_both_seqs = int(b)*3
            current_q_loc += nucs_in_both_seqs
            current_s_loc += nucs_in_both_seqs
        except ValueError:
            
            b_ordered = b
            if reversed_query:  # if subject animal is the one we are looking for its orthologs - read btop aa mismatch pair reversed
                b_ordered = b_ordered[::-1]    
            
            if b_ordered[0] == '-':
                current_s_loc += 3
            elif b_ordered[1] == '-':
                current_q_loc  += 3
            else:
                current_q_loc += 3
                current_s_loc += 3
        
        #find editing sites in query between the previous and current location in query mrna
        q_es_in_range = q_es[np.logical_and(previous_q_loc <= q_es['site_coding_loc_base0'], q_es['site_coding_loc_base0'] < current_q_loc)]
        
        #for each site in current coding mrna window, find corresponding location in subject mrna sequnce    
        for index, row in q_es_in_range.iterrows():
            if previous_s_loc==current_s_loc:  #if this condition is true and there are site in q_es_in_range, that means there is a gap in the subject sequence - thus no orthologs for those sites
                q_site_keys_to_subject_coding_location_dict.update({row['site_key']:'gap'})
            else:
                corresponding_coding_loc_in_subject = previous_s_loc + (row['site_coding_loc_base0'] - previous_q_loc)
                q_site_keys_to_subject_coding_location_dict.update({row['site_key']:corresponding_coding_loc_in_subject})

    return q_site_keys_to_subject_coding_location_dict


def connect_ortholog_sites_for_querie_editing_sites(row, queries_orthologs_dict, s_es, subject_animal):
    """
    This function adds for each row in qury editing sites (from query editing sites tbl) the data regarding the ortholog sequence,
    and if exists, the ortholog location in the subject coding sequence
    and if exists, the subject editing site key from the subject editing sites tbl (based on corresponding coding location in subject sequence)
    """
   
    queries_orthologs = queries_orthologs_dict[row['component']]
    
    if queries_orthologs == 'no_orthologs':   #no ortholog subject sequence for query sequence
        row[subject_animal+'_ortholog'] = 'no_ortholog'
        row[subject_animal+'_bitscore'] = None
        row[subject_animal+'_ortholog_coding_location'] = None
        row[subject_animal+'_ortholog_site_key'] = None
        
    elif queries_orthologs[2] == {}: #no editing sites in query table
        row[subject_animal+'_ortholog'] = queries_orthologs[0]
        row[subject_animal+'_bitscore'] = queries_orthologs[1]
        row[subject_animal+'_ortholog_coding_location'] = None
        row[subject_animal+'_ortholog_site_key'] = None
        print('unattended editing site ' + row['site_key'])
        
    else:
        ortholog_coding_location = queries_orthologs[2][row['site_key']]
        if type(ortholog_coding_location) == str: #ortholog subject sequence exists, but no subject ortholog location corresponsing to editing site query row
            row[subject_animal+'_ortholog'] = queries_orthologs[0]
            row[subject_animal+'_bitscore'] = queries_orthologs[1]
            row[subject_animal+'_ortholog_coding_location'] = ortholog_coding_location
            row[subject_animal+'_ortholog_site_key'] = None
        
        elif type(ortholog_coding_location) == np.int64 :#ortholog subject sequence exists, and corresponsing to editing site query row exists
            subject_component = queries_orthologs[0].split('|')[1]
            row[subject_animal+'_ortholog'] = queries_orthologs[0]
            row[subject_animal+'_bitscore'] = queries_orthologs[1]
            row[subject_animal+'_ortholog_coding_location'] = ortholog_coding_location
            s_es_row = s_es[np.logical_and(s_es['site_coding_loc_base0']==ortholog_coding_location, s_es['component']==subject_component)]   #subject editing events that correspond to the expected subject coding location calculated by the btop results parsing
            if len(s_es_row) == 0:  #subject corresponding coding location is not an editing event
                row[subject_animal+'_ortholog_site_key'] = None
            elif len(s_es_row) == 1: #subject corresponding coding location is an editing event
                s_es_row = s_es_row.squeeze()
                row[subject_animal+'_ortholog_site_key'] = s_es_row['site_key']
            else:   #this condition is not expected to be true as each query editing site should have maximum one corresponding coding location in subject sequence
                print('ambiguous ortholog site for ' + row['site_key'])
        
        else:
            raise Exception('unexpected type for ' + row['site_key'] + ': ' + str(type(ortholog_coding_location)))
    
    return row


def additional_match_data(row, query_animal, subject_animal):
    """
    calculate values for various potential filters
    """ 
    if row[query_animal+'_'+subject_animal+'_ortholog_site_key'] != None:
        row[query_animal+'_'+subject_animal+'_key'] = row[query_animal+'_site_key'].replace(';','_') + '_' + row[subject_animal+'_site_key'].replace(';','_')
        row[query_animal+'_'+subject_animal+'_ortholog_sequences'] = row[query_animal+'_component'] + '_' + row[subject_animal+'_component']
        row[query_animal+'_'+subject_animal+'_same_target_aa'] = row[query_animal+'_AA_after']==row[subject_animal+'_AA_after']
        row[query_animal+'_'+subject_animal+'_same_original_aa'] =  row[query_animal+'_AA_before']==row[subject_animal+'_AA_before'] 
        row[query_animal+'_'+subject_animal+'_same_position_in_codon'] = row[query_animal+'_site_loc_in_codon_base0']==row[subject_animal+'_site_loc_in_codon_base0'] 
    else:
        row[query_animal+'_'+subject_animal+'_key'] = None
        row[query_animal+'_'+subject_animal+'_ortholog_sequences'] = None
        row[query_animal+'_'+subject_animal+'_same_target_aa'] = None
        row[query_animal+'_'+subject_animal+'_same_original_aa'] =  None 
        row[query_animal+'_'+subject_animal+'_same_position_in_codon'] = None
    return row


def find_corresponding_unedited_nuc_data_in_subject_for(row, query_animal, subject_animal, subject_trinity_tbl, nucs_around_site):
    """
    For sites for which ortholog location exists, but not an editing site
    recalculate columns regarding ortholog nuc and sequence aroung this nuc
    (for sites that are orthologs, this data is available through the merge of the subject editing sites tbl to the query editing sites tbl)
    """
    subject_coding_location = row[query_animal+'_'+subject_animal+'_ortholog_coding_location']
    subject_ortholog_site = row[query_animal+'_'+subject_animal+'_ortholog_site_key']
    
    #if subject_coding_location is an integer it means that there exists a subject ortholog sequence 
    #also, if subject_ortholog_site is null it means that there was no corresponding editing site (to the subject_coding_location) in the subject ortholog sequence
    if pd.isnull(subject_ortholog_site) and type(subject_coding_location)==int:
        subject_coding_seq = subject_trinity_tbl.loc[row[query_animal+'_'+subject_animal+'_ortholog'].split('|')[-1], 'native_coding_sequence']
        row[subject_animal+'_nuc'] = subject_coding_seq[subject_coding_location]
        
        if subject_coding_location-nucs_around_site < 0:
            row[subject_animal+'_location_in_area'] = subject_coding_location
        else:
            row[subject_animal+'_location_in_area'] = nucs_around_site
        
        row[subject_animal+'_site_area'] = subject_coding_seq[max(0,subject_coding_location-nucs_around_site):min(len(subject_coding_seq),subject_coding_location+nucs_around_site)]
        
        site_pos_in_codon = subject_coding_location%3  
        if site_pos_in_codon == 0:
            codon = subject_coding_seq[subject_coding_location:subject_coding_location+3]
        elif site_pos_in_codon == 1:
            codon = subject_coding_seq[subject_coding_location-1:subject_coding_location+2]
        elif site_pos_in_codon == 2:
            codon = subject_coding_seq[subject_coding_location-2:subject_coding_location+1]
        
        row[subject_animal+'_site_coding_loc_base0'] = subject_coding_location
        row[subject_animal+'_site_loc_in_codon_base0'] = site_pos_in_codon
        row[subject_animal+'_original_codon'] = codon
        
    return row


def create_ortholog_locations_dict(query_animal, subject_animal, query_trinity_tbl, subject_trinity_tbl, all_best_hits, blastp, query_es_tbl, reversed_query = False):
    """
    This function uses find_corresponding_location_in_subject_coding_mrna to create the query to subject orthologs dictionary
    reversed_query - a parameter that defines for the sake of this calculation only, what will be the order in which btop mismatches pairs are read.
    (the blastp data assumes that one animal is alwayes the query and the other is the subject, so redefining the animals is necessary when calculating for 
    editing sites of the animal in which the blastp define as subject)
    """
    
    queries_orthologs_dict = {}
    for q_comp in set(list(query_trinity_tbl.index)):
        
        q_animal_comp = query_animal+'|'+q_comp
        s_animal_comp = None
        
        #serching the 2-ways best hit ortholog of the subject animal
        if q_animal_comp in list(all_best_hits['seq1']):
            ortholog_comps = all_best_hits[all_best_hits['seq1']==q_animal_comp]
            s_animal_comp = ortholog_comps[ortholog_comps['animal_2']==subject_animal]['seq2'].iloc[0]
        elif q_animal_comp in list(all_best_hits['seq2']):
            ortholog_comps = all_best_hits[all_best_hits['seq2']==q_animal_comp]
            s_animal_comp = ortholog_comps[ortholog_comps['animal_1']==subject_animal]['seq1'].iloc[0]
        
        if s_animal_comp is None:
#            print(q_animal_comp + ' has no best hits in any animal')
            queries_orthologs_dict.update({q_comp:'no_orthologs'})       
        elif type(s_animal_comp) is not str:
#            print(q_animal_comp  + ' has no ortholog in ' + subject_animal)
            queries_orthologs_dict.update({q_comp:'no_orthologs'})
        else:
#            print(q_animal_comp + ' ortholog is ' + s_animal_comp)
            q_es = query_es_tbl[query_es_tbl['component']==q_comp]
            q_trinity_comp = query_trinity_tbl.loc[q_comp,:]
            
#            if reversed_query:
#                blastp_data = blastp[np.logical_and(blastp['comp_query']==s_animal_comp,blastp['comp_subject']==q_animal_comp)]
#            else:
#                blastp_data = blastp[np.logical_and(blastp['comp_query']==q_animal_comp,blastp['comp_subject']==s_animal_comp)]
            pair = [q_animal_comp,s_animal_comp]
            blastp_data = blastp[np.logical_and(blastp['comp_query'].isin(pair),blastp['comp_subject'].isin(pair))]
            blastp_data = blastp_data.loc[blastp_data['bitscore'].idxmax(),:]
            
            if blastp_data['comp_query']==q_animal_comp:
                ortholog_locations = find_corresponding_location_in_subject_coding_mrna(blastp_data, q_trinity_comp, q_es)
            elif blastp_data['comp_query']==s_animal_comp:
                ortholog_locations = find_corresponding_location_in_subject_coding_mrna(blastp_data, q_trinity_comp, q_es, reversed_query=True)
            else:
                raise Exception('Unexpected query animal - ' + str(blastp_data['comp_query']))
            
            queries_orthologs_dict.update({q_comp:(s_animal_comp,blastp_data['bitscore'],ortholog_locations)})
     
    return queries_orthologs_dict


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Finding ortholog sites for 2 species based on orthomcl and blastp results')
    run_parser = parser.add_argument_group('generate random AG table')
    run_parser.add_argument('-parent_dir', dest='parent_dir', action='store', required = True, help='path to parent directory')
    run_parser.add_argument('-q', dest='query_animal', action='store', required = True, help='name of query animal')
    run_parser.add_argument('-s', dest='subject_animal', action='store', required = True, help='name of query animal')
    run_parser.add_argument('-nucs_around_site', dest='nucs_around_site', action='store', default = '30', help='sequence lengthe to keep down\up stream to sites')
    arguments = parser.parse_args()
    
    parent_dir = arguments.parent_dir
    query_animal = arguments.query_animal
    subject_animal = arguments.subject_animal
    nucs_around_site = int(arguments.nucs_around_site)
    
    blast_file = parent_dir + 'goodProteins_w_btop.blast'
    all_best_hits_file = parent_dir + 'all_best_hits.txt'
    editing_sites_path = parent_dir + 'editing_sites/'
    trinity_files_path = parent_dir + 'trinity_comps/'

#    parent_dir = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/'
#    blast_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/oct_sep_blast.txt'
#    all_best_hits_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/all_best_hits.txt'
#    editing_sites_path = 'E:/RNA_editing_Large_files/orthomcl/editing_sites/'
#    trinity_files_path = 'E:/RNA_editing_Large_files/trinity_comps/old/'    
#    nucs_around_site = 30
#    query_animal = 'oct'
#    subject_animal = 'sep'
    
    output_dir = parent_dir + 'results/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    triniy_comp_query = trinity_files_path + 'orfs_' + query_animal + '.fa'
    triniy_comp_subject = trinity_files_path + 'orfs_' + subject_animal + '.fa'
    editing_sites_query = editing_sites_path + query_animal + '_editing_site_info.txt'
    editing_sites_subject = editing_sites_path + subject_animal + '_editing_site_info.txt'

    print('\nFinding ortholog sites for ' + query_animal + ' ' + subject_animal + ':')
#    try:
    print('Reading 2-way best hits file')
    all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score'])
    all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
    all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
    all_best_hits = all_best_hits[np.logical_or(np.logical_and(all_best_hits['animal_1']==query_animal, all_best_hits['animal_2']==subject_animal), np.logical_and(all_best_hits['animal_1']==subject_animal, all_best_hits['animal_2']==query_animal))]
    
    #reading the blastp results, the transcriptome fasta files and the editing sites tables (all outputs of these functions are datafrmae)
    print('Reading BLASTP results')
    blastp = read_blast_results(blast_file, query_animal, subject_animal)
    
    print('Reading query (' + query_animal + ') data - fasta transcriptome and editing sites table')
    query_trinity_tbl = read_trinity_mrna_files(triniy_comp_query)
    query_es_tbl = read_editing_sites_tbl(editing_sites_query)
    query_trinity_tbl['native_coding_sequence'] = query_trinity_tbl.apply(lambda row: retrive_native_coding_sequence_in_orfs(row,query_es_tbl), axis = 1)
    query_es_tbl['trinity_exists'] = query_es_tbl.apply(lambda row: row['component'] in list(query_trinity_tbl.index), axis = 1)
    query_es_tbl = query_es_tbl[query_es_tbl['trinity_exists']]
    query_es_tbl = query_es_tbl.apply(lambda row: calculate_site_mrna_motif_data(row, query_trinity_tbl, nucs_around_site), axis = 1)
    
    print('Reading subject (' + subject_animal + ') data - fasta transcriptome and editing sites table')
    subject_trinity_tbl = read_trinity_mrna_files(triniy_comp_subject)
    subject_es_tbl = read_editing_sites_tbl(editing_sites_subject)
    subject_trinity_tbl['native_coding_sequence'] = subject_trinity_tbl.apply(lambda row: retrive_native_coding_sequence_in_orfs(row,subject_es_tbl), axis = 1)
    subject_es_tbl['trinity_exists'] = subject_es_tbl.apply(lambda row: row['component'] in list(subject_trinity_tbl.index), axis = 1)
    subject_es_tbl = subject_es_tbl[subject_es_tbl['trinity_exists']]
    subject_es_tbl = subject_es_tbl.apply(lambda row: calculate_site_mrna_motif_data(row, subject_trinity_tbl, nucs_around_site), axis = 1)
    
    
    print('Finding corresponding location in subject (' + subject_animal + ') coding ortholog sequence for each query (' + query_animal + ') sequence')
    queries_orthologs_dict = create_ortholog_locations_dict(query_animal, subject_animal, query_trinity_tbl, subject_trinity_tbl, all_best_hits, blastp, query_es_tbl)
    print('Connecting ' + subject_animal + ' ortholog sites for editing sites in ' + query_animal)
    query_es_tbl = query_es_tbl.apply(lambda row: connect_ortholog_sites_for_querie_editing_sites(row, queries_orthologs_dict, subject_es_tbl, subject_animal), axis = 1)    
    
    print('Finding corresponding location in query (' + query_animal + ') coding ortholog sequence for each subject (' + subject_animal + ') sequence')
    subject_ortholog_dict = create_ortholog_locations_dict(subject_animal, query_animal, subject_trinity_tbl, query_trinity_tbl, all_best_hits, blastp, subject_es_tbl)
    print('Connecting ' + query_animal + ' ortholog sites for editing sites in ' + subject_animal)
    subject_es_tbl = subject_es_tbl.apply(lambda row: connect_ortholog_sites_for_querie_editing_sites(row, subject_ortholog_dict, query_es_tbl, query_animal), axis = 1)    
    
    col_to_merge = ['component', 'protein', 'location', 'mm_type', 'DNA_A', 'DNA_T', 'DNA_G', 'DNA_C',
                    'RNA_A', 'RNA_T', 'RNA_G', 'RNA_C', 'Trinity', 'RNA_coverage', 'DNA_coverage', 
                    'p_val', 'AA_before', 'AA_after', 'type', 'protein_length', 'editing_level', 'strand',
                    'site_coding_loc_base0','site_loc_in_codon_base0','original_codon','nuc','location_in_area','site_area','site_key']
    
    print('Merging processed editing sites tables')
    for col in query_es_tbl.columns:
        new_col_name = query_animal + '_' + col
        query_es_tbl = query_es_tbl.rename(index = str, columns = {col:new_col_name}) 
    for col in subject_es_tbl.columns:
        new_col_name = subject_animal + '_' + col
        subject_es_tbl = subject_es_tbl.rename(index = str, columns = {col:new_col_name}) 
    col_to_merge_subject_to_query = [subject_animal+'_'+x for x in col_to_merge]
    col_to_merge_query_to_subject = [query_animal+'_'+x for x in col_to_merge]
    
    query_es_tbl_all_data = query_es_tbl.merge(subject_es_tbl[subject_es_tbl[subject_animal+'_site_key'].isin(query_es_tbl[query_animal+'_'+subject_animal+'_ortholog_site_key'])][col_to_merge_subject_to_query], left_on = query_animal+'_'+subject_animal+'_ortholog_site_key', right_on = subject_animal+'_site_key',how = 'outer')
    subject_es_tbl_all_data = subject_es_tbl.merge(query_es_tbl[query_es_tbl[query_animal+'_site_key'].isin(subject_es_tbl[subject_animal+'_'+query_animal+'_ortholog_site_key'])][col_to_merge_query_to_subject], left_on = subject_animal+'_'+query_animal+'_ortholog_site_key', right_on = query_animal+'_site_key',how = 'outer')
    
    print('Calculating additional orthologs data')
    query_es_tbl_all_data = query_es_tbl_all_data.apply(lambda row: additional_match_data(row, query_animal, subject_animal), axis = 1)
    subject_es_tbl_all_data = subject_es_tbl_all_data.apply(lambda row: additional_match_data(row, subject_animal, query_animal), axis = 1)
    query_es_tbl_all_data[query_animal+'_'+subject_animal+'_conserved'] = query_es_tbl_all_data.apply(lambda row: row[query_animal+'_'+subject_animal+'_same_target_aa'], axis = 1)
    subject_es_tbl_all_data[subject_animal+'_'+query_animal+'_conserved'] = subject_es_tbl_all_data.apply(lambda row: row[subject_animal+'_'+query_animal+'_same_target_aa'], axis = 1)
    
    print('Calculating nucleotides data where ortholog locations are not edited')
    query_es_tbl_all_data = query_es_tbl_all_data.apply(lambda row: find_corresponding_unedited_nuc_data_in_subject_for(row, query_animal, subject_animal, subject_trinity_tbl, nucs_around_site), axis = 1)
    subject_es_tbl_all_data = subject_es_tbl_all_data.apply(lambda row: find_corresponding_unedited_nuc_data_in_subject_for(row, subject_animal, query_animal, query_trinity_tbl, nucs_around_site), axis = 1)
    
    print('Writing processed editing sites table (with orthologs data)')
    query_es_tbl_all_data.to_csv(output_dir+query_animal+'_editing_sites_'+subject_animal+'_orthologs.txt',sep = '\t', index = False)
    subject_es_tbl_all_data.to_csv(output_dir+subject_animal+'_editing_sites_'+query_animal+'_orthologs.txt',sep = '\t', index = False)
    
    print('Writing conserved orthologs for ' + query_animal + '_' + subject_animal)
    conserved_orthologs = query_es_tbl_all_data[query_es_tbl_all_data[query_animal+'_'+subject_animal+'_conserved']==True]
    conserved_orthologs.to_csv(output_dir+query_animal+'_'+subject_animal+'_conserved_sites.txt', sep = '\t', index = False)
    print(query_animal + ' ' + subject_animal +' orthologs matchig finished')
    
    print(str(len(query_es_tbl_all_data[query_es_tbl_all_data[query_animal+'_'+subject_animal+'_conserved']==True])) + ' conserved orthologs ' + query_animal + '_' + subject_animal) 
    print(str(len(query_es_tbl_all_data[query_es_tbl_all_data[query_animal+'_'+subject_animal+'_conserved']==False])) + ' unconserved orthologs ' + query_animal + '_' + subject_animal)
    print(str(len(subject_es_tbl_all_data[subject_es_tbl_all_data[subject_animal+'_'+query_animal+'_conserved']==True])) + ' conserved orthologs ' + subject_animal + '_' + query_animal) 
    print(str(len(subject_es_tbl_all_data[subject_es_tbl_all_data[subject_animal+'_'+query_animal+'_conserved']==False])) + ' unconserved orthologs ' + subject_animal + '_' + query_animal)
        
#    except Exception as e:
#        print(query_animal + ' ' + subject_animal +' failed\n')
#        print(e)
#        
#    
#    query_animal = 'oct'
#    subject_animal = 'sep'
#
#    new_analysis_orthologs_sites = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/results/oct_sep_sites.txt'
#    squ_oct_brand_new = pd.read_csv(new_analysis_orthologs_sites, sep = '\t', index_col = False)
#    oct_seo_noa_new = ortholog_sites_df
#    oct_sep_brand_new = ortholog_sites_df
#    oct_sep_noa_new[query_animal+'_'+subject_animal+'_same_target_aa'] = oct_sep_noa_new.apply(lambda row: row['Species #1 Codon changes'][1]==row['Species #2 Codon changes'][1], axis = 1)
#    oct_sep_noa_new[query_animal+'_'+subject_animal+'_same_original_aa'] = oct_sep_noa_new.apply(lambda row: row['Species #1 Codon changes'][0]==row['Species #2 Codon changes'][0], axis = 1)
#    oct_sep_brand_new = oct_sep_brand_new[np.logical_and(oct_sep_brand_new[query_animal+'_'+subject_animal+'_same_target_aa'], oct_sep_brand_new[query_animal+'_'+subject_animal+'_same_original_aa'])].copy()
#
#    oct_sep_brand_new['site_in_noa_new'] = oct_sep_brand_new.apply(lambda row: row['oct_sep_key'] in list(oct_sep_noa_new['oct_sep_key']), axis = 1)
#    oct_sep_brand_new['orhtolog_in_noa_new'] = oct_sep_brand_new.apply(lambda row: row['oct_sep_key'] in list(oct_sep_noa_new['oct_sep_key']), axis = 1)
#    oct_sep_noa_new['site_in_brand_new'] = oct_sep_noa_new.apply(lambda row: row['oct_sep_key'] in list(oct_sep_brand_new['oct_sep_key']), axis = 1)
#    oct_sep_noa_new['ortholog_in_brand_new'] = oct_sep_noa_new.apply(lambda row: row['orthologs'] in list(oct_sep_brand_new[query_animal+'_'+subject_animal+'_ortholog_sequences']), axis = 1)
#
#   
#    
#
#    
#    new_analysis_orthologs_sites = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/oct_sep_orthologs_sites.txt'
#    oct_sep_brand_new = pd.read_csv(new_analysis_orthologs_sites, sep = '\t', index_col = False)
#    oct_sep_noa_new = ortholog_sites_df
#    oct_sep_noa_new[query_animal+'_'+subject_animal+'_same_target_aa'] = oct_sep_noa_new.apply(lambda row: row['Species #1 Codon changes'][1]==row['Species #2 Codon changes'][1], axis = 1)
#    oct_sep_noa_new[query_animal+'_'+subject_animal+'_same_original_aa'] = oct_sep_noa_new.apply(lambda row: row['Species #1 Codon changes'][0]==row['Species #2 Codon changes'][0], axis = 1)
#    oct_sep_brand_new = oct_sep_brand_new[np.logical_and(oct_sep_brand_new[query_animal+'_'+subject_animal+'_same_target_aa'], oct_sep_brand_new[query_animal+'_'+subject_animal+'_same_original_aa'])]
#    
#    oct_sep_brand_new['site_in_noa_new'] = oct_sep_brand_new.apply(lambda row: row['oct_sep_key'] in list(oct_sep_noa_new['oct_sep_key']), axis = 1)
#    oct_sep_brand_new['orhtolog_in_noa_new'] = oct_sep_brand_new.apply(lambda row: row['oct_sep_key'] in list(oct_sep_noa_new['oct_sep_key']), axis = 1)
#    oct_sep_noa_new['site_in_brand_new'] = oct_sep_noa_new.apply(lambda row: row['oct_sep_key'] in list(oct_sep_brand_new['oct_sep_key']), axis = 1)
#    oct_sep_noa_new['ortholog_in_brand_new'] = oct_sep_noa_new.apply(lambda row: row['orthologs'] in list(oct_sep_brand_new[query_animal+'_'+subject_animal+'_ortholog_sequences']), axis = 1)
#
#    blastp_data = blastp[np.logical_and(blastp['comp_query']=='oct|comp101167_c0_seq1',blastp['comp_subject']=='sep|comp241967_c0_seq1')]
#    blastp_data = blastp_data.loc[blastp_data['identity'].idxmax(),:]
#    a = queries_orthologs_dict['comp101167_c0_seq1']
#    q_comp = 'comp101167_c0_seq1'
#    s_comp = 'comp241967_c0_seq1'
#    q_es = query_es_tbl[query_es_tbl['component']==q_comp]
#    q_trinity_comp = query_trinity_tbl.loc[q_comp,:]
#    s_trinity_comp = subject_trinity_tbl.loc[s_comp,:]
#
#for k,v in subject_ortholog_dict.items():
#    if v[2] == {}:
#        print(k)
    
#
#for col in query_es_tbl.columns:
#    new_col_name = col.replace('oct_','')
#    query_es_tbl = query_es_tbl.rename(index = str, columns = {col:new_col_name})

