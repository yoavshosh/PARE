# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 18:51:27 2019

@author: shosh

This script iterate over multiple files containing msa results for super orthologs in many species.
It creates a table in which for each row represent a nucleotide, in which at least on of the animal had an editing event
and also contain data regarding all corresponding nucleotides in the other species (whether or not they are edited)
"""

import os
import re
import glob
import pandas as pd
import numpy as np
import itertools as it
import multiprocessing as mp
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
import argparse
import warnings


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def read_trinity_mrna_files(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = ['component','protein','orfs_start','orfs_end','strand','sequence']
    
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split('\t')
        if rec_data[0][-1] == ' ':  #some fastafiles have spaces after each id, so fixing it here.
            rec_data[0] = rec_data[0].replace(' ','')
        protein = rec_data[-1].split('|')[1]   #reading proteing from description assuming it was added to header using the transcriptome built pipeline we have for trinity
        rec_data = (rec_data[0],protein,int(rec_data[2]),int(rec_data[4]),rec_data[6],record.seq)
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

def create_super_orthologs_df(all_best_hits,animals):
    """
    from all_best_hits dataframe containing rows of 2-way best hits orthologs for pairs of animals
    return a dataframe in which each row contain sequences for multiple animals that are all the 2-way best hits of eachother.
    """

    super_orthologs = None
    for pair in it.combinations(animals,2):
        a1 = pair[0]
        a2 = pair[1]
        pair_best_hits = all_best_hits[np.logical_or(np.logical_and(all_best_hits['animal_1']==a1,all_best_hits['animal_2']==a2),np.logical_and(all_best_hits['animal_1']==a2,all_best_hits['animal_2']==a1))]

        #swapping a2 and a2 (if needed) to fit correct order of columns
        a1 = pair_best_hits['animal_1'].iloc[0]
        a2 = pair_best_hits['animal_2'].iloc[0]

        #keeping only relevant columns in new names
        pair_best_hits = pair_best_hits[['score','seq1','seq2']]
        pair_best_hits = pair_best_hits.rename(index = str, columns = {'score':a1+'_'+a2+'_score','seq1':a1,'seq2':a2}) 
        
        if super_orthologs is None:
            super_orthologs = pair_best_hits
        else:
            if all(a in super_orthologs.columns for a in [a1,a2]):
#                print('both animals are already in df')
                super_orthologs = super_orthologs.merge(pair_best_hits, on = [a1 ,a2], how = 'inner')
            elif a1 in super_orthologs.columns and a2 not in super_orthologs.columns:
#                print('animal_1 (' + a1 + ') is already in df')
                super_orthologs = super_orthologs.merge(pair_best_hits, on = a1 , how = 'inner')
            elif a2 in super_orthologs.columns and a1 not in super_orthologs.columns:
#                print('animal_2 (' + a2 + ') is already in df')
                super_orthologs = super_orthologs.merge(pair_best_hits, on = a2, how = 'inner')
            if all(a not in super_orthologs.columns for a in [a1,a2]):
                print('NO INTERSECTION OF ANIMALS')
    
    return super_orthologs


def calculate_site_mrna_motif_data(row, trinity_df):
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
    
    row['site_coding_loc_base0'] = site_coding_loc_base_0
    row['site_loc_in_codon_base0'] = site_pos_in_codon
    row['original_codon'] = codon
    row['site_key'] = row['component'] + ';' + str(row['location'])
    
    return row


def cluster_editing_sites(es_df_of_one_gene,c_range=50):
    """
    for each site in gene, cluster sites together if all are within maximum distance of c_range nucl 
    mark the strongest editing level site in each cluster
    """
    def cluster_sites(es_df_of_one_gene,c_range):
        
        clust_n = 1
        sites_clusters = []
        component = list(es_df_of_one_gene['component'])[0]
        sorted_site_locs_series = es_df_of_one_gene['site_coding_loc_base0'].sort_values()
        
        last_site_loc = sorted_site_locs_series.values[0]
        current_cluster = component+'_clust_'+str(clust_n)
        
        for i, loc in sorted_site_locs_series.iteritems():
            if loc-last_site_loc<=c_range:
                sites_clusters.append((current_cluster,loc))
            else:
                clust_n+=1
                current_cluster = component+'_clust_'+str(clust_n)
                sites_clusters.append((current_cluster,loc))
            last_site_loc = loc
                
        return pd.DataFrame(data=sites_clusters, columns=['range_'+str(c_range)+'_cluster', 'site_coding_loc_base0'])
                
    es_df_of_one_gene = es_df_of_one_gene.merge(cluster_sites(es_df_of_one_gene,c_range), on='site_coding_loc_base0', how='inner')
    es_df_of_one_gene['strongest_in_range_'+str(c_range)+'_cluster'] = es_df_of_one_gene.apply(lambda row: 1 if row.name==es_df_of_one_gene[es_df_of_one_gene['range_'+str(c_range)+'_cluster']==row['range_'+str(c_range)+'_cluster']]['editing_level'].idxmax() else 0, axis=1)
    
    return es_df_of_one_gene


def parse_msa_results(msa_num, alignment, relevant_editing_sites_dict,relevant_transcripts_dict,only_edited_position,proteins_msa,c_range,animals_sets_to_check_msa_identity, areas_for_msa_consensus_ratio=[7,10,16,22]):
    """
    This function parse the blast results, it iterate through the alignment matrix columns
    and for each column check if the nucleotides in it correspond to an editing site.
    if so, it mapps all the locations from all animals that correspond to that column
    and retrive the relevant data for it
    
    it returns a dataframe where each row is a location and columns are data regarding this location within each animal
    
    proteins_msa - a flag. if true then parsing amino acids msa and for each column scans all positions in codon.
    """
        
    def calc_msa_consensus_ratio_in_area(alignment,position,alignment_length,area_around_position,alignment_animals,animals_to_check):
        """
        Calculate precentage of msa positions in area_around_position that are in agreement for all animals
        """
        consensus_msa_position=0
        no_consensus_msa_position=0
        for i in range(max(0,position-area_around_position),min(position+area_around_position+1,alignment_length)):
            align_col = alignment[:,i]
            compare_list = []
            for j,a in enumerate(alignment_animals):
                if a in animals_to_check:
                    compare_list.append(align_col[j])
            if len(set(compare_list))==1 and '-' not in compare_list:
                consensus_msa_position+=1
            else:
                no_consensus_msa_position+=1
                
        return round(float(consensus_msa_position)/(consensus_msa_position+no_consensus_msa_position),3)
    
    def calc_distance_to_nearest_recoding_editing_site(nuc_coding_location,recoding_sites_coding_locations):
        if len(recoding_sites_coding_locations):
            return nuc_coding_location-recoding_sites_coding_locations[np.abs(recoding_sites_coding_locations-nuc_coding_location).idxmin()]
        else:    
            return '-'
           
            
    def calc_aa_change(nuc,mut_nuc,sequence,location,codon_loc):
        """
        This function checks how a given mutation changes the aa
        it returns the original aa, the target aa, and 0/1 for wether the swap is syn or non-syn
        """
        if codon_loc == 0:
            codon = nuc+sequence[location+1:location+3]
            mut_codon = mut_nuc+sequence[location+1:location+3]
        elif codon_loc == 1:
            codon = sequence[location-1]+nuc+sequence[location+1]
            mut_codon = sequence[location-1]+mut_nuc+sequence[location+1]
        elif codon_loc == 2:
            codon = sequence[location-2:location]+nuc
            mut_codon = sequence[location-2:location]+mut_nuc
        
        aa_before = str(Seq(codon, generic_dna).translate())
        aa_after = str(Seq(mut_codon, generic_dna).translate())
        if aa_before==aa_after:
            recoding=0
        else:
            recoding=1
        
        return aa_before, aa_after, recoding
        
    def calc_all_possible_nucs_changes(nuc,sequence,location,codon_loc):
        """
        run calc_aa_change for each missmatch type and return a dict of format mm:(recoding,aa_change)
        """
        aa_changes_dict = {}
        if nuc=='gap':
            for mm in all_mm:
                aa_changes_dict.update({mm:(0,'-')})
        else:    
            for mm in all_mm:
                if mm[0]==nuc:  
                    aa_before, aa_after, recoding = calc_aa_change(nuc,mm[1],sequence,location,codon_loc)
                    aa_changes_dict.update({mm:(recoding,aa_before+aa_after)})
                else:
                    aa_changes_dict.update({mm:(0,'-')})
        
        return aa_changes_dict
    
    
    alignment_components = [alignment[i].id for i in range(len(alignment))] #rows animals order in alignment matrix     
    alignment_animals = [comp[:3] for comp  in alignment_components]
    alignment_length = alignment.get_alignment_length()
    data = []

    columns = [a+'_'+col for a in alignment_animals for col in per_animals_columns_suffixes]
    columns += ['super_orthologs_id','alignment_length','alignment_pos_base_0']
    msa_type = 'aa' if proteins_msa else 'nucl'
    for animals_set in animals_sets_to_check_msa_identity:
        animals_str = '_'.join(animals_set)
        columns += [animals_str+'_consensus_ratio_'+str(l)+'_'+msa_type+'_range' for l in areas_for_msa_consensus_ratio]
    
    #initializing animals counter dict - to store current location (base0) for each animal's sequence in the MSA
    nucs_counters_dict = {}
    for a in alignment_animals:
        nucs_counters_dict.update({a:-1})
    
    if proteins_msa:
        nucs_per_column = 3
    else:
        nucs_per_column = 1
    
    for align_loc in range(alignment.get_alignment_length()):
        
        #calculate the ratio of msa positions that are all in consensus, withing some ranges of current align_loc and some sets of animals for which calculation should be preformed
        consensus_ratios = []
        for animals_set in animals_sets_to_check_msa_identity:
            for l in areas_for_msa_consensus_ratio:
                if (2*l+1)%3!=0 and proteins_msa==False:
                    warnings.warn("MSA identity tests: number of nucl in range is not a multiplet of 3. nucl around each position is "+str(l))
                consensus_ratios.append(calc_msa_consensus_ratio_in_area(alignment,align_loc,alignment_length,l,alignment_animals,animals_set))
                      
        for codon_position in range(nucs_per_column): #if proteins_msa, then iterating over all 3 nucleotides composing aa codon
        
            alingnment_nucl_position_base_0=align_loc*3+codon_position
#            print(alingnment_nucl_position_base_0)
            editing_sites_in_msa_loc = 0
            animals_current_loc = {}
            
            for i,a in enumerate(alignment_animals):
                if alignment[i,align_loc] != '-':  #if true then current animal location is not a gap
                    current_nuc_loc = nucs_counters_dict[a]+1  
                    nucs_counters_dict[a] = current_nuc_loc #updating current nuc location in animal sequence if not a gap
                    es_df_for_animal = relevant_editing_sites_dict[a]
                    if len(es_df_for_animal):
                        animals_es_row = es_df_for_animal[es_df_for_animal['site_coding_loc_base0']==current_nuc_loc].squeeze()
                        if len(animals_es_row): #if there is an editing site, animals_es_row will be a row series, otherwise, an empty dataframe (len = 0)
                            editing_sites_in_msa_loc += len(animals_es_row)    
                            animals_current_loc.update({a:animals_es_row})
                        else:
                            animals_current_loc.update({a:current_nuc_loc})
                    else:
                        animals_current_loc.update({a:current_nuc_loc})
                else:
                    animals_current_loc.update({a:'gap'})
            
            #only append data if at least one of the animals is edited in the current location (maybe will change this later)
            if editing_sites_in_msa_loc or not(only_edited_position):
                
                location_data = ()
                
                for a in alignment_animals:
    
                    loc = animals_current_loc[a]
                    if type(loc) == pd.core.series.Series:  #if this is True then loc is a row sereis from some editing sites dataframe of one of the animals
                        es_row = loc
                        animal_comp = alignment_components[alignment_animals.index(a)]
                        protein = relevant_transcripts_dict[a]['protein']
                        edited = 1
                        coding_location = es_row['site_coding_loc_base0']
                        nuc = relevant_transcripts_dict[a]['native_coding_sequence'][coding_location]
                        loc_inside_codon = es_row['site_loc_in_codon_base0']
                        site_key = es_row['site_key']
                        editing_level = es_row['editing_level']
                        
                        if es_row['AA_before']!=es_row['AA_after']: #the specigic loc examined int this block is in itself a recoding site
                            distance_to_nearest_recoding_site = 0
                        else:
                            editing_sites_for_animal = relevant_editing_sites_dict[a]
                            recoding_editing_sites_locations = editing_sites_for_animal[[r['AA_after']!=r['AA_before'] for i,r in editing_sites_for_animal.iterrows()]]
                            distance_to_nearest_recoding_site = calc_distance_to_nearest_recoding_editing_site(coding_location,recoding_editing_sites_locations['site_coding_loc_base0'])
                        strongest_in_gene = es_row['strongest_in_gene']
                        cluster = es_row['range_'+str(c_range)+'_cluster']
                        strongest_in_cluster = es_row['strongest_in_range_'+str(c_range)+'_cluster']
                        
                    elif type(loc)==int:
                        animal_comp = alignment_components[alignment_animals.index(a)]
                        protein = relevant_transcripts_dict[a]['protein']
                        edited = 0
                        coding_location = loc
                        nuc = relevant_transcripts_dict[a]['native_coding_sequence'][loc]
                        loc_inside_codon = loc%3
                        site_key = '-'
                        editing_level = 0.0
                
                        editing_sites_for_animal = relevant_editing_sites_dict[a]
                        if len(editing_sites_for_animal):
                            recoding_editing_sites_locations = editing_sites_for_animal[[r['AA_after']!=r['AA_before'] for i,r in editing_sites_for_animal.iterrows()]]
                            distance_to_nearest_recoding_site = calc_distance_to_nearest_recoding_editing_site(coding_location,recoding_editing_sites_locations['site_coding_loc_base0'])
                        else:
                            distance_to_nearest_recoding_site = '-'
                        strongest_in_gene = '-'
                        cluster = '-'
                        strongest_in_cluster = '-'
                        
                    elif loc == 'gap':
                        animal_comp = alignment_components[alignment_animals.index(a)]
                        protein = relevant_transcripts_dict[a]['protein']
                        edited = 0
                        coding_location = 'gap'
                        nuc = 'gap'
                        loc_inside_codon = 'gap'
                        site_key = 'gap'
                        editing_level = 0.0
                        distance_to_nearest_recoding_site = 'gap'
                        strongest_in_gene = '-'
                        cluster = '-'
                        strongest_in_cluster = '-'
                    
                    else: #this is not expected to be True, ever
                        print('type of data for animal ' + a +' is not identified')
                      
                    aa_changes_dict = calc_all_possible_nucs_changes(nuc,relevant_transcripts_dict[a]['native_coding_sequence'],coding_location,loc_inside_codon)
                    if nuc not in ['gap','A','C','G','T']:
                        warnings.warn('ilegal nucl type in msa number ' + str(msa_num) + ': ' + str(nuc))
                        
                    location_data += (animal_comp,protein,edited,coding_location,nuc,loc_inside_codon,site_key,editing_level,distance_to_nearest_recoding_site,strongest_in_gene,cluster,strongest_in_cluster)
                    
                    for mm in all_mm:
                        mm_data = aa_changes_dict[mm]
                        location_data+=(mm_data[0],mm_data[1],)
                
                location_data += (msa_num,alignment_length,alingnment_nucl_position_base_0)
                location_data += tuple(consensus_ratios)
                data.append(location_data)
    
    return pd.DataFrame(columns = columns, data = data)
    
                
def iterate_alignment_files_and_connect_sites(msa_path,super_ortholog_editing_sites_dict,super_ortholog_transcripts_dict, n_workers, only_edited_position, cluster_range, animals_sets_to_check_msa_identity, proteins_msa = False):
    """
    This function iterate over files in the msa_path (files are containing msa resutlts in fasta format)
    for each file it filter the editing sites from super_ortholog_editing_sites_dict relevent for the specific sequences
    and select relevant trinity data from super_ortholog_transcripts_dict
    if there are at least one editing site for the relevant sequences it calls parse_msa_results
    it does so in parallel using n_workers chiled peocesses 
    """
    pool = mp.Pool(n_workers)
    processes_results = []
    
    msa_files = natural_sort(glob.glob(os.path.join(msa_path,'msa_results*.fasta')))
    locations_dataframes = []

    for n, file in enumerate(msa_files):
        alignment = AlignIO.read(open(file), "fasta")
        alignment_components = [alignment[i].id for i in range(len(alignment))] #rows animals order in alignment matrix     
        alignment_animals = [comp[:3] for comp  in alignment_components]
        
        #filtering out only relevat compopnents in editing site and creating relevant_transcripts_dict
        total_es_in_animals = 0
        relevant_editing_sites_dict = {}
        relevant_transcripts_dict = {}
        for i,a in enumerate(alignment_animals):
            es_df = super_ortholog_editing_sites_dict[a]
            relevant_es_df = es_df[es_df['component']==alignment_components[i]].copy()
            if len(relevant_es_df):
                relevant_es_df['strongest_in_gene'] = relevant_es_df.apply(lambda row: 1 if row.name==relevant_es_df['editing_level'].idxmax() else 0, axis=1)
                relevant_es_df = cluster_editing_sites(relevant_es_df,c_range=cluster_range)
            relevant_editing_sites_dict.update({a:relevant_es_df})
            total_es_in_animals += len(relevant_editing_sites_dict[a])
            relevant_transcripts_dict.update({a:super_ortholog_transcripts_dict[a].loc[alignment_components[i],:]})
            
        if only_edited_position and total_es_in_animals==0:
            print('No editing sites for sequences in ' + file.split('/')[-1])
        else:
            print('Reading MSA resutls in ' + file.split('/')[-1])
#            locations_dataframes.append(parse_msa_results(alignment, relevant_editing_sites_dict, relevant_transcripts_dict))
            processes_results.append(pool.apply_async(parse_msa_results, (n+1, alignment, relevant_editing_sites_dict, relevant_transcripts_dict, only_edited_position, proteins_msa, cluster_range, animals_sets_to_check_msa_identity)))
            
    pool.close()
    pool.join()
    locations_dataframes = [r.get() for r in processes_results]
    
    return pd.concat(locations_dataframes,sort=False,ignore_index=True, axis=0)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Mapping of all editing sites in a bunch of animals supr ortholog genes to locations in all other animals')
    run_parser = parser.add_argument_group('Analize MSA results and construct a table for all editing events and corresponding locations in other animals')
    run_parser.add_argument('-parent_dir', dest='parent_dir', action='store', required = True, help='path to parent directory')
    run_parser.add_argument('-super_orthologs_msa_dir', dest='super_orthologs_msa_dir', action='store', required = True, help='sub dir within parent dir, in which super orthologs fastas and msa results')
    run_parser.add_argument('-animals', dest='animals', action='store', nargs = '+', default = ['oct','bim','squ','sep','bob','lin','nau'], help='subset of animals for finding super orthologs and perform MSA')
    run_parser.add_argument('-clust_range', dest='cluster_range', action='store', default = '100', help='editing sites clustering range - all sites in cluster have a neighbor in the cluster not further than cluster_range')
    run_parser.add_argument('-n_workers', dest='n_workers', action='store', default = '50', help='number of chiled processes for msa results analysis')
    run_parser.add_argument('-proteins_msa', dest='proteins_msa', action='store', default = 'False', help='read and parse muscle msa for proteins sequences')
    run_parser.add_argument('-only_edited_position', dest='only_edited_position', action='store', default = 'False', help='create table containing only msa columns (as rows) for which at least one animal is edited')
    run_parser.add_argument('-animals_order', dest='animals_order', action='store', nargs = '+', default = ['oct','bim','sep','squ','lin','bob','nau','apl'], help='animals in order for ordering columns after run')
    arguments = parser.parse_args()
    
    parent_dir = arguments.parent_dir
    super_orthologs_msa_dir = arguments.super_orthologs_msa_dir
    animals = arguments.animals
    cluster_range = int(arguments.cluster_range)
    n_workers = int(arguments.n_workers)
    proteins_msa = eval(arguments.proteins_msa)
    only_edited_position = eval(arguments.only_edited_position)
    animals_order = arguments.animals_order

#    parent_dir = 'C:/Users/shosh/OneDrive/Desktop/parent_test_dir/'
#    super_orthologs_msa_dir = 'super_ortholog_proteins_fasta_files_all8_'
#    animals = ['apl','oct','bim','squ','sep','bob','lin','nau']
#    cluster_range = 50
#    n_workers = 30
#    proteins_msa = True
#    only_edited_position = False

    per_animals_columns_suffixes = ('component','protein','edited','coding_location','nuc','loc_inside_codon','site_key','editing_level','distance_to_nearest_recoding_site','strongest_in_gene','range_'+str(cluster_range)+'_cluster','strongest_in_range_'+str(cluster_range)+'_cluster')
    all_mm = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG']
    aa_changes_columns = ()
    for mm in all_mm:
        aa_changes_columns+=('recoding_'+mm,mm+'_aa_change',)
    per_animals_columns_suffixes+=aa_changes_columns

    animals_sets_to_check_msa_identity = [tuple(sorted(tuple(animals)))]
    if 'nau' in animals:
        if 'apl' in animals:
            animals_sets_to_check_msa_identity.append(tuple(sorted([a for a in animals if a!='apl'])))
        animals_sets_to_check_msa_identity.append(tuple(sorted([a for a in animals if a not in ['apl','nau']])))
    minus_one_sets_combinations = []
    for l in animals_sets_to_check_msa_identity:
        minus_one_sets_combinations += [tuple(sorted(k)) for k in list(it.combinations(l,len(l)-1))]
    animals_sets_to_check_msa_identity+=minus_one_sets_combinations
    animals_sets_to_check_msa_identity = list(set(animals_sets_to_check_msa_identity))
    
    
    all_best_hits_file = parent_dir + 'all_best_hits.txt'
    editing_sites_path = parent_dir + 'editing_sites/'
    trinity_files_path = parent_dir + 'trinity_comps/'
    msa_results_path = parent_dir + super_orthologs_msa_dir+'/msa_results/'
    results_path = parent_dir+'all_nucs_matrix_from_'+super_orthologs_msa_dir+'/'
    if not os.path.exists(results_path):
        os.makedirs(results_path)
    else:
        warnings.warn('Output dir already exists and existing files will be overwritten')
    
    print('Reading 2-way best hits file for animals: ' + str(animals))
    all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score'])
    all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
    all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
    
    super_ortholog_genes = create_super_orthologs_df(all_best_hits,animals)
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
         
    #collencting fasta records for super ortholog sequences only
    print('Collecting super ortholog records from fasta transcriptomes')
    super_ortholog_transcripts_dict = {}    
    for a in animals:
        trinity_df = read_trinity_mrna_files(trinity_files_path+'orfs_'+a+'.fa')
        trinity_df['component'] = trinity_df.apply(lambda row: a+'|'+row['component'], axis=1)
        trinity_df.set_index('component', inplace = True)
        super_orthologs_trinity_df = trinity_df[trinity_df.index.isin(super_ortholog_genes[a])].copy()
        super_orthologs_trinity_df['native_coding_sequence'] = super_orthologs_trinity_df.apply(lambda row: retrive_native_coding_sequence_in_orfs(row,super_ortholog_editing_sites_dict[a]), axis = 1)
        super_ortholog_transcripts_dict.update({a:super_orthologs_trinity_df})
        
    print('Calculating additional data for editing sites')
    for a in animals:
        es_df = super_ortholog_editing_sites_dict[a]
        trinity_df = super_ortholog_transcripts_dict[a]
        es_df = es_df.apply(lambda row: calculate_site_mrna_motif_data(row, trinity_df), axis=1)
        super_ortholog_editing_sites_dict[a] = es_df
        
    print('Reading MSA results and mapping editing sites')
    grand_locations_table = iterate_alignment_files_and_connect_sites(msa_results_path,super_ortholog_editing_sites_dict,super_ortholog_transcripts_dict,n_workers,only_edited_position,cluster_range,animals_sets_to_check_msa_identity,proteins_msa=proteins_msa)
    print('Counting nucs types per row')
    grand_locations_table['edited_animals'] = grand_locations_table.apply(lambda row: sum([row[a+'_edited'] for a in animals]), axis=1)
    grand_locations_table['gapped_animals'] = grand_locations_table.apply(lambda row: sum([row[a+'_nuc']=='gap' for a in animals]), axis=1)
    grand_locations_table['A_animals'] = grand_locations_table.apply(lambda row: sum([row[a+'_nuc']=='A' for a in animals]), axis=1)
    grand_locations_table['C_animals'] = grand_locations_table.apply(lambda row: sum([row[a+'_nuc']=='C' for a in animals]), axis=1)
    grand_locations_table['G_animals'] = grand_locations_table.apply(lambda row: sum([row[a+'_nuc']=='G' for a in animals]), axis=1)
    grand_locations_table['T_animals'] = grand_locations_table.apply(lambda row: sum([row[a+'_nuc']=='T' for a in animals]), axis=1)
    
    #sortin columns to a specific order for competability of differnet runs outputs
    columns_ordered = ['super_orthologs_id','alignment_length','alignment_pos_base_0','edited_animals','gapped_animals','A_animals','C_animals','G_animals','T_animals']
    for a in animals_order:
        if a in animals:
            columns_ordered += [a+'_'+c for c in per_animals_columns_suffixes]
    for c in grand_locations_table.columns:
        if 'consensus_ratio' in c:
            columns_ordered += [c]
            
    print('Writing table')
    if only_edited_position:
        grand_locations_table[columns_ordered].to_csv(results_path+'sites_mappings.txt', sep='\t', index=False)
    else:
        grand_locations_table[columns_ordered].to_csv(results_path+'all_nucs_mappings.txt', sep='\t', index=False)
    print('Finished')
    
        

