# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 15:32:46 2020

@author: shosh

This script runs over paml4 rst output file, gets ancestors notations, and creates msa fasta file with ancestral sequences
"""

import re
import os
import  glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
import argparse

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='parsing paml4 rst file and creating msa fasta file with ancestors sequentialy for al genes')
    run_parser = parser.add_argument_group('Parse PAML4 results sequentialy for multiple genes')
    run_parser.add_argument('-msas_dirs_path', dest='msas_dirs_path', action='store', required = True, help='path to parent directory in which all codons-msas for all suprer orthologs and paml4 results, in seperate dirs')
    run_parser.add_argument('-tree', dest='tree_file', action='store', required = True, help='phylogenetic tree file including ancestors')
    arguments = parser.parse_args()
    
    msas_dirs_path = arguments.msas_dirs_path
    tree_file = arguments.tree_file
#    msas_dirs_path = 'C:/Users/shosh/OneDrive/Desktop/results/all8_ncbi/'
#    tree_file = 'C:/Users/shosh/OneDrive/Desktop/results/all8_ncbi/rooted_tree_all8_ncbi.txt'
    
    
    dirs = natural_sort(glob.glob(os.path.join(msas_dirs_path,'*/')))
    
    branches=[]
    distances = []
    question_marks = [] #for some reaone codeml decide to replace few codons (usually in alignment gaps). would like to quantify this phenomena
    for d in dirs:
        d=d.replace("\\",'/')
        dn=d.split('/')[-2]
        print('Parsing rst file in sub dir '+dn)
        msa_file = d+ 'codons_msa_for_super_orthologs_'+dn+'.fasta'
        rst_file = d+'rst'
        mlb_file = d+'mlb'
        
        msa_file_w_ancestors =  d+'/ancestors_and_leafs_'+msa_file.split('/')[-1]
        ancestors_nucl_confidence = d+'/ancestors_nucl_confidence_'+str(dn)
        
        leafs_sequences = []
        for record in SeqIO.parse(open(msa_file, "r"), "fasta"):
            leafs_sequences.append(record.id)
        
        ances_original_regex_str = ';(?<=;)(.*?)(?=[,)]|$)'
        animal_original_regex_str = '(?<=[(|,])(.*?)(?=[,|)])'
        animal_order_str = '(?<=[(,])(.*?)(?=_)'
        ances_original_regex = re.compile(ances_original_regex_str)
        animals_original_regex = re.compile(animal_original_regex_str)
        animal_order_regex = re.compile(animal_order_str)
        
        f_tree = open(tree_file, "r")
        original_tree = f_tree.readlines()[0].rstrip()
        f_tree.close()
        
        original_ance_in_tree_strs = [m.group() for m in ances_original_regex.finditer(original_tree)]
        original_animals_in_tree_strs = [m.group().replace('(','') for m in animals_original_regex.finditer(original_tree)]
        phylo_tree_no_ances = original_tree
        for ance in original_ance_in_tree_strs:
            phylo_tree_no_ances = phylo_tree_no_ances.replace(ance,'')
        original_ance_in_tree_strs = [ance.replace(';','') for ance in original_ance_in_tree_strs]
    
        f_rst = open(rst_file, "r")
        content = f_rst.readlines()
        f_rst.close()
        for i,l in enumerate(content):
            if l.rstrip()=="tree with node labels for Rod Page's TreeView":
                t=content[i+1].rstrip()
                t=t.replace(' ','')
                if t[-1]==';':
                    t=t[:-1]
                break
        
        animals_numbers = [m.group().replace('(','') for m in animal_order_regex.finditer(t)]
        original_animals_in_tree_strs_sorted = [e for _,e in sorted(zip(animals_numbers, original_animals_in_tree_strs))]
        
        ances_paml_regex_str = '(?<=[)])(.*?)(?=[,)]|$)'
        ances_paml_regex = re.compile(ances_paml_regex_str)
        paml_ance_in_tree_srts = [m.group() for m in ances_paml_regex.finditer(t)]
        
        #finding all sequences including ancestral reconstructed sequences and storing in a list
        for i,l in enumerate(content):
            if l.rstrip()=="List of extant and reconstructed sequences":
                break
        sequences=[]        
        for j,l in enumerate(content[i+4:]):
            sequences.append(l.rstrip()) 
            if j+1==len(leafs_sequences)+len(paml_ance_in_tree_srts):
                break
        
        #writing msa with ancestors into fasta file
        print('Writing full tree msa fasta file')
        writer =  FastaWriter(open(msa_file_w_ancestors, 'w'), wrap=None)
        writer.write_header()
        for s in sequences:
            sl=s.rsplit()
            if sl[0] in leafs_sequences:
                seq=''.join(sl[1:])
                question_marks.append(len(seq.split('?', -1))-1)
                seq=''.join(sl[1:]).replace('?','-')
                rec_id=sl[0]
                writer.write_record(SeqRecord(Seq(seq, generic_dna), id = rec_id, description=''))
            elif sl[0]=='node':
                seq=''.join(sl[2:])
                rec_id= original_ance_in_tree_strs[paml_ance_in_tree_srts.index(sl[1].replace('#',''))]+'|'+dn
                writer.write_record(SeqRecord(Seq(seq, generic_dna), id = rec_id, description=''))
            else:
                raise NameError(sl[0]+'is not a recognized part of sequence id')
        writer.write_footer()
        msa_length = len(seq)
        
        #finding probability of reconstructed nucl and creating a dataframe  
        ance_aa_prob_regex_str = '(?<=[(][A-Z]\s)(.*?)(?=[)])'
        ance_aa_prob_regex = re.compile(ance_aa_prob_regex_str)
        ance_codon_prob_regex_str = '(?<=\s[A-Z]\s)(.*?)(?=[\s][(])'
        ance_codon_prob_regex = re.compile(ance_codon_prob_regex_str)
        data = []
        sorted_original_ance = [x for _,x in sorted(zip([int(k) for k in paml_ance_in_tree_srts],original_ance_in_tree_strs))]        
        columns = ['position_base_0','frequency']+[a+'_codon_prob' for a in sorted_original_ance]+[a+'_aa_prob' for a in sorted_original_ance]
        for i,l in enumerate(content):
            if l.rstrip()=="Prob of best state at each node, listed by site":
                break
        for j,l in enumerate(content[i+4:]):
            line = l.rstrip().lstrip()
            aa_position = int(line.split(' ')[0])
            for k in range(3):    
                position_base_0 = (aa_position-1)*3+k
                position_data = ()    
                while '  ' in line:
                    line = line.replace('  ',' ')
                probs = line.split(':')[-1].rstrip().lstrip()
                position_data+=(position_base_0,int(line.split(' ')[1]))
                for codon_prob in ance_codon_prob_regex.findall(probs):
                    position_data+=(float(codon_prob),)
                for aa_prob in ance_aa_prob_regex.findall(probs):
                    position_data+=(float(aa_prob),)
                data.append(position_data)
            if j*3+3==msa_length:
                break
        
        print('Writing ancestors nucl confidence table')
        nucl_prob_df = pd.DataFrame(data=data, columns=columns)
        nucl_prob_df.to_csv(ancestors_nucl_confidence,sep='\t',index=False)
        
        #finiding branches length
        f_mlb = open(mlb_file, "r")
        content = f_mlb.readlines()
        f_mlb.close()
        for i, l in enumerate(content):
            if "lnL(ntime:" in l.rstrip():
                break
        branches.append(content[i+1].rstrip().lstrip())
        distances.append([float(x) for x in content[i+2].rstrip().lstrip().split(' ')])
    
    
    if len(set(branches))>1:
        print('WARNING!! Branches are not in identical order for all mlb files!')
    branches=branches[0]
    while '  ' in branches:
        branches=branches.replace('  ',' ')
    branches = branches.split(' ')
    
    translated_branches=[]
    all_nods = original_animals_in_tree_strs_sorted+sorted_original_ance
    
    distances = np.array(distances)
    print(branches)
    for b in branches:
        b=b.split('..')
        translated_branches.append(all_nods[int(b[0])-1]+'->'+all_nods[int(b[1])-1])
    
    print(translated_branches)
    print(np.mean(distances,axis=0)[:len(translated_branches)])   
    print(str(sum(question_marks))+' ? characters in all sequences together')
    print('Finished')

#9..1     9..10   10..2    10..11   11..12   12..3    12..4    11..13   13..14   14..5    14..6    13..15   15..7    15..8        
#2.021301 1.955118 0.442719 0.282412 0.415505 0.007813 0.017725 0.192185 0.000004 0.063779 0.080807 0.004737 0.130658 0.071163 2.160284 0.331490 0.189419      0.728584 0.340918 0.428405        
#(apl: 2.021301, (nau: 0.442719, ((oct: 0.007813, bim: 0.017725): 0.415505, ((sep: 0.063779, squ: 0.080807): 0.000004, (bob: 0.130658, lin: 0.071163): 0.004737): 0.192185): 0.282412): 1.955118);
        
        
        
            
            
    
    
    
    
    
    