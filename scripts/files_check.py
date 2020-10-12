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
from Bio import AlignIO
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
#    run_parser.add_argument('-tree', dest='tree_file', action='store', required = True, help='phylogenetic tree file including ancestors')
    arguments = parser.parse_args()
    
    msas_dirs_path = arguments.msas_dirs_path
#    tree_file = arguments.tree_file

#    msas_dirs_path = 'C:/Users/shosh/OneDrive/Desktop/parent_test_dir/super_ortholog_proteins_fasta_files_all8/codons_msa_results/'
#    tree_file = 'C:/Users/shosh/OneDrive/Desktop/parent_test_dir/super_ortholog_proteins_fasta_files_all8/codons_msa_results/full_tree.txt'
    
    dirs = natural_sort(glob.glob(os.path.join(msas_dirs_path,'*/')))
    
    msa_fmt='clustal'
    branches=[]
    distances = []
    question_marks = [] #for some reaone codeml decide to replace few codons (usually in alignment gaps). would like to quantify this phenomena
    
    numbers = [int(d.split('/')[-2]) for d in dirs]
    
    print(str(len(dirs))+' from '+str(min(numbers))+' to '+str(max(numbers)))
    
    for d in dirs:
        
        try:
            dn=d.split('/')[-2]
            codons_msa_file_clustal = 'codons_msa_for_super_orthologs_'+str(dn)+'.aln'
            alignment = AlignIO.read(open(d+codons_msa_file_clustal), msa_fmt)
        except ValueError:
            print('problem reading '+codons_msa_file_clustal)
        