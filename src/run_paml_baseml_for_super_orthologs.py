# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 14:28:13 2020

@author: shosh

This script creates tree file and control file for each super orthologs codons msa
Then executes PAML4 for ASR and other stats

if execution of PAML4 does not work,
you can run this cmd from msas_dirs_path after control file and tree were created for each msa:

for d in $(ls | grep -v txt); do cd $d; ctl=$(ls | grep ctl); echo "nohup ~/PAML4/paml4.8/bin/baseml "$ctl" > nohup_baseml_log" | sh; cd ..; done &

or maybe rewrite the script using PAML for python https://biopython.org/wiki/PAML
"""

import os
import re
import glob
import subprocess
from Bio import AlignIO
import argparse
import time

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def create_tree_for_super_orthologs(phylo_tree, ids, leaf_regex):
    
    tree = phylo_tree
    leafs_in_tree_strs = [m.group() for m in leaf_regex.finditer(tree)]
    for a in leafs_in_tree_strs:
        seq_id = [i for i in ids if i.split('|')[0]==a][0]
        tree = tree.replace(a,seq_id)
    return tree
        
  
"""
"tree with node labels for Rod Page's TreeView"
"(1_apl|comp99602_c0_seq1, (2_nau|comp200214_c5_seq3, ((3_oct|comp146819_c0_seq1, 4_bim|comp136326_c0_seq1) 12 , ((5_lin|comp768090_c0_seq4, 6_bob|comp377873_c8_seq1) 14 , (7_sep|comp181912_c0_seq1, 8_squ|comp135029_c1_seq2) 15 ) 13 ) 11 ) 10 ) 9 ;"
"""

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Runnig paml4 for each super ortholog protein msa')
    run_parser = parser.add_argument_group('Run PAML4 for a list of super orthologs proteins msa given a single tree for all')
    run_parser.add_argument('-msas_dirs_path', dest='msas_dirs_path', action='store', required = True, help='path to parent directory in which all codons-msas for all suprer orthologs, in seperate dirs')
    run_parser.add_argument('-tree', dest='tree_file', action='store', required = True, help='phylogenetic tree file including ancestors')
    run_parser.add_argument('-msa_fmt', dest='msa_fmt', action='store', default='clustal', help='msa format')
    arguments = parser.parse_args()
    
    msas_dirs_path = arguments.msas_dirs_path
    tree_file = arguments.tree_file
    msa_fmt = arguments.msa_fmt
    
    phylo_tree_w_ance = open(tree_file, "r").readlines()[0].rstrip()
    print('Phylogenetic tree with ancestors is: \n'+phylo_tree_w_ance)
    
    #remoning ancestors from phylo tree as PAML taked neweek tree with no ancestors. just the leafs. and ';' at the end.
    ances_regex_str = ';(?<=;)(.*?)(?=[,)]|$)'
    ances_regex = re.compile(ances_regex_str)
    ance_in_tree_srts = [m.group() for m in ances_regex.finditer(phylo_tree_w_ance)]
    phylo_tree_no_ances=phylo_tree_w_ance
    for m in ance_in_tree_srts:
        phylo_tree_no_ances = phylo_tree_no_ances.replace(m,'')
    phylo_tree_no_ances = phylo_tree_no_ances+';'
    print('\nPhylogenetic tree with no ancestors is (for PAML4): \n'+phylo_tree_no_ances)
    
    
    leaf_regex_str = '[a-z](?<=[a-z])([a-z])(?=[a-z])[a-z]'
    leaf_regex = re.compile(leaf_regex_str)
    paml4_processes = []
    dirs = natural_sort(glob.glob(os.path.join(msas_dirs_path,'*/')))
    print('\nCreating trees, control files, and running PAML4 for all MSAs')
    for d in dirs:
        dn = d.split('/')[-2]
        codons_msa_file_clustal = 'codons_msa_for_super_orthologs_'+str(dn)+'.aln'
        codons_msa_file_fasta = 'codons_msa_for_super_orthologs_'+str(dn)+'.fasta'
        ids=[]
        alignment = AlignIO.read(open(d+codons_msa_file_clustal), msa_fmt)
        for row in range(len(alignment)):
            ids.append(alignment[row].id)
        AlignIO.write(alignment, open(d+codons_msa_file_fasta, "w"), 'fasta')
        
        #writing tree for codons msa
        t_file = str(dn)+'.tree'
        t=open(d+t_file, "w")
        t.write(create_tree_for_super_orthologs(phylo_tree_no_ances,ids,leaf_regex))
        t.close()
        
        #writing baseml control file for codons msa
        ctl_file = 'baseml_'+str(dn)+'.ctl'
        ctl = open(d+ctl_file, "w")
        ctl.write('seqfile = '+codons_msa_file_fasta+'\n')
        ctl.write('treefile = '+t_file+'\n')
        ctl.write('outfile = mlb'+'\n')
        ctl.write('noisy = 2'+'\n')
        ctl.write('verbose = 0'+'\n')
        ctl.write('runmode = 0'+'\n')
        ctl.write('model = 7'+'\n')
        ctl.write('Mgene = 0'+'\n')
        ctl.write('*ndata = 100'+'\n')
        ctl.write('clock = 0'+'\n')
        ctl.write('fix_kappa = 0'+'\n')
        ctl.write('kappa = 5'+'\n')
        ctl.write('fix_alpha = 0'+'\n')
        ctl.write('alpha = 0.5'+'\n')
        ctl.write('Malpha = 0'+'\n')
        ctl.write('ncatG = 5'+'\n')
        ctl.write('nparK = 0'+'\n')
        ctl.write('nhomo = 0'+'\n')
        ctl.write('getSE = 0'+'\n')
        ctl.write('RateAncestor = 1'+'\n')
        ctl.write('Small_Diff = 7e-6'+'\n')
        ctl.write('cleandata = 0'+'\n')
        ctl.write('*icode = 0'+'\n')
        ctl.write('*fix_blength = -1'+'\n')
        ctl.write('method = 0'+'\n')
        ctl.close()
        
        #run PAML4
        try:
            paml_cmd = 'nohup ~/PAML4/paml4.8/bin/baseml ./' + ctl_file + ' > ' + d+'nohup_baseml_log'
            p = subprocess.Popen(paml_cmd, shell = True, universal_newlines = True, cwd=d)
            paml4_processes.append(p)
            while p.poll() is None:
                time.sleep(1)
        except subprocess.CalledProcessError as e:
            print('Could not Run PAML for super orthologs ' + str(dn))
            print(e.output.decode())
            
    print('Finished')
            
    
            