# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 18:05:35 2020

@author: shosh

This script uses pal2nal for each proteins super orthologs msa
It first iterates over all super orthologs proteins fasta's and re-construct the mrna sequnces used for protein in silico translation
then, for each super ortholog gene it cals pal2nal with the mrna sequnces fasta file and the proteins msa result

if log file for this run shows an error for certain files, use this cmd to fined files
for f in *fasta; do echo $f; grep "<an_id_from_problematic_file>" $f; done | grep -b1 "<an_id_from_problematic_file>"

"""

import os
import re
import glob
import subprocess
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
import argparse
import time

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Runnig pal2nal for each super ortholog protein msa')
    run_parser = parser.add_argument_group('Run pal2nal for a list of super orthologs proteins msa')
    run_parser.add_argument('-native_coding_mrnas', dest='native_coding_mrnas', action='store', required = True, help='path to native mrna sequences directory resides')
    run_parser.add_argument('-proteins_super_orthologs', dest='proteins_super_orthologs', action='store', required = True, help='path to proteins super orthologs dir in which msa_results direcotry with all msa files for all super orthologs genes')
    run_parser.add_argument('-o', dest='output_pal2nal_codon_msa', action='store', default = 'codons_msa_results', help='name of output subdirectory in super_orthologs_proteins_msa')
    run_parser.add_argument('-outfmt', dest='outfmt', action='store', default = 'clustal', help='output format of pal2nal')
    run_parser.add_argument('-animals', dest='animals', action='store', nargs = '+', default = ['oct','bim','squ','sep','bob','lin','nau'], help='subset of animals for finding super orthologs and perform MSA')
    arguments = parser.parse_args()
    
    native_coding_mrnas = arguments.native_coding_mrnas
    proteins_super_orthologs = arguments.proteins_super_orthologs
    animals = arguments.animals
    outfmt = arguments.outfmt
    output_pal2nal_codon_msa = proteins_super_orthologs+arguments.output_pal2nal_codon_msa+'/'
    


    #creating nested dict. each dict contains all coding native mrnas for animal
    print('Reading native coding mrna files')
    coding_mrnas_dict = {}
    for a in animals:
        native_coding_mrna_file = native_coding_mrnas+'native_coding_mrna_orfs_'+a+'.fa'
        animal_coding_mrnas_dict = {}
        for record in SeqIO.parse(open(native_coding_mrna_file, "r"), "fasta"):
            animal_coding_mrnas_dict.update({record.id:record.seq})
        coding_mrnas_dict.update({a:animal_coding_mrnas_dict})
        
    #creating native coding mrna sequences file for each protein super orthologs file
    print('Creating native coding mrna file for each proteins super orthologs in ' + proteins_super_orthologs)
    super_orthologs_protins_files = natural_sort(glob.glob(os.path.join(proteins_super_orthologs,'super_orthologs*.fasta')))
    for f in super_orthologs_protins_files:        
        mrnas_file = f.split('/')[-1].split('.')[0]+'_coding_mrnas.fasta'
        writer =  FastaWriter(open(proteins_super_orthologs+mrnas_file, 'w'), wrap=None)
        writer.write_header()
        for record in SeqIO.parse(open(f, "r"), "fasta"):
            animal = record.id.split('|')[0]
            native_sequence = coding_mrnas_dict[animal][record.id]
            writer.write_record(SeqRecord(native_sequence, id = record.id, description=''))
        writer.write_footer()
    print('Finished Writeing native coding mrna files')
    time.sleep(10)
        
        
    if not os.path.exists(output_pal2nal_codon_msa):
        os.makedirs(output_pal2nal_codon_msa)
    else:
        import warnings
        warnings.warn('Output dir already exists and existing files will be overwritten')    
    
    print('Creating codonds msa using pal2nal')
    pal2nal_processes = []
    for f in super_orthologs_protins_files:
        
        fn = [int(s) for s in f.replace('_',' ').replace('.',' ').split() if s.isdigit()][-1]
        gene_path = output_pal2nal_codon_msa+str(fn)+'/'
        if not os.path.exists(gene_path):
            os.makedirs(gene_path)
        else:
            import warnings
            warnings.warn(str(fn)+': dir already exists and existing files will be overwritten')   
        
        msa = '/'.join(f.split('/')[:-1])+'/msa_results/'+'msa_results_for_'+f.split('/')[-1].split('.')[0]+'.aln'
        native_coding_mrnas = '/'.join(f.split('/')[:-1])+'/'+f.split('/')[-1].split('.')[0]+'_coding_mrnas.fasta'
        pal2nal_cmd = '~/PAL2NAL/pal2nal.v14/pal2nal.pl ' + msa + ' ' + native_coding_mrnas + ' -output ' + outfmt + ' > ' + gene_path + 'codons_msa_for_'+f.split('/')[-1].split('.')[0]+'.aln'
        try:
            p = subprocess.Popen(pal2nal_cmd, shell = True, universal_newlines = True)
            pal2nal_processes.append(p)
            while p.poll() is None:
                time.sleep(1)
        except subprocess.CalledProcessError as e:
            print('file was not processed: ' + f)
            print(e.output.decode())
            
    print('Finished')
    
