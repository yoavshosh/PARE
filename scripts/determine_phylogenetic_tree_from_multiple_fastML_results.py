# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 11:01:06 2019

@author: shosh

Just a small script that iterates through multiple fastML results and calculate prevalences of predicted newick trees

"""
import os
import glob
import subprocess
import argparse
from Bio import Phylo, AlignIO
from Bio.Phylo.Consensus import _BitString
import operator
import time

try:
    from cStringIO import StringIO ## for Python 2
except ImportError:
    from io import BytesIO ## for Python 3


fastml_perliminary_cmd = 'for i in {10..24}; do wget -nd -r -e robots=off -P ./S047_Pediatric_Brain_Cancer/"$i"CBTTC_PBT_Proteome_HMS_20180626/ -A raw "https://cptc-xfer.uis.georgetown.edu/publicData/External/S047_Pediatric_Brain_Cancer/""$i""CBTTC_PBT_Proteome_HMS_20180626/""$i""CBTTC_PBT_Proteome_HMS_20180626_raw/"; done &'

examp_tree1 = '(sep|comp176695_c0_seq1:0.251443,(bob|comp399435_c0_seq1:0.307890,lin|comp740065_c0_seq1:0.216070)N2:0.126805,(squ|comp140345_c0_seq1:0.250529,(nau|comp193952_c1_seq2:2.282868,(oct|comp140610_c0_seq1:0.034362,bim|comp146893_c0_seq2:0.014226)N5:1.053954)N4:0.391607)N3:0.113700)N1;'
examp_tree2 = '(lin|comp771355_c1_seq75:0.085036,(squ|comp138476_c0_seq1:0.026776,(nau|comp200952_c0_seq6:0.261511,(oct|comp169726_c0_seq3:0.003481,bim|comp148733_c0_seq2:0.005260)N4:0.141647)N3:0.132393)N2:0.019505,(bob|comp374464_c6_seq11:0.119102,sep|comp263529_c0_seq9:0.067529)N5:0.005050)N1;'

animal_id_deliminator = '|'

def getSubTree(newick_tree,animals=('sep','squ','bob','lin')):
    """
    recursivly serch tree clades and keep the subtree contaning animals
    """
    if type(newick_tree) is str:
        if newick_tree in animals:
            return newick_tree
    else:
        temp_subtree = tuple([x for x in [getSubTree(t,animals=animals) for t in newick_tree] if x])
        if len(temp_subtree)==1:
            return temp_subtree[0]
        else:
            return temp_subtree
            

def countUniqSubtrees(sorted_distinct_trees,animals=('sep','squ','bob','lin')):
    
    distinct_subtrees = {}
    for t,n in sorted_distinct_trees:
        t_str=t
        for a in animals_in_msa:
            t_str=t_str.replace(a,'"'+a+'"')
        subtree_str=str(getSubTree(eval(t_str),animals=animals))
        try:
            handle = StringIO(subtree_str) #for python2
        except NameError:
            handle = BytesIO(subtree_str) #for python3
        subtree_obj=Phylo.read(handle, "newick")
        subtree_obj.name = subtree_str
        
        tree_exist=False
        for k in distinct_subtrees.keys():
            if compare(subtree_obj,k):
                distinct_subtrees[k]+=n
                tree_exist = True
                break
        if not tree_exist:
            distinct_subtrees.update({subtree_obj:n})
            
    return distinct_subtrees



def run_perliminary_fastML_analysis(msa_path, fastML_proc, seq_type, files_prefix):
    
    """
    This function execute fastML for a list of msa files without any postulated tree
    this is for perliminary analysis in order to determine the most frequent phyogenetic tree among all msa files
    """
    
    fastML_logs_dir = msa_path+'perliminary_fastML_LOGS/'
    if not os.path.exists(fastML_logs_dir):
        os.makedirs(fastML_logs_dir)
        
    processes = []
    files_list = glob.glob(os.path.join(msa_path,files_prefix+'*.aln'))
    pennding_processes = True
    current_proc = 0    
    while pennding_processes:
        if sum([p.poll() is None for p in processes]) < fastML_proc:
            file_path = files_list[current_proc]
            file_name = file_path.split('/')[-1]
            out_dir_str = '/'.join(file_path.split('/')[:-1])+'/'+file_name.replace('.','_')+'_perliminary_fastML_analysis/'
            print('starting peliminary fastML analysis for ' + file_name)
            fastML_cmd = 'nohup perl ~/FastML/FastML.v3.11/www/fastml/FastML_Wrapper.pl --outDir '+out_dir_str + ' --MSA_File '+file_name + ' --seqType '+seq_type + ' > '+fastML_logs_dir+'peliminary_fastML_LOG_for_' + file_name + '.txt'
            processes.append(subprocess.Popen(fastML_cmd, shell = True, universal_newlines = True))
            current_proc+=1
        if len(files_list)==current_proc:
            pennding_processes = False
    
    #wati until all processes are finished
    while sum([p.poll() is None for p in processes]):
        time.sleep(1)
        


def simplify_newick_tree(newick_tree,animal_id_deliminator):
    """
    This funciton simplifies a full newick tree (with distances) and sequences names
    it returns the tree with no distances and just the animals names 
    animals names have to be a part of the sequence id and seperated by animal_id_deliminator from the rest of the seq id
    """
    
    tree_list = newick_tree.split(':')
    tree_list_no_dist = []
    for tree_fragment in tree_list:
        splited_tree_fragment = tree_fragment.split(',')
        for frag in splited_tree_fragment:
            if any([k in frag for k in [animal_id_deliminator,'(',')']]):
                tree_list_no_dist.append(frag)
                
    for i,tree_fragment in enumerate(tree_list_no_dist):
        if ')' in tree_fragment:
            tree_list_no_dist[i]=')'
        elif animal_id_deliminator in tree_fragment:
            tree_list_no_dist[i]=tree_fragment.split(animal_id_deliminator)[0]
            
    tree_no_dist_str = tree_list_no_dist[0]
    for i in tree_list_no_dist[1:]:
        if i == ')':
            tree_no_dist_str+=i
        else:
            tree_no_dist_str+= ','+i
    
    return tree_no_dist_str     

def sum_distinct_subtrees(simplified_distinct_trees_cnt,animals=('apl','oct','sep','squ')):
    
    simplified_trees_list_phylo_obj = []
    for t,n in enumerate(simplified_distinct_trees_cnt):
        try:
            handle = StringIO(t) #for python2
        except NameError:
            handle = BytesIO(t) #for python3
        tree_obj = Phylo.read(handle, "newick")
        tree_obj.name = t
        simplified_trees_list_phylo_obj.append(tree_obj)
    
    distinct_trees_dict = {simplified_trees_list_phylo_obj[0]:1}
    for t in simplified_trees_list_phylo_obj[1:]:
        tree_exist = False
        for k in distinct_trees_dict.keys():
            if compare(t,k):
                distinct_trees_dict[k]+=1
                tree_exist = True
                break
        if not tree_exist:
            distinct_trees_dict.update({t:1})
            
    return distinct_trees_dict
      

def collect_trees_structures(msa_path,tree_subpath_str,animal_id_deliminator):
    """
    iterate through a list of newick trees, simlify them and return a list of all simplified trees
    """
    trees_list = []
#    files_list = glob.glob(os.path.join(msa_path,'*_perliminary_fastML_analysis/tree.newick.txt'))
    files_list = glob.glob(os.path.join(msa_path,tree_subpath_str))
#    print(files_list)
    for file in files_list:
#        print('Reading ' + file)
        tree = open(file, "r").read().replace('\n','')
        trees_list.append(simplify_newick_tree(tree,animal_id_deliminator))
            
    return trees_list


def _bitstrs(tree):
    """
    from biopython tutorial
    https://biopython.org/wiki/Phylo
    """ 
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs

def compare(tree1, tree2):
    """
    from biopython tutorial
    https://biopython.org/wiki/Phylo
    """
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False      

def distinct_trees_count(simplified_trees_list):
    """
    This function recives a list of newick trees strings 
    and returns a dictionary where keys are Bio.Phylo objcets and values are number of times the tree apeared in the list
    """
    simplified_trees_list_phylo_obj = []
    for i,t in enumerate(simplified_trees_list):
        try:
            handle = StringIO(t) #for python2
        except NameError:
            handle = BytesIO(t) #for python3
        tree_obj = Phylo.read(handle, "newick")
        tree_obj.name = t
        simplified_trees_list_phylo_obj.append(tree_obj)
    
    distinct_trees_dict = {simplified_trees_list_phylo_obj[0]:1}
    for t in simplified_trees_list_phylo_obj[1:]:
        tree_exist = False
        for k in distinct_trees_dict.keys():
            if compare(t,k):
                distinct_trees_dict[k]+=1
                tree_exist = True
                break
        if not tree_exist:
            distinct_trees_dict.update({t:1})
            
    return distinct_trees_dict
  
      
def concatMSAfiles(animals,msas_path,files_prefix):
    joined_msas_dict={}
    for a in animals_in_msa:
        joined_msas_dict.update({a:''})
    for f in glob.glob(os.path.join(msa_path,files_prefix+'*.aln')):
        alignment = AlignIO.read(open(f), "clustal")
        for i in range(len(animals)):
            current_id = alignment[i].id
            current_animal = current_id.split('|')[0]
            current_seq = str(alignment[i].seq)
            joined_msas_dict[current_animal]+=current_seq 
    
    with open(msas_path+'concatinated_msas.fasta',"w") as handle:
        for k,v in joined_msas_dict.items():
            handle.write('>'+k+'\n')
            handle.write(v+'\n')
            


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Run fastML on a list of msa files without a postulated phylogenetic tree')
    run_parser = parser.add_argument_group('Run perliminary fastML')
    run_parser.add_argument('-msa_path', dest='msa_path', action='store', required = True, help='path to directory of all msa results files')
    run_parser.add_argument('-animal_id_deliminator', dest='animal_id_deliminator', action='store', default = '|', help='deliminator of animal name and rest of sequence name in sequences ids')
    run_parser.add_argument('-files_prefix', dest='files_prefix', action='store', default = 'codons_msa_for_super_orthologs_', help='')
    run_parser.add_argument('-seq_type', dest='seq_type', action='store', default = 'aa', help='types of sequences in msa files aa/nuc/codon')
    run_parser.add_argument('-run_fastML', dest='run_fastML', action='store', default = 'False', help='run fastML to reconstruct the trees - needed if trees were not already calculated')
    run_parser.add_argument('-fastML_proc', dest='fastML_proc', action='store', default = '30', help='nuber of fastML processes to run simultaniously')
    run_parser.add_argument('-tree_subpath_str', dest='tree_subpath_str', action='store', default = 'msa_results_for_super_orthologs_*/tree.newick.txt', help='sub path of trees within msa path')
    run_parser.add_argument('-animals_in_subtree', dest='animals_in_subtree', action='store', default = "('sep','squ','bob','lin')", help='sub path of trees within msa path')
    run_parser.add_argument('-only_subtrees_from_trees', dest='only_subtrees_from_trees', action='store', default = False)
    run_parser.add_argument('-uniq_trees_count', dest='uniq_trees_count', action='store', default = 'uniq_trees')
    run_parser.add_argument('-animals_in_msa', dest='animals_in_msa', action='store', default = ['apl','nau','oct', 'bim', 'sep', 'squ', 'lin', 'bob'])
    run_parser.add_argument('-concat_msas', dest='concat_msas', action='store', default = True, help='run a single process based on a concatination of all msas construct tree')
    
    
    arguments = parser.parse_args()
    msa_path = arguments.msa_path
    files_prefix = arguments.files_prefix
    animal_id_deliminator = arguments.animal_id_deliminator
    seq_type = arguments.seq_type
    run_fastML = eval(arguments.run_fastML)
    fastML_proc = int(arguments.fastML_proc)
    tree_subpath_str = arguments.tree_subpath_str
    animals_in_subtree = eval(arguments.animals_in_subtree)
    only_subtrees_from_trees = arguments.only_subtrees_from_trees
    uniq_trees_count_file = msa_path+arguments.uniq_trees_count
    animals_in_msa = arguments.animals_in_msa
    concat_msas = arguments.concat_msas
    
    if concat_msas:
        
        concatMSAfiles(animals_in_msa,msa_path,files_prefix)
    
        if run_fastML:
            fastML_logs_dir = msa_path+'perliminary_fastML_LOGS/'
            if not os.path.exists(fastML_logs_dir):
                os.makedirs(fastML_logs_dir)
            
            out_dir_str = msa_path+'concatinated_msa_results/'
            fastML_cmd = 'nohup perl ~/FastML/FastML.v3.11/www/fastml/FastML_Wrapper.pl --outDir '+out_dir_str + ' --MSA_File '+msa_path+'concatinated_msas.fasta' + ' --seqType '+seq_type + ' > '+fastML_logs_dir+'fastML_LOG_for_concatinated_msas.txt'
            subprocess.Popen(fastML_cmd, shell = True, universal_newlines = True)
        
    
    else:
        if run_fastML:
            run_perliminary_fastML_analysis(msa_path, fastML_proc, seq_type, files_prefix)
            only_subtrees_from_trees=False
            
            
        if not only_subtrees_from_trees:
            simplified_trees_list = collect_trees_structures(msa_path,tree_subpath_str,animal_id_deliminator)
            distinct_trees_dict = distinct_trees_count(simplified_trees_list)
            sorted_distinct_trees = sorted(distinct_trees_dict.items(), key=operator.itemgetter(1), reverse = True)
            total_distinct_trees = 0
            total_trees = 0
            with open(uniq_trees_count_file,"w") as handle:
                for t, n in sorted_distinct_trees:
                    handle.write(t.name + ' : ' + str(n)+'\n')
                    total_distinct_trees += 1
                    total_trees += n
            #        Phylo.draw(t)
            handle.close()
            print('Total trees: ' + str(total_trees))
            print('Total distinct trees: ' + str(total_distinct_trees))
            print('Subtrees of '+str(animals_in_subtree))
        
    
        sorted_distinct_trees=[(k[0],int(k[1])) for k in [l.split(' : ') for l in open(uniq_trees_count_file,"r").readlines()]]
        uniq_subtrees=countUniqSubtrees(sorted_distinct_trees,animals=animals_in_subtree)
        sorted_uniq_subtrees = sorted(uniq_subtrees.items(), key=operator.itemgetter(1), reverse = True)
        total_distinct_subtrees = 0
        total_subtrees = 0
        with open(msa_path+'uniq_subtrees',"w") as handle:
            for t, n in sorted_uniq_subtrees:
                handle.write(t.name + ' : ' + str(n)+'\n')
                total_distinct_subtrees += 1
                total_subtrees += n
        #        Phylo.draw(t)
        handle.close()    
    
        print('Total trees: ' + str(total_subtrees))
        print('Total distinct trees: ' + str(total_distinct_subtrees))