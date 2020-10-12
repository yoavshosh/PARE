# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 12:06:10 2020

@author: shosh
"""

import os
import sys
import pandas as pd
import numpy as np
from functools import reduce
from scipy import stats
from Bio import Phylo


try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
    


class Hypothesis:
    
    filter_types_for_suffix_dict = {'id':[str,tuple],
                            'length':int,
                            'pos_base_0':[int,tuple,list],
                            'animals':[int,tuple,list],
                            'protein':[str,tuple],
                            'edited':[int,tuple,list],
                            'coding_location':[int,tuple,list],
                            'nuc':[str,tuple],
                            'site_key':[str,tuple],
                            'editing_level':list,
                            'distance_to_nearest_recoding_site':list,
                            'strongest_in_gene':int,
                            'range_100_cluster':[str,tuple],
                            'strongest_in_range_100_cluster':int,
                            'original_aa':[str,tuple],
                            'recoding':int,
                            'target_aa':[str,tuple],
                            'nuc_prob':list,
                            'nucl_range':list}
    
    
    def __init__(self, nucl_mat, filters_dict={},filter_now=True,tree_str=None):
        self.tree_str=tree_str
        self.tree = None
        self.leaves_to_ancestors_dict = {}
        self.ancestors_to_leaves_dict = {}
        self.ancestors_to_downstream_ancestors_dict = {}
        self.ancestors_to_upstream_ancestors_dict = {}
        self.filters_dict=filters_dict
        self.original_filters_dict = filters_dict
        self.nucl_mat=nucl_mat
        self.original_full_matrix=nucl_mat
        self.edited = None
        self.mutated_editing_sites = None
        self.edited_leaves = None
        self.editing_level_method=None
        self.ancestor = None
        self.intermediate = None
        self.intermediate_nucl = None
        self.leaf = None
        self.leaf_nucl = None
        self.rates = {}
        self.sites_in_intermediate_nodes={}
        self.strict_data = {}
        self.liberal_data = {}
        self.no_depletion_data = {}
        self.mutations_count = {}
        self.editing_levels_histogram = {}
        self.mutated_sites_editing_levels_histogram = {}
        self.editing_sites_by_types = {}
        self.editing_types_rates = {}
        self.editing_levels_distribution_by_type = {}
        self.sites_substitutions = {}
        self.editing_rates_by_ancestral_state = {}
        self.editing_levels_distribution_by_ancestral_state = {}
        self.editing_level_method = None
        self.General_model = self.create_General_model()
        self.HPM = self.create_Harm_permitting_model()
        
        if tree_str is not None:
            self.create_tree_relations_dicts()
        
        if filter_now:
            self.filter_nucl_matrix()
            
    def create_tree_relations_dicts(self):
    
        handle = StringIO(self.tree_str)
        tree = Phylo.read(handle, "newick")
        self.tree = tree
        root_name=tree.root.name
        
        self.leaves_to_ancestors_dict.clear()
        leaves=tree.get_terminals()
        for l in leaves:
            ancestors = (root_name,)
            for n in tree.get_path(l.name):
                if n.name != l.name:
                    ancestors+=(n.name,)
            self.leaves_to_ancestors_dict.update({l.name:ancestors})
        
        self.ancestors_to_leaves_dict.clear()
        ancestors = tree.get_nonterminals()
        for a in ancestors:
            self.ancestors_to_leaves_dict.update({a.name:[l.name for l in a.get_terminals()]})
            
        self.ancestors_to_downstream_ancestors_dict.clear()
        ancestors = tree.get_nonterminals()
        for a in ancestors:
            self.ancestors_to_downstream_ancestors_dict.update({a.name:[i.name for i in a.get_nonterminals() if i.name!=a.name]})
        
        self.ancestors_to_upstream_ancestors_dict.clear()
        ancestors = tree.get_nonterminals()
        for a in ancestors:
            if a.name==root_name:
                self.ancestors_to_upstream_ancestors_dict.update({a.name:[]})
            else:
                upstream_ancestors = (root_name,)
                for n in tree.get_path(a.name):
                    if n.name!=a.name:
                        upstream_ancestors+=(n.name,)
                self.ancestors_to_upstream_ancestors_dict.update({a.name:upstream_ancestors})
                
    def find_recent_common_ancestor(self, leavs):
        assert (self.tree is not None), "phylogenetic tree is not defined"
        return self.tree.common_ancestor(leavs)
        
        
    def reset_all(self):
        self.filters_dict= self.original_filters_dict
        self.nucl_mat=self.original_full_matrix
        self.edited = None
        self.mutated_editing_sites = None
        self.edited_leaves = None
        self.editing_level_method=None
        self.ancestor = None
        self.intermediate = None
        self.intermediate_nucl = None
        self.leaf = None
        self.leaf_nucl = None
        self.rates = {}
        self.sites_in_intermediate_nodes={}
        self.strict_data = {}
        self.liberal_data = {}
        self.no_depletion_data = {}
        self.mutations_count = {}
        self.editing_levels_histogram = {}
        self.mutated_sites_editing_levels_histogram = {}
        self.editing_sites_by_types = {}
        self.editing_types_rates = {}
        self.editing_levels_distribution_by_type = {}
        self.sites_substitutions = {}
        self.editing_rates_by_ancestral_state = {}
        self.editing_levels_distribution_by_ancestral_state = {}
        self.editing_level_method = None
        self.General_model = self.create_General_model()
        self.HPM = self.create_Harm_permitting_model()
      
        
    def calc_editing_levels_histogram(self, matrix, bin_size=0.05):
        sites=[]
        lower=0.0
        for upper in np.linspace(bin_size,1,1/bin_size):
            sites.append((lower,upper,len(matrix[np.logical_and(matrix['combined_editing_level']>lower, matrix['combined_editing_level']<=upper)])))
            lower =  upper
        return pd.DataFrame(data = sites, columns=('lower','upper','sites'))

    
    
    def remove_filters(self):
        self.nucl_mat=self.original_full_matrix
      
        
    def filter_nucl_matrix(self, filters_dict=None):
        """
        This function apply each filter in filters_dict to slice self.nucl_mat
        """
        if filters_dict is None:
            pass
        else:
            self.filters_dict = filters_dict
            
        self.remove_filters()
        filtered_nucl_mat = self.nucl_mat
        for k,v in self.filters_dict.items():
            
            key_suffix = ''
            for d in self.filter_types_for_suffix_dict.keys():
                if k.endswith(d):
                    key_suffix=d
                    break
            assert (key_suffix!=''), "key suffix was not found in filter_types_for_suffix_dict"
            
            if type(self.filter_types_for_suffix_dict[key_suffix])!=type:
                assert (type(v) in self.filter_types_for_suffix_dict[key_suffix]),"Ilegal filter type for "+k+" got "+type(v).__name__+", expects "+str(self.filter_types_for_suffix_dict[key_suffix])
            else:
                assert (type(v)==self.filter_types_for_suffix_dict[key_suffix]),"Ilegal filter type for "+k+" got "+type(v).__name__+", expects "+self.filter_types_for_suffix_dict[key_suffix].__name__
            
            if type(v)==str or type(v)==int:
                filtered_nucl_mat = filtered_nucl_mat[filtered_nucl_mat[k]==v]
            elif type(v)==tuple:
                filtered_nucl_mat = filtered_nucl_mat[filtered_nucl_mat[k].isin(v)]
            elif type(v)==list:
                filtered_nucl_mat = filtered_nucl_mat[filtered_nucl_mat[k]>v[0]]
                filtered_nucl_mat = filtered_nucl_mat[filtered_nucl_mat[k]<=v[-1]]
        
        self.nucl_mat=filtered_nucl_mat
        #filteting editeding sites
        if self.edited is not None and self.edited_leaves is not None:
            self.collect_editing_sites(self.edited_leaves)
        #filtering mutated editing sites
        if self.mutated_editing_sites is not None:
            self.collect_mutated_editing_sites(self.leaf,self.leaf_nucl)
            
        #recalculating mutation rates
        for k in self.rates.keys():
            internal_nod=k.split('_')[0]
            end_nod=k.split('_')[-1]
            for mm in self.rates[k].keys():
                self.calc_mutation_rate(internal_nod,mm[0],end_nod,mm[1])
        
                
    def define_nodes(self,ancestral,intermediate,leaf,intermediate_nucl='A',leaf_nucl='G'):
        self.ancestor = ancestral
        self.intermediate = intermediate
        self.leaf = leaf
        self.calc_mutation_rate(intermediate,intermediate_nucl,leaf,leaf_nucl)
                
        
    def calc_mutation_rate(self,ancestral_nod,ancestral_nod_nucl,end_nod,mutated_end_nod_nucl):
        """
        this function updates self.rates with a new claculation acording to branch and desired mutation along brnach
        branch is determined be ancestral_nod to end_nod
        desired mutation is ancestral_nod_nucl to mutated_end_nod_nucl
        """
        exist=False
        branch=ancestral_nod+'_to_'+end_nod
        mm=ancestral_nod_nucl+mutated_end_nod_nucl
        if branch in self.rates.keys():
            if mm in self.rates[branch].keys():
                exist=True
                
        if not exist:
            mat = self.nucl_mat
            mat = mat[mat[ancestral_nod+'_nuc']==ancestral_nod_nucl]
            total_syn = mat[mat[ancestral_nod+'_'+ancestral_nod_nucl+mutated_end_nod_nucl+'_recoding']==0]          
            total_nonsyn = mat[mat[ancestral_nod+'_'+ancestral_nod_nucl+mutated_end_nod_nucl+'_recoding']==1]
            mutated_syn = total_syn[total_syn[end_nod+'_nuc']==mutated_end_nod_nucl]
            mutated_non_syn = total_nonsyn[total_nonsyn[end_nod+'_nuc']==mutated_end_nod_nucl]         
            if branch in self.rates:
                if mm in self.rates[branch]:
                    self.rates[branch][mm]=(('syn', float(len(mutated_syn)), float(len(mutated_syn))/float(len(total_syn))), ('nonsyn', float(len(mutated_non_syn)), float(len(mutated_non_syn))/float(len(total_nonsyn))))
                else:
                    self.rates[branch].update({mm:(('syn', float(len(mutated_syn)), float(len(mutated_syn))/float(len(total_syn))), ('nonsyn', float(len(mutated_non_syn)), float(len(mutated_non_syn))/float(len(total_nonsyn))))})
            else:
                self.rates.update({branch:{mm:(('syn', float(len(mutated_syn)), float(len(mutated_syn))/float(len(total_syn))), ('nonsyn', float(len(mutated_non_syn)), float(len(mutated_non_syn))/float(len(total_nonsyn))))}})
    
    
    def get_groups_of_edited_and_unedited_sites(self):
        
        groups_for_product=[]
        other_intermediate_leaves=[l for l in self.ancestors_to_leaves_dict[self.intermediate] if l!=self.leaf and all([k not in set.intersection(set(self.ancestors_to_downstream_ancestors_dict[self.intermediate]),set(self.leaves_to_ancestors_dict[self.leaf])) for k in self.leaves_to_ancestors_dict[l]])]
        downstream_to_intermediate_ancestors=self.ancestors_to_downstream_ancestors_dict[self.intermediate]
        upstream_to_intermediate_ancestor=[ance for ance in self.ancestors_to_upstream_ancestors_dict[self.intermediate] if ance not in self.ancestors_to_upstream_ancestors_dict[self.ancestor]]
        for ance in downstream_to_intermediate_ancestors:
            if ance not in self.leaves_to_ancestors_dict[self.leaf]:
                groups_for_product.append(self.ancestors_to_leaves_dict[ance])
        outside_intermediate_leafs=[]
        for ance in upstream_to_intermediate_ancestor:
            outside_intermediate_leafs+=[l for l in self.ancestors_to_leaves_dict[ance] if l not in self.ancestors_to_leaves_dict[self.intermediate]]
        outside_intermediate_leafs=list(set(outside_intermediate_leafs))
        if len(outside_intermediate_leafs):
            groups_for_product.append(outside_intermediate_leafs)
        leaves_for_product=[l for group in groups_for_product for l in group]
        edited_leaves=[]
        for l in other_intermediate_leaves:
            for m in leaves_for_product:
                if l!=m and not len(set.intersection(set(self.leaves_to_ancestors_dict[l]),set(self.leaves_to_ancestors_dict[m]),set(self.ancestors_to_downstream_ancestors_dict[self.intermediate]))):
                    edited_leaves.append(tuple(sorted([l,m])))
        edited_leaves = list(set(edited_leaves))
        self.edited_leaves = [x for x in edited_leaves if x]
        
        non_edited_leaves = []
        all_ancestors_leaves=self.ancestors_to_leaves_dict[self.ancestor]
        all_root_leaves=self.ancestors_to_leaves_dict[self.tree.root.name]
        for l in all_root_leaves:
            if l not in all_ancestors_leaves:
                non_edited_leaves.append(l)
        non_edited_leaves = non_edited_leaves
        self.non_edited_leaves = non_edited_leaves
    
    
    def collect_editing_sites(self, edited_leaves = None, non_edited_leaves=None, editing_level_method='average', bin_size=0.05):
        """
        edited_leaves - a list of tuples. each tuple is a group of animals in which all should be edited
        returns a data frame in which each row comply to at least one of the tuples in edited_leaves
        """
        def in_unification_of_intersections(row, edited_leaves,non_edited_leaves):
            in_unified_intersections = 0
            for group in edited_leaves:
                if len(group)==sum([row[a+'_edited'] for a in group]):
                    in_unified_intersections=1
            for leaf in non_edited_leaves:
                if row[leaf+'_edited']==1:
                    in_unified_intersections=0
                    break
            return in_unified_intersections
        
        if edited_leaves is None:
            self.get_groups_of_edited_and_unedited_sites()
            edited_leaves=self.edited_leaves
            if non_edited_leaves is None:
                non_edited_leaves=self.non_edited_leaves
            else:
                self.non_edited_leaves=non_edited_leaves
        elif non_edited_leaves is None:
            self.get_groups_of_edited_and_unedited_sites()
            non_edited_leaves=self.non_edited_leaves
            self.edited_leaves = edited_leaves
        else:
            self.edited_leaves = edited_leaves
            self.non_edited_leaves = non_edited_leaves
        
        mat = self.nucl_mat.copy()
        mat['in_unified_intersections'] = mat.apply(lambda row: in_unification_of_intersections(row, edited_leaves, non_edited_leaves), axis=1)
        mat = mat[mat['in_unified_intersections']==1]  
        self.editing_level_method=editing_level_method
        self.edited = mat.drop(labels='in_unified_intersections',axis=1)
        self.editing_level_method = editing_level_method
        self.calc_combined_editing_levels(editing_level_method=editing_level_method)
                
    def collect_mutated_editing_sites(self, leaf, leaf_nucl, bin_size=0.05):
        """
        collect editing sites from self.edited for which given leaf was mutated to leaf_nucl
        """
        self.leaf_nucl = leaf_nucl
        self.leaf = leaf
        assert (self.edited is not None), "You should define edited sites slice from nucl mat first"
        mat = self.edited
        mat = mat[np.logical_and(mat[leaf+'_nuc']==leaf_nucl, mat[leaf+'_edited']==0)]
        self.mutated_editing_sites = mat
        
    
    def calc_combined_editing_levels(self,table='edited',editing_level_method='average'):
        """
        prints number of sites in table for each editing level range
        combined editing level of edited_animals detemined by editing_level (averag/minimal/maximal)
        """
        def calc_combined_editing_levels_for_method(row,editing_level_method,edited_animals):
            assert (editing_level_method in ['average','minimal','maximal']),str(editing_level_method)+ "editing level type is not defined, please pass average/minimal/maximal"
            if editing_level_method=='average':
                return float(sum([row[a+'_editing_level'] for a in edited_animals]))/float(sum([row[a+'_edited'] for a in edited_animals]))
            elif editing_level_method=='minimal':
                return min([row[a+'_editing_level'] for a in edited_animals if row[a+'_edited']])
            elif editing_level_method=='maximal':
                return max([row[a+'_editing_level'] for a in edited_animals])

        editing_level_method = self.editing_level_method
        edited_animals = list(set([a for group in self.edited_leaves for a in group]))
        if table=='edited':
            mat = self.edited
            mat['combined_editing_level']=mat.apply(lambda row: calc_combined_editing_levels_for_method(row,editing_level_method,edited_animals), axis=1)
            self.edited=mat.copy()
        elif table=='mutated':
            mat = self.mutated_editing_sites
            mat['combined_editing_level']=mat.apply(lambda row: calc_combined_editing_levels_for_method(row,editing_level_method,edited_animals), axis=1)
            self.mutated_editing_sites=mat.copy()
            
    
    def write_data(self, path, file_name='hypotheses_analysis', data_to_write = ['strict','liberal','no_depletion'], file_type='xlsx', sub_name='sheet1'):
        
        def write(file_type, data_frame, path, file_name, sub_name='sheet1'):
            
            if file_type=='xlsx':
                if os.path.isfile(path+file_name):
                    from openpyxl import load_workbook
                    book = load_workbook(path+file_name)
                    writer = pd.ExcelWriter(path+file_name, engine = 'openpyxl')
                    writer.book = book
                else:
                    writer = pd.ExcelWriter(path+file_name+'.xlsx', engine='xlsxwriter')
                data_frame.to_excel(writer, sheet_name = sub_name)
                writer.save()
                writer.close()
            elif file_type=='csv':
                name=file_name+'_'+sub_name
                while name[-1]=='_':
                    name=name[:-1]
                data_frame.to_csv(path+file_name+'_'+sub_name, sep='\t')
             
        for data in data_to_write:
            if data=='strict':
                data_frame=pd.concat([s for s in self.strict_data.values()], axis=1, sort=False)    
                write(file_type, data_frame, path, file_name, sub_name='strict')
            elif data=='liberal':
                data_frame=pd.concat([s for s in self.liberal_data.values()], axis=1, sort=False)
                write(file_type, data_frame, path, file_name, sub_name='liberal') 
            elif data=='no_depletion':
                data_frame=pd.concat([s for s in self.no_depletion_data.values()], axis=1, sort=False)
                write(file_type, data_frame, path, file_name, sub_name='no_depletion') 
            elif data=='editing_types_count':
                data_frame=pd.concat([s for s in self.editing_sites_by_types.values()], axis=1, sort=False).transpose()
                write(file_type, data_frame, path, file_name, sub_name=sub_name)
            elif data=='editing_types_rates':
                data_frame=pd.concat([s for s in self.editing_types_rates.values()], axis=1, sort=False).transpose()
                write(file_type, data_frame, path, file_name, sub_name=sub_name)
            elif data=='mutations_count':
                data_frame=pd.concat([s for s in self.mutations_count.values()], axis=1, sort=False).transpose()
                write(file_type, data_frame, path, file_name, sub_name=sub_name)
            elif data=='editing_ancestral_rates':
                data_frame=pd.concat([s for s in self.editing_rates_by_ancestral_state.values()], axis=1, sort=False).transpose()
                write(file_type, data_frame, path, file_name, sub_name=sub_name)
            
            elif data=='sites_substitutions':
                data_frame=pd.concat([s for s in self.sites_substitutions.values()], axis=1, sort=False).transpose()
                write(file_type, data_frame, path, file_name, sub_name=sub_name)
            elif data=='editing_levels_by_types':
                if file_type=='xlsx':
                    data_frame=pd.concat([s for s in self.editing_levels_distribution_by_type.values()], axis=1, sort=False).transpose()
                    write(file_type, data_frame, path, file_name, sub_name=sub_name)
                elif file_type=='csv':
                    for k,s in self.editing_levels_distribution_by_type.items():
                        s.to_csv(path+file_name+'_'+k, sep='\t')
            elif data=='editing_levels_by_ancestral_state':
                if file_type=='xlsx':
                    data_frame=pd.concat([s for s in self.editing_levels_distribution_by_ancestral_state.values()], axis=1, sort=False).transpose()
                    write(file_type, data_frame, path, file_name, sub_name=sub_name)
                elif file_type=='csv':
                    for k,s in self.editing_levels_distribution_by_ancestral_state.items():
                        s.to_csv(path+file_name+'_'+k, sep='\t')
            
            
    def collect_sites_rates_in_intermediate_nodes(self, ancestors,intermediate,editing_level_method='average',editing_levels=[0,0.05],original_nucl='A',target_nucl='G'):
        """
        This function collect all edited sites at intermediate node (generated from ancestor to that that node)
        and devide to syn/nonsyn sites
        """
        mat = self.nucl_mat
        syn_a = len(mat[mat[intermediate+'_'+original_nucl+target_nucl+'_recoding']==0])
        non_syn_a = len(mat[mat[intermediate+'_'+original_nucl+target_nucl+'_recoding']==1])
        intermediate_decendants = self.ancestors_to_leaves_dict[intermediate]
        ancestors_decendants = self.ancestors_to_leaves_dict[ancestors]
        outside_intermediate_decendants = [l for l in ancestors_decendants if l not in intermediate_decendants]
        edited_leaves = list(set.intersection(set(outside_intermediate_decendants),set(intermediate_decendants)))
        self.collect_editing_sites(edited_leaves=edited_leaves,non_edited_leaves=[],editing_level_method=editing_level_method)
        edited=self.edited.copy()
        edited=edited[edited[intermediate+'_nuc']==original_nucl]
        edited=edited[np.logical_and(edited['combined_editing_level']>editing_levels[0],edited['combined_editing_level']<=editing_levels[1])]
        syn_s=len(edited[edited[intermediate+'_'+original_nucl+target_nucl+'_recoding']==0])
        non_syn_s=len(edited[edited[intermediate+'_'+original_nucl+target_nucl+'_recoding']==1])
        self.sites_in_intermediate_nodes.update({intermediate,(float(syn_s)/syn_a,float(non_syn_s)/non_syn_a)})
    
    
    def create_Harm_permitting_model(self):
        return Hypothesis.HPM(self)
    class HPM:
        def __init__(self, hypothesis):
            self.hpm = hypothesis
            
            
        def determine_editing_type(self, row, ancestors, leaf, leaf_original_nucl='A', leaf_target_nucl='G', exclude_ancestry=['N0']):
            """
            This function determine for each nonsyn mutation wrt to leaf_original_nucl->leaf_target_nucl
            whether it is restorative (mutated leaf target aa is present in at least one ancestor in ancerstors for the examined location) 
            or diversifying (mutated leaf target aa is not present at any of the ancerstors for the examined location)
            """
            
            ancestors_to_check=[a for a in ancestors if a not in exclude_ancestry]
            if row[leaf+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0:
                return 'syn'
            else:
                target_aa = row[leaf+'_'+leaf_original_nucl+leaf_target_nucl+'_target_aa']
                if target_aa in [row[a+'_original_aa'] for a in ancestors_to_check]:
                    return 'res'
                else:
                    return 'div'
        
    
        def determine_editing_type_for_multiple_leaves(self, row, leaves, leaf_original_nucl='A', leaf_target_nucl='G', exclude_ancestry=['N0'], ancestor_for_type=None):
            
            """
            for sites shared by multiple leaves, editing type is determined wet to the ancestry of the recent common ancestor of all leaves
            """
            def determine_type_from_ancestry(row,target_aa,original_aa,ancestry):
                if target_aa==original_aa:
                    return 'syn'
                else:
                    if target_aa in [row[a+'_original_aa'] for a in ancestry]:
                        return 'res'
                    elif len(ancestry)==0:
                        return 'nonsyn_unknown_type'
                    else:
                        return 'div'
            
            if ancestor_for_type is None:
                ancestor_for_type = self.hpm.tree.common_ancestor(leaves).name
            ancestry = self.hpm.ancestors_to_upstream_ancestors_dict[ancestor_for_type]
            ancestors_to_check = [a for a in ancestry if a not in exclude_ancestry]
            original_aa = row[ancestor_for_type+'_original_aa']
            target_aa = row[ancestor_for_type+'_'+leaf_original_nucl+leaf_target_nucl+'_target_aa']
            editing_type=determine_type_from_ancestry(row,target_aa,original_aa,ancestors_to_check)            
            return editing_type


        def find_substituations_sites(self, write_source_table=True, path=None, ancestor_for_type=None, examined_animals=None, leaves_groups=None, non_edited_leaves=[], editing_levels=[0,1], leaf_original_nucl='A', leaf_target_nucl='G', editing_level_method='average', count_multiple_subs_per_sites=False):
            
            def calc_p_val_fisher_exact(hits,pop,control_hits,control_pop):
                if pop!=0:
                    return stats.fisher_exact([[hits,pop],[control_hits,control_pop]])[1]
                else:
                    return '-'
                
            
            assert(type(leaves_groups)==list or type(leaves_groups)==tuple), "leaves_groups must be a nested list"
            assert(write_source_table and path is not None), "must provide a valid path if write_source_table=True"
            
            mat=self.hpm.nucl_mat

            if examined_animals is None:
                examined_animals = list(set([a for group in leaves_groups for a in group]))
            mat['editing_type'] = mat.apply(lambda row: self.determine_editing_type_for_multiple_leaves(row,examined_animals,ancestor_for_type=ancestor_for_type,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl), axis=1)
            
            self.hpm.collect_editing_sites(edited_leaves=leaves_groups,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method)
            editing_level_lower_bound = editing_levels[0]
            editing_level_upper_bound = editing_levels[1]
            
            edited=self.hpm.edited.copy()    
            edited = edited[np.logical_and(edited['combined_editing_level']>editing_level_lower_bound,edited['combined_editing_level']<=editing_level_upper_bound)]
            
            if not count_multiple_subs_per_sites:
                edited['has_a_genomic'+leaf_target_nucl] = edited.apply(lambda row: 1 if any([row[a+'_nuc']==leaf_target_nucl for a in examined_animals]) else 0, axis=1)
            else:
                edited['has_a_genomic'+leaf_target_nucl] = edited.apply(lambda row: sum([1 for a in examined_animals if row[a+'_nuc']==leaf_target_nucl]), axis=1)
            
            edited_subs = edited[edited['has_a_genomic'+leaf_target_nucl]>0]
                   
            name='conserved_in'
            for g in leaves_groups:
                name=name+'_'+'-'.join(g)
            
            if not count_multiple_subs_per_sites:
                name=name+'_one_'+leaf_target_nucl+'sub_per_sites_in_'+'-'.join(examined_animals)
            else:
                name=name+'_multiple_'+leaf_target_nucl+'subs_per_site_in_'+'-'.join(examined_animals)
            name=name+'_'+str(editing_levels[0])+'_'+str(editing_levels[1])
            
            if write_source_table:
                edited.to_csv(path+'edited_rows_'+name,sep='\t')
                                
            div = len(edited[edited['editing_type']=='div'])
            res = len(edited[edited['editing_type']=='res'])
            syn = len(edited[edited['editing_type']=='syn'])
            nonsyn_unknown_type = len(edited[edited['editing_type']=='nonsyn_unknown_type'])
            tot_nonsyn = nonsyn_unknown_type+res+div
            
            div_subs = edited_subs[edited_subs['editing_type']=='div']['has_a_genomic'+leaf_target_nucl].sum()
            res_subs = edited_subs[edited_subs['editing_type']=='res']['has_a_genomic'+leaf_target_nucl].sum()
            syn_subs = edited_subs[edited_subs['editing_type']=='syn']['has_a_genomic'+leaf_target_nucl].sum()        
            nonsyn_unknown_type_subs = edited_subs[edited_subs['editing_type']=='nonsyn_unknown_type']['has_a_genomic'+leaf_target_nucl].sum()
            tot_nonsyn_subs = nonsyn_unknown_type_subs+res_subs+div_subs
            
            div_p = calc_p_val_fisher_exact(div_subs,div,syn_subs,syn)
            res_p = calc_p_val_fisher_exact(res_subs,res,syn_subs,syn)
            nonsyn_unknown_type_p = calc_p_val_fisher_exact(nonsyn_unknown_type_subs,nonsyn_unknown_type,syn_subs,syn)
            tot_nonsyn_p = calc_p_val_fisher_exact(tot_nonsyn_subs,tot_nonsyn,syn_subs,syn)
            
            self.hpm.sites_substitutions.update({name:pd.Series(data=(ancestor_for_type,leaves_groups,non_edited_leaves,examined_animals,editing_level_lower_bound,editing_level_upper_bound,count_multiple_subs_per_sites,syn_subs,syn,res_subs,res,res_p,div_subs,div,div_p,nonsyn_unknown_type_subs,nonsyn_unknown_type,nonsyn_unknown_type_p,tot_nonsyn_subs,tot_nonsyn,tot_nonsyn_p),
                                                                index=('common_ancestor_for_sites_type','edited_groups','non_edited_leaves','subs_in_any_of','editing_level_lower_bound','editing_level_upper_bound','multiple_subs_per_sites','syn_subs','syn','res_subs','res','res_p','div_subs','div','div_p','nonsyn_unknown_type_subs','nonsyn_unknown_type','nonsyn_unknown_type_p','total_nonsyn_subs','total_nonsyn','tot_nonsyn_p'), name=name)})
            
            
            
        def collect_editing_levels_distributions_by_types(self, leaves, leaf_original_nucl='A', leaf_target_nucl='G', ancestors=None, editing_level_method='average'):
            
            if ancestors is None:
                if type(leaves) is str:
                    ancestors=self.hpm.leaves_to_ancestors_dict[leaves]
                if type(leaves) in [list,tuple]:
                    ancestors=list(reduce(set.intersection, [set(ances) for l,ances in self.hpm.leaves_to_ancestors_dict.items() if l in leaves]))
            
            if type(leaves) is str:
                name=leaves
                self.hpm.collect_editing_sites(edited_leaves=[(leaves,)],non_edited_leaves=[],editing_level_method=editing_level_method)
                mat=self.hpm.edited.copy()
                mat[leaves+'_editing_type']=mat.apply(lambda row: self.determine_editing_type(row,ancestors,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl),axis=1)
                self.hpm.editing_levels_distribution_by_type.update({name+'_syn':mat[mat[leaves+'_editing_type']=='syn'][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_res':mat[mat[leaves+'_editing_type']=='res'][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_div':mat[mat[leaves+'_editing_type']=='div'][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_specific_syn':mat[np.logical_and(mat[leaves+'_editing_type']=='syn',mat['edited_animals']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_specific_res':mat[np.logical_and(mat[leaves+'_editing_type']=='res',mat['edited_animals']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_specific_div':mat[np.logical_and(mat[leaves+'_editing_type']=='div',mat['edited_animals']==1)][leaves+'_editing_level']})
                
            elif type(leaves) in [list,tuple]: 
                name='_'.join(sorted(leaves))
                self.hpm.collect_editing_sites(edited_leaves=[leaves],non_edited_leaves=[],editing_level_method=editing_level_method)
                mat=self.hpm.edited.copy()
                mat['editing_type'] = mat.apply(lambda row: self.determine_editing_type_for_multiple_leaves(row,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl), axis=1)
                self.hpm.editing_levels_distribution_by_type.update({name+'_syn':mat[mat['editing_type']=='syn']['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_res':mat[mat['editing_type']=='res']['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_div':mat[mat['editing_type']=='div']['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_specific_syn':mat[np.logical_and(mat['editing_type']=='syn',mat['edited_animals']==len(leaves))]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_specific_res':mat[np.logical_and(mat['editing_type']=='res',mat['edited_animals']==len(leaves))]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_type.update({name+'_specific_div':mat[np.logical_and(mat['editing_type']=='div',mat['edited_animals']==len(leaves))]['combined_editing_level']})       
            
            
        def calc_editing_types_rates(self, leaves, leaf_original_nucl='A', leaf_target_nucl='G', ancestors=None, editing_levels=[0,1]):
            
            if ancestors is None:
                if type(leaves) is str:
                    ancestors=self.hpm.leaves_to_ancestors_dict[leaves]
                if type(leaves) in [list,tuple]:
                    ancestors=list(reduce(set.intersection, [set(ances) for l,ances in self.hpm.leaves_to_ancestors_dict.items() if l in leaves]))
                           
            mat=self.hpm.nucl_mat
            
            if type(leaves) is str:
            
                name=leaves+'_editing_level_'+str(editing_levels[0])+'-'+str(editing_levels[1])
                mat = mat[mat[leaves+'_nuc']==leaf_original_nucl].copy()
                mat[leaves+'_editing_type'] = mat.apply(lambda row: self.determine_editing_type(row,ancestors,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl),axis=1)
                edited=mat[mat[leaves+'_edited']==1].copy()
                edited = edited[np.logical_and(edited[leaves+'_editing_level']>editing_levels[0],edited[leaves+'_editing_level']<=editing_levels[1])]
                
                syn_edited = len(edited[edited[leaves+'_editing_type']=='syn'])
                specie_specific_syn_edited = len(edited[np.logical_and(edited[leaves+'_editing_type']=='syn', edited['edited_animals']==1)])
                syn = len(mat[mat[leaves+'_editing_type']=='syn'])
                res_edited = len(edited[edited[leaves+'_editing_type']=='res'])
                specie_specific_res_edited = len(edited[np.logical_and(edited[leaves+'_editing_type']=='res', edited['edited_animals']==1)])
                res = len(mat[mat[leaves+'_editing_type']=='res'])
                div_edited = len(edited[edited[leaves+'_editing_type']=='div'])
                specie_specific_div_edited = len(edited[np.logical_and(edited[leaves+'_editing_type']=='div', edited['edited_animals']==1)])
                div = len(mat[mat[leaves+'_editing_type']=='div'])
                
            if type(leaves) in [list,tuple]: 
                
                common_ances = self.hpm.tree.common_ancestor(leaves).name
                mat = mat[mat[common_ances+'_nuc']==leaf_original_nucl].copy()
                name='_'.join(sorted(leaves))+'_editing_level_'+str(editing_levels[0])+'-'+str(editing_levels[1])
                mat['editing_type'] = mat.apply(lambda row: self.determine_editing_type_for_multiple_leaves(row,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl), axis=1)
                edited=mat
                for l in leaves:
                    edited=edited[edited[l+'_edited']==1].copy()
                edited['combied_editing_level']=edited.apply(lambda row: sum([row[l+'_editing_level'] for l in leaves])/float(len(leaves)) , axis=1)
                edited = edited[np.logical_and(edited['combied_editing_level']>editing_levels[0],edited['combied_editing_level']<=editing_levels[1])]
                
                syn_edited = len(edited[edited['editing_type']=='syn'])
                specie_specific_syn_edited = len(edited[np.logical_and(edited['editing_type']=='syn', edited['edited_animals']==len(leaves))])
                syn = len(mat[mat['editing_type']=='syn'])
                res_edited = len(edited[edited['editing_type']=='res'])
                specie_specific_res_edited = len(edited[np.logical_and(edited['editing_type']=='res', edited['edited_animals']==len(leaves))])
                res = len(mat[mat['editing_type']=='res'])
                div_edited = len(edited[edited['editing_type']=='div'])
                specie_specific_div_edited = len(edited[np.logical_and(edited['editing_type']=='div', edited['edited_animals']==len(leaves))])
                div = len(mat[mat['editing_type']=='div'])
                
            self.hpm.editing_types_rates.update({name:pd.Series(data=(syn_edited,specie_specific_syn_edited,syn,res_edited,specie_specific_res_edited,res,div_edited,specie_specific_div_edited,div),
                                                                      index=('syn_edited','specie_specific_syn_edited','syn','res_edited','specie_specific_res_edited','res','div_edited','specie_specific_div_edited','div'), name=name)})
            
                 
        
        def count_editing_types(self, leaves, leaf_original_nucl='A', leaf_target_nucl='G', ancestors=None):
                
            if ancestors is None:
                if type(leaves) is str:
                    ancestors=self.hpm.leaves_to_ancestors_dict[leaves]
                if type(leaves) in [list,tuple]:
                    ancestors=list(reduce(set.intersection, [set(ances) for l,ances in self.hpm.leaves_to_ancestors_dict.items() if l in leaves]))
                             
            mat=self.hpm.nucl_mat
            if type(leaves) is str:
                name=leaves
                mat=mat[mat[leaves+'_edited']==1].copy()
                mat[leaves+'_editing_type']=mat.apply(lambda row: self.determine_editing_type(row,ancestors,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl),axis=1)
                syn_n = len(mat[mat[leaves+'_editing_type']=='syn'])
                res_n = len(mat[mat[leaves+'_editing_type']=='res'])
                div_n = len(mat[mat[leaves+'_editing_type']=='div'])
                specie_specific_syn_n = len(mat[np.logical_and(mat[leaves+'_editing_type']=='syn',mat['edited_animals']==1)])
                specie_specific_res_n = len(mat[np.logical_and(mat[leaves+'_editing_type']=='res',mat['edited_animals']==1)])
                specie_specific_div_n = len(mat[np.logical_and(mat[leaves+'_editing_type']=='div',mat['edited_animals']==1)])
                
            elif type(leaves) in [list,tuple]: 
                name='_'.join(sorted(leaves))
                for l in leaves:
                    mat=mat[mat[l+'_edited']==1].copy()
                mat['editing_type'] = mat.apply(lambda row: self.determine_editing_type_for_multiple_leaves(row,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl), axis=1)
                syn_n = len(mat[mat['editing_type']=='syn'])
                res_n = len(mat[mat['editing_type']=='res'])
                div_n = len(mat[mat['editing_type']=='div'])
                specie_specific_syn_n=len(mat[np.logical_and(mat['editing_type']=='syn',mat['edited_animals']==len(leaves))])
                specie_specific_res_n=len(mat[np.logical_and(mat['editing_type']=='res',mat['edited_animals']==len(leaves))])
                specie_specific_div_n=len(mat[np.logical_and(mat['editing_type']=='div',mat['edited_animals']==len(leaves))])

            
            self.hpm.editing_sites_by_types.update({name:pd.Series(data=(syn_n,specie_specific_syn_n,res_n,specie_specific_res_n,div_n,specie_specific_div_n),
                                                                       index=('synonymous','synonymous_specie_specific','restorative','restorative_specie_specific','diversifying','diversifying_specie_specific'), name=name)})
                
        
        
        
        def determine_ancestral_state(self, row, ancestors, leaf, leaf_original_nucl='A', leaf_target_nucl='G', exclude_ancestry=['N0']):
            
            #count ancestral C\T
            ancestors_to_check = [a for a in ancestors if a not in exclude_ancestry]
            if leaf_target_nucl in [row[a+'_nuc'] for a in ancestors_to_check]:
                return 'res'
            else:
                return 'div'
            
        def determine_ancestral_state_for_multiple_leaves(self, row, leaves, leaf_original_nucl='A', leaf_target_nucl='G', exclude_ancestry=['N0'], ancestor_for_type=None):
            
            """
            for sites shared by multiple leaves, editing type is determined wet to the ancestry of the recent common ancestor of all leaves
            """
            def determine_state_from_ancestry(row,ancestry,target_nucl):

                if target_nucl in [row[a+'_nuc'] for a in ancestry]:
                    return 'res'
                elif len(ancestry)==0:
                    return 'unknown_type'
                else:
                    return 'div'
        
            if ancestor_for_type is None:
                ancestor_for_type = self.hpm.tree.common_ancestor(leaves).name
            ancestry = self.hpm.ancestors_to_upstream_ancestors_dict[ancestor_for_type]
            ancestors_to_check = [a for a in ancestry if a not in exclude_ancestry]
            editing_type=determine_state_from_ancestry(row,ancestors_to_check,target_nucl=leaf_target_nucl)            
            return editing_type

        
        def calc_acestral_states_rates(self, leaves, leaf_original_nucl='A', leaf_target_nucl='G', ancestors=None, editing_levels=[0,1]):
            
            if ancestors is None:
                if type(leaves) is str:
                    ancestors=self.hpm.leaves_to_ancestors_dict[leaves]
                if type(leaves) in [list,tuple]:
                    ancestors=list(reduce(set.intersection, [set(ances) for l,ances in self.hpm.leaves_to_ancestors_dict.items() if l in leaves]))
                           
            mat=self.hpm.nucl_mat
            
            if type(leaves) is str:
            
                name=leaves+'_editing_level_'+str(editing_levels[0])+'-'+str(editing_levels[1])
                mat = mat[mat[leaves+'_nuc']==leaf_original_nucl].copy()
                mat[leaves+'_editing_type'] = mat.apply(lambda row: self.determine_ancestral_state(row,ancestors,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl),axis=1)
                edited=mat[mat[leaves+'_edited']==1].copy()
                edited = edited[np.logical_and(edited[leaves+'_editing_level']>editing_levels[0],edited[leaves+'_editing_level']<=editing_levels[1])]
                
                syn_edited = edited[edited[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0]
                specie_specific_syn_edited = edited[np.logical_and(edited[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0, edited['edited_animals']==1)]
                syn = mat[mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0]
                
                nonsyn_edited = edited[edited[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1]
                specie_specific_nonsyn_edited = edited[np.logical_and(edited[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1, edited['edited_animals']==1)]
                nonsyn = mat[mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1]
                
                syn_div = len(syn[syn[leaves+'_editing_type']=='div'])
                syn_div_edited = len(syn_edited[syn_edited[leaves+'_editing_type']=='div'])
                syn_div_species_specific_edited = len(specie_specific_syn_edited[specie_specific_syn_edited[leaves+'_editing_type']=='div'])
                syn_res = len(syn[syn[leaves+'_editing_type']=='res'])
                syn_res_edited = len(syn_edited[syn_edited[leaves+'_editing_type']=='res'])
                syn_res_species_specific_edited = len(specie_specific_syn_edited[specie_specific_syn_edited[leaves+'_editing_type']=='res'])
                
                nonsyn_div = len(nonsyn[nonsyn[leaves+'_editing_type']=='div'])
                nonsyn_div_edited = len(nonsyn_edited[nonsyn_edited[leaves+'_editing_type']=='div'])
                nonsyn_div_species_specific_edited = len(specie_specific_nonsyn_edited[specie_specific_nonsyn_edited[leaves+'_editing_type']=='div'])
                nonsyn_res = len(nonsyn[nonsyn[leaves+'_editing_type']=='res'])
                nonsyn_res_edited = len(nonsyn_edited[nonsyn_edited[leaves+'_editing_type']=='res'])
                nonsyn_res_species_specific_edited = len(specie_specific_nonsyn_edited[specie_specific_nonsyn_edited[leaves+'_editing_type']=='res'])
                

            if type(leaves) in [list,tuple]: 
                
                common_ances = self.hpm.tree.common_ancestor(leaves).name
                mat = mat[mat[common_ances+'_nuc']==leaf_original_nucl].copy()
                name='_'.join(sorted(leaves))+'_editing_level_'+str(editing_levels[0])+'-'+str(editing_levels[1])
                mat['editing_type'] = mat.apply(lambda row: self.determine_ancestral_state_for_multiple_leaves(row,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl), axis=1)
                edited=mat
                for l in leaves:
                    edited=edited[edited[l+'_edited']==1].copy()
                edited['combied_editing_level']=edited.apply(lambda row: sum([row[l+'_editing_level'] for l in leaves])/float(len(leaves)) , axis=1)
                edited = edited[np.logical_and(edited['combied_editing_level']>editing_levels[0],edited['combied_editing_level']<=editing_levels[1])]
                
                common_ancestor = self.hpm.tree.common_ancestor(leaves).name
            
                syn_edited = edited[edited[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0]
                specie_specific_syn_edited = edited[np.logical_and(edited[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0, edited['edited_animals']==len(leaves))]
                syn = mat[mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0]
                
                nonsyn_edited = edited[edited[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1]
                specie_specific_nonsyn_edited = edited[np.logical_and(edited[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1, edited['edited_animals']==len(leaves))]
                nonsyn = mat[mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1]
                
            
                syn_div = len(syn[syn['editing_type']=='div'])
                syn_div_edited = len(syn_edited[syn_edited['editing_type']=='div'])
                syn_div_species_specific_edited = len(specie_specific_syn_edited[specie_specific_syn_edited['editing_type']=='div'])
                syn_res = len(syn[syn['editing_type']=='res'])
                syn_res_edited = len(syn_edited[syn_edited['editing_type']=='res'])
                syn_res_species_specific_edited = len(specie_specific_syn_edited[specie_specific_syn_edited['editing_type']=='res'])
                
                nonsyn_div = len(nonsyn[nonsyn['editing_type']=='div'])
                nonsyn_div_edited = len(nonsyn_edited[nonsyn_edited['editing_type']=='div'])
                nonsyn_div_species_specific_edited = len(specie_specific_nonsyn_edited[specie_specific_nonsyn_edited['editing_type']=='div'])
                nonsyn_res = len(nonsyn[nonsyn['editing_type']=='res'])
                nonsyn_res_edited = len(nonsyn_edited[nonsyn_edited['editing_type']=='res'])
                nonsyn_res_species_specific_edited = len(specie_specific_nonsyn_edited[specie_specific_nonsyn_edited['editing_type']=='res'])
                
            data=(syn_div,syn_div_edited,syn_div_species_specific_edited,syn_res,syn_res_edited,syn_res_species_specific_edited,nonsyn_div,nonsyn_div_edited,nonsyn_div_species_specific_edited,nonsyn_res,nonsyn_res_edited,nonsyn_res_species_specific_edited)
            index=('syn_div','syn_div_edited','syn_div_species_specific_edited','syn_res','syn_res_edited','syn_res_species_specific_edited','nonsyn_div','nonsyn_div_edited','nonsyn_div_species_specific_edited','nonsyn_res','nonsyn_res_edited','nonsyn_res_species_specific_edited')
            self.hpm.editing_rates_by_ancestral_state.update({name:pd.Series(data=data,index=index,name=name)})
        
            
        def collect_editing_levels_distributions_by_ancestral_state(self, leaves, leaf_original_nucl='A', leaf_target_nucl='G', ancestors=None, editing_level_method='average'):
            
            if ancestors is None:
                if type(leaves) is str:
                    ancestors=self.hpm.leaves_to_ancestors_dict[leaves]
                if type(leaves) in [list,tuple]:
                    ancestors=list(reduce(set.intersection, [set(ances) for l,ances in self.hpm.leaves_to_ancestors_dict.items() if l in leaves]))
            
            if type(leaves) is str:
                name=leaves
                self.hpm.collect_editing_sites(edited_leaves=[(leaves,)],non_edited_leaves=[],editing_level_method=editing_level_method)
                mat=self.hpm.edited.copy()
                mat[leaves+'_editing_type']=mat.apply(lambda row: self.determine_ancestral_state(row,ancestors,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl),axis=1)
                
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_res':mat[np.logical_and(mat[leaves+'_editing_type']=='res',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_div':mat[np.logical_and(mat[leaves+'_editing_type']=='div',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_res':mat[np.logical_and(mat[leaves+'_editing_type']=='res',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_div':mat[np.logical_and(mat[leaves+'_editing_type']=='div',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_res_species_specific':mat[np.logical_and(np.logical_and(mat[leaves+'_editing_type']=='res',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0),mat['edited_animals']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_div_species_specific':mat[np.logical_and(np.logical_and(mat[leaves+'_editing_type']=='div',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0),mat['edited_animals']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_res_species_specific':mat[np.logical_and(np.logical_and(mat[leaves+'_editing_type']=='res',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1),mat['edited_animals']==1)][leaves+'_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_div_species_specific':mat[np.logical_and(np.logical_and(mat[leaves+'_editing_type']=='div',mat[leaves+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1),mat['edited_animals']==1)][leaves+'_editing_level']})
                
                
            elif type(leaves) in [list,tuple]: 
                name='_'.join(sorted(leaves))
                self.hpm.collect_editing_sites(edited_leaves=[leaves],non_edited_leaves=[],editing_level_method=editing_level_method)
                mat=self.hpm.edited.copy()
                mat['editing_type'] = mat.apply(lambda row: self.determine_ancestral_state_for_multiple_leaves(row,leaves,leaf_original_nucl=leaf_original_nucl,leaf_target_nucl=leaf_target_nucl), axis=1)
                
                common_ancestor = self.hpm.tree.common_ancestor(leaves).name
            
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_res':mat[np.logical_and(mat['editing_type']=='res',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0)]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_div':mat[np.logical_and(mat['editing_type']=='div',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0)]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_res':mat[np.logical_and(mat['editing_type']=='res',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1)]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_div':mat[np.logical_and(mat['editing_type']=='div',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1)]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_res_species_specific':mat[np.logical_and(np.logical_and(mat['editing_type']=='res',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0),mat['edited_animals']==len(leaves))]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_syn_div_species_specific':mat[np.logical_and(np.logical_and(mat['editing_type']=='div',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==0),mat['edited_animals']==len(leaves))]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_res_species_specific':mat[np.logical_and(np.logical_and(mat['editing_type']=='res',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1),mat['edited_animals']==len(leaves))]['combined_editing_level']})
                self.hpm.editing_levels_distribution_by_ancestral_state.update({name+'_nonsyn_div_species_specific':mat[np.logical_and(np.logical_and(mat['editing_type']=='div',mat[common_ancestor+'_'+leaf_original_nucl+leaf_target_nucl+'_recoding']==1),mat['edited_animals']==len(leaves))]['combined_editing_level']})
                
    
        
    def create_General_model(self):
        return Hypothesis.General_model(self)
    class General_model:
        def __init__(self, hypothoesis):
            self.general_model = hypothoesis
            
            
             
        def strict(self, ancestor, intermediate, leaf, leaf_nucl, edited_leaves=None, non_edited_leaves=None, editing_level_method='average', strong_levels=[0.1,1], confidence_level=0.95, n_random=1000000, sites_recalculation = True):
        
            def calc_strict_model_distribuation(syn_a,non_syn_a,syn_ag_mut,strong_syn_sites,strong_non_syn_sites,confidence_level,n):
    
                try:
                    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
                    p=float(strong_syn_sites)/syn_a
                    z = 1-(1-confidence_level)/2
                    
                    expected_non_syn_sites_dist = []
                    for i in range(n):
                        p_rand = np.random.normal(p,z*np.sqrt(p*(1-p)/syn_a))
                        while p_rand<0:
                            p_rand = np.random.normal(p,z*np.sqrt(p*(1-p)/syn_a))
                        lam = non_syn_a*p_rand
                        expected_non_syn_sites_dist.append(np.random.poisson(lam))
                    
                    expected_non_syn_eg_mut_dist = []
                    for i in range(n):
                        expected_eg_mut = (strong_non_syn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
                        expected_non_syn_eg_mut_dist.append(np.random.poisson(expected_eg_mut))    
                    
                    return expected_non_syn_eg_mut_dist
            
                except Exception as e:
                    return str(e)
            
            intermediate_nucl='A'
            self.general_model.define_nodes(ancestor,intermediate,leaf,intermediate_nucl=intermediate_nucl,leaf_nucl=leaf_nucl)
            if sites_recalculation or (self.general_model.edited is None or self.general_model.mutated_editing_sites is None):
                self.general_model.collect_editing_sites(edited_leaves=edited_leaves,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method)
                self.general_model.collect_mutated_editing_sites(leaf, leaf_nucl)
            
            intermediate_nucl_mat = self.general_model.nucl_mat[self.general_model.nucl_mat[intermediate+'_nuc']==intermediate_nucl]
            edited_mat = self.general_model.edited[self.general_model.edited[intermediate+'_nuc']==intermediate_nucl]
            mutations_mat = self.general_model.mutated_editing_sites[self.general_model.mutated_editing_sites[intermediate+'_nuc']==intermediate_nucl]
            
            syn_a = len(intermediate_nucl_mat[intermediate_nucl_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0])
            non_syn_a = len(intermediate_nucl_mat[intermediate_nucl_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1])
            syn_ag_mut = int(self.general_model.rates[intermediate+'_to_'+leaf][intermediate_nucl+leaf_nucl][0][1])
            non_syn_ag_mut = int(self.general_model.rates[intermediate+'_to_'+leaf][intermediate_nucl+leaf_nucl][1][1])
            strong_syn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0, np.logical_and(edited_mat['combined_editing_level']>strong_levels[0], edited_mat['combined_editing_level']<=strong_levels[1]))])
            strong_nonsyn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1, np.logical_and(edited_mat['combined_editing_level']>strong_levels[0], edited_mat['combined_editing_level']<=strong_levels[1]))])
            actual_strong_nonsyn_eg_mutations = len(mutations_mat[np.logical_and(mutations_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1,np.logical_and(mutations_mat['combined_editing_level']>strong_levels[0], mutations_mat['combined_editing_level']<=strong_levels[1]))])
            
            syn_sites_creation_rate = float(strong_syn_sites)/float(syn_a)
            expected_strong_nonsyn_sites = syn_sites_creation_rate*float(non_syn_a)
            strong_nonsyn_sites_excess = strong_nonsyn_sites-expected_strong_nonsyn_sites
            expected_nonsyn_EG_mutations = strong_nonsyn_sites_excess*self.general_model.rates[intermediate+'_to_'+leaf]['A'+leaf_nucl][0][1]/float(syn_a)
            excess_of_unmutated_sites = expected_nonsyn_EG_mutations-actual_strong_nonsyn_eg_mutations
                      
            strict_expected_non_syn_eg_mut_dist = calc_strict_model_distribuation(syn_a, non_syn_a, syn_ag_mut, strong_syn_sites, strong_nonsyn_sites, confidence_level, n_random)
            if type(strict_expected_non_syn_eg_mut_dist)==str:
                strict_std=strict_expected_non_syn_eg_mut_dist
                strict_p=strict_expected_non_syn_eg_mut_dist
            else:
                strict_std =np.std(strict_expected_non_syn_eg_mut_dist)
                strict_p = float(sum([1 for j in strict_expected_non_syn_eg_mut_dist if j<actual_strong_nonsyn_eg_mutations]))/float(n_random)
     
            self.general_model.strict_data.update({ancestor+'_to_'+intermediate+'_to_'+leaf+'_'+editing_level_method+'_'+str(strong_levels[0]):pd.Series(data = [ancestor,intermediate,leaf,editing_level_method,strong_levels[0],strong_levels[1],syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,float(syn_ag_mut)/float(syn_a),float(non_syn_ag_mut)/float(non_syn_a),self.general_model.edited_leaves,self.general_model.non_edited_leaves,strong_syn_sites,strong_nonsyn_sites,expected_strong_nonsyn_sites,strong_nonsyn_sites_excess,expected_nonsyn_EG_mutations,strict_std,actual_strong_nonsyn_eg_mutations,excess_of_unmutated_sites,strict_p],
                                                                                      index = ['ancestor','intermediate','leaf','combined_editing_level_method','strong_levels_lower_limit','strong_levels_upper_limit','intermediate_syn_a','intermediate_non_syn_a','syn_ag_mut_in_leaf','non_syn_ag_mut_in_leaf','syn_ag_mut_rate','non_syn_ag_mut_rate','groups_of_edited_leaves','non_edited_leaves','strong_syn_sites','strong_non_syn_sites','expected_strong_nonsyn_sites','strong_nonsyn_sites_excess','expected_nonusn_EG_mutations','std','actual_nonsyn_eg_mutations','excess_of_unmutated_sites','p'], name=leaf)})
        
        def liberal(self, ancestor, intermediate, leaf, leaf_nucl, edited_leaves=None, non_edited_leaves=None, editing_level_method='average', weak_levels=[0,0.02], strong_levels=[0.1,1], confidence_level=0.95, n_random=1000000, sites_recalculation=True):
            
            def calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,weak_syn_sites,weak_non_syn_sites,strong_syn_sites,strong_non_syn_sites,confidence_level,n):
                
                try:
                    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
                    non_syn_ag_mut_rate = float(non_syn_ag_mut)/non_syn_a
                    z = 1-(1-confidence_level)/2
                    p_weak_nss=float(weak_non_syn_sites)/non_syn_a
                    p_weak_ss=float(weak_syn_sites)/syn_a
                    p_strong_ss=float(strong_syn_sites)/syn_a
                    
                    expected_non_syn_sites_dist = []
                    for i in range(n):
                        p_weak_nss_rand = np.random.normal(p_weak_nss,z*np.sqrt(p_weak_nss*(1-p_weak_nss)/non_syn_a))
                        while p_weak_nss_rand<0:
                            p_weak_nss_rand = np.random.normal(p_weak_nss,z*np.sqrt(p_weak_nss*(1-p_weak_nss)/non_syn_a))
                        p_weak_ss_rand = np.random.normal(p_weak_ss,z*np.sqrt(p_weak_ss*(1-p_weak_ss)/syn_a))
                        while p_weak_ss_rand<0:
                            p_weak_ss_rand = np.random.normal(p_weak_ss,z*np.sqrt(p_weak_ss*(1-p_weak_ss)/syn_a))
                        p_strong_ss_rand = np.random.normal(p_strong_ss,z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a))
                        while p_strong_ss_rand<0:
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
                
                except Exception as e:
                    return str(e)

            
            intermediate_nucl='A'
            self.general_model.define_nodes(ancestor,intermediate,leaf,intermediate_nucl=intermediate_nucl,leaf_nucl=leaf_nucl)
            if sites_recalculation or (self.general_model.edited is None or self.general_model.mutated_editing_sites is None):
                self.general_model.collect_editing_sites(edited_leaves=edited_leaves,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method)
                self.general_model.collect_mutated_editing_sites(leaf, leaf_nucl)
                        
            intermediate_nucl_mat = self.general_model.nucl_mat[self.general_model.nucl_mat[intermediate+'_nuc']==intermediate_nucl]
            edited_mat = self.general_model.edited[self.general_model.edited[intermediate+'_nuc']==intermediate_nucl]
            mutations_mat = self.general_model.mutated_editing_sites[self.general_model.mutated_editing_sites[intermediate+'_nuc']==intermediate_nucl]
            
            syn_a = len(intermediate_nucl_mat[intermediate_nucl_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0])
            non_syn_a = len(intermediate_nucl_mat[intermediate_nucl_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1])
            syn_ag_mut = int(self.general_model.rates[intermediate+'_to_'+leaf][intermediate_nucl+leaf_nucl][0][1])
            non_syn_ag_mut = int(self.general_model.rates[intermediate+'_to_'+leaf][intermediate_nucl+leaf_nucl][1][1])
            weak_syn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0, np.logical_and(edited_mat['combined_editing_level']>weak_levels[0], edited_mat['combined_editing_level']<=weak_levels[1]))])
            weak_nonsyn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1, np.logical_and(edited_mat['combined_editing_level']>weak_levels[0], edited_mat['combined_editing_level']<=weak_levels[1]))])
            strong_syn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0, np.logical_and(edited_mat['combined_editing_level']>strong_levels[0], edited_mat['combined_editing_level']<=strong_levels[1]))])
            strong_nonsyn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1, np.logical_and(edited_mat['combined_editing_level']>strong_levels[0], edited_mat['combined_editing_level']<=strong_levels[1]))])
            actual_strong_nonsyn_eg_mutations = len(mutations_mat[np.logical_and(mutations_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1,np.logical_and(mutations_mat['combined_editing_level']>strong_levels[0], mutations_mat['combined_editing_level']<=strong_levels[1]))])
            
            syn_sites_creation_rate = float(strong_syn_sites)/float(syn_a)
            non_syn_over_syn_a = float(non_syn_a)/float(syn_a)
            non_syn_over_syn_weak_sites = float(weak_nonsyn_sites)/float(weak_syn_sites)
            nonsyn_sites_depletion_factor = non_syn_over_syn_weak_sites/non_syn_over_syn_a
            expected_strong_nonsyn_sites = nonsyn_sites_depletion_factor*syn_sites_creation_rate*float(non_syn_a)
            strong_nonsyn_sites_excess = strong_nonsyn_sites-expected_strong_nonsyn_sites
            expected_nonsyn_EG_mutations_from_excess = strong_nonsyn_sites_excess*self.general_model.rates[intermediate+'_to_'+leaf]['A'+leaf_nucl][0][1]/float(syn_a)
            expected_nonsyn_EG_mutations_from_expected = expected_strong_nonsyn_sites*self.general_model.rates[intermediate+'_to_'+leaf]['A'+leaf_nucl][1][1]/float(non_syn_a)
            total_expected_EG_mutations = expected_nonsyn_EG_mutations_from_excess+expected_nonsyn_EG_mutations_from_expected
            excess_of_unmutated_sites = total_expected_EG_mutations-actual_strong_nonsyn_eg_mutations
                      
            liberal_expected_non_syn_eg_mut_dist = calc_liberal_model_distribuation(syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,weak_syn_sites,weak_nonsyn_sites,strong_syn_sites,strong_nonsyn_sites,confidence_level,n_random)
            if type(liberal_expected_non_syn_eg_mut_dist)==str:
                liberal_std=liberal_expected_non_syn_eg_mut_dist
                liberal_p=liberal_expected_non_syn_eg_mut_dist
            else:
                liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
                liberal_p = float(sum([1 for j in liberal_expected_non_syn_eg_mut_dist if j<actual_strong_nonsyn_eg_mutations]))/float(n_random)
            
            self.general_model.liberal_data.update({ancestor+'_to_'+intermediate+'_to_'+leaf+'_'+editing_level_method+'_'+str(weak_levels[1])+'_'+str(strong_levels[0]):pd.Series(data = [ancestor,intermediate,leaf,editing_level_method,strong_levels[0],strong_levels[1],weak_levels[0],weak_levels[1],syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,float(syn_ag_mut)/float(syn_a),float(non_syn_ag_mut)/float(non_syn_a),self.general_model.edited_leaves,self.general_model.non_edited_leaves,weak_syn_sites,weak_nonsyn_sites,nonsyn_sites_depletion_factor,strong_syn_sites,strong_nonsyn_sites,expected_strong_nonsyn_sites,strong_nonsyn_sites_excess,expected_nonsyn_EG_mutations_from_excess,expected_nonsyn_EG_mutations_from_expected,total_expected_EG_mutations,liberal_std,actual_strong_nonsyn_eg_mutations,excess_of_unmutated_sites,liberal_p],
                                                                                       index = ['ancestor','intermediate','leaf','combined_editing_level_method','strong_levels_lower_limit','strong_levels_upper_limit','weak_levels_lower_limit','weak_levels_upper_limit','intermediate_syn_a','intermediate_non_syn_a','syn_ag_mut_in_leaf','non_syn_ag_mut_in_leaf','syn_ag_mut_rate','non_syn_ag_mut_rate','groups_of_edited_leaves','non_edited_leaves','weak_syn_sites','weak_nonsyn_sites','nonsyn_sites_depletion_factor','strong_syn_sites','strong_non_syn_sites','expected_strong_nonsyn_sites','strong_nonsyn_sites_excess','expected_nonsyn_EG_mutations_from_excess','expected_nonsyn_EG_mutations_from_expected','total_expected_EG_mutations','std','actual_nonsyn_eg_mutations','excess_of_unmutated_sites','p'], name=leaf)})
    
            
        def no_depletion_factor(self, ancestor, intermediate, leaf, leaf_nucl, edited_leaves=None, non_edited_leaves=None,  strong_levels=[0.1,1], editing_level_method='average', confidence_level=0.95, n_random=1000000, sites_recalculation=True):
            
            def calc_no_depletion_factor_model_distribuation(syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,strong_syn_sites,strong_non_syn_sites,confidence_level,n):
                
                try:
                    syn_ag_mut_rate = float(syn_ag_mut)/syn_a
                    non_syn_ag_mut_rate = float(non_syn_ag_mut)/non_syn_a
                    z = 1-(1-confidence_level)/2
                    p_strong_ss=float(strong_syn_sites)/syn_a
                    
                    expected_non_syn_sites_dist = []
                    for i in range(n):
                        p_strong_ss_rand = np.random.normal(p_strong_ss,z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a))
                        while p_strong_ss_rand<0:
                            p_strong_ss_rand = np.random.normal(p_strong_ss,z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a))
                        lam = p_strong_ss_rand*non_syn_a
                        expected_non_syn_sites_dist.append(np.random.poisson(lam))
                        
                    expected_non_syn_eg_mut_dist = []
                    for i in range(n):
                        expected_eg_mut_from_excess = (strong_non_syn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
                        p_non_syn_ag_mut_rate_rand = np.random.normal(non_syn_ag_mut_rate,z*np.sqrt(non_syn_ag_mut_rate*(1-non_syn_ag_mut_rate)/non_syn_a))
                        expected_eg_mut_from_expected_nss = expected_non_syn_sites_dist[i]*p_non_syn_ag_mut_rate_rand
                        expected_non_syn_eg_mut_dist.append(np.random.poisson(expected_eg_mut_from_excess)+np.random.poisson(expected_eg_mut_from_expected_nss))
                        
                    return expected_non_syn_eg_mut_dist
                
                except Exception as e:
                    return str(e)

            
            intermediate_nucl='A'
            self.general_model.define_nodes(ancestor,intermediate,leaf,intermediate_nucl=intermediate_nucl,leaf_nucl=leaf_nucl)
            if sites_recalculation or (self.general_model.edited is None or self.general_model.mutated_editing_sites is None):
                self.general_model.collect_editing_sites(edited_leaves=edited_leaves,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method)
                self.general_model.collect_mutated_editing_sites(leaf, leaf_nucl)
                        
            intermediate_nucl_mat = self.general_model.nucl_mat[self.general_model.nucl_mat[intermediate+'_nuc']==intermediate_nucl]
            edited_mat = self.general_model.edited[self.general_model.edited[intermediate+'_nuc']==intermediate_nucl]
            mutations_mat = self.general_model.mutated_editing_sites[self.general_model.mutated_editing_sites[intermediate+'_nuc']==intermediate_nucl]
            
            syn_a = len(intermediate_nucl_mat[intermediate_nucl_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0])
            non_syn_a = len(intermediate_nucl_mat[intermediate_nucl_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1])
            syn_ag_mut = int(self.general_model.rates[intermediate+'_to_'+leaf][intermediate_nucl+leaf_nucl][0][1])
            non_syn_ag_mut = int(self.general_model.rates[intermediate+'_to_'+leaf][intermediate_nucl+leaf_nucl][1][1])
            strong_syn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0, np.logical_and(edited_mat['combined_editing_level']>strong_levels[0], edited_mat['combined_editing_level']<=strong_levels[1]))])
            strong_nonsyn_sites = len(edited_mat[np.logical_and(edited_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1, np.logical_and(edited_mat['combined_editing_level']>strong_levels[0], edited_mat['combined_editing_level']<=strong_levels[1]))])
            
            actual_strong_nonsyn_eg_mutations = len(mutations_mat[np.logical_and(mutations_mat[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1,np.logical_and(mutations_mat['combined_editing_level']>strong_levels[0], mutations_mat['combined_editing_level']<=strong_levels[1]))])
            
            syn_sites_creation_rate = float(strong_syn_sites)/float(syn_a)
            expected_strong_nonsyn_sites = syn_sites_creation_rate*float(non_syn_a)
            strong_nonsyn_sites_excess = strong_nonsyn_sites-expected_strong_nonsyn_sites
            expected_nonsyn_EG_mutations_from_excess = strong_nonsyn_sites_excess*self.general_model.rates[intermediate+'_to_'+leaf]['A'+leaf_nucl][0][1]/float(syn_a)
            expected_nonsyn_EG_mutations_from_expected = expected_strong_nonsyn_sites*self.general_model.rates[intermediate+'_to_'+leaf]['A'+leaf_nucl][1][1]/float(non_syn_a)
            total_expected_EG_mutations = expected_nonsyn_EG_mutations_from_excess+expected_nonsyn_EG_mutations_from_expected
            excess_of_unmutated_sites = total_expected_EG_mutations-actual_strong_nonsyn_eg_mutations
                      
            liberal_expected_non_syn_eg_mut_dist = calc_no_depletion_factor_model_distribuation(syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,strong_syn_sites,strong_nonsyn_sites,confidence_level,n_random)
            if type(liberal_expected_non_syn_eg_mut_dist)==str:
                liberal_std=liberal_expected_non_syn_eg_mut_dist
                liberal_p=liberal_expected_non_syn_eg_mut_dist
            else:
                liberal_std =np.std(liberal_expected_non_syn_eg_mut_dist)
                liberal_p = float(sum([1 for j in liberal_expected_non_syn_eg_mut_dist if j<actual_strong_nonsyn_eg_mutations]))/float(n_random)
            
            self.general_model.no_depletion_data.update({ancestor+'_to_'+intermediate+'_to_'+leaf+'_'+editing_level_method+'_'+str(strong_levels[0]):pd.Series(data = [ancestor,intermediate,leaf,editing_level_method,strong_levels[0],strong_levels[1],syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,float(syn_ag_mut)/float(syn_a),float(non_syn_ag_mut)/float(non_syn_a),self.general_model.edited_leaves,self.general_model.non_edited_leaves,strong_syn_sites,strong_nonsyn_sites,expected_strong_nonsyn_sites,strong_nonsyn_sites_excess,expected_nonsyn_EG_mutations_from_excess,expected_nonsyn_EG_mutations_from_expected,total_expected_EG_mutations,liberal_std,actual_strong_nonsyn_eg_mutations,excess_of_unmutated_sites,liberal_p],
                                                                                       index = ['ancestor','intermediate','leaf','combined_editing_level_method','strong_levels_lower_limit','strong_levels_upper_limit','intermediate_syn_a','intermediate_non_syn_a','syn_ag_mut_in_leaf','non_syn_ag_mut_in_leaf','syn_ag_mut_rate','non_syn_ag_mut_rate','groups_of_edited_leaves','non_edited_leaves','strong_syn_sites','strong_non_syn_sites','expected_strong_nonsyn_sites','strong_nonsyn_sites_excess','expected_nonsyn_EG_mutations_from_excess','expected_nonsyn_EG_mutations_from_expected','total_expected_EG_mutations','std','actual_nonsyn_eg_mutations','excess_of_unmutated_sites','p'], name=leaf)})
       
        
        def compare_edited_and_unedited_substitution(self, ancestor, intermediate, leaf, intermediate_nucl='A', leaf_nucl='G', syn=True, edited_leaves=None, non_edited_leaves=[],  levels=[0.1,1], editing_level_method='average', sites_recalculation=True):
            
           
            self.general_model.define_nodes(ancestor,intermediate,leaf,intermediate_nucl=intermediate_nucl,leaf_nucl=leaf_nucl)
            if sites_recalculation or (self.general_model.edited is None):
                self.general_model.collect_editing_sites(edited_leaves=edited_leaves,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method)
                        
            intermediate_mat = self.general_model.nucl_mat[self.general_model.nucl_mat[intermediate+'_nuc']==intermediate_nucl].copy()
            edited = self.general_model.edited[self.general_model.edited[intermediate+'_nuc']==intermediate_nucl].copy()
            unedited = intermediate_mat[~intermediate_mat.index.isin(list(edited.index))]
            edited = edited[np.logical_and(edited['combined_editing_level']>levels[0], edited['combined_editing_level']<=levels[1])]
            
            
            if syn:
                name = intermediate+'_'+leaf+'_syn_'+intermediate_nucl+leaf_nucl+'_'+editing_level_method+'_levels_'+str(levels[0])+'-'+str(levels[1])
                unedited=unedited[unedited[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0]
                edited = edited[edited[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==0]
                recoding = 'syn'
            
            else:
                name = intermediate+'_'+leaf+'_nonsyn_'+intermediate_nucl+leaf_nucl+'_'+editing_level_method+'_levels_'+str(levels[0])+'-'+str(levels[1])
                unedited=unedited[unedited[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1]
                edited=edited[edited[intermediate+'_'+intermediate_nucl+leaf_nucl+'_recoding']==1]
                recoding = 'nonsyn'
            
            edited_mutations = edited[edited[leaf+'_nuc']==leaf_nucl]
            unedited_mutations = unedited[unedited[leaf+'_nuc']==leaf_nucl]
            
            edited_n = len(edited)
            edited_mutations_n = len(edited_mutations)
            unedited_n = len(unedited)
            unedited_mutations_n = len(unedited_mutations)
            mismatch=intermediate_nucl+leaf_nucl
            
            
            
            self.general_model.mutations_count.update({name:pd.Series(data = [ancestor,intermediate,leaf,mismatch,recoding,editing_level_method,levels[0],levels[1],self.general_model.edited_leaves,self.general_model.non_edited_leaves,edited_n,edited_mutations_n,unedited_n,unedited_mutations_n],
                                                                      index = ['ancestor','intermediate','leaf','mismatch','type','editing_level_method','editing_level_lower_bound','editing_level_upper_bound','edited_leaves','unedited_leaves','edited','edited_mutations','unedited','unedited_mutations'], name=leaf)})
       
        




if __name__=='__main__':
    
    nucl_mat_file = sys.argv[1]
    tree=sys.argv[2]
    ancestor = sys.argv[3]
    intermediate = sys.argv[4]
    leaf = sys.argv[5]
    # leaf_mutated_nucl = sys.argv[6]
    # ancestor_nucl = sys.argv[7]
    ancestor_nucl = None
    # leaf_mutated_nucl = 'G'
    only_N1=True
    if ancestor_nucl=='None':
        ancestor_nucl=None
    
    trees={'coleoids_rooted':"((oct,bim)O,((sep,squ)S,(bob,lin)B)D)C",
           'coleoids_unrooted':"((oct,bim)O,(sep,squ)S,(bob,lin)B)C",
           'all8_rooted':"(apl,(nau,((oct,bim)O,((sep,squ)S,(bob,lin)B)D)C)N1)N0",
           'all8_unrooted':"(apl,nau,((oct,bim)O,((sep,squ)S,(bob,lin)B)D)C)N0",
           'no_boblin_rooted':"(apl,(nau,((oct,bim)O,(sep,squ)S)C)N1)N0",
           'no_boblin_unrooted':"(apl,nau,((oct,bim)O,(sep,squ)S)C)N0",
           'all8_rooted_ncbi':"(apl,(nau,((oct,bim)O,(squ,(bob,(sep,lin)S1)S0)D)C)N1)N0",   
           'coleoids_rooted_ncbi':"((oct,bim)O,(squ,(bob,(sep,lin)S1)S0)D)C"}
    
    newick_tree_str=trees[tree]
    
    if 'all8' in tree:
        animals=['apl','nau','oct','bim','sep','squ','bob','lin']
    elif 'coleoids' in tree:
        animals=['oct','bim','sep','squ','bob','lin']
    elif 'no_boblin' in tree:
        animals=['apl','nau','oct','bim','sep','squ']
        
    
# =============================================================================
#     outpath = '/'.join(nucl_mat_file.split('/')[0:-1])+'/'
#     print('Reading matrix')
#     nucl_mat=pd.read_csv(nucl_mat_file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')
#     nucl_mat=nucl_mat.apply(pd.to_numeric, errors='ignore')
#     print(str(len(nucl_mat)) + ' rows in nucl matrix')
# =============================================================================
    
    
# =============================================================================
#     animals=['oct','bim','sep','squ','bob','lin']
#     columns=['anmimal','syn_sites','syn_a','nonsyn_sites','nonsyn_a']
#     data=[]
#     for a in animals:
#         n(nucl_mat[np.logical_and(nucl_mat[a+'_edited']==1,nucl_mat[a+'_AG_recoding']==0)])
#         syn_a = len(nu
#         syn_sites = lecl_mat[np.logical_and(nucl_mat[a+'_nuc']=="A",nucl_mat[a+'_AG_recoding']==0)])
#         nonsyn_sites = len(nucl_mat[np.logical_and(nucl_mat[a+'_edited']==1,nucl_mat[a+'_AG_recoding']==1)])
#         nonsyn_a = len(nucl_mat[np.logical_and(nucl_mat[a+'_nuc']=="A",nucl_mat[a+'_AG_recoding']==1)])
#     
#         data.append((a,syn_sites,syn_a,nonsyn_sites,nonsyn_a))
#     
#     df = pd.DataFrame(data=data,columns=columns)
#     df.to_csv(outpath+'syn_non_syn_incedence',sep='\t',index=False)
# =============================================================================
    
    
# =============================================================================
#     outpath = '/'.join(nucl_mat_file.split('/')[0:-1])+'/'
#     print('Reading matrix')
#     nucl_mat=pd.read_csv(nucl_mat_file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')
#     nucl_mat=nucl_mat.apply(pd.to_numeric, errors='ignore')
#     print(str(len(nucl_mat)) + ' rows in nucl matrix')
#     
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
# 
#     hpm.calc_editing_types_rates('sep')
#     hpm.calc_editing_types_rates('squ')
#     hpm.calc_editing_types_rates('oct')
#     hpm.calc_editing_types_rates('bim')
#     hpm.calc_editing_types_rates('bob')
#     hpm.calc_editing_types_rates('lin')
#     hpm.calc_editing_types_rates(['bim','oct'])
#     hpm.calc_editing_types_rates(['sep','squ'])
#     hpm.calc_editing_types_rates(['bob','lin'])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin'])
#     hpm.calc_editing_types_rates(['sep','squ','oct','bim'])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin','oct','bim'])
#     hpm.hpm.write_data(outpath,file_name='editing_types_rates_strong',data_to_write=['editing_types_rates'],file_type='csv',sub_name='')
# 
#     del(hpm)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_editing_types_rates('sep',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('squ',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('oct',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('bim',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('bob',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('lin',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['bim','oct'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','squ'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['bob','lin'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','squ','oct','bim'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0.1,1])
#     hpm.hpm.write_data(outpath,file_name='editing_types_rates_strong',data_to_write=['editing_types_rates'],file_type='csv',sub_name='')
#     
#     del(hpm)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_editing_types_rates('sep',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('squ',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('oct',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('bim',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('bob',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('lin',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['bim','oct'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','squ'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['bob','lin'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','squ','oct','bim'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0,0.1])
#     hpm.hpm.write_data(outpath,file_name='editing_types_rates_weak',data_to_write=['editing_types_rates'],file_type='csv',sub_name='')
# #    
# 
#     del(hpm)
#     hpm=Hypothesis.HPM(hyp) 
#     hpm.collect_editing_levels_distributions_by_types('sep')
#     hpm.collect_editing_levels_distributions_by_types('squ')
#     hpm.collect_editing_levels_distributions_by_types('oct')
#     hpm.collect_editing_levels_distributions_by_types('bim')
#     hpm.collect_editing_levels_distributions_by_types('bob')
#     hpm.collect_editing_levels_distributions_by_types('lin')
#     hpm.collect_editing_levels_distributions_by_types(['bim','oct'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','squ'])
#     hpm.collect_editing_levels_distributions_by_types(['bob','lin'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','squ','bob','lin'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','squ','oct','bim'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','squ','bob','lin','oct','bim'])
#     hpm.hpm.write_data(outpath,file_name='editing_levels_by_types',data_to_write=['editing_levels_by_types'],file_type='csv',sub_name='')
#     
#     # hpm.calc_acestral_states_rates('sep')
#     # hpm.calc_acestral_states_rates('squ')
#     # hpm.calc_acestral_states_rates('oct')
#     # hpm.calc_acestral_states_rates('bim')
#     # hpm.calc_acestral_states_rates('bob')
#     # hpm.calc_acestral_states_rates('lin')
#     # hpm.calc_acestral_states_rates(['bim','oct'])
#     # hpm.calc_acestral_states_rates(['sep','squ'])
#     # hpm.calc_acestral_states_rates(['bob','lin'])
#     # hpm.calc_acestral_states_rates(['sep','squ','bob','lin'])
#     # hpm.calc_acestral_states_rates(['sep','squ','oct','bim'])
#     # hpm.calc_acestral_states_rates(['sep','squ','bob','lin','oct','bim'])
#     # hpm.hpm.write_data(outpath,file_name='ancestral_editing_types_rates',data_to_write=['editing_ancestral_rates'],file_type='csv',sub_name='')
# 
#     # del(hpm)
#     # hpm=Hypothesis.HPM(hyp)
#     # hpm.calc_acestral_states_rates('sep',editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates('squ',editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates('oct',editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates('bim',editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates('bob',editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates('lin',editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates(['bim','oct'],editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates(['sep','squ'],editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates(['bob','lin'],editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates(['sep','squ','bob','lin'],editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates(['sep','squ','oct','bim'],editing_levels=[0.1,1])
#     # hpm.calc_acestral_states_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0.1,1])
#     # hpm.hpm.write_data(outpath,file_name='ancestral_editing_types_rates_strong',data_to_write=['editing_ancestral_rates'],file_type='csv',sub_name='')
#     
#     # del(hpm)
#     # hpm=Hypothesis.HPM(hyp)
#     # hpm.calc_acestral_states_rates('sep',editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates('squ',editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates('oct',editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates('bim',editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates('bob',editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates('lin',editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates(['bim','oct'],editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates(['sep','squ'],editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates(['bob','lin'],editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates(['sep','squ','bob','lin'],editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates(['sep','squ','oct','bim'],editing_levels=[0,0.1])
#     # hpm.calc_acestral_states_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0,0.1])
#     # hpm.hpm.write_data(outpath,file_name='ancestral_editing_types_rates_weak',data_to_write=['editing_ancestral_rates'],file_type='csv',sub_name='')
#     
#     # del(hpm)
#     # hpm=Hypothesis.HPM(hyp)
#     # hpm.collect_editing_levels_distributions_by_ancestral_state('sep')
#     # hpm.collect_editing_levels_distributions_by_ancestral_state('squ')
#     # hpm.collect_editing_levels_distributions_by_ancestral_state('oct')
#     # hpm.collect_editing_levels_distributions_by_ancestral_state('bim')
#     # hpm.collect_editing_levels_distributions_by_ancestral_state('bob')
#     # hpm.collect_editing_levels_distributions_by_ancestral_state('lin')
#     # hpm.collect_editing_levels_distributions_by_ancestral_state(['bim','oct'])
#     # hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','squ'])
#     # hpm.collect_editing_levels_distributions_by_ancestral_state(['bob','lin'])
#     # hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','squ','bob','lin'])
#     # hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','squ','oct','bim'])
#     # hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','squ','bob','lin','oct','bim'])
#     # hpm.hpm.write_data(outpath,file_name='ancestral_editing_leveils_by_types',data_to_write=['editing_levels_by_ancestral_state'],file_type='csv',sub_name='')
# =============================================================================
    
    

# =============================================================================
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
# 
#     hpm.calc_editing_types_rates('sep')
#     hpm.calc_editing_types_rates('squ')
#     hpm.calc_editing_types_rates('oct')
#     hpm.calc_editing_types_rates('bim')
#     hpm.calc_editing_types_rates('bob')
#     hpm.calc_editing_types_rates('lin')
#     hpm.calc_editing_types_rates(['bim','oct'])
#     hpm.calc_editing_types_rates(['sep','lin'])
#     hpm.calc_editing_types_rates(['sep','bob','lin'])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin'])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin','oct','bim'])
#     hpm.hpm.write_data(outpath,file_name='editing_types_rates_strong',data_to_write=['editing_types_rates'],file_type='csv',sub_name='')
# 
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_editing_types_rates('sep',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('squ',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('oct',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('bim',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('bob',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates('lin',editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['bim','oct'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','lin'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','bob','lin'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin'],editing_levels=[0.1,1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0.1,1])
#     hpm.hpm.write_data(outpath,file_name='editing_types_rates_strong',data_to_write=['editing_types_rates'],file_type='csv',sub_name='')
#     
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_editing_types_rates('sep',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('squ',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('oct',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('bim',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('bob',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates('lin',editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['bim','oct'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','lin'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','bob','lin'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin'],editing_levels=[0,0.1])
#     hpm.calc_editing_types_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0,0.1])
#     hpm.hpm.write_data(outpath,file_name='editing_types_rates_weak',data_to_write=['editing_types_rates'],file_type='csv',sub_name='')
# 
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.collect_editing_levels_distributions_by_types('sep')
#     hpm.collect_editing_levels_distributions_by_types('squ')
#     hpm.collect_editing_levels_distributions_by_types('oct')
#     hpm.collect_editing_levels_distributions_by_types('bim')
#     hpm.collect_editing_levels_distributions_by_types('bob')
#     hpm.collect_editing_levels_distributions_by_types('lin')
#     hpm.collect_editing_levels_distributions_by_types(['bim','oct'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','lin'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','bob','lin'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','squ','bob','lin'])
#     hpm.collect_editing_levels_distributions_by_types(['sep','squ','bob','lin','oct','bim'])
#     hpm.hpm.write_data(outpath,file_name='editing_levels_by_types',data_to_write=['editing_levels_by_types'],file_type='csv',sub_name='')
# 
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_acestral_states_rates('sep')
#     hpm.calc_acestral_states_rates('squ')
#     hpm.calc_acestral_states_rates('oct')
#     hpm.calc_acestral_states_rates('bim')
#     hpm.calc_acestral_states_rates('bob')
#     hpm.calc_acestral_states_rates('lin')
#     hpm.calc_acestral_states_rates(['bim','oct'])
#     hpm.calc_acestral_states_rates(['sep','lin'])
#     hpm.calc_acestral_states_rates(['sep','bob','lin'])
#     hpm.calc_acestral_states_rates(['sep','squ','bob','lin'])
#     hpm.calc_acestral_states_rates(['sep','squ','bob','lin','oct','bim'])
#     hpm.hpm.write_data(outpath,file_name='ancestral_editing_types_rates',data_to_write=['editing_ancestral_rates'],file_type='csv',sub_name='')
# 
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_acestral_states_rates('sep',editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates('squ',editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates('oct',editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates('bim',editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates('bob',editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates('lin',editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates(['bim','oct'],editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates(['sep','lin'],editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates(['sep','bob','lin'],editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates(['sep','squ','bob','lin'],editing_levels=[0.1,1])
#     hpm.calc_acestral_states_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0.1,1])
#     hpm.hpm.write_data(outpath,file_name='ancestral_editing_types_rates_strong',data_to_write=['editing_ancestral_rates'],file_type='csv',sub_name='')
#     
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.calc_acestral_states_rates('sep',editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates('squ',editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates('oct',editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates('bim',editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates('bob',editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates('lin',editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates(['bim','oct'],editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates(['sep','lin'],editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates(['sep','bob','lin'],editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates(['sep','squ','bob','lin'],editing_levels=[0,0.1])
#     hpm.calc_acestral_states_rates(['sep','squ','bob','lin','oct','bim'],editing_levels=[0,0.1])
#     hpm.hpm.write_data(outpath,file_name='ancestral_editing_types_rates_weak',data_to_write=['editing_ancestral_rates'],file_type='csv',sub_name='')
#     
#     del(hpm)
#     del(hyp)
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     hpm.collect_editing_levels_distributions_by_ancestral_state('sep')
#     hpm.collect_editing_levels_distributions_by_ancestral_state('squ')
#     hpm.collect_editing_levels_distributions_by_ancestral_state('oct')
#     hpm.collect_editing_levels_distributions_by_ancestral_state('bim')
#     hpm.collect_editing_levels_distributions_by_ancestral_state('bob')
#     hpm.collect_editing_levels_distributions_by_ancestral_state('lin')
#     hpm.collect_editing_levels_distributions_by_ancestral_state(['bim','oct'])
#     hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','lin'])
#     hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','bob','lin'])
#     hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','squ','bob','lin'])
#     hpm.collect_editing_levels_distributions_by_ancestral_state(['sep','squ','bob','lin','oct','bim'])
#     hpm.hpm.write_data(outpath,file_name='ancestral_editing_leveils_by_types',data_to_write=['editing_levels_by_ancestral_state'],file_type='csv',sub_name='')
# =============================================================================
    
# =============================================================================
#     outpath = '/'.join(nucl_mat_file.split('/')[0:-1])+'/'
#     print('Reading matrix')
#     nucl_mat=pd.read_csv(nucl_mat_file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')
#     nucl_mat=nucl_mat.apply(pd.to_numeric, errors='ignore')
#     print(str(len(nucl_mat)) + ' rows in nucl matrix')
#     
#     filtered_nucl_mat = nucl_mat[nucl_mat['gapped_animals']==0]
# 
#     hyp=Hypothesis(filtered_nucl_mat.copy(),tree_str=newick_tree_str)
#     hpm=Hypothesis.HPM(hyp)
#     
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,leaves_groups=[('sep','bim','bob'),('sep','bim','lin'),('sep','oct','bob'),('sep','oct','lin'),('squ','bim','bob'),('squ','bim','lin'),('squ','oct','bob'),('squ','oct','lin')])
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,leaves_groups=[('sep','bim','bob'),('sep','bim','lin'),('sep','oct','bob'),('sep','oct','lin'),('squ','bim','bob'),('squ','bim','lin'),('squ','oct','bob'),('squ','oct','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,leaves_groups=[('sep','squ','bob'),('sep','squ','lin'),('sep','bob','lin'),('squ','bob','lin')])
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,leaves_groups=[('sep','squ','bob'),('sep','squ','lin'),('sep','bob','lin'),('squ','bob','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin','oct','bim'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')])
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin','oct','bim'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim'),('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')])
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim'),('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')])
#     hpm.find_substituations_sites(ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')],count_multiple_subs_per_sites=True)
# 
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,leaves_groups=[('sep','bim','bob'),('sep','bim','lin'),('sep','oct','bob'),('sep','oct','lin'),('squ','bim','bob'),('squ','bim','lin'),('squ','oct','bob'),('squ','oct','lin')])
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,leaves_groups=[('sep','bim','bob'),('sep','bim','lin'),('sep','oct','bob'),('sep','oct','lin'),('squ','bim','bob'),('squ','bim','lin'),('squ','oct','bob'),('squ','oct','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,leaves_groups=[('sep','squ','bob'),('sep','squ','lin'),('sep','bob','lin'),('squ','bob','lin')])
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,leaves_groups=[('sep','squ','bob'),('sep','squ','lin'),('sep','bob','lin'),('squ','bob','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin','oct','bim'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')])
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin','oct','bim'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim'),('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')])
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim'),('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')])
#     hpm.find_substituations_sites(editing_levels=[0.1,1],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')],count_multiple_subs_per_sites=True)
# 
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,leaves_groups=[('sep','bim','bob'),('sep','bim','lin'),('sep','oct','bob'),('sep','oct','lin'),('squ','bim','bob'),('squ','bim','lin'),('squ','oct','bob'),('squ','oct','lin')])
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,leaves_groups=[('sep','bim','bob'),('sep','bim','lin'),('sep','oct','bob'),('sep','oct','lin'),('squ','bim','bob'),('squ','bim','lin'),('squ','oct','bob'),('squ','oct','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,leaves_groups=[('sep','squ','bob'),('sep','squ','lin'),('sep','bob','lin'),('squ','bob','lin')])
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,leaves_groups=[('sep','squ','bob'),('sep','squ','lin'),('sep','bob','lin'),('squ','bob','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin','oct','bim'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')])
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin','oct','bim'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim'),('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')])
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim'),('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')],count_multiple_subs_per_sites=True)
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')])
#     hpm.find_substituations_sites(editing_levels=[0,0.99999999999999994],ancestor_for_type=intermediate,examined_animals=['sep','squ','bob','lin'],leaves_groups=[('sep','bob'),('sep','lin'),('squ','bob'),('squ','lin')],count_multiple_subs_per_sites=True)
#     
#     hpm.hpm.write_data(outpath,file_name='subs_count_from_'+intermediate,data_to_write=['sites_substitutions'],file_type='csv',sub_name='')
# =============================================================================
    

# =============================================================================
#     outpath = '/'.join(nucl_mat_file.split('/')[0:-1])+'/'
#     print('Reading matrix')
#     nucl_mat=pd.read_csv(nucl_mat_file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')
#     nucl_mat=nucl_mat.apply(pd.to_numeric, errors='ignore')
#     print(str(len(nucl_mat)) + ' rows in nucl matrix')
#     
#     non_edited_leaves = []
#     nucl_mat = nucl_mat[nucl_mat[intermediate+'_nuc']=="A"].copy()
#     
#     if ancestor_nucl is not None:
#         if 'unrooted' in tree:
#             nucl_mat = nucl_mat[nucl_mat['N0_nuc']==ancestor_nucl].copy()
#         elif 'rooted' in tree:
#             if only_N1:
#                 nucl_mat = nucl_mat[nucl_mat['N1_nuc']==ancestor_nucl].copy()
#             else:
#                 nucl_mat = nucl_mat[np.logical_or(nucl_mat['N1_nuc']==ancestor_nucl,nucl_mat['N0_nuc']==ancestor_nucl)].copy()
#         
#     hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#     gm = Hypothesis.General_model(hyp)
#     file_name='ances'+str(ancestor_nucl)+'_'+ancestor+'_'+intermediate+'_'+leaf+'_no_msa_filter'
#     for editing_level_method in ['average','maximal','minimal']:
#         for strong_levels in [[0.1,1],[0.15,1],[0.2,1],[0.3,1]]:
#             print('Strict - editing levels method: '+editing_level_method+' strong levels:'+str(strong_levels))
# #            gm.strict(ancestor,intermediate,leaf,leaf_mutated_nucl,strong_levels=strong_levels,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method,sites_recalculation=False)
#             for weak_levels in [[0,0.01],[0,0.02],[0,0.05]]:
#                 print('liberal - editing levels method: '+editing_level_method+' strong levels:'+str(strong_levels)+'_weak_levels:'+str(weak_levels))
#                 gm.liberal(ancestor,intermediate,leaf,leaf_mutated_nucl,weak_levels=weak_levels,strong_levels=strong_levels,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method,sites_recalculation=False)
# #    gm.general_model.write_data(outpath,file_name=file_name,file_type='csv',data_to_write = ['strict','liberal'])
#     gm.general_model.write_data(outpath,file_name=file_name,file_type='csv',data_to_write = ['liberal'])
#     
#     
#     iden_filter_cols = [col for col in nucl_mat.columns if ('aa_range' in col and leaf not in col and all([a in col for a in animals if a!=leaf ]))]
#     print('Creating general models from C')
#     for r in [5,10,15,20]:
#         filter_col = [c for c in iden_filter_cols if str(r) in c][0] 
#         for iden in [0.3,0.5,0.7]:
#             print('filtering '+filter_col+' for '+str(iden)+' and above')
#             filtered_nucl_mat = nucl_mat[nucl_mat[filter_col]>=iden]
#             hyp=Hypothesis(filtered_nucl_mat.copy(),tree_str=newick_tree_str)
#             gm = Hypothesis.General_model(hyp)
#             file_name='ances'+str(ancestor_nucl)+'_'+ancestor+'_'+intermediate+'_'+leaf+'_iden'+str(iden)+'_in_range'+str(r)
#             for editing_level_method in ['average','maximal','minimal']:
#                 for strong_levels in [[0.1,1],[0.15,1],[0.2,1],[0.3,1]]:
#                     print('Strict - editing levels method: '+editing_level_method+' strong levels:'+str(strong_levels))
# #                    gm.strict(ancestor,intermediate,leaf,leaf_mutated_nucl,strong_levels=strong_levels,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method,sites_recalculation=False)
#                     for weak_levels in [[0,0.01],[0,0.02],[0,0.05]]:
#                         print('liberal - editing levels method: '+editing_level_method+' strong levels:'+str(strong_levels)+'_weak_levels:'+str(weak_levels))
#                         gm.liberal(ancestor,intermediate,leaf,leaf_mutated_nucl,weak_levels=weak_levels,strong_levels=strong_levels,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method,sites_recalculation=False)
# #            gm.general_model.write_data(outpath,file_name=file_name,file_type='csv',data_to_write = ['strict','liberal'])
#             gm.general_model.write_data(outpath,file_name=file_name,file_type='csv',data_to_write = ['liberal'])
# =============================================================================
    
    
    outpath = '/'.join(nucl_mat_file.split('/')[0:-1])+'/'
    print('Reading matrix')
    nucl_mat=pd.read_csv(nucl_mat_file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')
    nucl_mat=nucl_mat.apply(pd.to_numeric, errors='ignore')
    print(str(len(nucl_mat)) + ' rows in nucl matrix')
    
    if ancestor_nucl is not None:
        if 'unrooted' in tree:
            nucl_mat = nucl_mat[nucl_mat['N0_nuc']==ancestor_nucl].copy()
        elif 'rooted' in tree:
            if only_N1:
                nucl_mat = nucl_mat[nucl_mat['N1_nuc']==ancestor_nucl].copy()
            else:
                nucl_mat = nucl_mat[np.logical_or(nucl_mat['N1_nuc']==ancestor_nucl,nucl_mat['N0_nuc']==ancestor_nucl)].copy()
    
    file_name='mutations_count_in_edited_unedited_sites_'+ancestor+'_'+intermediate+'_'+leaf
    iden_filter_cols = [col for col in nucl_mat.columns if ('aa_range' in col and leaf not in col and all([a in col for a in animals if a!=leaf ]))]
    filter_col = [c for c in iden_filter_cols if str(10) in c][0] 
    filtered_nucl_mat = nucl_mat[nucl_mat[filter_col]>=0.3].copy()
    hyp=Hypothesis(filtered_nucl_mat,tree_str=newick_tree_str)
    gm = Hypothesis.General_model(hyp)
    gm.compare_edited_and_unedited_substitution(ancestor, intermediate, leaf, intermediate_nucl='A', leaf_nucl='G', syn=True)
    gm.compare_edited_and_unedited_substitution(ancestor, intermediate, leaf, intermediate_nucl='A', leaf_nucl='C', syn=False)
    gm.compare_edited_and_unedited_substitution(ancestor, intermediate, leaf, intermediate_nucl='A', leaf_nucl='T', syn=False)
    gm.general_model.write_data(outpath,file_name=file_name,file_type='csv',data_to_write = ['mutations_count'])
    
    
    
# =============================================================================
#     arguments_groups=[('C','D','sep',[('lin','bob','oct','bim')]),
#                       ('C','D','squ',[('lin','bob','oct','bim')]),
#                       ('C','D','bob',[('sep','squ','oct','bim')]),
#                       ('C','D','lin',[('sep','squ','oct','bim')]),
#                       ('C','B','bob',[('lin','sep','squ','oct','bim')]),
#                       ('C','B','lin',[('bob','sep','squ','oct','bim')]),
#                       ('C','S','sep',[('lin','bob','squ','oct','bim')]),
#                       ('C','S','squ',[('lin','bob','sep','oct','bim')])]
#     
#     for args in arguments_groups:
#         
#         non_edited_leaves = []
#         ancestor_nucl=None
#         ancestor = args[0]
#         intermediate = args[1]
#         leaf = args[2]
#         edited_leaves = args[3]
#         
#         print('\n'+ancestor+'_'+intermediate+'_'+leaf+':')
#         outpath = nucl_mat_file
#         print('Reading matrix')
#         file = nucl_mat_file+intermediate+'nucl_is_A'
#         nucl_mat=pd.read_csv(file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')
#         nucl_mat=nucl_mat.apply(pd.to_numeric, errors='ignore')
#         print(str(len(nucl_mat)) + ' rows in nucl matrix')
#         
#         
#         non_edited_leaves = []
#         nucl_mat = nucl_mat[nucl_mat[intermediate+'_nuc']=="A"].copy()
#         
#         if ancestor_nucl is not None:
#             if 'unrooted' in tree:
#                 nucl_mat = nucl_mat[nucl_mat['N0_nuc']==ancestor_nucl].copy()
#             elif 'rooted' in tree:
#                 nucl_mat = nucl_mat[np.logical_or(nucl_mat['N1_nuc']==ancestor_nucl,nucl_mat['N0_nuc']==ancestor_nucl)].copy()
#             
#         hyp=Hypothesis(nucl_mat.copy(),tree_str=newick_tree_str)
#         gm = Hypothesis.General_model(hyp)
#         file_name='conserved_sites'+ancestor+'_'+intermediate+'_'+leaf+'_no_msa_filter'
#         for editing_level_method in ['average','maximal','minimal']:
#             for strong_levels in [[0.1,1],[0.15,1],[0.2,1],[0.3,1]]:
#                 print('liberal - editing levels method: '+editing_level_method+' strong levels:'+str(strong_levels))
#                 gm.no_depletion_factor(ancestor,intermediate,leaf,leaf_mutated_nucl,edited_leaves=edited_leaves,strong_levels=strong_levels,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method,sites_recalculation=False)
#         gm.general_model.write_data(outpath,file_name=file_name,data_to_write = ['no_depletion'],file_type='csv')
#     
#         iden_filter_cols = [col for col in nucl_mat.columns if ('aa_range' in col and leaf not in col and all([a in col for a in animals if a!=leaf ]))]
#         print('Creating general models from C')
#         for r in [5,10,15,20]:
#             filter_col = [c for c in iden_filter_cols if str(r) in c][0] 
#             for iden in [0.3,0.5,0.7]:
#                 print('filtering '+filter_col+' for '+str(iden)+' and above')
#                 filtered_nucl_mat = nucl_mat[nucl_mat[filter_col]>=iden]
#                 hyp=Hypothesis(filtered_nucl_mat.copy(),tree_str=newick_tree_str)
#                 gm = Hypothesis.General_model(hyp)
#                 file_name='conserved_sites'+ancestor+'_'+intermediate+'_'+leaf+'_iden'+str(iden)+'_in_range'+str(r)
#                 for editing_level_method in ['average','maximal','minimal']:
#                     for strong_levels in [[0.1,1],[0.15,1],[0.2,1],[0.3,1]]:
#                         print('liberal - editing levels method: '+editing_level_method+' strong levels:'+str(strong_levels))
#                         gm.no_depletion_factor(ancestor,intermediate,leaf,leaf_mutated_nucl,strong_levels=strong_levels,edited_leaves=edited_leaves,non_edited_leaves=non_edited_leaves,editing_level_method=editing_level_method,sites_recalculation=False)
#                 gm.general_model.write_data(outpath,file_name=file_name,data_to_write = ['no_depletion'],file_type='csv')
# =============================================================================
