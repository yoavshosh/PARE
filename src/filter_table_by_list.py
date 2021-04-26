# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 01:50:27 2020

@author: shosh
"""

import pandas as pd
import sys

nucl_mat_file = sys.argv[1]
l_file = sys.argv[2]
on_col = sys.argv[3]

nucl_mat=pd.read_csv(nucl_mat_file,sep='\t',error_bad_lines=False, index_col=False, dtype='unicode')

l = open(l_file, "r")
genes_subset = l.readlines()
l.close()
genes_subset = [g.rstrip().lstrip() for g in genes_subset]
filtered_nucl_mat = nucl_mat[nucl_mat[on_col].isin(genes_subset)]
filtered_nucl_mat.to_csv(nucl_mat_file+'_from_'+l_file.split('/')[-1], sep='\t', index=False)