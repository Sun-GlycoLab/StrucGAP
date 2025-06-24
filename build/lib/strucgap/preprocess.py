# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import copy
import re
import os
from ast import literal_eval
from scipy import stats
from scipy.stats import kstest
import functools
import operator
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
import random
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import chi2
# pd.set_option('display.max_columns',7)
from gprofiler import GProfiler
import gseapy as gp
import ast
from tqdm import tqdm
from sklearn.impute import KNNImputer
from sklearn.preprocessing import RobustScaler
from sklearn import preprocessing
import statistics        
import requests
from itertools import chain
import math
from math import ceil
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.utils import ImageReader
from reportlab.lib.utils import simpleSplit
from pyecharts.charts import Radar
from pyecharts import options as opts
from pyecharts_snapshot.main import make_a_snapshot
from snapshot_phantomjs import snapshot
from pyecharts.render import make_snapshot
from pyecharts.globals import RenderType
from svglib.svglib import svg2rlg
import cairosvg
import fitz
from scipy.stats import spearmanr
from reportlab.graphics import renderPDF
from pyecharts.charts import Polar
from pyecharts.charts import Funnel
from pyecharts.charts import Parallel
from pyecharts.charts import Pie
from pyecharts.charts import Sankey
from pyecharts.charts import Sunburst
import palettable.colorbrewer.qualitative as brewer_qualitative
import palettable.cartocolors.qualitative as carto_qualitative
from pyecharts.charts import Boxplot
from pyecharts.charts import Bar
from pyecharts.commons.utils import JsCode
import scipy.cluster.hierarchy as sch
from pyecharts.charts import HeatMap
import PyComplexHeatmap as pch
import matplotlib.pylab as plt
from pyecharts.charts import Line
from pyecharts.charts import Scatter
from pyecharts.charts import Tree
import matplotlib.pyplot as plt
from venn import venn
from itertools import product
from scipy.cluster.hierarchy import linkage
import seaborn as sns
from upsetplot import UpSet, generate_counts
from pyecharts.charts import Page
from pyecharts.charts import Grid
from PyComplexHeatmap import *
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import plotly.io as pio
import matplotlib.colors as mcolors
from matplotlib.patches import Ellipse
import requests
import matplotlib.colors as mcolors
from matplotlib.patches import Wedge
import networkx as nx
import pickle
import types
import werkzeug.local
from difflib import get_close_matches
from datetime import datetime
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.platypus import Image
from typing import Dict, List
from reportlab.lib.utils import ImageReader
from PIL import Image, ImageChops
import matplotlib
from reportlab.lib.colors import HexColor
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
## 数据质控模块--7
class StrucGAP_Preprocess:
    def __init__(self, data_manager, search_engine = 'StrucGP',
                 data_dir=None, sample_group_data_dir=None, data_sheet_name=None, branch_list_dir=None):
        """        
        Parameters:
            search_engine: Support multiple search_engine in ['StrucGP','MSFragger-Glyco','pGlyco3','Glyco-Decipher'].
        """
        if data_dir == None:
            data_dir = input(f"Please enter your data file path (such as: 'D:\\doctor\\analysisys\\data\\mouse uterus.xlsx'): ")
            
        data = pd.read_excel(data_dir, sheet_name=data_sheet_name or 0)
        
        if sample_group_data_dir == None:
            sample_group_data_dir = input(f"Please enter your sample group data file path (such as: 'D:\\doctor\\analysisys\\data\\sample_group.xlsx'): ")
        sample_group = pd.read_excel(sample_group_data_dir)
        sample_group = sample_group.set_index('sample',drop=True)
        self.sample_group = sample_group
 
        if branch_list_dir is not None:
            branch_list = pd.read_excel(branch_list_dir)
            self.branch_list = branch_list['Structure coding']
        self.branch_list_dir = branch_list_dir
        # data = data[~data['PeptideSequence+structure_coding+ProteinID'].duplicated()]
        self.data = data
        self.data_fdr_filtered = None
        self.data_peptide_fdr_data = None
        self.data_glycan_fdr_data = None
        self.data_outliers_filtered = None
        self.data_cv_filtered = None
        self.data_psm_filtered = None
        #######################################################################
        if search_engine not in ['StrucGP','MSFragger-Glyco','pGlyco3','Glyco-Decipher']:
            print("Select search_engine in ['StrucGP','MSFragger-Glyco','pGlyco3','Glyco-Decipher']")
            search_engine = 'StrucGP'
            if search_engine is not 'StrucGP':
                print('You can only use a samll subset of the functions in StrucGAP_GlycanStructure and StrucGAP_GlycoSite module!')
        self.search_engine = search_engine
        #######################################################################
        self.data_manager = data_manager
        self.data_manager.register_module('StrucGAP_Preprocess', self, {'input_data': data_dir, 'sample_group': sample_group_data_dir})
        self.data_manager.log_params('StrucGAP_Preprocess', 'input_data', {'input_data': data_dir, 'sample_group': sample_group_data_dir})
        
    def median_cheng(self, data):
        """An auxiliary function called by other functions to calculates the median."""
        filtered_data = [x for x in data if not np.isnan(x)]
        filtered_data.sort()
        half = len(filtered_data) // 2
        if not filtered_data:
            return np.nan 
        if len(filtered_data) % 2 == 0:
            return (filtered_data[half - 1] * filtered_data[half]) ** 0.5
        else:
            return filtered_data[half]
    
    def data_cleaning(self, data_type=None):
        """
        Data cleaning and glycan substructure features extraction.
        
        Parameters:
            data_type: output data of StrucGP with TMT or label-free type. Select data_type in ['tmt', 'label free'].
        
        Returns:
            self.data (cleaned data). 
            
        Return type:
            dataframe
        
        """
        if self.search_engine == 'StrucGP':
            if data_type not in ['tmt', 'label free']:
                input("Please select data_type in ['tmt', 'label free']")
                data_type='tmt'
            
            self.data['structure_coding'] = self.data['structure_coding'].str.replace(r'\+Ammonium\(\+17\)', '', regex=True)
            self.data['Glycosite_Position'] = self.data['Glycosite_Position'].astype(str)
            self.data['GlycanComposition'] = self.data['GlycanComposition'].str.replace(r'\+Ammonium\(\+17\)', '', regex=True)
            self.data['PeptideSequence+structure_coding+ProteinID'] = self.data['PeptideSequence'] + '+' + self.data['structure_coding'] + '+' + self.data['ProteinID']
            self.data['structure_coding'] = self.data['structure_coding'].replace(np.nan, np.nan)
            self.data['GlycanComposition'] = self.data['GlycanComposition'].str.replace(r'\+Ammonium\(\+17\)', '', regex=True)
            # branches
            if self.branch_list_dir is not None:
                branches = []
                for i in list(self.data['structure_coding']):
                    branch = []
                    if pd.notna(i):
                        for j in self.branch_list:
                            if j in i:
                                branch.append(j)
                    branches.append(str(branch))
                self.data['Branches'] = branches
            
            def parse_branches(branch_str):
                try:
                    return ast.literal_eval(branch_str)
                except (ValueError, SyntaxError):
                    return []
            self.data['Branches'] = self.data['Branches'].apply(parse_branches)
            def expand_branches(row):
                structure = row['structure_coding']
                branches = row['Branches']
                expanded_branches = []
                for branch in branches:
                    count = structure.count(branch)  
                    expanded_branches.extend([branch] * count)  
                return expanded_branches
            self.data['Branches'] = self.data.apply(expand_branches, axis=1)
            def format_list_as_string(lst):
                return str(lst)
            self.data['Branches'] = self.data['Branches'].apply(format_list_as_string)
            
            # glycan type
            glycantype = []
            for i, j, k, l in zip(list(self.data['GlycanComposition']), 
                                  list(self.data['structure_coding']), 
                                  list(self.data['Bisection']), 
                                  list(self.data['BranchNumber'])):
                if pd.notnull(i) and pd.notnull(j):
                    if 'N2' in i:
                        glycantype.append('Oligo mannose')
                    elif 'N3' in i:
                        if k == 0:
                            glycantype.append('Hybrid')
                        else:
                            glycantype.append('Oligo mannose')
                    elif 'N4' in i:
                        # if 'D1d' in j:
                        #     glycantype.append('Complex')
                        # else:
                        if k == 0:
                            if l == 1:
                                if 'D1d' in j:
                                    glycantype.append('Complex')
                                else:
                                    glycantype.append('Hybrid')
                            elif l == 2:
                                if 'E1' in j:
                                    glycantype.append('Hybrid')
                                else:
                                    glycantype.append('Complex')
                        else:
                            glycantype.append('Hybrid')
                    elif 'N5' in i:
                        if k == 0:
                            if 'E1' in j:
                                glycantype.append('Hybrid')
                            else:
                                glycantype.append('Complex')
                        else:
                            if l == 1:
                                glycantype.append('Hybrid')
                            elif l == 2:
                                if 'E1' in j:
                                    glycantype.append('Hybrid')
                                else:
                                    glycantype.append('Complex')
                    elif any(n in i for n in ['N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12', 'N13', 'N14', 'N15']):
                        if 'E1' in j:
                            glycantype.append('Hybrid')
                        else:
                            glycantype.append('Complex')
                else:
                    glycantype.append(np.nan)
                
            self.data['Glycan_type'] = glycantype
            
            # branch number
            # self.data['BranchNumber'] = self.data['structure_coding'].apply(lambda x: x.count('E'))
            self.data['BranchNumber'] = self.data['structure_coding'].apply(
                lambda x: np.nan if pd.isnull(x) else x.count('E')
            )
            
            # core structure
            temp_data = self.data
            core_structure_list = []
            for i, row in temp_data.iterrows():
                if (row['GlycanComposition'] == 'N2H2')|(row['GlycanComposition'] == 'N2H2F1')|('A2B2C1D2d' in row['structure_coding']):
                    core_structure_list.append(np.nan)
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD1')
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD1dcbB5')
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD2dD1')
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD2dD1dcbB5')
                else:
                    core_structure_list.append(np.nan)
            self.data['core_structure'] = core_structure_list
            
            # lacdinac
            def extract_key_strings(s):
                if pd.isnull(s):  
                    return ' '
                matches = re.findall(r'(E2F2.*?fe)', s)
                return ', '.join(matches) if matches else ' '
            self.data['lacdinac'] = self.data['structure_coding'].apply(extract_key_strings)
            
            # Ac Gc
            acgc = []
            for value in list(self.data['structure_coding']):
                if pd.notnull(value):
                    contains_3 = "3" in value
                    contains_4 = "4" in value
                    if contains_3 and contains_4:
                        acgc.append("dual")
                    elif contains_3:
                        acgc.append("Ac")
                    elif contains_4:
                        acgc.append("Gc")
                    else:
                        acgc.append(' ')
                else:
                    acgc.append(' ')
                
            self.data['Ac/Gc'] = acgc
            
            # core antenna fucosylated
            def fucosylated_type(row):
                if pd.isnull(row):
                    return ' '
                if '5' in row:  
                    if 'B5' in row:  
                        if any(x in row for x in ['E5', 'F5', 'G5', 'H5']):  
                            return 'dual' 
                        else:
                            return 'core fucosylated'  
                    else:
                        return 'antenna fucosylated'  
                else:
                    return ' ' 
            self.data['fucosylated type'] = self.data['structure_coding'].apply(fucosylated_type)
            
            # F S G
            fsg = []
            for value in list(self.data['GlycanComposition']):
                if 'F' in value and 'S' not in value and 'G' not in value:
                    fsg.append('F')
                elif 'F' not in value and 'S' in value and 'G' not in value:
                    fsg.append('S')
                elif 'F' not in value and 'S' not in value and 'G' in value:
                    fsg.append('G')
                elif 'F' in value and 'S' in value and 'G' not in value:
                    fsg.append('F + S')
                elif 'F' in value and 'S' not in value and 'G' in value:
                    fsg.append('F + G')
                elif 'F' not in value and 'S' in value and 'G' in value:
                    fsg.append('S + G')
                elif 'F' in value and 'S' in value and 'G' in value:
                    fsg.append('F + S + G')
                elif 'F' not in value and 'S' not in value and 'G' not in value:
                    fsg.append('Others')
            self.data['FSG'] = fsg
            
            if data_type == 'label free':
                quantnum = pd.DataFrame(self.data['PeptideSequence+structure_coding+ProteinID'].value_counts())
                quantnum.columns=['quantnum']
                self.data = pd.merge(self.data.set_index('PeptideSequence+structure_coding+ProteinID',drop=False), quantnum, left_index=True, right_index=True, how='left')
                self.data = self.data[~self.data['PeptideSequence+structure_coding+ProteinID'].duplicated()]

        elif self.search_engine == 'MSFragger-Glyco':
            composition_map = {
                'HexNAc': 'N',
                'Hex': 'H',
                'Fuc': 'F',
                'NeuAc': 'S',
                'NeuGc': 'G'
            }
            def extract_glycan_composition(mod_string):
                if not isinstance(mod_string, str):
                    return np.nan
                matches = re.findall(r'([A-Za-z:]+)\((\d+)\)', mod_string)
                counts = {letter: 0 for letter in composition_map.values()}
                for name, num in matches:
                    if name in composition_map:
                        counts[composition_map[name]] += int(num)
                    elif name not in composition_map:
                        return np.nan  
                order = ['N', 'H', 'F', 'S', 'G']
                result = ''.join(f"{k}{counts[k]}" for k in order if counts[k] > 0)
                return result if result else np.nan
            self.data['GlycanComposition'] = self.data['Observed Modifications'].apply(extract_glycan_composition)
            self.data = self.data.dropna(subset=['GlycanComposition'])
            self.data = self.data.rename(columns={'Peptide': 'PeptideSequence',
                                                  'Protein ID': 'ProteinID',
                                                  'Gene': 'GeneName'})
            self.data['PeptideSequence+structure_coding+ProteinID'] = self.data['Spectrum']
            
        elif self.search_engine == 'pGlyco3':
            self.data = self.data.rename(columns={'Peptide': 'PeptideSequence',
                                                  'ProSites': 'Glycosite_Position',
                                                  'Proteins': 'ProteinID',
                                                  'Genes': 'GeneName'})
            self.data['GeneName'] = self.data['GeneName'].replace(r'^;+$', 'null', regex=True)
            self.data['GeneName'] = self.data['GeneName'].replace(np.nan, 'null', regex=True)
            self.data['Glycosite_Position'] = self.data['Glycosite_Position'].astype(str)
            # protein id
            def extract_protein_ids(text):
                matches = re.findall(r'sp\|([^|]+)\|', text)
                return ';'.join(matches)
            self.data['ProteinID'] = self.data['ProteinID'].apply(extract_protein_ids)
            # glycan composition
            mapping = {'H': 'H', 'N': 'N', 'A': 'S', 'F': 'F', 'G': 'G'}
            order = ['N', 'H', 'F', 'S', 'G']
            def transform_composition(comp):
                parts = re.findall(r'([A-Za-z]+)\((\d+)\)', comp)
                if not parts:
                    return np.nan
                counts = {}
                for elem, num in parts:
                    if elem not in mapping or len(elem) != 1:
                        return np.nan
                    mapped_elem = mapping[elem]
                    counts[mapped_elem] = counts.get(mapped_elem, 0) + int(num)
                result = ''
                for key in order:
                    if key in counts:
                        result += f'{key}{counts[key]}'
                return result
            self.data['GlycanComposition'] = self.data['GlycanComposition'].apply(transform_composition)
            self.data = self.data.dropna(subset=['GlycanComposition'])
            self.data = self.data[~self.data['GlycanComposition'].str.startswith('N1', na=False)]
            # structure coding
            data = self.data['PlausibleStruct']
            data = data.values
            data = data.tolist()
            results = []
            for each in data:
                level = 0
                st_s = ""
                st = ord('A') - 1
                st_x = ord('a') - 1
                for x in each:
                    if x == '(':
                        st += 1
                        level += 1
                    elif x == ')':
                        st_s += chr(st_x + level)
                        level -= 1
                        st -= 1
                    else:
                        st_s += chr(st)
                        if x == 'N':
                            st_s += '2'
                        elif x == 'H':
                            st_s += '1'
                        elif x == 'F':
                            st_s += '5'
                        elif x == 'A':
                            st_s += '3'
                if "B5" in st_s:
                    st_s = st_s.replace("B5b", "")
                    index_a = st_s.rfind('a')
                    st_s = st_s[:index_a] + "B5b" + st_s[index_a:]
                if "D2d" in st_s:
                    d2d_index = st_s.find("D2d")
                    st_s = st_s[:d2d_index] + st_s[d2d_index + 3:]
                    next_d_index = st_s.find('d', d2d_index)
                    if next_d_index != -1:
                        st_s = st_s[:next_d_index + 1] + "D2d" + st_s[next_d_index + 1:]
                results.append(st_s)
            self.data['structure_coding'] = results
            # branches
            def get_branch(coding):
                branches = []
                start = 0
                def get_mannose_branch(s):
                    for char in s:
                        if char.isdigit() and char != '1':
                            return False
                    return True 
                while True:
                    e_start = coding.find("E", start)
                    if e_start == -1:
                        break
                    e_end = coding.find("e", e_start + 1)
                    if e_end == -1:
                        break
                    branch = coding[e_start: e_end + 1]  
                    branches.append(branch)
                    start = e_end + 1
                branches = [b for b in branches if not get_mannose_branch(b)]
                return branches
            self.data['Branches'] = self.data['structure_coding'].apply(get_branch)
            # branch number
            def get_branch_number(branches):
                return len(branches)
            self.data['BranchNumber'] = self.data['Branches'].apply(get_branch_number)
            # bisection
            def get_besic(coding):
                bisection=0
                if 'D2' in coding:
                    bisection=1
                return bisection
            self.data['Bisection'] = self.data['structure_coding'].apply(get_besic)
            # glycan type
            glycantype = []
            for i, j, k, l in zip(list(self.data['GlycanComposition']), 
                                  list(self.data['structure_coding']), 
                                  list(self.data['Bisection']), 
                                  list(self.data['BranchNumber'])):
                if pd.notnull(i) and pd.notnull(j):
                    if 'N2' in i:
                        glycantype.append('Oligo mannose')
                    elif 'N3' in i:
                        if k == 0:
                            glycantype.append('Hybrid')
                        else:
                            glycantype.append('Oligo mannose')
                    elif 'N4' in i:
                        if k == 0:
                            if l == 1:
                                if 'D1d' in j:
                                    glycantype.append('Complex')
                                else:
                                    glycantype.append('Hybrid')
                            elif l != 1:
                                if 'E1' in j:
                                    glycantype.append('Hybrid')
                                else:
                                    glycantype.append('Complex')
                        else:
                            glycantype.append('Hybrid')
                    elif 'N5' in i:
                        if k == 0:
                            if 'E1' in j:
                                glycantype.append('Hybrid')
                            else:
                                glycantype.append('Complex')
                        else:
                            if l == 1:
                                glycantype.append('Hybrid')
                            elif l != 1:
                                if 'E1' in j:
                                    glycantype.append('Hybrid')
                                else:
                                    glycantype.append('Complex')
                    elif any(n in i for n in ['N6', 'N7', 'N8', 'N9', 'N10', 'N11', 'N12', 'N13', 'N14', 'N15']):
                        if 'E1' in j:
                            glycantype.append('Hybrid')
                        else:
                            glycantype.append('Complex')
                else:
                    glycantype.append(np.nan)
            self.data['Glycan_type'] = glycantype
            # branch number
            self.data['BranchNumber'] = self.data['structure_coding'].apply(
                lambda x: np.nan if pd.isnull(x) else x.count('E')
            )
            # branches
            if self.branch_list_dir is not None:
                branches = []
                for i in list(self.data['structure_coding']):
                    branch = []
                    if pd.notna(i):
                        for j in self.branch_list:
                            if j in i:
                                branch.append(j)
                    branches.append(str(branch))
                self.data['Branches'] = branches
            def parse_branches(branch_str):
                try:
                    return ast.literal_eval(branch_str)
                except (ValueError, SyntaxError):
                    return []
            self.data['Branches'] = self.data['Branches'].apply(parse_branches)
            def expand_branches(row):
                structure = row['structure_coding']
                branches = row['Branches']
                expanded_branches = []
                for branch in branches:
                    count = structure.count(branch)  
                    expanded_branches.extend([branch] * count)  
                return expanded_branches
            self.data['Branches'] = self.data.apply(expand_branches, axis=1)
            def format_list_as_string(lst):
                return str(lst)
            self.data['Branches'] = self.data['Branches'].apply(format_list_as_string)
            # core structure
            temp_data = self.data
            core_structure_list = []
            for i, row in temp_data.iterrows():
                if (row['GlycanComposition'] == 'N2H2')|(row['GlycanComposition'] == 'N2H2F1')|('A2B2C1D2d' in row['structure_coding']):
                    core_structure_list.append(np.nan)
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD1')
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD1dcbB5')
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD2dD1')
                elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                    core_structure_list.append('A2B2C1D1dD2dD1dcbB5')
                else:
                    core_structure_list.append(np.nan)
            self.data['core_structure'] = core_structure_list
            # lacdinac
            def extract_key_strings(s):
                if pd.isnull(s):  
                    return ' '
                matches = re.findall(r'(E2F2.*?fe)', s)
                return ', '.join(matches) if matches else ' '
            self.data['lacdinac'] = self.data['structure_coding'].apply(extract_key_strings)
            # Ac Gc
            acgc = []
            for value in list(self.data['structure_coding']):
                if pd.notnull(value):
                    contains_3 = "3" in value
                    contains_4 = "4" in value
                    if contains_3 and contains_4:
                        acgc.append("dual")
                    elif contains_3:
                        acgc.append("Ac")
                    elif contains_4:
                        acgc.append("Gc")
                    else:
                        acgc.append(' ')
                else:
                    acgc.append(' ')
            self.data['Ac/Gc'] = acgc
            # core antenna fucosylated
            def fucosylated_type(row):
                if pd.isnull(row):
                    return ' '
                if '5' in row:  
                    if 'B5' in row:  
                        if any(x in row for x in ['E5', 'F5', 'G5', 'H5']):  
                            return 'dual' 
                        else:
                            return 'core fucosylated'  
                    else:
                        return 'antenna fucosylated'  
                else:
                    return ' ' 
            self.data['fucosylated type'] = self.data['structure_coding'].apply(fucosylated_type)
            # F S G
            fsg = []
            for value in list(self.data['GlycanComposition']):
                if 'F' in value and 'S' not in value and 'G' not in value:
                    fsg.append('F')
                elif 'F' not in value and 'S' in value and 'G' not in value:
                    fsg.append('S')
                elif 'F' not in value and 'S' not in value and 'G' in value:
                    fsg.append('G')
                elif 'F' in value and 'S' in value and 'G' not in value:
                    fsg.append('F + S')
                elif 'F' in value and 'S' not in value and 'G' in value:
                    fsg.append('F + G')
                elif 'F' not in value and 'S' in value and 'G' in value:
                    fsg.append('S + G')
                elif 'F' in value and 'S' in value and 'G' in value:
                    fsg.append('F + S + G')
                elif 'F' not in value and 'S' not in value and 'G' not in value:
                    fsg.append('Others')
            self.data['FSG'] = fsg
            # core fucose
            def core_fucose(coding):
                bisection=0
                if 'B5' in coding:
                    bisection=1
                return bisection
            self.data['Corefucose'] = self.data['structure_coding'].apply(core_fucose)
            self.data['PeptideSequence+structure_coding+ProteinID'] = self.data['PeptideSequence'] + '+' + self.data['structure_coding'] + '+' + self.data['ProteinID']
        
        elif self.search_engine == 'Glyco-Decipher':
            self.data = self.data.rename(columns={'Peptide': 'PeptideSequence',
                                                  'GlycoSite': 'Glycosite_Position',
                                                  'Protein': 'ProteinID'})
            # protein id genename
            def extract_info(text):
                entries = text.split(';')
                accessions = []
                genes = []
                for entry in entries:
                    acc_match = re.search(r'sp\|([^|]+)\|', entry)
                    if acc_match:
                        accession = acc_match.group(1)
                        accessions.append(accession)
                    else:
                        accession = None
                    gene_match = re.search(r'\|([^|_]+)_', entry)
                    if gene_match:
                        gene = gene_match.group(1)
                        gene = gene.capitalize()
                        genes.append(gene)
                    else:
                        genes.append('')
                return ';'.join(accessions), ';'.join(genes)
            self.data[['ProteinID', 'GeneName']] = self.data['ProteinID'].apply(lambda x: pd.Series(extract_info(x)))
            # glycosite
            self.data['Glycosite_Position'] = self.data['Glycosite_Position'].str.rstrip(';')
            self.data['Glycosite_Position'] = self.data['Glycosite_Position'].astype(str)
            # glycan composition
            mapping = {
                'HexNAc': 'N',
                'Hex': 'H',
                'Fuc': 'F',
                'NeuAc': 'S',
                'NeuGc': 'G'
            }
            order = ['N', 'H', 'F', 'S', 'G']
            def transform_glycan(comp):
                parts = re.findall(r'([A-Za-z]+)\((\d+)\)', comp)
                if not parts:
                    return np.nan
                counts = {}
                for elem, num in parts:
                    if elem not in mapping:
                        return np.nan
                    mapped = mapping[elem]
                    counts[mapped] = counts.get(mapped, 0) + int(num)
                result = ''
                for key in order:
                    if key in counts:
                        result += f'{key}{counts[key]}'
                return result
            self.data['GlycanComposition'] = self.data['GlycanComposition'].apply(transform_glycan)
            self.data['PeptideSequence+structure_coding+ProteinID'] = self.data['Title']

        self.data_manager.log_params('StrucGAP_Preprocess', 'data_cleaning', {})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_cleaned', self.data)
        
        return self
        
    def fdr(self, feature_type = None):
        """
        Supports four filtering levels—no control, peptide-level, glycan-level, and both—providing customizable stringency for data confidence.
        
        Parameters:
            feature_type: fdr control level from ['peptide','glycan','both' or 'no'].
        
        Returns:
            self.data_peptide_fdr_data (peptide-level fdr filtered data). 
            self.data_glycan_fdr_data (glycan-level fdr filtered data). 
            self.data_fdr_data (both-level fdr filtered data). 
            
        Return type:
            dataframe
        
        """
        if self.search_engine == 'StrucGP':
            if feature_type == None:
                feature_type = input("Please enter fdr level from ['peptide','glycan','both' or 'no']: ")
                expected_options = ['peptide','glycan','both', 'no']
                matches = get_close_matches(feature_type, expected_options, n=1, cutoff=0.5)
                if matches:
                    feature_type = matches[0]
                    print(f"Using '{feature_type}' as the input.")
                else:
                    print("No close match found. Using 'no' as the input.")
                    feature_type = 'no'
                
            if feature_type == 'peptide':
                self.data_fdr_filtered = self.data[(self.data['PeptideScore']>24.22)&(self.data['Status']=='NORMAL')]
                print('peptide')
                self.data_peptide_fdr_data = self.data_fdr_filtered
                self.data_glycan_fdr_data = self.data_fdr_filtered
            if feature_type == 'glycan':
                self.data_fdr_filtered = self.data[((self.data['PMZShift']==0)&(self.data['GlycanScore']>53.80))|
                     ((self.data['PMZShift']==1)&(self.data['GlycanScore']>128.88))|
                     ((self.data['PMZShift']==2)&(self.data['GlycanScore']>175.28))]
                print('glycan')
                self.data_peptide_fdr_data = self.data
                self.data_glycan_fdr_data = self.data_fdr_filtered
            if feature_type == 'both':
                if 'Status' not in self.data.columns:
                    self.data_peptide_fdr_data = self.data[((self.data['PeptideScore']>24.22))]
                    self.data_glycan_fdr_data = self.data_peptide_fdr_data[
                                          (((self.data_peptide_fdr_data['PMZShift']==0)&(self.data_peptide_fdr_data['GlycanScore']>53.80))|
                                           ((self.data_peptide_fdr_data['PMZShift']==1)&(self.data_peptide_fdr_data['GlycanScore']>128.88))|
                                           ((self.data_peptide_fdr_data['PMZShift']==2)&(self.data_peptide_fdr_data['GlycanScore']>175.28)))]
                    self.data_fdr_filtered = self.data_glycan_fdr_data
                    # self.data_fdr_filtered = self.data[((self.data['PeptideScore']>24.22))&
                    #                       (((self.data['PMZShift']==0)&(self.data['GlycanScore']>53.80))|
                    #                        ((self.data['PMZShift']==1)&(self.data['GlycanScore']>128.88))|
                    #                        ((self.data['PMZShift']==2)&(self.data['GlycanScore']>175.28)))]
                else:  
                    self.data_peptide_fdr_data = self.data[((self.data['PeptideScore']>24.22)&(self.data['Status']=='NORMAL'))]
                    self.data_glycan_fdr_data = self.data_peptide_fdr_data[
                                          (((self.data_peptide_fdr_data['PMZShift']==0)&(self.data_peptide_fdr_data['GlycanScore']>53.80))|
                                           ((self.data_peptide_fdr_data['PMZShift']==1)&(self.data_peptide_fdr_data['GlycanScore']>128.88))|
                                           ((self.data_peptide_fdr_data['PMZShift']==2)&(self.data_peptide_fdr_data['GlycanScore']>175.28)))]
                    self.data_fdr_filtered = self.data_glycan_fdr_data
                    # self.data_fdr_filtered = self.data[((self.data['PeptideScore']>24.22)&(self.data['Status']=='NORMAL'))&
                    #                       (((self.data['PMZShift']==0)&(self.data['GlycanScore']>53.80))|
                    #                        ((self.data['PMZShift']==1)&(self.data['GlycanScore']>128.88))|
                    #                        ((self.data['PMZShift']==2)&(self.data['GlycanScore']>175.28)))]
                print('both')
            if feature_type == 'no':
                self.data_fdr_filtered = self.data
                self.data_peptide_fdr_data = self.data
                self.data_glycan_fdr_data = self.data
                
            self.data_fdr_filtered = self.data_fdr_filtered.set_index('PeptideSequence+structure_coding+ProteinID',drop=False)

        elif self.search_engine != 'StrucGP':
            self.data_fdr_filtered = self.data
            self.data_peptide_fdr_data = self.data
            self.data_glycan_fdr_data = self.data
            self.data_fdr_filtered = self.data_fdr_filtered.set_index('PeptideSequence+structure_coding+ProteinID',drop=False)

        self.data_manager.log_params('StrucGAP_Preprocess', 'fdr', {'feature_type': feature_type})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_peptide_fdr_data', self.data_peptide_fdr_data)
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_glycan_fdr_data', self.data_glycan_fdr_data)
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_fdr_filtered', self.data_fdr_filtered)
        
        return self
    
    def outliers(self, abundance_ratio=None):
        """
        Corrects and normalizes TMT quantification data from matched reporter ions for each IGP.
        
        Parameters:
            abundance_ratio: normalized factors from global proteomics data.
        
        Returns:
            self.data_outliers_filtered (corrected and normalized TMT quantification data). 
            
        Return type:
            dataframe
        
        """
        if self.search_engine == 'StrucGP':
            Pep0 = pd.DataFrame(self.data_fdr_filtered['PeptideSequence+structure_coding+ProteinID'])
            Pep1 = pd.DataFrame(self.data_fdr_filtered['Matched_Reporter_Ions'])
            
            if abundance_ratio is None:
                abundance_ratio = input("Please enter the abundance ratio as a comma-separated list (e.g., 1, 1.019488, 1.740756, 1.554661, 2.674981,1.071297, 1.145732, 1.082529, 1.733719, 1.850238), when you done, please click 'Enter' in your keyboard:  ")
                abundance_ratio = [float(x) for x in abundance_ratio.split(',')]
                print('Your data has been successfully entered, please waiting ... ')
            self.abundance_ratio = abundance_ratio
            #    
            glycopep= []
            ions = []
            for i in Pep0.values.tolist():
                glycopep.append(i[0])
            z = 0
            #
            for i in Pep1.values.tolist(): 
                tmp0 = str(i[0]).split('),')  
                tmp1 = [glycopep[z]]
                for j in tmp0:
                    try:
                        tmp1.append(float(j.split(',')[-1]))
                    except:
                        tmp1.append(float(j.split(',')[-1][:-2]))
                z+=1
                tmp2 = [tmp1[0]]
                for j in range(len(tmp1[1:])):
                    try:
                        rat = abundance_ratio[j]
                        abundance = tmp1[j+1]
                        value = (abundance*rat)
                        tmp2.append(value)
                    except:
                        continue
                ions.append(tmp2)
            #
            header = ['PeptideSequence+structure_coding+ProteinID', *map(str, self.sample_group.index), 'psm']
            result = [header]
            result1 = [header]
            # 转 ions 为 DataFrame 便于分组处理
            ion_df = pd.DataFrame(ions)
            ion_df.rename(columns={0: 'peptide'}, inplace=True)
            # Melt所有的定量值（后续使用groupby）
            value_cols = ion_df.columns[1:]
            ion_melted = ion_df.melt(id_vars='peptide', value_vars=value_cols, var_name='channel', value_name='abundance')
            # Group by peptide → 列表化每个 channel 对应的 abundance 值
            grouped = ion_melted.groupby(['peptide', 'channel'])['abundance'].apply(list).unstack(fill_value=[]).reset_index()
            # 开始逐个处理 unique 的 glycopeptides
            for _, row in grouped.iterrows():
                pep_id = row['peptide']
                values_per_channel = row[1:].tolist()
                # 拆分前/后半部分（即：control/sample）
                half = len(values_per_channel) // 2
                data_c = pd.DataFrame(values_per_channel[:half]).replace(0, np.nan)
                data_s = pd.DataFrame(values_per_channel[half:]).replace(0, np.nan)
                # 归一化（使用 module1.median_cheng）
                for col in data_c.columns:
                    median_val = self.median_cheng(data_c[col].tolist())
                    if not np.isnan(median_val) and median_val != 0:
                        data_c[col] = data_c[col] / median_val
                        data_s[col] = data_s[col] / median_val
                # 汇总统计 → 中位数输出行
                tmp3 = [pep_id]
                # 对每一行（通道）进行归一化并取中位数
                for l in range(data_c.shape[0]):
                    row_values = data_c.loc[l].tolist()
                    tmp3.append(self.median_cheng(row_values))
                for l in range(data_s.shape[0]):
                    row_values = data_s.loc[l].tolist()
                    tmp3.append(self.median_cheng(row_values))
                tmp3.append(data_s.shape[1])  # psm count
                result1.append(tmp3)
            #
            self.data_outliers_filtered = pd.DataFrame(result1)
            self.data_outliers_filtered.columns = self.data_outliers_filtered.iloc[0]
            self.data_outliers_filtered = self.data_outliers_filtered.drop(self.data_outliers_filtered.index[0])
            self.data_outliers_filtered = self.data_outliers_filtered.set_index('PeptideSequence+structure_coding+ProteinID',drop=False)
            
            data_outliers_filtered_cleaned = self.data_outliers_filtered.dropna(axis=1, how='all')
            dropped_columns = self.data_outliers_filtered.columns.difference(data_outliers_filtered_cleaned.columns).tolist()
            dropped_columns = [float(x) for x in dropped_columns]
            self.sample_group = self.sample_group.drop(index=dropped_columns)
            self.data_outliers_filtered = data_outliers_filtered_cleaned
            
            # self.data_outliers_filtered = pd.concat([self.data_fdr_filtered, self.data_outliers_filtered],axis=1,join='inner')
            self.data_outliers_filtered = pd.merge(self.data_fdr_filtered, self.data_outliers_filtered, left_index=True, right_index=True, how='left')
            self.data_outliers_filtered = self.data_outliers_filtered.rename(columns={'PeptideSequence+structure_coding+ProteinID_x':'PeptideSequence+structure_coding+ProteinID',
                                                                                      'PeptideSequence+structure_coding+ProteinID_y':'PeptideSequence+structure_coding+ProteinID'})
            
            self.data_outliers_filtered = self.data_outliers_filtered[~self.data_outliers_filtered.index.duplicated()]

        elif self.search_engine != 'StrucGP':
            self.data_outliers_filtered = self.data_fdr_filtered
            self.data_outliers_filtered = self.data_outliers_filtered[~self.data_outliers_filtered.index.duplicated()]

        self.data_manager.log_params('StrucGAP_Preprocess', 'outliers', {'abundance_ratio': abundance_ratio})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_outliers_filtered', self.data_outliers_filtered)
        
        return self
    
    def cv(self, threshold = None):
        """
        Enables optional coefficient-of-variation filtering based on user-defined sample groupings.
        
        Parameters:
            threshold: cv filter threshold (e.g. 0.3).
        
        Returns:
            self.data_cv_filtered (cv filtered data). 
            
        Return type:
            dataframe
        
        """
        if self.search_engine == 'StrucGP':
            self.data_cv_filtered = copy.deepcopy(self.data_outliers_filtered)      
            self.data_cv_filtered['cv_control'] = ( self.data_cv_filtered[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[0]].index)]].std(axis=1, ddof=1) / self.data_cv_filtered[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[0]].index)]].mean(axis=1) )
            self.data_cv_filtered['cv_sample'] = ( self.data_cv_filtered[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[1]].index)]].std(axis=1, ddof=1) / self.data_cv_filtered[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[1]].index)]].mean(axis=1) )
            #
            if threshold == None:
                threshold = input("Please enter a threshold for CV filtering (e.g., 0.3), or 'no' to skip: ")
    
            if threshold != 'no':  
                try:
                    threshold = float(threshold)
                    self.data_cv_filtered = self.data_cv_filtered[
                        (self.data_cv_filtered['cv_control'] < threshold) & 
                        (self.data_cv_filtered['cv_sample'] < threshold)
                    ]
                except ValueError:
                    print("Invalid input. Skipping CV filtering.")
            else:
                print("Skipping CV filtering.")

        elif self.search_engine != 'StrucGP':
            self.data_cv_filtered = self.data_outliers_filtered.copy()

        self.data_manager.log_params('StrucGAP_Preprocess', 'cv', {'threshold': threshold})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_cv_filtered', self.data_cv_filtered)
        
        return self
    
    def psm(self, psm_number = None):
        """
        Filters IGPs by the minimum number of supporting PSMs.
        
        Parameters:
            psm_filter: psm filter threshold (e.g. 3).
        
        Returns:
            self.data_psm_filtered (psm filtered data). 
            
        Return type:
            dataframe
        
        """
        if self.search_engine == 'StrucGP':
            if psm_number is None:
                psm_number = input("Please enter a PSM number for filtering (e.g., 3), or 'no' to skip: ")
            if psm_number.lower() != 'no':
                try:
                    psm_number = int(psm_number)
                    self.data_psm_filtered = self.data_cv_filtered[self.data_cv_filtered['psm']>=psm_number]
                except ValueError:
                    print("Invalid input. Skipping PSM filtering.")
            else:
                print("Skipping PSM filtering.")
                self.data_psm_filtered = self.data_cv_filtered
            
            # imoputation
            sample_size = int(self.sample_group.shape[0] / 2)
            control_data = self.data_psm_filtered[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[0]].index)]]
            experiment_data = self.data_psm_filtered[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[1]].index)]]
            
            knn_imputer = KNNImputer(n_neighbors=sample_size)
            
            control_filled = knn_imputer.fit_transform(control_data)
            experiment_filled = knn_imputer.fit_transform(experiment_data)
            
            control_filled_df = pd.DataFrame(control_filled, columns=control_data.columns, index=control_data.index)
            experiment_filled_df = pd.DataFrame(experiment_filled, columns=experiment_data.columns, index=control_data.index)
            no_missing_value_data = pd.concat([control_filled_df, experiment_filled_df], axis=1)
                    
            self.data_psm_filtered.loc[:, self.sample_group.index.astype(str)] = no_missing_value_data

        elif self.search_engine != 'StrucGP':
            self.data_psm_filtered = self.data_cv_filtered.copy()
            psm_number = None

        self.data_manager.log_params('StrucGAP_Preprocess', 'psm', {'psm_number': psm_number})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_psm_filtered', self.data_psm_filtered)
        
        return self
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        
        with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_Preprocess.xlsx'), engine='xlsxwriter') as writer:
            if self.data is not None and not self.data.empty:
                self.data.to_excel(writer, sheet_name='cleaned_data')
            if self.data_peptide_fdr_data is not None and not self.data_peptide_fdr_data.empty:
                self.data_peptide_fdr_data.to_excel(writer, sheet_name='data_peptide_fdr_data')
            if self.data_glycan_fdr_data is not None and not self.data_glycan_fdr_data.empty:
                self.data_glycan_fdr_data.to_excel(writer, sheet_name='data_glycan_fdr_data')
            if self.data_fdr_filtered is not None and not self.data_fdr_filtered.empty:
                self.data_fdr_filtered.to_excel(writer, sheet_name='data_fdr_filtered')
            if self.data_outliers_filtered is not None and not self.data_outliers_filtered.empty:
                self.data_outliers_filtered.to_excel(writer, sheet_name='data_outliers_filtered')
            if self.data_cv_filtered is not None and not self.data_cv_filtered.empty:
                self.data_cv_filtered.to_excel(writer, sheet_name='data_cv_filtered')
            if self.data_psm_filtered is not None and not self.data_psm_filtered.empty:
                self.data_psm_filtered.to_excel(writer, sheet_name='data_psm_filtered')
            


