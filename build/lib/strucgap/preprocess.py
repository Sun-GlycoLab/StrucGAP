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
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from itertools import combinations 
from statsmodels.stats.multitest import multipletests
import urllib.parse
import time
from urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
from collections import OrderedDict
from typing import Any, Optional
from collections import defaultdict
from collections import Counter
from collections import defaultdict, deque
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
## 数据质控模块--7
class StrucGAP_Preprocess:
    """
    Parameters:
        data_manager: Data manager instance, such as 'data_manager' if data_manager = StrucGAP_InsightTracker().
        search_engine: Support multiple search_engine in ['StrucGP','MSFragger-Glyco','pGlyco3','Glyco-Decipher','Byonic','GlycanFinder'].  
        data_dir: Path to search results data files.
        sample_group_data_dir: Path to sample grouping information.
        data_sheet_name: Sheet name of input Excel file.
        branch_list_dir: Path to branch list file.
    
    Both file templates were provided in Github repository (https://github.com/Sun-GlycoLab/StrucGAP/tree/main/tests).
    
    """
    def __init__(self, data_manager, search_engine = 'StrucGP',
                 data_dir=None, sample_group_data_dir=None, data_sheet_name=None, branch_list_dir=None):
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
        if search_engine not in ['StrucGP','MSFragger-Glyco','pGlyco3','Glyco-Decipher','Byonic','GlycanFinder']:
            print("Select search_engine in ['StrucGP','MSFragger-Glyco','pGlyco3','Glyco-Decipher','Byonic','GlycanFinder']")
            search_engine = 'StrucGP'
            if search_engine is not 'StrucGP':
                print('You can only use a samll subset of the functions in StrucGAP_GlycanStructure and StrucGAP_GlycoSite module!')
        self.search_engine = search_engine
        #######################################################################
        self.data_manager = data_manager
        self.data_manager.register_module('StrucGAP_Preprocess', self, {'input_data': data_dir, 'sample_group': sample_group_data_dir, 'search_engine': search_engine})
        self.data_manager.log_params('StrucGAP_Preprocess', 'input_data', {'input_data': data_dir, 'sample_group': sample_group_data_dir, 'search_engine': search_engine})
        
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
        
    def fc_recommendation(self, data):
        """
        StrucGAP recommends a data-driven fold change (FC) thresholding pipeline that integrates (i) CV-based theoretical minimal detectable fold change (MDFC), (ii) empirical null distributions from control–control comparisons FC, and (iii) empirical P values derived from these nulls. This strategy defines FC cut-offs that are specific to dataset variability, ensuring that at least 95% of control–control comparisons are non-significant, and combines them with empirical significance testing to balance effect-size and statistical rigor.

        """
        # Optional: BH-FDR from statsmodels if available; otherwise implement a tiny BH helper
        def bh_fdr(pvals):
            _, q, _, _ = multipletests(pvals, method="fdr_bh")
            return q
        # -----------------------------
        # Utility functions
        # -----------------------------
        def cv(arr):
            """Coefficient of variation on linear scale (sample SD / mean)."""
            a = np.asarray(arr, dtype=float)
            a = a[~np.isnan(a)]
            if a.size == 0:
                return np.nan
            m = np.mean(a)
            if m == 0:
                return np.nan
            return np.std(a, ddof=1) / m
        def sigma_ln_from_cv(cv_val):
            """Map linear CV -> natural-log SD using sigma_ln^2 = ln(1 + CV^2)."""
            if np.isnan(cv_val):
                return np.nan
            return np.sqrt(np.log(1.0 + cv_val**2))
        def mdfc_from_cv(cv_val, n_reps, z=1.96):
            """
            Theoretical minimal detectable fold-change (MDFC) at confidence z,
            derived from CV. Returns a dict with:
              - group_mean_log: Δ_ln^(n) for group means
              - group_mean_fc: MDFC^(n) in linear space
              - single_pair_log: Δ_ln^(single) for single-channel vs single-channel
              - single_pair_fc: MDFC^(single) in linear space
            """
            sig_ln = sigma_ln_from_cv(cv_val)
            if np.isnan(sig_ln):
                return dict(group_mean_log=np.nan, group_mean_fc=np.nan,
                            single_pair_log=np.nan, single_pair_fc=np.nan)
            # two groups, each n_reps; variance of group means difference in ln-space
            delta_ln_group = z * np.sqrt(2 * (sig_ln**2) / n_reps)
            # single-channel (pairwise) conservative bound
            delta_ln_single = z * np.sqrt(2) * sig_ln
            return dict(
                group_mean_log=delta_ln_group,
                group_mean_fc=float(np.exp(delta_ln_group)),
                single_pair_log=delta_ln_single,
                single_pair_fc=float(np.exp(delta_ln_single)),
            )
        def empirical_null_ctrl_vs_ctrl(df_ctrl_log2):
            """
            Build empirical null of |log2 differences| by enumerating all control-channel pairs
            for all features. Returns the pooled null distribution as a 1D numpy array.
            """
            k = df_ctrl_log2.shape[1]
            diffs = []
            for i, j in combinations(range(k), 2):
                d = np.abs(df_ctrl_log2.iloc[:, i] - df_ctrl_log2.iloc[:, j])
                diffs.append(d.values)
            return np.concatenate(diffs) if diffs else np.array([])
        def empirical_p_from_null(observed_abs_log2fc, null_abs_log2_diffs):
            """
            Empirical right-tail p-value for each observed |log2FC| against the null distribution.
            p = (count(null >= obs) + 1) / (len(null) + 1)
            """
            null_sorted = np.sort(null_abs_log2_diffs)
            N = null_sorted.size
            if N == 0:
                return np.full_like(observed_abs_log2fc, np.nan, dtype=float)
            # vectorized: for each obs, count how many null >= obs
            idx = np.searchsorted(null_sorted, observed_abs_log2fc, side="left")
            ge_counts = N - idx
            p = (ge_counts + 1.0) / (N + 1.0)
            return p
        ####
        cols = data.columns.tolist()
        ctrl_cols = cols[:len(cols)//2]
        exp_cols = cols[len(cols)//2:]
        # ctrl_cols = ['126.1277', '127.1248', '127.1311', '128.1281', '128.1344']
        # exp_cols  = ['129.1315', '129.1378', '130.1348', '130.1411', '131.1382']
        ####
        # -----------------------------
        # A) CV-based theoretical MDFC（自然对数）
        # -----------------------------
        data["CV_ctrl_raw"] = data[ctrl_cols].apply(lambda r: cv(r.values), axis=1)
        cv_median = float(np.nanmedian(data["CV_ctrl_raw"].values))
        mdfc_stats = mdfc_from_cv(cv_median, n_reps=5, z=1.96)  # 95%置信
        # 方便展示
        mdfc_table = pd.DataFrame({
            "Metric": ["Median CV (controls)",
                       "Δ_ln (group means, n=5)",
                       "MDFC (group means, linear)",
                       "Δ_ln (single pair)",
                       "MDFC (single pair, linear)"],
            "Value": [cv_median,
                      mdfc_stats["group_mean_log"],
                      mdfc_stats["group_mean_fc"],
                      mdfc_stats["single_pair_log"],
                      mdfc_stats["single_pair_fc"]]
        })
        ## fc
        # -----------------------------
        # B1) 经验空分布（Ctrl–Ctrl）→ 单一 FC cut-off
        # -----------------------------
        log2_ctrl = np.log2(data[ctrl_cols])
        null_abs_log2 = empirical_null_ctrl_vs_ctrl(log2_ctrl)
        tau_log2_95 = float(np.quantile(null_abs_log2[~np.isnan(null_abs_log2)], 0.95))
        tau_fc_95 = float(2**tau_log2_95)
        #
        tau_log2_99 = float(np.quantile(null_abs_log2[~np.isnan(null_abs_log2)], 0.99))
        tau_fc_99 = float(2**tau_log2_99)
        # -----------------------------
        # C) 5v5 真实差异的经验 p 值 + BH-FDR + 复合判定
        # -----------------------------
        log2_mean_ctrl = log2_ctrl.mean(axis=1)
        log2_mean_exp = np.log2(data[exp_cols]).mean(axis=1)
        data["log2FC"] = log2_mean_exp - log2_mean_ctrl
        data["abs_log2FC"] = np.abs(data["log2FC"])
        # 经验 p 值 & FDR
        p_emp = empirical_p_from_null(data["abs_log2FC"].values, null_abs_log2)
        q_emp = bh_fdr(p_emp)
        data["p_empirical"] = p_emp
        # data["q_empirical"] = q_emp
        data["emp_cutoff_log2_95"] = tau_log2_95
        data["emp_cutoff_FC_95"]   = tau_fc_95
        #
        data["emp_cutoff_log2_99"] = tau_log2_99
        data["emp_cutoff_FC_99"]   = tau_fc_99
        data["pass_FC_cutoff"]  = data["abs_log2FC"] >= tau_log2_95
        # data["pass_FDR_0.05"]   = data["q_empirical"] < 0.05
        # data["significant"]     = data["pass_FC_cutoff"] & data["pass_FDR_0.05"]
        # -----------------------------
        # 展示结果
        # -----------------------------
        # from ace_tools import display_dataframe_to_user
        # summary = data[ctrl_cols + exp_cols + [
        #     "CV_ctrl_raw","log2FC","abs_log2FC","emp_cutoff_log2","emp_cutoff_FC",
        #     "p_empirical","q_empirical","pass_FC_cutoff","pass_FDR_0.05","significant"
        # ]].round(5)
        summary = data[ctrl_cols + exp_cols + [
            "log2FC","abs_log2FC","emp_cutoff_log2_95","emp_cutoff_FC_95","emp_cutoff_log2_99","emp_cutoff_FC_99",
            "p_empirical","pass_FC_cutoff",
        ]].round(5)
        #        
        return mdfc_table, summary
    
    def data_cleaning(self, quantification_from_no_strucgp = False, data_type=None, 
                      quantification_data_dir = None, sheet_name = None, quant_cols = None):
        """
        Data cleaning and glycan substructure features extraction.
        
        Parameters:
            data_type: output data of StrucGP with TMT or label-free type. Select data_type in ['tmt', 'label free'].
            quantification_from_no_strucgp: whether execute quantification for other search engine input (except StrucGP).
            quantification_data_dir: quantification data path.
            sheet_name: sheet name of quantification data excel.
            quant_cols: quantification values columns (e.g. ["Young-SN-1","Young-SN-2","Young-SN-3","PD-SN-1","PD-SN-2","PD-SN-3"]).
        
        Returns:
            self.data (cleaned data). 
            
        Return type:
            dataframe
        
        """
        self.quantification_from_no_strucgp = quantification_from_no_strucgp
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
            #    
            self.data['Glycan_type'] = glycantype
            #
            # branch number
            # self.data['BranchNumber'] = self.data['structure_coding'].apply(lambda x: x.count('E'))
            self.data['BranchNumber'] = self.data['structure_coding'].apply(
                lambda x: np.nan if pd.isnull(x) else x.count('E')
            )
            #
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
            #
            # lacdinac
            def extract_key_strings(s):
                if pd.isnull(s):  
                    return ' '
                matches = re.findall(r'(E2F2.*?fe)', s)
                return ', '.join(matches) if matches else ' '
            self.data['lacdinac'] = self.data['structure_coding'].apply(extract_key_strings)
            #
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
            #
            self.data['Ac/Gc'] = acgc
            #
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
            #
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
            #
            if data_type == 'label free':
                quantnum = pd.DataFrame(self.data['PeptideSequence+structure_coding+ProteinID'].value_counts())
                quantnum.columns=['quantnum']
                self.data = pd.merge(self.data.set_index('PeptideSequence+structure_coding+ProteinID',drop=False), quantnum, left_index=True, right_index=True, how='left')
                self.data = self.data[~self.data['PeptideSequence+structure_coding+ProteinID'].duplicated()]
        #
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
        #
        elif self.search_engine == 'pGlyco3':
            # quant
            if quantification_from_no_strucgp:
                quant_data = pd.read_excel(quantification_data_dir, sheet_name = sheet_name)
                # 你要使用的定量列（按你数据中的实际列名顺序填写）
                quant_cols = quant_cols
                # 若是 n=10，就把 10 列名按顺序填满 quant_cols
                channels_10 = ["126.1277","127.1248","127.1311","128.1281","128.1344",
                               "129.1315","129.1378","130.1348","130.1411","131.1382"]
                channels_6 = ["126.1277","127.1311","128.1344","129.1378","130.1411","131.1382"]
                channel_meta = {
                    "126.1277": (126.127725, 126.12754821777344),
                    "127.1248": (127.12476,  0.0),
                    "127.1311": (127.131079, 127.1312484741211),
                    "128.1281": (128.128114, 0.0),
                    "128.1344": (128.134433, 128.1345977783203),
                    "129.1315": (129.131468, 0.0),
                    "129.1378": (129.137787, 129.13775634765625),
                    "130.1348": (130.134822, 0.0),
                    "130.1411": (130.141141, 130.1413116455078),
                    "131.1382": (131.138176, 131.1381072998047),
                }
                def build_matched_dict(row, qcols, use_channels):
                    # 准备：已使用通道 -> 强度
                    used = {ch: float(row[qcols[i]]) for i, ch in enumerate(use_channels)}
                    # 按 channels_10 的固定顺序构造
                    od = OrderedDict()
                    for ch in channels_10:
                        theo, obs = channel_meta[ch]
                        intensity = used.get(ch, 0.0)   # 未使用的通道置 0
                        # 对未使用通道，obs 也按规则置 0.0
                        od[ch] = (theo, obs if ch in used else 0.0, intensity if ch in used else 0.0)
                    return str(dict(od))
                # ===== 核心逻辑 =====
                n = len(quant_cols)
                if n == 6:
                    quant_data = quant_data.rename(columns=dict(zip(quant_cols, channels_6)))
                    use_channels = channels_6
                elif n == 10:
                    quant_data = quant_data.rename(columns=dict(zip(quant_cols, channels_10)))
                    use_channels = channels_10
                else:
                    raise ValueError("只支持 6 或 10 个定量列")
                # 丢掉 use_channels 中任一列为 NaN 或 0 的行
                quant_data = quant_data.dropna(subset=use_channels)
                # 删掉所有 use_channels 列都为 0 的行
                quant_data = quant_data[~(quant_data[use_channels] == 0).all(axis=1)]
                # 把 0 当缺失
                quant_data[use_channels] = quant_data[use_channels].replace(0, np.nan)
                # 从 log2 还原到线性空间
                quant_data[use_channels] = np.power(2.0, quant_data[use_channels])
                # 每列中位数归一化（除以各列中位数）
                col_medians = quant_data[use_channels].median(axis=0, skipna=True)
                quant_data[use_channels] = quant_data[use_channels].divide(col_medians, axis=1)
                # （可选）乘回全局中位数，保持整体量级
                global_median = np.nanmedian(quant_data[use_channels].values)
                quant_data[use_channels] = quant_data[use_channels] * global_median
                # 缺失值填补：用 1% 分位数的一半作为填充值
                baseline = np.nanpercentile(quant_data[use_channels].values, 1)
                fill_value = baseline * 0.5
                quant_data[use_channels] = quant_data[use_channels].fillna(fill_value)
                # 生成字符串列（顺序严格按 channels_10）
                quant_data["Matched_Reporter_Ions"] = quant_data.apply(lambda r: build_matched_dict(r, use_channels, use_channels), axis=1)
                # combine quant and struc
                # -------- 参数（可按需调整） --------
                DEFAULT_PPM_TOL = 10.0
                FDR_MAX = 0.01
                ALLOW_FDR_RELAX = False
                FDR_RELAX_MAX = 0.05
                REQUIRE_NON_DECOY = True
                USE_RT_TIEBREAKER = True
                # -------- 规范化工具 --------
                def normalize_multivalue_field(val: Any) -> str:
                    if pd.isna(val):
                        return ''
                    s = str(val)
                    parts = re.split(r'[;,|\s]+', s)
                    parts = [p.strip() for p in parts if p.strip()]
                    parts = sorted(set(parts))
                    return ';'.join(parts)
                def normalize_mod_field(val: Any) -> str:
                    if pd.isna(val):
                        return ''
                    s = str(val).strip().lower()
                    # 可在此扩展修饰命名映射
                    s = s.replace('oxidation(m)', 'm(ox)')
                    tokens = re.split(r'[;,|\s]+', s)
                    tokens = [t for t in tokens if t]
                    tokens = sorted(set(tokens))
                    return ';'.join(tokens)
                def _counts_to_key(counts: dict) -> str:
                    """固定顺序输出 H/N/A/G/F"""
                    return f"H{counts['H']}N{counts['N']}A{counts['A']}G{counts['G']}F{counts['F']}"
                def _parse_tuple_to_comp(val: Any) -> str | None:
                    """
                    解析 'Glycan(H,N,A,G,F)' 类列，如 '4,2,0,0,0' 或 '(4,2,0,0,0)'，返回 H#N#A#G#F# 字符串。
                    """
                    if pd.isna(val):
                        return None
                    s = str(val)
                    nums = re.findall(r'-?\d+', s)
                    if len(nums) >= 5:
                        h, n, a, g, f = map(int, nums[:5])
                        return f"H{h}N{n}A{a}G{g}F{f}"
                    return None
                def normalize_glycan_composition(val: Any) -> str:
                    """
                    解析多种写法并规约为固定顺序 'H#N#A#G#F#'：
                    - H(4)N(2)、H4N2、H:4;N:2、H=4 N=2
                    - HEX(4)HEXNAC(2)、NEUAC/NEUGC/FUC/SIA 同义词
                    - 支持全角括号（（ ））
                    """
                    counts = {'H': 0, 'N': 0, 'A': 0, 'G': 0, 'F': 0}
                    if pd.isna(val):
                        return _counts_to_key(counts)
                    s = str(val).strip()
                    if not s:
                        return _counts_to_key(counts)
                    # 统一大小写/括号/空白
                    s = s.upper()
                    s = s.replace('（', '(').replace('）', ')')
                    s = s.replace(' ', '')
                    # 同义词映射（先长词后短词，避免 HEXNAC 被先替换成 H）
                    synonyms = [
                        ('HEXNAC', 'N'),
                        ('HEX', 'H'),
                        ('NEUAC', 'A'),
                        ('NEUGC', 'G'),
                        ('FUC', 'F'),
                        ('SIA', 'A'),   # 若你的数据里 SIA 意味 NeuAc
                    ]
                    for old, new in synonyms:
                        s = re.sub(old, new, s)
                    # 匹配模式：
                    # 1) L(123)  2) L=123 / L:123 / L-123  3) L123
                    pattern = re.compile(r'([HNAGF])(?:\((\d+)\)|[:=\-]?(\d+))')
                    for m in pattern.finditer(s):
                        L = m.group(1)
                        num = m.group(2) or m.group(3)
                        if num is not None:
                            counts[L] += int(num)
                    return _counts_to_key(counts)
                def to_int(val: Any):
                    try:
                        if pd.isna(val):
                            return None
                        return int(val)
                    except Exception:
                        try:
                            return int(float(val))
                        except Exception:
                            return None
                # -------- 列名标准化 --------
                def standardize_columns(df: pd.DataFrame) -> pd.DataFrame:
                    """
                    常见同义列名标准化：
                    - Protein -> Proteins
                    - ProSites/ProSite -> Prosite
                    - Peptide -> peptide
                    - Modification -> Mod
                    - 若缺少 GlycanComposition 而有 CorrectedComposition，则回退
                    """
                    col_map = {}
                    cols = set(df.columns)
                    if 'Proteins' not in cols and 'Protein' in cols:
                        col_map['Protein'] = 'Proteins'
                    if 'Prosite' not in cols:
                        if 'ProSites' in cols:
                            col_map['ProSites'] = 'Prosite'
                        elif 'ProSite' in cols:
                            col_map['ProSite'] = 'Prosite'
                    if 'peptide' not in cols and 'Peptide' in cols:
                        col_map['Peptide'] = 'peptide'
                    if 'Mod' not in cols and 'Modification' in cols:
                        col_map['Modification'] = 'Mod'
                    if 'GlycanComposition' not in cols and 'CorrectedComposition' in cols:
                        col_map['CorrectedComposition'] = 'GlycanComposition'
                    if col_map:
                        df = df.rename(columns=col_map)
                    return df
                # -------- 基于母离子的聚类 --------
                def cluster_by_precursor_mh(masses: np.ndarray, ppm_tol: float) -> np.ndarray:
                    if len(masses) == 0:
                        return np.array([], dtype=int)
                    order = np.argsort(masses)
                    clusters = np.zeros(len(masses), dtype=int)
                    cid = 0
                    clusters[order[0]] = cid
                    for i in range(1, len(order)):
                        prev_idx = order[i-1]
                        cur_idx = order[i]
                        prev_m = masses[prev_idx]
                        cur_m = masses[cur_idx]
                        ref_m = 0.5*(prev_m + cur_m)
                        tol_da = ref_m * ppm_tol * 1e-6
                        if abs(cur_m - prev_m) <= tol_da:
                            clusters[cur_idx] = cid
                        else:
                            cid += 1
                            clusters[cur_idx] = cid
                    return clusters
                # -------- 候选排序 --------
                def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
                    def safe_col(col, default):
                        return df[col] if col in df.columns else default
                    core = safe_col('CoreMatched', False)
                    if core.dtype != bool:
                        core = core.astype(str).str.lower().isin(['true','1','yes','y'])
                    core_sort = (~core).astype(int)  # True 优先（True->0, False->1）
                    abs_ppm = df['PPM'].abs().fillna(np.inf) if 'PPM' in df.columns else pd.Series(np.inf, index=df.index)
                    if USE_RT_TIEBREAKER and 'RT' in df.columns:
                        med_rt = df['RT'].median(skipna=True)
                        rt_dev = (df['RT'] - med_rt).abs()
                        rt_dev = rt_dev.fillna(rt_dev.max() if rt_dev.notna().any() else np.inf)
                    else:
                        rt_dev = pd.Series(np.inf, index=df.index)
                    glycan_fdr = df['GlycanFDR'] if 'GlycanFDR' in df.columns else pd.Series(np.inf, index=df.index)
                    glyscore   = df['GlyScore']  if 'GlyScore'  in df.columns else pd.Series(-np.inf, index=df.index)
                    totalscore = df['TotalScore'] if 'TotalScore' in df.columns else pd.Series(-np.inf, index=df.index)
                    sort_keys = []
                    sort_keys.append(glycan_fdr.fillna(np.inf))      # 升
                    sort_keys.append(-glyscore.fillna(-np.inf))      # 降
                    sort_keys.append(core_sort.values)               # True 优先
                    sort_keys.append(abs_ppm.values)                 # 升
                    sort_keys.append(-totalscore.fillna(-np.inf))    # 降
                    sort_keys.append(rt_dev.values)                  # 升
                    order = np.lexsort(sort_keys[::-1])
                    return df.iloc[order]
                # -------- 主流程 --------
                def assign_structures(quant_df: pd.DataFrame, qual_df: pd.DataFrame,
                                      ppm_tol: float = DEFAULT_PPM_TOL,
                                      fdr_max: float = FDR_MAX,
                                      allow_relax: bool = ALLOW_FDR_RELAX,
                                      fdr_relax_max: float = FDR_RELAX_MAX,
                                      require_non_decoy: bool = REQUIRE_NON_DECOY) -> pd.DataFrame:
                    qn = standardize_columns(quant_df.copy())
                    qa = standardize_columns(qual_df.copy())
                    required_keys = ['Proteins', 'Prosite', 'peptide', 'GlycanComposition', 'Charge', 'Mod']
                    for col in required_keys:
                        if col not in qn.columns:
                            raise ValueError(f"Quant table missing required key column: {col}")
                        if col not in qa.columns:
                            raise ValueError(f"Qual table missing required key column: {col}")
                    # --- 生成规范匹配键（加入 GlycanComposition 解析失败时的回退策略） ---
                    def row_comp_key(row: pd.Series) -> str:
                        # 1) 主解析：GlycanComposition
                        primary = normalize_glycan_composition(row.get('GlycanComposition', ''))
                        # 2) 若解析为全 0，尝试几种回退
                        if primary == 'H0N0A0G0F0':
                            # 2a) CorrectedComposition
                            cc = row.get('CorrectedComposition', None)
                            if cc is not None and str(cc).strip():
                                alt = normalize_glycan_composition(cc)
                                if alt != 'H0N0A0G0F0':
                                    return alt
                            # 2b) Glycan(H,N,A,G,F)
                            gh = row.get('Glycan(H,N,A,G,F)', None)
                            if gh is not None and str(gh).strip():
                                alt = _parse_tuple_to_comp(gh)
                                if alt:
                                    return alt
                            # 2c) CorrectedGlycan(H,N,A,G,F)
                            cgh = row.get('CorrectedGlycan(H,N,A,G,F)', None)
                            if cgh is not None and str(cgh).strip():
                                alt = _parse_tuple_to_comp(cgh)
                                if alt:
                                    return alt
                        return primary
                    def canonize(df: pd.DataFrame) -> pd.DataFrame:
                        df = df.copy()
                        df['__key_proteins'] = df['Proteins'].apply(normalize_multivalue_field)
                        df['__key_prosite']  = df['Prosite'].apply(normalize_multivalue_field)
                        df['__key_peptide']  = df['peptide'].astype(str)
                        df['__key_mod']      = df['Mod'].apply(normalize_mod_field)
                        # 组合键来自行级函数（带回退）
                        df['__key_gcomp']    = df.apply(row_comp_key, axis=1)
                        df['__key_charge']   = df['Charge'].apply(to_int)
                        df['__key'] = (
                            df['__key_proteins'] + '|' +
                            df['__key_prosite']  + '|' +
                            df['__key_peptide']  + '|' +
                            df['__key_gcomp']    + '|' +
                            df['__key_mod']      + '|' +
                            df['__key_charge'].astype(str)
                        )
                        df['__row_id'] = np.arange(len(df))
                        return df
                    qn = canonize(qn)
                    qa = canonize(qa)
                    # FDR & 诱饵过滤
                    def fdr_filter(df: pd.DataFrame, fdr: float) -> pd.Series:
                        ok = pd.Series(True, index=df.index)
                        for col in ['GlycanFDR', 'PeptideFDR', 'TotalFDR']:
                            if col in df.columns:
                                ok &= (df[col] <= fdr)
                        return ok
                    ok_main = fdr_filter(qa, fdr_max)
                    if require_non_decoy:
                        for decoy_col in ['GlycanDecoy', 'PepDecoy']:
                            if decoy_col in qa.columns:
                                qa[decoy_col] = qa[decoy_col].astype(str).str.lower().isin(['false','0','no','n','f'])
                        for decoy_col in ['GlycanDecoy', 'PepDecoy']:
                            if decoy_col in qa.columns:
                                ok_main &= qa[decoy_col]
                    qa_main = qa[ok_main].copy()
                    if len(qa_main) == 0 and allow_relax:
                        ok_relax = fdr_filter(qa, fdr_relax_max)
                        if require_non_decoy:
                            for decoy_col in ['GlycanDecoy', 'PepDecoy']:
                                if decoy_col in qa.columns:
                                    ok_relax &= qa[decoy_col]
                        qa_main = qa[ok_relax].copy()
                    qa_groups = qa_main.groupby('__key')
                    # —— 输出列（含 AssignedGenes / AssignedGlySpec） ——
                    qn['PlausibleStruct']   = pd.Series([pd.NA]*len(qn), dtype="object")
                    qn['AssignedGlyID']             = pd.Series([pd.NA]*len(qn), dtype="object")
                    qn['AssignedCorrectedComposition'] = pd.Series([pd.NA]*len(qn), dtype="object")
                    qn['ChosenGlycanFDR']           = np.nan
                    qn['ChosenGlyScore']            = np.nan
                    qn['ChosenTotalScore']          = np.nan
                    qn['ChosenPPM']                 = np.nan
                    qn['ChosenRT']                  = np.nan
                    qn['QualRowIndex']              = pd.Series([-1]*len(qn), dtype="Int64")
                    qn['CandidateCount']            = pd.Series([0]*len(qn), dtype="Int64")
                    qn['ClusterSize']               = pd.Series([0]*len(qn), dtype="Int64")
                    qn['IsomerAmbiguous']           = pd.Series([False]*len(qn), dtype="boolean")
                    qn['AssignReason']              = pd.Series([""]*len(qn), dtype="string")
                    qn['Genes']                     = pd.Series([pd.NA]*len(qn), dtype="object")   # <---
                    qn['GlySpec']                   = pd.Series([pd.NA]*len(qn), dtype="object")   # <---
                    # 逐行指派
                    for idx, row in qn.iterrows():
                        k = row['__key']
                        if k not in qa_groups.groups:
                            qn.at[idx, 'AssignReason'] = 'NoQualCandidate'
                            continue
                        cand = qa_main.loc[qa_groups.groups[k]].copy()
                        qn.at[idx, 'CandidateCount'] = int(len(cand))
                        # 母离子聚类（定性候选内部）
                        if 'PrecursorMH' not in cand.columns:
                            cand['__cluster'] = 0
                        else:
                            masses = cand['PrecursorMH'].astype(float).values
                            clusters = cluster_by_precursor_mh(masses, ppm_tol=ppm_tol)
                            cand['__cluster'] = clusters
                        # 选最大簇；并列时看 |PPM| 中位数
                        cluster_sizes = cand.groupby('__cluster').size()
                        top_size = cluster_sizes.max()
                        top_clusters = cluster_sizes[cluster_sizes == top_size].index.tolist()
                        if len(top_clusters) > 1 and 'PPM' in cand.columns:
                            best_c = None
                            best_med = np.inf
                            for c_id in top_clusters:
                                med = cand.loc[cand['__cluster']==c_id, 'PPM'].abs().median()
                                if med < best_med:
                                    best_med = med
                                    best_c = c_id
                            chosen_cluster = best_c if best_c is not None else top_clusters[0]
                        else:
                            chosen_cluster = top_clusters[0]
                        cand_cluster = cand[cand['__cluster'] == chosen_cluster].copy()
                        qn.at[idx, 'ClusterSize'] = int(len(cand_cluster))
                        ranked = rank_candidates(cand_cluster)
                        best = ranked.iloc[0]
                        is_ambiguous = (len(ranked) > 1 and
                                        np.isfinite(best.get('GlycanFDR', np.inf)) and
                                        (best.get('GlycanFDR', np.inf) == ranked.iloc[1].get('GlycanFDR', np.inf)))
                        qn.at[idx, 'IsomerAmbiguous'] = bool(is_ambiguous)
                        # 写回结果（含 Genes / GlySpec）
                        qn.at[idx, 'PlausibleStruct'] = best.get('PlausibleStruct', pd.NA)
                        qn.at[idx, 'AssignedGlyID'] = best.get('GlyID', pd.NA)
                        qn.at[idx, 'AssignedCorrectedComposition'] = best.get('CorrectedComposition', pd.NA)
                        qn.at[idx, 'ChosenGlycanFDR'] = best.get('GlycanFDR', np.nan)
                        qn.at[idx, 'ChosenGlyScore'] = best.get('GlyScore', np.nan)
                        qn.at[idx, 'ChosenTotalScore'] = best.get('TotalScore', np.nan)
                        qn.at[idx, 'ChosenPPM'] = best.get('PPM', np.nan)
                        qn.at[idx, 'ChosenRT'] = best.get('RT', np.nan)
                        qn.at[idx, 'QualRowIndex'] = int(best.get('__row_id', -1))
                        qn.at[idx, 'Genes'] = best.get('Genes', pd.NA)          # <---
                        qn.at[idx, 'GlySpec'] = best.get('GlySpec', pd.NA)      # <---
                        reason = []
                        if 'GlycanFDR' in ranked.columns:
                            reason.append(f"GlycanFDR={best.get('GlycanFDR', np.nan)} (min)")
                        if 'GlyScore' in ranked.columns:
                            reason.append(f"GlyScore={best.get('GlyScore', np.nan)} (max)")
                        if 'PPM' in ranked.columns:
                            reason.append(f"PPM={best.get('PPM', np.nan)} (min|abs)")
                        if 'CoreMatched' in ranked.columns:
                            reason.append(f"CoreMatched={best.get('CoreMatched', np.nan)}")
                        if 'TotalScore' in ranked.columns:
                            reason.append(f"TotalScore={best.get('TotalScore', np.nan)}")
                        reason.append(f"ClusterSize={len(cand_cluster)}/{len(cand)}")
                        qn.at[idx, 'AssignReason'] = '; '.join(map(str, reason))
                    return qn
                # ----------- 在这里改为你的路径 -----------
                # 读取 Excel（默认第一个工作表；如需指定请加 sheet_name='Sheet1'）
                quant_df = quant_data
                qual_df  = self.data.copy()
                # 执行分配
                out_df = assign_structures(
                    quant_df, qual_df,
                    ppm_tol=DEFAULT_PPM_TOL,
                    fdr_max=FDR_MAX,
                    allow_relax=ALLOW_FDR_RELAX,
                    fdr_relax_max=FDR_RELAX_MAX,
                    require_non_decoy=REQUIRE_NON_DECOY
                )
                # 保存为 Excel
                self.data = out_df
            #
            self.data = self.data.rename(columns={'Peptide': 'PeptideSequence',
                                                  'peptide': 'PeptideSequence',
                                                  'ProSites': 'Glycosite_Position',
                                                  'Prosite': 'Glycosite_Position',
                                                  'Proteins': 'ProteinID',
                                                  'Genes': 'GeneName',
                                                  'GlySpec': 'MS2Scan',
                                                  })
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
            self.data = self.data[~self.data['PlausibleStruct'].isnull()]
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
        #
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
                        pass  
                return ';'.join(accessions), ';'.join([g for g in genes if g])
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
        #
        elif self.search_engine == 'Byonic':
            def extract_glycosite_positions(peptide_full, start_pos):
                start_pos = int(start_pos)
                parts = peptide_full.split(".")
                pep = parts[1] if len(parts) == 3 else peptide_full
                clean_seq = re.sub(r'\[\+?\d+\.?\d*\]', '', pep)
                gly_sites = []
                for match in re.finditer(r'([A-Z])\[(\+?\d+\.?\d*)\]', pep):
                    aa = match.group(1)
                    mass_shift = float(match.group(2))
                    prefix = re.sub(r'\[\+?\d+\.?\d*\]', '', pep[:match.start()])
                    idx_in_pep = len(prefix)  # 0-based
                    abs_pos = start_pos + idx_in_pep
                    if aa == "N" and mass_shift > 500:
                        gly_sites.append(str(abs_pos))
                return ";".join(gly_sites) if gly_sites else None
            self.data["Glycosite_Position"] = self.data.apply(
                lambda row: extract_glycosite_positions(row["Peptide\n< ProteinMetrics Confidential >"], row["Starting\nposition"]),
                axis=1
            )
            self.data = self.data.rename(columns={'Peptide\n< ProteinMetrics Confidential >': 'PeptideSequence',
                                                  'Glycans\nNHFAGNa': 'GlycanComposition',
                                                  'Protein Name': 'ProteinID'})
            # protein id genename
            self.data[['ProteinID_new', 'GeneName']] = self.data['ProteinID'].str.extract(r'\|([A-Z0-9]+)\|.*GN=([A-Za-z0-9]+)')
            self.data['GeneName'] = self.data['GeneName'].str.capitalize()
            self.data = self.data.drop(columns=['ProteinID'])
            self.data = self.data.rename(columns={'ProteinID_new': 'ProteinID'})
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
                if not isinstance(comp, str):
                    return np.nan
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
            # peptidesequence
            def extract_peptide(seq):
                if not isinstance(seq, str):
                    return ''
                seq = re.sub(r'^[A-Z]\.|[.][A-Z]$', '', seq)
                seq = re.sub(r'\[[^\]]*\]', '', seq)
                seq = ''.join(re.findall(r'[A-Z]', seq))
                return seq
            self.data['PeptideSequence'] = self.data['PeptideSequence'].apply(extract_peptide)
            self.data['PeptideSequence+structure_coding+ProteinID'] = self.data['Query #:z']
        #
        elif self.search_engine == 'GlycanFinder':
            self.data = self.data.rename(columns={'Peptide': 'PeptideSequence',
                                                  'Glycosylation position': 'Glycosite_Position',
                                                  'Glycan': 'GlycanComposition',
                                                  'Accession': 'ProteinID'})
            # proteinid
            def split_and_extract(s):
                if pd.isnull(s):
                    return '', ''
                pairs = s.split(';')
                proteins = []
                genenames = []
                for pair in pairs:
                    items = pair.split('|')
                    if len(items) == 2:
                        proteins.append(items[0])
                        genename = items[1].split('_')[0].capitalize()
                        genenames.append(genename)
                return ';'.join(proteins), ';'.join(genenames)
            self.data[['ProteinID', 'GeneName']] = self.data['ProteinID'].apply(lambda x: pd.Series(split_and_extract(x)))
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
                if not isinstance(comp, str):
                    return np.nan
                parts = re.findall(r'\((\w+)\)(\d+)', comp)
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
            # peptidesequence
            def clean_peptide(seq):
                if not isinstance(seq, str):
                    return ''
                cleaned = re.sub(r'\([^)]+\)', '', seq)
                return cleaned.replace(' ', '')
            self.data['PeptideSequence'] = self.data['PeptideSequence'].apply(clean_peptide)
            self.data['PeptideSequence+structure_coding+ProteinID'] = self.data['Scan'].astype(str) + '+' + self.data['m/z'].astype(str)
        #
        self.data_manager.log_params('StrucGAP_Preprocess', 'data_cleaning', {})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_cleaned', self.data)
        #
        return self
    
    def cv_raw(self, threshold = None, fc_recommendation = True):
        """
        Raw Data CV Assessment: For each glycopeptide, prior to any normalization or aggregation, we will calculate the CV of TMT channel intensities (across biological replicates or relevant channels) directly from the raw data. This will provide an objective measure of data reproducibility and enable users to evaluate overall dataset quality prior to downstream analysis.
        
        Parameters:
            threshold: cv(raw) filter threshold (e.g. 0.3).
            fc_recommendation: whether execute fc recommendation pipeline (True or False).
        
        Returns:
            self.data_with_cv_raw (cleaned data with cv(raw) information). 
            self.data (cleaned and cv(raw) filtered data). 
            
        Return type:
            dataframe
        
        """
        def calculate_cv(values):
            # 过滤掉无效值（如0.0）
            valid_values = [v for v in values if v != 0.0 and not np.isnan(v)]
            # 计算CV值，只有有效值才参与计算
            if len(valid_values) > 0:
                mean = np.mean(valid_values)
                std = np.std(valid_values)
                cv = std / mean if mean != 0 else np.nan  # 避免除以0
            else:
                cv = np.nan  # 如果没有有效值，返回NaN
            return cv
        #
        def calculate_group_cvs(df, column_name='Matched_Reporter_Ions'):
            # 初始化两个空列表，用来存储实验组和对照组的CV
            exp_group_cvs = []
            control_group_cvs = []
            # 遍历数据框的每一行
            for _, row in df.iterrows():
                # 获取Matched_Reporter_Ions列的字符串，转换为字典
                reporter_ions_str = row[column_name]
                reporter_ions = ast.literal_eval(reporter_ions_str)  # 将字符串解析为字典
                # 提取实验组和对照组的通道值
                exp_group_values = [reporter_ions[channel][2] for channel in list(reporter_ions.keys())[:5]]  # 前五个通道是实验组
                control_group_values = [reporter_ions[channel][2] for channel in list(reporter_ions.keys())[5:]]  # 后五个通道是对照组
                # 分别计算实验组和对照组的CV
                exp_group_cv = calculate_cv(exp_group_values)
                control_group_cv = calculate_cv(control_group_values)
                # 将CV添加到相应的列表中
                exp_group_cvs.append(exp_group_cv)
                control_group_cvs.append(control_group_cv)
            # 将新的CV列添加到数据框中
            df['Raw_Group1_CV'] = exp_group_cvs
            df['Raw_Group2_CV'] = control_group_cvs
            return df
        #
        if self.search_engine == 'StrucGP' or (self.search_engine == 'pGlyco3' and self.quantification_from_no_strucgp):
            df = self.data.copy()
            df = calculate_group_cvs(df)
            self.data_with_cv_raw = df.copy()
            #
            if threshold == None:
                threshold = input("Please enter a threshold for raw CV filtering (e.g., 0.3), or 'no' to skip: ")
            if threshold != 'no':  
                try:
                    threshold = float(threshold)
                    self.data = self.data_with_cv_raw[
                        (self.data_with_cv_raw['Raw_Group1_CV'] < threshold) & 
                        (self.data_with_cv_raw['Raw_Group2_CV'] < threshold)
                    ]
                except ValueError:
                    print("Invalid input. Skipping raw CV filtering.")
            else:
                print("Skipping raw CV filtering.")
        #
        else:
            self.data = self.data.copy()
        #
        if fc_recommendation:
            df = self.data.copy()
            df["Matched_Reporter_Ions_dict"] = df["Matched_Reporter_Ions"].apply(ast.literal_eval)
            # 提取强度值（每个字典的 value 第3个元素，即 index=2）
            for channel in ["126.1277","127.1248","127.1311","128.1281","128.1344",
                            "129.1315","129.1378","130.1348","130.1411","131.1382"]:
                df[channel] = df["Matched_Reporter_Ions_dict"].apply(lambda d: d.get(channel, (None,None,None))[2])
            # 如果不需要中间的字典列，可以删除
            df = df[["126.1277","127.1248","127.1311","128.1281","128.1344",
                    "129.1315","129.1378","130.1348","130.1411","131.1382"]]
            # 计算每列 0 值的比例
            zero_ratio = (df == 0).sum() / len(df)
            df = df.loc[:, zero_ratio <= 0.5]
            self.fc_recommendation_based_on_raw_cv, self.fc_recommendation_based_on_raw_intensity = self.fc_recommendation(data = df)
            #
            self.data_manager.log_params('StrucGAP_Preprocess', 'data_cleaning', {'cv_threshold': threshold, 'fc_recommendation':fc_recommendation})
            self.data_manager.log_output('StrucGAP_Preprocess', 'data_with_cv_raw', self.data_with_cv_raw)
            self.data_manager.log_output('StrucGAP_Preprocess', 'data_cleaned_cv_filtered', self.data)
            self.data_manager.log_output('StrucGAP_Preprocess', 'fc_recommendation_based_on_raw_cv', self.fc_recommendation_based_on_raw_cv)
            self.data_manager.log_output('StrucGAP_Preprocess', 'fc_recommendation_based_on_raw_intensity', self.fc_recommendation_based_on_raw_intensity)
        else:
            self.data_with_cv_raw = pd.DataFrame()
            self.fc_recommendation_based_on_raw_cv = pd.DataFrame()
            self.fc_recommendation_based_on_raw_intensity = pd.DataFrame()
        #
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
        #
        elif self.search_engine != 'StrucGP':
            self.data_fdr_filtered = self.data
            self.data_peptide_fdr_data = self.data
            self.data_glycan_fdr_data = self.data
            self.data_fdr_filtered = self.data_fdr_filtered.set_index('PeptideSequence+structure_coding+ProteinID',drop=False)
        #
        self.data_manager.log_params('StrucGAP_Preprocess', 'fdr', {'feature_type': feature_type})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_peptide_fdr_data', self.data_peptide_fdr_data)
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_glycan_fdr_data', self.data_glycan_fdr_data)
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_fdr_filtered', self.data_fdr_filtered)
        #
        return self
    
    def outliers(self, samplewise_normalization = False, abundance_ratio=None):
        """
        Corrects and normalizes TMT quantification data from matched reporter ions for each IGP.
        
        Parameters:
            abundance_ratio: normalized factors from global proteomics data.
            samplewise_normalization: sample-wise median normalization.
        
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
            if samplewise_normalization is True:
                records = []
                # 假设 grouped 是原始 DataFrame，第一列是 peptide 名，其它列是样本
                for idx, row in grouped.iterrows():
                    glycopeptide = row.peptide
                    for sample in grouped.columns[1:]:
                        intensities = row[sample]
                        for val in intensities:
                            records.append((glycopeptide, sample, val))
                # 一次性创建 DataFrame，极大加快速度
                df_flattened = pd.DataFrame(records, columns=['Glycopeptide', 'Sample', 'Intensity']) 
                # 计算每个样本自己的中位数
                sample_medians = df_flattened.groupby('Sample')['Intensity'].median()
                # 每个 intensity 除以所属样本的中位数
                df_flattened['Normalized_Intensity'] = df_flattened.apply(
                    lambda row: row['Intensity'] / sample_medians[row['Sample']],
                    axis=1
                )
                grouped = df_flattened.groupby(['Glycopeptide', 'Sample'])['Normalized_Intensity'].apply(list).unstack()
                grouped.reset_index(inplace=True)
                grouped.rename(columns={'Glycopeptide': 'peptide'}, inplace=True)
                def replace_nan_in_list(x):
                    if isinstance(x, list):
                        return [0.0 if (pd.isna(i) or (isinstance(i, float) and math.isnan(i))) else i for i in x]
                    return x
                # 仅对样本列进行替换（排除 peptide 列）
                sample_cols = [col for col in grouped.columns if col != 'peptide']
                grouped[sample_cols] = grouped[sample_cols].applymap(replace_nan_in_list)
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
                    # if not np.isnan(median_val) and median_val != 0:
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
            #
            self.data_outliers_filtered = self.data_outliers_filtered[~self.data_outliers_filtered.index.duplicated()]
        #
        elif self.search_engine != 'StrucGP':
            self.data_outliers_filtered = self.data_fdr_filtered
            self.data_outliers_filtered = self.data_outliers_filtered[~self.data_outliers_filtered.index.duplicated()]
            self.abundance_ratio = abundance_ratio
            expected = self.sample_group.index.astype(str)
            actual   = self.data_outliers_filtered.columns.astype(str)
            to_drop  = expected.difference(actual)  # 这是一个 Index
            mask = self.sample_group.index.astype(str).isin(to_drop)
            self.sample_group = self.sample_group.loc[~mask]
        #
        self.data_manager.log_params('StrucGAP_Preprocess', 'outliers', {'abundance_ratio': abundance_ratio, 'samplewise_normalization': samplewise_normalization})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_outliers_filtered', self.data_outliers_filtered)
        #
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
        if self.search_engine == 'StrucGP' or (self.search_engine == 'pGlyco3' and self.quantification_from_no_strucgp):
            self.data_cv_filtered = copy.deepcopy(self.data_outliers_filtered)      
            # 取两个组各自需要的列名（转成字符串以和列名类型对齐）
            g0 = self.sample_group.loc[
                self.sample_group['group'] == self.sample_group['group'].unique()[0]
            ].index.astype(str)
            g1 = self.sample_group.loc[
                self.sample_group['group'] == self.sample_group['group'].unique()[1]
            ].index.astype(str)
            sub0 = self.data_cv_filtered.loc[:, self.data_cv_filtered.columns.isin(g0)]
            sub1 = self.data_cv_filtered.loc[:, self.data_cv_filtered.columns.isin(g1)]
            # 计算 CV（均值为 0 时设为 NaN 以避免除零）
            mean0, std0 = sub0.mean(axis=1), sub0.std(axis=1, ddof=1)
            mean1, std1 = sub1.mean(axis=1), sub1.std(axis=1, ddof=1)
            self.data_cv_filtered['cv_control'] = std0 / mean0.replace(0, np.nan)
            self.data_cv_filtered['cv_sample']  = std1 / mean1.replace(0, np.nan)
            #
            if threshold == None:
                threshold = input("Please enter a threshold for CV filtering (e.g., 0.3), or 'no' to skip: ")
            #
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
        #
        else:
            self.data_cv_filtered = self.data_outliers_filtered.copy()
        #
        self.data_manager.log_params('StrucGAP_Preprocess', 'cv', {'threshold': threshold})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_cv_filtered', self.data_cv_filtered)
        #
        return self
    
    def psm(self, psm_number = None, fc_recommendation = True):
        """
        Filters IGPs by the minimum number of supporting PSMs.
        
        Parameters:
            psm_filter: psm filter threshold (e.g. 3).
            fc_recommendation: whether execute fc recommendation pipeline (True or False).
        
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
        #
        elif self.search_engine != 'StrucGP':
            self.data_psm_filtered = self.data_cv_filtered.copy()
            psm_number = None
        # 
        if fc_recommendation:
            df = self.data_psm_filtered.copy()
            cols = ["126.1277","127.1248","127.1311","128.1281","128.1344",
                    "129.1315","129.1378","130.1348","130.1411","131.1382"]
            df = df[df.columns.intersection(cols)]
            # 计算每列 0 值的比例
            zero_ratio = (df == 0).sum() / len(df)
            df = df.loc[:, zero_ratio <= 0.5]
            self.fc_recommendation_based_on_preprocessed_cv, self.fc_recommendation_based_on_preprocessed_intensity = self.fc_recommendation(data = df)
            self.data_manager.log_output('StrucGAP_Preprocess', 'fc_recommendation_based_on_preprocessed_cv', self.fc_recommendation_based_on_preprocessed_cv)
            self.data_manager.log_output('StrucGAP_Preprocess', 'fc_recommendation_based_on_preprocessed_intensity', self.fc_recommendation_based_on_preprocessed_intensity)
        else:
            self.fc_recommendation_based_on_preprocessed_cv = pd.DataFrame()
            self.fc_recommendation_based_on_preprocessed_intensity = pd.DataFrame()
        #
        self.data_manager.log_params('StrucGAP_Preprocess', 'psm', {'psm_number': psm_number})
        self.data_manager.log_output('StrucGAP_Preprocess', 'data_psm_filtered', self.data_psm_filtered)
        #
        return self
    
    def annotation(self, glytoucan = True, glytoucan_structure = False, glytoucan_wurcs_file = None,
                   biosynthetic_pathways = False, glycobiology_filter = True):
        """
        Glycan plausibility annotation, whereby each glycan composition is mapped to GlyTouCan identifiers, cross-checked against KEGG biosynthetic rules, and further evaluated by rule-based filters derived from biosynthetic logic and chemical constraints, with user-selectable exclusion of implausible glycans.
        
        Parameters:
            glytoucan: whether the glycan composition has been annotated by GlyTouCan.
            glytoucan_structure: whether the plausible glycan structure has been annotated by GlyTouCan.
            biosynthetic_pathways: whether the glycan composition has been annotated by KEGG biosynthetic pathways.
            glycobiology_filter: whether the glycan composition has been annotated by rule-based filters.
        
        Returns:
            self.data_psm_filtered columns ['Glytoucan id', 'in_biosynthetic_pathways', 'RuleFlags']. 
            
        Return type:
            dataframe
        
        """
        ## Glytoucan
        API = "https://api.glycosmos.org/glycancompositionconverter/1.0.0/composition2wurcs"
        # N→hexnac，H→hex，F→dhex，S→neu5ac，G→neu5gc
        LETTER_TO_KEY = {
            "N": "hexnac",
            "H": "hex",
            "F": "dhex",
            "S": "neu5ac",
            "G": "neu5gc",
        }
        def parse_comp(s: str):
            s = (s or "").strip()
            counts = { "hex":"0","hexnac":"0","dhex":"0","neu5ac":"0","neu5gc":"0","P":"0","S":"0","Ac":"0" }
            for m in re.finditer(r'([A-Za-z]+)(\d+)', s):
                letter, num = m.group(1).upper(), m.group(2)
                if letter in LETTER_TO_KEY:
                    counts[LETTER_TO_KEY[letter]] = str(int(num))
                elif letter in ("P","S","AC"):
                    key = "Ac" if letter == "AC" else letter
                    counts[key] = str(int(num))
            return counts
        def build_payload(values):
            return [parse_comp(x) for x in values]
        def make_session(total_retries=5, backoff_factor=0.8):
            retry = Retry(
                total=total_retries,
                read=total_retries,
                connect=total_retries,
                status=total_retries,
                backoff_factor=backoff_factor,
                status_forcelist=(429, 500, 502, 503, 504),
                allowed_methods=frozenset(["POST","GET"]),
                raise_on_status=False,
            )
            s = requests.Session()
            s.headers.update({"User-Agent": "glycomp-pipeline/1.0"})
            adapter = HTTPAdapter(max_retries=retry, pool_connections=10, pool_maxsize=10)
            s.mount("https://", adapter)
            s.mount("http://", adapter)
            return s
        def post_chunk(session, chunk, timeout=(5, 120)):
            """对一个payload列表调用API；返回data列表"""
            r = session.post(API, json=chunk, timeout=timeout)
            r.raise_for_status()
            j = r.json()
            return j.get("data", [])
        def post_chunk_adaptive(session, chunk, timeout=(5, 120), min_size=1):
            """
            自适应批量：如果超时/网络错误，递归二分该批次直到成功或降到单条。
            返回与chunk等长的结果列表（顺序一致）。
            """
            try:
                data = post_chunk(session, chunk, timeout=timeout)
                if len(data) != len(chunk):
                    # 保护：返回条数不匹配时也按失败处理走拆分
                    raise RuntimeError(f"mismatched length: got {len(data)} for {len(chunk)}")
                return data
            except Exception as e:
                if len(chunk) <= min_size:
                    # 单条仍失败：最后再做几次指数退避重试
                    for i in range(3):
                        try:
                            time.sleep(1.5 * (2 ** i))
                            return post_chunk(session, chunk, timeout=(5, 180))
                        except Exception:
                            continue
                    # 实在不行，回一个“空结果”，占位
                    return [{"id": None, "wurcs": None}] * len(chunk)
                # 二分拆小再试
                mid = len(chunk) // 2
                left = post_chunk_adaptive(session, chunk[:mid], timeout=timeout, min_size=min_size)
                right = post_chunk_adaptive(session, chunk[mid:], timeout=timeout, min_size=min_size)
                return left + right
        def add_glytoucan_ids(df: pd.DataFrame, col_in: str, col_out: str = "glytoucan id",
                              batch_start=100, timeout=(5,120)) -> pd.DataFrame:
            """
            更稳健的批量转换：去重→分批→自适应批量→回填。
            """
            s = make_session()
            # 只对去重后的值调用一次
            uniq = pd.Index(df[col_in].astype(str).fillna("").unique())
            mapping = {}  # composition -> GlyTouCan ID (或None)
            # 按起始批量切分
            payload_all = build_payload(uniq)
            for i in range(0, len(payload_all), batch_start):
                chunk = payload_all[i:i+batch_start]
                data = post_chunk_adaptive(s, chunk, timeout=timeout)
                # 映射回 composition 字符串
                for comp_str, row in zip(uniq[i:i+batch_start], data):
                    gid = (row or {}).get("id") or None  # 有些没有ID，会返回空/None
                    mapping[comp_str] = gid
            # 回填到原DataFrame
            df[col_out] = df[col_in].astype(str).map(mapping)
            return df
        ## Glytoucan structure
        UP = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        LOW = "abcdefghijklmnopqrstuvwxyz"
        def residue_to_digit(res_desc: str) -> int:
            """
            残基描述 → 数字编码
            优先级：5（deoxyhexose/Fuc）> 4（*NCCO/3=O）> 3（Aad+*NCC/3=O）> 2（*NCC/3=O）> 1（默认 Hex）
            """
            s = res_desc or ""
            # 5：Fuc 等去氧糖（常见 a1221m-...）
            if "1221m" in s:
                return 5
            # 4：HexA 等（你给的 G4 场景：*NCCO/3=O）
            if "*NCCO/3=O" in s:
                return 4
            # 3：唾液酸（示例：Aad... 且 *NCC/3=O）
            if "Aad" in s and "*NCC/3=O" in s:
                return 3
            # 2：HexNAc（含 *NCC/3=O）
            if "*NCC/3=O" in s:
                return 2
            # 1：默认 Hex
            return 1

        def parse_wurcs(wurcs: str):
            """
            返回:
              - node_order: ['a','b',...]
              - node_to_resid: {'a':'1', ...}
              - resid_to_desc: {'1':'[...]', ...}
              - edges: 有序列表 [(parent, child), ...]  支持 a?-b1、g2-f3|f6、e4~n 等
            """
            # 1) 拆出连接片段（最后一个 / 后面）
            prefix, topology = wurcs.rsplit("/", 1)
            # 2) 拆出节点序列（倒数第二个 / 后面）
            prefix2, node_seq = prefix.rsplit("/", 1)
            # 残基列表（逐一编号 1,2,3...）
            resid_list = re.findall(r"\[([^\]]+)\]", prefix2)
            resid_to_desc = {str(i+1): desc for i, desc in enumerate(resid_list)}
            # 节点序列（a,b,c,... 对应 1-2-3-...）
            node_res_seq = node_seq.split("-")
            node_order = [chr(ord('a') + i) for i in range(len(node_res_seq))]
            node_to_resid = {node: node_res_seq[i] for i, node in enumerate(node_order)}
            # 3) 解析连接：保持“有序列表”，不要用 set
            edges = []
            for token in topology.split("_"):
                token = token.strip()
                if not token:
                    continue
                # 形如：
                #   a4-b1
                #   a?-b1
                #   g2-f3|f6
                #   e4~n   （这种没有 '-'，要跳过）
                m = re.match(r'^([a-z])[0-9\?]*-(.+)$', token)
                if not m:
                    # 例如 e4~n 没有 '-'，不是父子连边，忽略
                    continue
                parent = m.group(1)
                right = m.group(2)
                # 多候选子节点用 '|' 分割；每一段只取首个子字母
                for part in right.split("|"):
                    cm = re.match(r'([a-z])[0-9\?]*', part)
                    if cm:
                        child = cm.group(1)
                        edges.append((parent, child))
            return node_order, node_to_resid, resid_to_desc, edges

        def build_tree(edges, start='a'):
            """
            无向生成树：
            - 用 edges 构造无向图；
            - 以 start ('a') 为根（若不在图内则取最小字母）做 BFS，得到父子关系；
            - children 按字母序排列（后续还会在 DFS 时把 num==5 的孩子放到最后）。
            返回: root, children, reached_nodes
            """
            # 无向邻接表
            adj = defaultdict(list)
            nodes = set()
            for u, v in edges:
                adj[u].append(v); adj[v].append(u)
                nodes.add(u); nodes.add(v)
            if not nodes:
                return None, {}, set()
            root = start if start in nodes else min(nodes)
            # BFS 构造生成树
            parent = {root: None}
            order = [root]
            q = deque([root])
            while q:
                u = q.popleft()
                # 相邻节点按字母序，保证输出稳定
                for v in sorted(adj[u]):
                    if v not in parent:
                        parent[v] = u
                        order.append(v)
                        q.append(v)
            # 构建 children
            children = defaultdict(list)
            for v, p in parent.items():
                if p is not None:
                    children[p].append(v)
            for k in children:
                children[k] = sorted(children[k])
            return root, children, set(order)

        # def wurcs_to_custom(wurcs: str) -> str | None:
        from typing import Optional
        def wurcs_to_custom(wurcs: str) -> Optional[str]:
            try:
                node_order, node_to_resid, resid_to_desc, edges = parse_wurcs(wurcs)
                # 用无向生成树来还原父子关系（以 'a' 为首选根）
                root, children, nodes = build_tree(edges, start='a')
                if root is None:
                    return None
                # 节点数字映射
                node_digit = {}
                for node in nodes:
                    resid = node_to_resid.get(node)
                    if resid is None:
                        return None
                    desc = resid_to_desc.get(resid, "")
                    node_digit[node] = residue_to_digit(desc)
                # DFS 输出；把 num==5 的孩子排到最后，避免 F5fF1
                out = []
                def dfs(node: str, depth: int):
                    if depth >= len(UP):
                        raise ValueError("Tree depth exceeds supported levels.")
                    out.append(f"{UP[depth]}{node_digit[node]}")
                    ordered = sorted(children.get(node, []), key=lambda ch: (node_digit[ch] == 5, ch))
                    for ch in ordered:
                        dfs(ch, depth + 1)
                    out.append(LOW[depth])
                dfs(root, 0)
                return "".join(out)
            except Exception:
                return None
        ## kegg
        # ============== 全局配置（可按需调） ==============
        TIMEOUT = 20                 # 单次请求超时（秒）
        MAX_RETRY = 4                # 最大重试次数（含指数退避）
        BASE_PAUSE = 0.35            # 基础节流（秒），KEGG 建议 ≤ 3 req/s
        JITTER = 0.15                # 抖动（秒），避免“节律性”触发风控
        REST_BASES = ["https://rest.kegg.jp", "http://rest.kegg.jp"]
        DBGET_BASES = ["https://www.genome.jp"]     # 需要时可追加 "http://www.genome.jp"
        USE_PROXY = False            # 若内网必须走代理，设为 True 并填 PROXIES
        PROXIES = {"http": "http://your.proxy:port", "https": "http://your.proxy:port"}
        TARGET_PATHWAYS_PRIMARY = {
            "map00510": "N-glycan biosynthesis",
            "map00512": "Mucin type O-glycan biosynthesis",
            "map00513": "Various types of N-glycan biosynthesis",
        }
        # ============== 会话与重试 ==============
        def make_session(name="StrucGAP/1.2"):
            s = requests.Session()
            s.headers.update({"User-Agent": name})
            retry = Retry(
                total=MAX_RETRY,
                connect=MAX_RETRY,
                read=MAX_RETRY,
                backoff_factor=0.5,          # 指数退避 (0.5, 1, 2, 4 ...)
                status_forcelist=(429, 500, 502, 503, 504),
                allowed_methods=frozenset(["GET"]),
                raise_on_status=False,
            )
            adapter = HTTPAdapter(max_retries=retry, pool_connections=50, pool_maxsize=50)
            s.mount("https://", adapter)
            s.mount("http://", adapter)
            if USE_PROXY:
                s.proxies.update(PROXIES)
            return s
        SESS_DBGET = make_session("StrucGAP-DBGET/1.2")
        SESS_REST  = make_session("StrucGAP-KEGG/1.2")
        def _sleep():
            time.sleep(BASE_PAUSE + random.uniform(0, JITTER))
        # ============== 解析 structure_coding ==============
        def parse_structure_coding(code: str):
            counts = {"GlcNAc":0, "Man":0, "Gal":0, "Neu5Ac":0, "Neu5Gc":0, "Fuc":0}
            last_letter = None
            for ch in str(code):
                if ch.isalpha():
                    last_letter = ch
                elif ch.isdigit():
                    d = int(ch)
                    if d == 2: counts["GlcNAc"] += 1
                    elif d == 3: counts["Neu5Ac"] += 1
                    elif d == 4: counts["Neu5Gc"] += 1
                    elif d == 5: counts["Fuc"] += 1
                    elif d == 1:
                        if last_letter is not None and last_letter.upper() >= 'F':
                            counts["Gal"] += 1
                        else:
                            counts["Man"] += 1
            return counts
        # ============== 关键词（空格分隔；含/不含 Asn） ==============
        def build_bfind_queries_space(comp: dict, include_asn=True):
            parts = []
            if comp.get("GlcNAc"): parts.append(f"(GlcNAc){comp['GlcNAc']}")
            if comp.get("Man"):    parts.append(f"(Man){comp['Man']}")
            if comp.get("Gal"):    parts.append(f"(Gal){comp['Gal']}")
            if comp.get("Fuc"):    parts.append(f"(Fuc){comp['Fuc']}")
            if comp.get("Neu5Ac"): parts.append(f"(Neu5Ac){comp['Neu5Ac']}")
            if comp.get("Neu5Gc"): parts.append(f"(Neu5Gc){comp['Neu5Gc']}")
            q_main = " ".join(parts).strip()
            queries = []
            if q_main:
                queries.append(q_main)
                if "(Neu5Ac)" in q_main:
                    queries.append(q_main.replace("(Neu5Ac)", "(NeuAc)"))
                if include_asn:
                    queries.append(q_main + " (Asn)1")
                    if "(Neu5Ac)" in q_main:
                        queries.append(q_main.replace("(Neu5Ac)", "(NeuAc)") + " (Asn)1")
            # 去重
            out, seen = [], set()
            for q in queries:
                if q and q not in seen:
                    seen.add(q); out.append(q)
            return out
        # ============== 小缓存 ==============
        _cache_bfind = {}     # query(str) -> [G...]
        _cache_get   = {}     # G -> {comp dict}
        _cache_link  = {}     # G -> [path ids]
        # ============== DBGET bfind（稳健：多基址 + 重试 + 降级） ==============
        def bfind_glycan_gids(query: str, max_hit=200):
            if query in _cache_bfind:
                return _cache_bfind[query]
            params = {"dbkey":"glycan", "keywords":query, "mode":"bfind", "max_hit":str(max_hit)}
            qs = urllib.parse.urlencode(params, quote_via=urllib.parse.quote_plus)
            last_err = None
            for base in DBGET_BASES:
                url = f"{base}/dbget-bin/www_bfind_sub?{qs}"
                try:
                    resp = SESS_DBGET.get(url, timeout=TIMEOUT)
                    if resp.status_code == 200:
                        html = resp.text
                        gids = re.findall(r'G\d{5}', html)
                        # 去重保序
                        seen = set(); uniq = []
                        for g in gids:
                            if g not in seen:
                                seen.add(g); uniq.append(g)
                        _cache_bfind[query] = uniq
                        _sleep()
                        return uniq
                    last_err = f"HTTP {resp.status_code}"
                except requests.exceptions.RequestException as e:
                    last_err = str(e)
                    continue
            # —— 降级：尝试 REST /find/glycan/（召回差，但比完全失败好）
            try:
                # 把空格转成 '+'，同时去掉括号，用于 REST 粗搜
                rough = query.replace("(", "").replace(")", "").replace("  ", " ")
                rough = rough.replace(" ", "+")
                url = f"{REST_BASES[0]}/find/glycan/{rough}"
                resp = SESS_REST.get(url, timeout=TIMEOUT)
                if resp.status_code == 200:
                    gids = []
                    for ln in resp.text.strip().splitlines():
                        if "\t" not in ln: 
                            continue
                        gid = ln.split("\t",1)[0].strip()
                        if gid.startswith("G"):
                            gids.append(gid)
                    # 去重保序
                    seen = set(); uniq = []
                    for g in gids:
                        if g not in seen:
                            seen.add(g); uniq.append(g)
                    _cache_bfind[query] = uniq
                    _sleep()
                    return uniq
            except requests.exceptions.RequestException:
                pass
            _cache_bfind[query] = []
            return []
        # ============== REST get / link（带缓存与回退） ==============
        def rest_get_composition(gid: str):
            if gid in _cache_get: return _cache_get[gid]
            for base in REST_BASES:
                try:
                    url = f"{base}/get/gl:{gid}"
                    r = SESS_REST.get(url, timeout=TIMEOUT)
                    if r.status_code == 200:
                        comp_line = ""
                        for ln in r.text.splitlines():
                            if ln.startswith("COMPOSITION") or ln.startswith("Composition"):
                                comp_line = ln.split(None,1)[-1].strip(); break
                        comp = {name:int(cnt) for name,cnt in re.findall(r"\(([A-Za-z0-9/]+)\)\s*(\d+)", comp_line)}
                        _cache_get[gid] = comp
                        _sleep()
                        return comp
                except requests.exceptions.RequestException:
                    continue
            _cache_get[gid] = {}
            return {}
        def rest_link_pathways(gid: str):
            if gid in _cache_link: return _cache_link[gid]
            for base in REST_BASES:
                try:
                    url = f"{base}/link/pathway/gl:{gid}"
                    r = SESS_REST.get(url, timeout=TIMEOUT)
                    if r.status_code == 200:
                        lines = r.text.strip().splitlines()
                        paths = [ln.split("\t")[1].split(":")[1].strip() for ln in lines if "\t" in ln]
                        _cache_link[gid] = paths
                        _sleep()
                        return paths
                except requests.exceptions.RequestException:
                    continue
            _cache_link[gid] = []
            return []
        # ============== 组成严格比对 ==============
        def _normname(x):
            x = x.strip().lower()
            if x in ("neuac","neu5ac"): return "Neu5Ac"
            if x == "neu5gc": return "Neu5Gc"
            if x == "glcnac": return "GlcNAc"
            if x == "man": return "Man"
            if x == "gal": return "Gal"
            if x == "fuc": return "Fuc"
            if x == "asn": return "Asn"
            return x
        def comp_match(ours: dict, kegg_comp: dict, include_asn: bool):
            kc = {_normname(k): v for k, v in kegg_comp.items()}
            keys = ["GlcNAc","Man","Gal","Fuc","Neu5Ac","Neu5Gc"]
            if include_asn: keys.append("Asn")
            return all(ours.get(k,0) == kc.get(k,0) for k in keys)
        # ============== 核心：DataFrame 批量注释（批量去重、降请求量） ==============
        def annotate_df(df: pd.DataFrame, col="structure_coding"):
            # 先解析所有行，构造关键词（含/不含 Asn），并做“关键词级去重”
            parsed = []
            for idx, row in df.iterrows():
                # 行ID：优先用 df['id']，否则用索引标签 idx
                rid = row['id'] if ('id' in df.columns and pd.notna(row.get('id'))) else idx
                code = str(row.get(col, "") or "")
                comp = parse_structure_coding(code)
                queries = build_bfind_queries_space(comp, include_asn=True)
                parsed.append((rid, code, comp, queries))
            # 关键词去重后统一查询（对上万行非常关键）
            all_queries = []
            for _, _, _, qs in parsed:
                all_queries.extend(qs)
            seen_q, uniq_queries = set(), []
            for q in all_queries:
                if q and q not in seen_q:
                    seen_q.add(q); uniq_queries.append(q)
            # 批量 bfind（唯一关键词）
            query2gids = {}
            for q in uniq_queries:
                query2gids[q] = bfind_glycan_gids(q, max_hit=200)
            # 汇总每行候选 G（并集 → 去重），并选最佳
            rows = []
            for rid, code, comp, qs in parsed:
                cand_gids = []
                for q in qs:
                    cand_gids.extend(query2gids.get(q, []))
                seen, gids_all = set(), []
                for g in cand_gids:
                    if g not in seen:
                        seen.add(g); gids_all.append(g)
                best_gid, best_paths, mode = "", [], "no_candidate"
                # 1) 精确匹配（忽略 Asn）
                for g in gids_all:
                    kc = rest_get_composition(g)
                    if kc and comp_match(comp, kc, include_asn=False):
                        best_gid, best_paths, mode = g, rest_link_pathways(g), "exact_no_asn"; break
                # 2) 精确匹配（含 Asn）
                if not best_gid:
                    for g in gids_all:
                        kc = rest_get_composition(g)
                        if kc and comp_match(comp, kc, include_asn=True):
                            best_gid, best_paths, mode = g, rest_link_pathways(g), "exact_with_asn"; break
                # 3) 命中 00510/512/513 最多者
                if not best_gid and gids_all:
                    top_g, top_paths, score = None, [], -1
                    for g in gids_all:
                        paths = rest_link_pathways(g)
                        s = sum(1 for p in paths if p in TARGET_PATHWAYS_PRIMARY)
                        if s > score:
                            top_g, top_paths, score = g, paths, s
                    if top_g:
                        best_gid, best_paths, mode = top_g, top_paths, "biosyn_score"
                biosyn = [p for p in best_paths if p in TARGET_PATHWAYS_PRIMARY]
                rows.append({
                    "id": rid,
                    # "structure_coding": code,
                    # "composition_detected": " ".join([f"({k}){v}" for k,v in comp.items() if v]) or "(empty)",
                    # "bfind_queries": " ; ".join(qs),
                    # "candidate_gids": ",".join(gids_all),
                    # "best_gid": best_gid,
                    # "pathways_all": ",".join(best_paths),
                    "in_biosynthetic_pathways": ",".join(biosyn),
                    # "match_mode": mode,
                })
            return pd.DataFrame(rows)
        ## rule base
        def annotate_glycans(df):
            results = []
            for _, row in df.iterrows():
                tags = []
                comp = row['GlycanComposition']
                struct = row['structure_coding'][0]
                gtype = row['Glycan_type']
                # 从 GlycanComposition 提取数量
                def get_count(letter):
                    m = re.search(fr"{letter}(\d+)", comp)
                    return int(m.group(1)) if m else 0
                N = get_count("N")
                H = get_count("H")
                F = get_count("F")
                S = get_count("S")
                # 分支数
                antenna_est = struct.count("E")
                # 外部岩藻糖判断
                external_fuc = bool(re.search(r'(?!B)[A-Z]5', struct))  # B5 是核心，其它字母+5 是外部
                # 岩藻糖总个数（structure_coding 中数字5出现次数）
                total_fuc_count = struct.count("5")
                # 规则1 核心结构不完整
                if N < 2 or H < 3:
                    tags.append("CoreFail")
                # # 规则2 高甘露糖 + 外部岩藻糖
                # if gtype == "Oligo mannose" and F > 1 and external_fuc:
                #     tags.append("HighMan_ExtFuc")
                # 规则2 高甘露糖 + 外部或核心岩藻糖
                core_fuc = "B5" in struct 
                if gtype == "Oligo mannose" and (external_fuc or core_fuc):
                    tags.append("HighMan_Fuc")
                # 规则3 唾液酸超量
                if S > 2 * antenna_est:
                    tags.append("ExcessSia")
                # 规则4a 岩藻糖超量（严格）
                if F > (1 + antenna_est + 2):
                    tags.append("ExcessFuc")
                # 规则4b 岩藻糖超量（罕见）
                if total_fuc_count >= 5:
                    tags.append("ExcessFuc")
                # 规则5 天线数异常
                if antenna_est >= 4:
                    tags.append("ExcessAntenna")
                # 规则6 高甘露糖体积过大
                if gtype == "Oligo mannose" and H >= 10:
                    tags.append("HighMan_Huge")
                # 规则7 少分支且多 LacNAc 串联
                if antenna_est in [2, 3]:
                    # 查找连续多个 LacNAc 模式（X2Y1X2Y1...）
                    if re.search(r'(?:[E-Z]\d+F1){2,}', struct) or re.search(r'(?:[E-Z]2[A-Z]1){2,}', struct):
                        tags.append("MultiLacNAc_LowAnt")
                # 规则8 高甘露糖 + Bisect (N3H6-10)
                if gtype == "Oligo mannose" and N == 3 and 6 <= H <= 10:
                    tags.append("HighMan_Bisect")
                results.append(";".join(tags) if tags else np.nan)
            df['RuleFlags'] = results
            return df
        #
        if glytoucan:
            self.data_psm_filtered = add_glytoucan_ids(self.data_psm_filtered, col_in="GlycanComposition", col_out="Glytoucan id")
        #
        if glytoucan_structure:
            # df = pd.read_csv(glytoucan_wurcs_file)
            # df["structure_coding"] = df["WURCS"].apply(wurcs_to_custom)
            df = pd.read_csv(glytoucan_wurcs_file)
            df_grouped = (
                df.groupby("structure_coding")["GlyTouCan ID"]
                  .apply(lambda x: ",".join(x))
                  .reset_index()
            )
            df_grouped = df_grouped.rename(columns={"GlyTouCan ID": "GlyTouCan ID", 
                                                    "structure_coding": "structure_coding"})
            def extract_branches(code: str):
                """
                从结构编码中提取“分支字符串”多重集（忽略顺序）：
                1) 去掉前缀 A?B?C?（如 A2B2C1）；若不存在就不去掉；
                2) 用小写字母块分割；保留非空段（如 D1E2F1G4）。
                """
                if not isinstance(code, str):
                    return []
                # 去掉前缀 A[1-5]B[1-5]C[1-5]
                s = re.sub(r'^A[1-5]B[1-5]C[1-5]', '', code)
                # 按小写字母分割，得到每个分支的“顺序串”
                segs = [seg for seg in re.split(r'[a-z]+', s) if seg]
                return segs

            def structure_signature(code: str) -> str:
                """
                组合签名：
                - part1: 大写+数字对 以及 小写字母 的种类与计数（全局计数）
                - part2: 分支多重集（去掉 A?B?C? 前缀后，按小写切分得到的分支串的多重集），忽略分支顺序
                返回可比较的稳定字符串签名。
                """
                if not isinstance(code, str):
                    return None
                # 全局计数（大写+数字对 + 小写字母）
                caps_pairs = re.findall(r'[A-Z][1-5]', code)
                lowers = re.findall(r'[a-z]', code)
                cnt_all = Counter(caps_pairs + lowers)
                part1 = "|".join(f"{k}:{cnt_all[k]}" for k in sorted(cnt_all.keys()))
                # 分支多重集（忽略顺序）
                branches = extract_branches(code)
                cnt_br = Counter(branches)
                part2 = "|".join(f"{k}:{cnt_br[k]}" for k in sorted(cnt_br.keys()))
                return f"{part1}||BR||{part2}"
            
            left = self.data_psm_filtered.copy()
            # left = left.assign(structure_coding=left.index)
            left["__sig__"] = left["structure_coding"].map(structure_signature)
            right = df_grouped.copy()
            # right = right.assign(structure_coding=right.index)
            right["__sig__"] = right["structure_coding"].map(structure_signature)
            def _agg_glytoucan(series):
                ids = []
                for s in series.dropna().astype(str):
                    ids.extend([x.strip() for x in s.split(",") if x.strip()])
                ids = sorted(set(ids))
                return ",".join(ids)
            right_agg = (
                right.groupby("__sig__", dropna=False)
                     .agg({"GlyTouCan ID": _agg_glytoucan})
                     .reset_index()
                     .rename(columns={"GlyTouCan ID": "GlyTouCan structure"})
            )
            left["__sig__"] = left["__sig__"].astype("object")
            right_agg["__sig__"] = right_agg["__sig__"].astype("object")
            ra = right_agg.set_index("__sig__")
            ra.index = ra.index.astype("object")
            result_both = left.join(ra["GlyTouCan structure"], on="__sig__")
            result_both = result_both.drop(columns="__sig__")
            # result_both = left.join(
            #     right_agg.set_index("__sig__")["GlyTouCan structure"],
            #     on="__sig__"
            # )
            # result_both = result_both.drop(columns="__sig__")
            self.data_psm_filtered = result_both
            if glytoucan:
                mask = self.data_psm_filtered["Glytoucan id"].isna() & self.data_psm_filtered["GlyTouCan structure"].notna()
                self.data_psm_filtered.loc[mask, "Glytoucan id"] = "struc:" + self.data_psm_filtered.loc[mask, "GlyTouCan structure"]
        #
        if biosynthetic_pathways:
            df = self.data_psm_filtered.copy()
            df = annotate_df(df, col="structure_coding")
            df = df.set_index('id', drop=True)
            self.data_psm_filtered = pd.concat([self.data_psm_filtered, df], axis=1)
        #
        if glycobiology_filter:
            self.data_psm_filtered = annotate_glycans(self.data_psm_filtered)
        #
        # self.data_psm_filtered['structure_coding'] = self.data_psm_filtered['structure_coding'].str.extract(r'\+(.*?)\+')
        return self
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        
        with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_Preprocess.xlsx'), engine='xlsxwriter') as writer:
            if self.data is not None and not self.data.empty:
                self.data.to_excel(writer, sheet_name='cleaned_data')
            if self.data_with_cv_raw is not None and not self.data_with_cv_raw.empty:
                self.data_with_cv_raw.to_excel(writer, sheet_name='data_with_cv_raw')
            if self.fc_recommendation_based_on_raw_cv is not None and not self.fc_recommendation_based_on_raw_cv.empty:
                self.fc_recommendation_based_on_raw_cv.to_excel(writer, sheet_name='fc_recomm_raw_fc')
            if self.fc_recommendation_based_on_raw_intensity is not None and not self.fc_recommendation_based_on_raw_intensity.empty:
                self.fc_recommendation_based_on_raw_intensity.to_excel(writer, sheet_name='fc_recomm_raw_intens')
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
            if self.fc_recommendation_based_on_preprocessed_cv is not None and not self.fc_recommendation_based_on_preprocessed_cv.empty:
                self.fc_recommendation_based_on_preprocessed_cv.to_excel(writer, sheet_name='fc_recomm_preprocessed_fc')
            if self.fc_recommendation_based_on_preprocessed_intensity is not None and not self.fc_recommendation_based_on_preprocessed_intensity.empty:
                self.fc_recommendation_based_on_preprocessed_intensity.to_excel(writer, sheet_name='fc_recomm_preprocessed_intens')
            


