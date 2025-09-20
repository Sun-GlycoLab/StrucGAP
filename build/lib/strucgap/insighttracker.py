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
from scipy.stats import pearsonr
from strucgap.functionannotation import StrucGAP_FunctionAnnotation
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
## 数据分析管理模块
class StrucGAP_InsightTracker:
    def __init__(self):
        self.module_records = {}
        self.analysis_params = {}
        self.outputs = {}
        
    def register_module(self, module_name, module_instance, params):
        """
        Registers a module instance, records the initial parameters, and output locations.

        Parameters:
            module_name (str): The name of the module being registered.
            module_instance (object): The instance of the module to be registered.
            params (dict): The initial parameters used to configure the module.
    
        Returns:
            None
            
        """
        self.module_records[module_name] = {'instance': module_instance, 'params': params, 'outputs': {}}
        self.analysis_params[module_name] = {}

    def log_params(self, module_name, function_name, params):
        """
        Logs the parameters for each function call.
    
        Parameters:
            module_name (str): The name of the module.
            function_name (str): The name of the function being logged.
            params (dict): The parameters passed to the function.
    
        Returns:
            None
            
        """
        if module_name not in self.analysis_params:
            self.analysis_params[module_name] = {}
        self.analysis_params[module_name][function_name] = params

    def log_output(self, module_name, output_name, output_data):
        """
        Logs the output data of a module.
    
        Parameters:
            module_name (str): The name of the module.
            output_name (str): The name of the output being logged.
            output_data (object): The output data produced by the module.
    
        Returns:
            None
            
        """
        if module_name in self.module_records:
            self.module_records[module_name]['outputs'][output_name] = output_data
            self.outputs[output_name] = output_data

    def retrieve_data(self, module_name, output_name):
        """
        Retrieves output data from a specified module for use in other modules.
    
        Parameters:
            module_name (str): The name of the module.
            output_name (str): The name of the output to retrieve.
    
        Returns:
            object or None: The output data if found, otherwise None.
            
        """
        return self.module_records[module_name]['outputs'].get(output_name, None)

    def show_params(self, module_name):
        """
        Displays the analysis parameters for a specified module.

        Parameters:
            module_name (str): The name of the module whose parameters are to be displayed.
    
        Returns:
            dict or None: The analysis parameters for the specified module, or None if not found.
        """
        return self.analysis_params.get(module_name, None)

    def get_all_data(self):
        """
        Returns all registered modules' parameters and output data.
    
        Parameters:
            None.
    
        Returns:
            dict: A dictionary containing all registered module parameters and output data.
            
        """
        return self.module_records
    
    def output_analysis_params(self, output_dir='./analysis_result', output_file='GAP_analysis_params.xlsx'):
        """
        Outputs the analysis parameters to an Excel file.
    
        Parameters:
            output_dir (str): The directory where the output file will be saved (default is './analysis_result').
            output_file (str): The name of the output Excel file (default is 'GAP_analysis_params.xlsx').
    
        Returns:
            None
            
        """
        os.makedirs(output_dir, exist_ok=True)  # 创建输出目录
        output_path = os.path.join(output_dir, output_file)

        # 使用ExcelWriter管理多个Sheet
        with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
            for module_name, module_params in self.analysis_params.items():
                rows = []
                for function_name, params in module_params.items():
                    if isinstance(params, dict) and len(params) != 0:  
                        for sub_key, sub_value in params.items():
                            rows.append([function_name, sub_key, sub_value])
                    else:  # 如果不是字典，直接存储
                        rows.append([function_name, np.nan, np.nan])
                # 创建一个DataFrame用于存储结果
                df = pd.DataFrame(rows, columns=["Function", "Parameter", "Value"])
                df.to_excel(writer, sheet_name=module_name, index=False)
                
    def key_information_extraction(self, module, target_percentage = 0.08, min_value = 1):
        """
        Mines outputs from both StrucGAP_GlycoPeptideQuant, StrucGAP_FunctionAnnotation and StrucGAP_GlycoNetwork to identify structurally and functionally relevant glycan substructural features.
        
        Parameters:
            module:  StrucGAP_GlycoPeptideQuant, StrucGAP_FunctionAnnotation or StrucGAP_GlycoNetwork.
            target_percentage: when extract key information from StrucGAP_FunctionAnnotation, target_percentage means the key information ratio from both analysis results.
            min_value: when extract key information from StrucGAP_FunctionAnnotation, target_percentage represents the minimum number of substructural features contained in a term.
        
        Returns:
            None (output key information directly as a table).
        
        """
        if module not in ['StrucGAP_GlycoPeptideQuant', 'StrucGAP_FunctionAnnotation', 'StrucGAP_GlycoNetwork']:
            print("Please select module in ['StrucGAP_GlycoPeptideQuant', 'StrucGAP_FunctionAnnotation' or 'StrucGAP_GlycoNetwork']!")
            return
        
        if module == 'StrucGAP_GlycoPeptideQuant' and 'StrucGAP_GlycoPeptideQuant' in self.module_records.keys():
            StrucGAP_GlycoPeptideQuant_key_set = {}
            print('StrucGAP_GlycoPeptideQuant module has been executed, key information extracting ...')
            for i in ['core_structure','branches_structure','glycan_type','branches_count',
                      'glycan_composition','lacdinac','fucosylated_type','acgc']: # i='glycan_type'
                ratio_sheetname = f"result_{i}_ratio"
                up_sheetname = f"result_{i}_up"
                up_ratio_sheetname = f"result_{i}_up_ratio"
                down_sheetname = f"result_{i}_down"
                down_ratio_sheetname = f"result_{i}_down_ratio"
                #
                data_ratio = getattr(self.module_records['StrucGAP_GlycoPeptideQuant']['instance'], ratio_sheetname)
                data_ratio.iloc[:,1:] = data_ratio.iloc[:,1:].replace(0, np.nan)
                data_up = getattr(self.module_records['StrucGAP_GlycoPeptideQuant']['instance'], up_sheetname)
                data_up.iloc[:,1:] = data_up.iloc[:,1:].replace(0, np.nan)
                data_up_ratio = getattr(self.module_records['StrucGAP_GlycoPeptideQuant']['instance'], up_ratio_sheetname)
                data_up_ratio.iloc[:,1:] = data_up_ratio.iloc[:,1:].replace(0, np.nan)
                data_down = getattr(self.module_records['StrucGAP_GlycoPeptideQuant']['instance'], down_sheetname)
                data_down.iloc[:,1:] = data_down.iloc[:,1:].replace(0, np.nan)
                data_down_ratio = getattr(self.module_records['StrucGAP_GlycoPeptideQuant']['instance'], down_ratio_sheetname)
                data_down_ratio.iloc[:,1:] = data_down_ratio.iloc[:,1:].replace(0, np.nan)
                
                def is_monotonic(row):
                    values = row[1:].astype(float)  
                    increasing = all(values[i] <= values[i+1] for i in range(len(values)-1))
                    decreasing = all(values[i] >= values[i+1] for i in range(len(values)-1))
                    return increasing or decreasing
                
                output = pd.DataFrame()
                
                def stack_with_empty_row(final_df, new_df):
                    if not new_df.empty:
                        column_names = new_df.columns.tolist()
                        new_df.columns = [0,1,2,3,4,5]
                        new_df.loc[-1] = column_names  
                        new_df.index = new_df.index + 1 
                        new_df = new_df.sort_index()
                        # 堆叠数据
                        final_df = pd.concat([final_df, new_df], ignore_index=True)
                        # 插入空行
                        empty_row = pd.DataFrame({col: np.nan for col in new_df.columns}, index=[0])
                        final_df = pd.concat([final_df, empty_row], ignore_index=True)
                    return final_df
                
                # ratio
                df = pd.DataFrame(data_ratio)
                for col in df.columns[1:]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df_filtered = df.dropna()
                result = df_filtered[df_filtered.apply(is_monotonic, axis=1)]
                if not result.empty:
                    output = stack_with_empty_row(output, result)
                    
                # up
                df = pd.DataFrame(data_up_ratio)
                for col in df.columns[1:]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df_filtered = df.dropna()
                result = df_filtered[df_filtered.apply(is_monotonic, axis=1)]
                if not result.empty:
                    result1 = data_up[data_up[data_up.columns[0]].isin(result[data_up.columns[0]])]
                    result1['min_value'] = result1.iloc[:, 2:].min(axis=1)
                    result1 = result1[result1['min_value'] >= 3]
                    result1 = result1.drop(columns=['min_value'])
                    output = stack_with_empty_row(output, result1)
                    if not result1.empty:
                        result = result[result[result.columns[0]].isin(result1[result1.columns[0]])]
                        output = stack_with_empty_row(output, result)
                        
                # down
                df = pd.DataFrame(data_down_ratio)
                for col in df.columns[1:]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df_filtered = df.dropna()
                result = df_filtered[df_filtered.apply(is_monotonic, axis=1)]
                if not result.empty:
                    result1 = data_down[data_down[data_down.columns[0]].isin(result[data_down.columns[0]])]
                    result1['min_value'] = result1.iloc[:, 2:].min(axis=1)
                    result1 = result1[result1['min_value'] >= 3]
                    result1 = result1.drop(columns=['min_value'])
                    output = stack_with_empty_row(output, result1)
                    if not result1.empty:
                        result = result[result[result.columns[0]].isin(result1[result1.columns[0]])]
                        output = stack_with_empty_row(output, result)
                        
                StrucGAP_GlycoPeptideQuant_key_set[i] = output
                
                if i != 'glycan_composition':
                    ratio_sheetname = f"differential_analysis_{i}"
                    data_da = getattr(self.module_records['StrucGAP_GlycoPeptideQuant']['instance'], ratio_sheetname)
                    data_da.iloc[:,1:] = data_da.iloc[:,1:].replace(0, np.nan)
                    max_up_index = data_da.loc[
                        (data_da.iloc[:, 1] == data_da.iloc[:, 1].max()) & 
                        (data_da.iloc[:, 3] == data_da.iloc[:, 3].max())].index
                    max_down_index = data_da.loc[
                        (data_da.iloc[:, 2] == data_da.iloc[:, 2].max()) & 
                        (data_da.iloc[:, 4] == data_da.iloc[:, 4].max())].index
                    final_index = max_up_index.union(max_down_index)
                    data_da = data_da.loc[final_index].copy()
                    for idx in data_da.index:
                        if idx in max_up_index:
                            data_da.iloc[data_da.index.get_loc(idx), 2] = np.nan 
                            data_da.iloc[data_da.index.get_loc(idx), 4] = np.nan  
                        elif idx in max_down_index:
                            data_da.iloc[data_da.index.get_loc(idx), 1] = np.nan 
                            data_da.iloc[data_da.index.get_loc(idx), 3] = np.nan 
                    StrucGAP_GlycoPeptideQuant_key_set[f'da_{i}'] = data_da
                
                output_dir = './analysis_result'
                os.makedirs(output_dir, exist_ok=True)
                with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_GlycoPeptideQuant_key_information.xlsx'), engine='xlsxwriter') as writer:
                    for sheet_name, df in StrucGAP_GlycoPeptideQuant_key_set.items():
                        if 'da' in sheet_name:
                            df.to_excel(writer, sheet_name=sheet_name, index=True)
                        else:
                            df.iloc[:-1,:].to_excel(writer, sheet_name=sheet_name, index=True)
                        
        if module == 'StrucGAP_FunctionAnnotation' and 'StrucGAP_FunctionAnnotation' in self.module_records.keys():
            StrucGAP_FunctionAnnotation_key_set = {}
            function_data_type = self.analysis_params['StrucGAP_FunctionAnnotation']['go_function_structure']['function_data']
            database = self.analysis_params['StrucGAP_FunctionAnnotation']['ora']['terms']
            if 'GO' in database[0]:
                if len(database) == 3:
                    database = 'GO'
                else:
                    database = database
            if 'KEGG' in database[0]:
                database = 'KEGG'
            
            if database == 'GO':
                term_list = ['bp_core_structure','mf_core_structure','cc_core_structure',
                      'bp_glycan_type','mf_glycan_type','cc_glycan_type',
                      'bp_branches_structure','mf_branches_structure','cc_branches_structure',
                      'bp_branches_count','mf_branches_count','cc_branches_count',
                      'bp_sialicacid_count','mf_sialicacid_count','cc_sialicacid_count',
                      'bp_fucose_count','mf_fucose_count','cc_fucose_count',
                      'bp_sialicacid_structure','mf_sialicacid_structure','cc_sialicacid_structure',
                      'bp_fucose_structure','mf_fucose_structure','cc_fucose_structure',
                      'bp_lacdinac','mf_lacdinac','cc_lacdinac',
                      'bp_structurecoding','mf_structurecoding','cc_structurecoding',
                      'bp_fucosylated_type','mf_fucosylated_type','cc_fucosylated_type',
                      'bp_acgc','mf_acgc','cc_acgc',
                      ]
            elif database == 'KEGG':
                term_list = ['kegg_core_structure',
                      'kegg_glycan_type',
                      'kegg_branches_structure',
                      'kegg_branches_count',
                      'kegg_sialicacid_count',
                      'kegg_fucose_count',
                      'kegg_sialicacid_structure',
                      'kegg_fucose_structure',
                      'kegg_lacdinac',
                      'kegg_structurecoding',
                      'kegg_fucosylated_type',
                      'kegg_acgc',
                      ]
            
            for i in term_list: # i='mf_branches_count'
                
                print(i)
                count_sheetname = f"{i}_count"
                ratio_sheetname = f"{i}"
                
                data_count = getattr(self.module_records['StrucGAP_FunctionAnnotation']['instance'], count_sheetname)
                data_count.iloc[:,1:] = data_count.iloc[:,1:].replace(0, np.nan)
                data_ratio = getattr(self.module_records['StrucGAP_FunctionAnnotation']['instance'], ratio_sheetname)
                data_ratio.iloc[:,1:] = data_ratio.iloc[:,1:].replace(0, np.nan)
                
                def extract_columns_with_dynamic_threshold(df, target_percentage=target_percentage):
                    num_columns = len(df.columns)-1
                    if num_columns <= 10:
                        target_percentage=0.5
                    target_columns = int(num_columns * target_percentage)
                    first_column = df.columns[0]
                    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                    
                    best_threshold = None
                    best_column_count = float('inf')
                    
                    for threshold in thresholds:
                        selected_columns = [col for col in df.columns[1:] if df[col].max() > threshold]
                        if abs(len(selected_columns) - target_columns) < abs(best_column_count - target_columns):
                            best_threshold = threshold
                            best_column_count = len(selected_columns)
                    
                    selected_columns = [first_column] + [col for col in df.columns[1:] if df[col].max() > best_threshold]
                    
                    print(f"Using threshold {best_threshold}, selected {len(selected_columns)} columns out of {num_columns}")
                    return df[selected_columns]
                
                def filter_columns_by_min_value(df, min_value=1):
                    first_column = df.columns[0]
                    selected_columns = [first_column] + [col for col in df.columns[1:] if df[col].min() > min_value]
                    return df[selected_columns]
                
                output = pd.DataFrame()
                
                def stack_with_empty_row(final_df, new_df):
                    if not new_df.empty:
                        column_names = new_df.columns.tolist()
                        new_df.columns = list(range(0,len(new_df.columns)))
                        new_df.loc[-1] = column_names  
                        new_df.index = new_df.index + 1 
                        new_df = new_df.sort_index()
                        # 堆叠数据
                        final_df = pd.concat([final_df, new_df], ignore_index=True)
                        # 插入空行
                        empty_row = pd.DataFrame({col: np.nan for col in new_df.columns}, index=[0])
                        final_df = pd.concat([final_df, empty_row], ignore_index=True)
                    return final_df
                
                df = pd.DataFrame(data_ratio)
                columns_to_keep = df.iloc[:, 1:].columns[~df.iloc[:, 1:].isin([1]).any()]
                df = df[[df.columns[0]] + list(columns_to_keep)]
                result = extract_columns_with_dynamic_threshold(df, target_percentage = target_percentage)
                if not result.iloc[:,1:].empty:
                    result1 = data_count[result.columns]
                    if result.iloc[:,1:].shape[1]>10:
                        result1 = filter_columns_by_min_value(result1, min_value = min_value)
                    output = stack_with_empty_row(output, result1)
                    if not result1.empty:
                        result = result[list(result1.loc[0])]
                        output = stack_with_empty_row(output, result)
                
                if not output.iloc[:,1:].empty:
                    StrucGAP_FunctionAnnotation_key_set[i] = output
                output_dir = './analysis_result'
                os.makedirs(output_dir, exist_ok=True)
                with pd.ExcelWriter(os.path.join(output_dir, f'StrucGAP_FunctionAnnotation_{database}_{function_data_type}_key_information.xlsx'), engine='xlsxwriter') as writer:
                    for sheet_name, df in StrucGAP_FunctionAnnotation_key_set.items():
                        df.iloc[:-1,:].to_excel(writer, sheet_name=sheet_name, index=True)

        if module == 'StrucGAP_GlycoNetwork' and 'StrucGAP_GlycoNetwork' in self.module_records.keys():
            StrucGAP_GlycoNetwork_key_set = {}
            
            upregulation_only_in_glycopeptide = getattr(self.module_records['StrucGAP_GlycoNetwork']['instance'], 'protein_no_glyco_up').copy() 
            downregulation_only_in_glycopeptide = getattr(self.module_records['StrucGAP_GlycoNetwork']['instance'], 'protein_no_glyco_down').copy()
            
            def annotate_glyco_feature(row):
                features = []
                # 1. LacNAc
                if 'E2F1fe' in str(row['structure_coding']):
                    features.append('LacNAc')
                # 2. GlcNAc
                if 'E2e' in str(row['structure_coding']):
                    features.append('GlcNAc')
                # 3. Lewis
                structure = str(row['structure_coding'])
                if (
                    'E2F1G5gfF5fe' in structure or
                    'E2F1fF5fe' in structure or
                    'E2F1G3gfF5fe' in structure
                ):
                    features.append('Lewis')
                # 4. LacdiNAc (正则)
                if row['structure_coding'] and re.findall(r'(E2F2.*?fe)', str(row['structure_coding'])):
                    features.append('LacdiNAc')
                # 5. Sialylation
                if pd.notnull(row['Ac/Gc']) and str(row['Ac/Gc']).strip() != ' ':
                    features.append('Sialylation')
                # 6. Neu5Ac sialylation
                if 'Ac' in str(row['Ac/Gc']):
                    features.append('Neu5Ac sialylation')
                # 7. Neu5Gc sialylation
                if 'Gc' in str(row['Ac/Gc']):
                    features.append('Neu5Gc sialylation')
                # 8. Fucosylation
                if pd.notnull(row['fucosylated type']) and str(row['fucosylated type']).strip() != ' ':
                    features.append('Fucosylation')
                # 9. Core fucosylation
                if str(row['fucosylated type']) == 'core fucosylated':
                    features.append('Core fucosylation')
                # 10. Antenna fucosylation
                if str(row['fucosylated type']) == 'antenna fucosylated':
                    features.append('Antenna fucosylation')
                # combine
                features = list(dict.fromkeys(features))
                return ';'.join(features)
            
            upregulation_only_in_glycopeptide['GlycoFeature'] = upregulation_only_in_glycopeptide.apply(annotate_glyco_feature, axis=1)
            downregulation_only_in_glycopeptide['GlycoFeature'] = downregulation_only_in_glycopeptide.apply(annotate_glyco_feature, axis=1)
            
            def judge_lacnac(description):
                exclude = r"(?i)(keratan\s*sulfate|keratan-sulfate|keratan\s*polymer|within\s*(glycan\s*)?chains|keratan$)"
                strong_pattern = r"(?i)(involved in the synthesis of terminal.*LacNAc|synthesis of terminal.*LacNAc|formation of terminal.*LacNAc|catalyzes.*transfer of Gal.*to the non-reducing terminal N-acetyl[ -]?glucosamine|fucosylates.*LacNAc|modifies.*LacNAc unit)"
                if re.search(exclude, description):
                    return False
                elif re.search(strong_pattern, description):
                    return True
                else:
                    return False
                
            def judge_glcnac(description):
                strong_glcnac_pattern = r"""(?i)(
                    transfer(s|ed|ing)?\s.*\s(to|on)\s(terminal\s)?(beta-)?N-acetylglucosamine(\s*\(GlcNAc\))? |
                    (beta-)?N-acetylglucosamine(\s*\(GlcNAc\))?\s(residue|unit|structure|acceptor) |
                    N-acetylglucosaminyltransferase |
                    alpha-D-GlcNAc |
                    GlcNAc(\(2\))? |
                    N-glycosylation |
                    lactosaminide\s*\(Gal-beta-1,4-GlcNAc-R\) |
                    (GlcNAc|N-acetylglucosamine).*(acceptor|donor|unit|residue)
                )"""
                exclude_pattern = r"""(?i)(
                    may\salso\sshow.*activity\stoward.*GlcNAc.*(unsure|in vitro|synthetic\s+substrate|not.*major) |
                    O-mannosyl\s*glycosylation |
                    O-glycosylation |
                    O-linked |
                    O-GlcNAc |
                    O-Ser/Thr |
                    O-fucosylation |
                    O-glycosidic |
                    in vivo.*unsure |
                    synthetic substrate
                )"""
                if re.search(exclude_pattern, description, re.VERBOSE):
                    return False
                elif re.search(strong_glcnac_pattern, description, re.VERBOSE):
                    return True
                else:
                    return False
                
            def judge_lewis(description):
                exclude = r"(?i)(serve as (scaffold|backbone|epitope) for.*lewis|may be involved in.*lewis|selectively involved in.*lewis|can also.*lewis|required for.*lewis|backbone for.*lewis|scaffold for.*lewis|participating in (the )?biosynthesis of.*lewis|can lead to.*lewis|(forming|may contribute to the generation of|leads to the formation of).{0,40}(lewis\s*[axy]|lex|ley|sialyl[\s-]?lewis\s*[axy]|slex))"
                lewis_pattern = r"(?i)(catalyzes|fucosylates|biosynthesi[sz]es?|synthesi[sz]es?|generates|forms|transfer[s]?).{0,100}(lewis\s*[xy]|lex|ley|sialyl[\s-]?lewis\s*x|slex)"
                if re.search(exclude, description):
                    return False
                elif re.search(lewis_pattern, description):
                    return True
                else:
                    return False
                
            def judge_lacdinac(description):
                lacdi_pattern = r"""(?ix)
                    lacdinac
                    | n[,']?n'-?diacetyllactosediamine
                    | galnac[-\sβ]*1[-\s,>]*4[-\s]*glcnac
                    | n-acetylgalactosamine\s*\(galnac\)[^\.\n]{0,80}?beta-?1,?4-?linkage[^\.\n]{0,40}?glcnac
                    | galnac[^\.\n]{0,80}?beta-?1,?4-?linkage[^\.\n]{0,40}?glcnac
                    | galnac[^\.\n]{0,80}?beta-?1,?4[^\.\n]{0,40}?glcnac
                """
                if re.search(lacdi_pattern, description):
                    return True
                else:
                    return False

            def judge_sialylation(description):
                sialylation_pattern = r"""(?ix)
                (
                    sialyltransferase |
                    CMP-?sial(ic|yl)[ -]?acid |
                    catalyzes[^.]{0,50}transfer[^.]{0,20}sial(ic|yl)[ -]?acid |
                    transfer[sred]?[^.]{0,60}sial(ic|yl)[ -]?acid |
                    catalyzes[^.]{0,50}formation[^.]{0,40}NeuAc
                )
                """
                O_glycan_pattern = r"""(?ix)
                (
                    core 1 O-glycan |
                    O-?glycan(?!c) |
                    O-?linked |
                    Ser/?Thr |
                    Gal(?:beta|β)-(1->3)-GalNAc |
                    Gal(?:beta|β)?[^A-Za-z]{0,4}1[^A-Za-z]{0,4}3[^A-Za-z]{0,4}GalNAc[^A-Za-z]{0,5}O-Ser(?:/Thr)? |
                    mucin
                )
                """
                N_glycan_pattern = r"""(?ix)
                (
                    N-?glycan |
                    N-?linked |
                    N-?oligosaccharide |
                    type II glycan |
                    neolactoside |
                    terminal galactose residue |
                    Gal(?:beta|β)-1-4-GlcNAc |
                    NeuAc[^.;,()]{0,20}?Gal[^.;,()]{0,20}?GlcNAc
                )
                """
                if not re.search(sialylation_pattern, description):
                    return False
                else:
                    if re.search(O_glycan_pattern, description):
                        if re.search(N_glycan_pattern, description):
                            return True
                        else:
                            return False
                    else:
                        return True
                    
            def judge_neu5ac_sialylation(description):
                sialylation_pattern = r"(?i)(sialyltransferase|CMP-?sial(ic|yl)[ -]?acid|catalyzes[^.]{0,50}transfer[^.]{0,20}sial(ic|yl)[ -]?acid|transfer[sred]?[^.]{0,60}sial(ic|yl)[ -]?acid)"
                Neu5Ac_pattern = r"(?i)(Neu5Ac|N-acetyl-neuraminic acid|N-acetylneuraminic acid|CMP-?Neu5Ac|CMP-?N-acetyl-neuraminic acid|CMP-?N-acetyl-beta-neuraminate)"
                N_glycan_pattern = r"(?i)(N-?glycan|N-?linked|N-?oligosaccharide|type II glycan|Galbeta1-4GlcNAc|terminal galactose residue|type II lactosamine)"
                glycolipid_pattern = r"(?i)(glycolipid|ganglioside|globo series|ganglio series|GM1a?|GA1|GD1b|SSEA3|SSEA4)"
                O_glycan_pattern = r"(?i)(core 1 O-glycan|O-?glycan(?!c)|O-?linked|Ser/Thr|Gal(beta|β)-(1->3)-GalNAc|mucin)"
                if not re.search(sialylation_pattern, description):
                    return False
                elif not re.search(Neu5Ac_pattern, description):
                    return False
                elif re.search(N_glycan_pattern, description):
                    return True
                elif re.search(glycolipid_pattern, description):
                    return False
                elif re.search(O_glycan_pattern, description):
                    return False
                else:
                    return False
                
            def judge_neu5gc_sialylation(description):
                return False
                
            def judge_fucosylation(description):
                fucosylation_pattern = re.compile(
                    r"(fucosyltransferase|catalyzes[^.]{0,50}(transfer|addition)[^.]{0,20}fucos|"
                    r"transfer[sred]?[^.]{0,60}fucos|addition of fucose|fucosylates|fucosylation|"
                    r"GDP-?fucose|GDP-?β?-?L-?fucose|L-?fucose)",
                    re.I
                )
                N_glycan_pattern = re.compile(
                    r"("
                        r"N-?\s?linked|N-?\s?glycan|N-?\s?oligosaccharide|N-\s*or\s*O-linked|"
                        r"both\s+O-\s*and\s+N-linked|O-\s*and\s*N-linked|"
                        r"Asn[^.;,]{0,20}N[^P][ST]|N[^P][ST][^.;,]{0,20}Asn"
                        r"|"
                        r"(?:type\s*2\s*(?:lacto)?samine|LacNAc|polylactosamine|Galβ?1-4GlcNAc|"
                        r"type\s*II\s*(?:glycan|lactosamine)|neolactoside|terminal galactose residue)"
                        r"(?![^.]{0,120}\b(?:lactose|milk oligosaccharide)\b)"
                    r")",
                    re.I | re.S
                )
                N_exclusion = re.compile(r"\b(lactose|milk oligosaccharide)\b", re.I)
                if not fucosylation_pattern.search(description):
                    return False
                elif N_glycan_pattern.search(description) and not N_exclusion.search(description):
                    return True
                else:
                    return False

            def judge_core_fucosylation(description):
                core_fucosylation_pattern = (
                    r"(?i)(core fucosylation"
                    r"|fucose in (alpha|α)[ -]?1[,->]*6 linkage"
                    r"|fucosylation (at|of|on|to) (the )?(core|first) (GlcNAc|N-acetylglucosamine)"
                    r"|FUT8"
                    r"|catalyzes[^.]{0,50}(addition|transfer)[^.]{0,20}(fucose|fucosyl)[^.]{0,30}(alpha|α)[ -]?1[,->]*6"
                    r")"
                )
                N_glycan_pattern = r"(?i)(N-?glycan|N-?linked|N-?oligosaccharide|glycoprotein|N-acetylglucosamine|GlcNAc)"
                if re.search(core_fucosylation_pattern, description) and re.search(N_glycan_pattern, description):
                    return True
                else:
                    return False
                
            def judge_antenna_fucosylation(description):
                antenna_fucosylation_pattern = (
                    r"(?i)("
                    r"fucosylation (at|of|on|to) (the )?(terminal|distal|outer|non-reducing|antenna|arm|branch) (GlcNAc|N-acetylglucosamine)"
                    r"|catalyzes[^.]{0,80}(addition|transfer)[^.]{0,30}(fucose|fucosyl)[^.]{0,30}(alpha|α)[ -]?(1[,->]*3|1[,->]*4)"
                    r"|Lewis[ -]?(X|Y|x|y|A|a|B|b)"
                    r"|sialyl[ -]?Lewis[ -]?(X|x|A|a)"
                    r"|Lex"
                    r"|Ley"
                    r"|selectin ligand"
                    r")"
                )
                N_glycan_pattern = r"(?i)(N-?glycan|N-?linked|N-?oligosaccharide|glycoprotein|N-acetylglucosamine|GlcNAc)"
                if re.search(antenna_fucosylation_pattern, description) and re.search(N_glycan_pattern, description):
                    return True
                else:
                    return False

            # 其它特征也类似，可按需定制
            glyco_feature_judge_dict = {
                'LacNAc': judge_lacnac,
                'GlcNAc': judge_glcnac,
                'Lewis': judge_lewis,
                'LacdiNAc': judge_lacdinac,
                'Sialylation': judge_sialylation,
                'Neu5Ac sialylation': judge_neu5ac_sialylation,
                'Neu5Gc sialylation': judge_neu5gc_sialylation,
                'Fucosylation': judge_fucosylation, 
                'Core fucosylation': judge_core_fucosylation, 
                'Antenna fucosylation': judge_antenna_fucosylation, 
            }

            glyco_feature_to_enzyme_tables = {
                'LacNAc': ['glycosyltransferases'],
                'GlcNAc': ['glycosyltransferases'],
                'Lewis': ['glycosyltransferases'],
                'LacdiNAc': ['glycosyltransferases'],
                'Sialylation': ['glycosyltransferases', 'sialyltransferases'],
                'Neu5Ac sialylation': ['glycosyltransferases', 'sialyltransferases'],
                'Neu5Gc sialylation': ['glycosyltransferases', 'sialyltransferases'],
                'Fucosylation': ['glycosyltransferases', 'fucosyltransferase'],
                'Core fucosylation': ['glycosyltransferases', 'fucosyltransferase'],
                'Antenna fucosylation': ['glycosyltransferases', 'fucosyltransferase'],
            }

            def judge_lacnac_gbp(description):
                exclude_pattern = r"(?i)(does not bind lactose|does not bind (carbohydrate|galactoside|galactose))"
                lacunac_pattern = (
                    r"(?i)\b(beta-?galactoside|galactose-specific|lactose|LacNAc|N-acetyllactosamine|"
                    r"galectin|cell-surface glycans|glycan recognition)\b"
                )
                if re.search(exclude_pattern, description):
                    return False
                elif re.search(lacunac_pattern, description):
                    return True
                else:
                    return False
                
            def judge_glcnac_gbp(description):
                glcnac_pattern = (
                    r"(?i)\b("
                    r"binds?([^.;]*?)n-acetylglucosamine|"
                    r"binds?([^.;]*?)glcnac|"
                    r"n-acetylglucosamine-?binding|"
                    r"glcnac-?binding|"
                    r"catalyz(es|ing) ([^.;]*?)n-acetylglucosamine|"
                    r"catalyz(es|ing) ([^.;]*?)glcnac|"
                    r"transfer of n-acetylglucosamine|"
                    r"transfer of glcnac|"
                    r"addition of n-acetylglucosamine|"
                    r"addition of glcnac|"
                    r"n-acetylglucosaminyltransferase|"
                    r"glcnac-?transferase|"
                    r"glcnac\([^\)]*\)man|"
                    r"galnac[^.;]*?glcnac"
                    r")"
                )
                exclude_o_glycan_pattern = r"(?i)(O-?linked|O-?mannosyl|O-?glycosyl|O-?mannosylation|O-?glycan|O-?Ser/Thr)"
                exclude_undp_glcnac_pattern = (
                    r"(?i)("
                    r"UDP-GlcNAc|"
                    r"UDP-GalNAc|"
                    r"uridine diphosphate[- ]?(N-)?acetyl(glucosamine|galactosamine)|"
                    r"sugar nucleotide|"
                    r"hexosamine biosynthetic pathway|"
                    r"biosynthesis of (UDP-)?n-acetylglucosamine|"
                    r"converts [^.;]* to (UDP-)?n-acetylglucosamine|"
                    r"catalyzing the formation of (UDP-)?n-acetylglucosamine|"
                    r"pyrophosphorylation of (UDP-)?n-acetylglucosamine|"
                    r"precursor for [^.;]*glycosylation|"
                    r"GlcNAc[- ]?1[- ]?phosphate"
                    r")"
                )
                if (
                    re.search(exclude_o_glycan_pattern, description)
                    or re.search(exclude_undp_glcnac_pattern, description)
                ):
                    return False
                elif re.search(glcnac_pattern, description):
                    return True
                else:
                    return False
                
            def judge_lewis_gbp(description):
                lewis_pattern = (
                    r"(?i)\b("
                    r"sialyl[\s-]?lewis[\s-]?(x|a|y)|"
                    r"lewis[\s-]?(x|a|y)|"
                    r"\b(slex|sl[a-z]x|lex|lea|ley)\b"
                    r")"
                )
                if re.search(lewis_pattern, description):
                    return True
                else:
                    return False
                
            def judge_lacdinac_gbp(description):
                lacdiNAc_pattern = (
                    r"(?i)\b("
                    r"lacdiNac|"
                    r"galnac[- ]?β?1-4[- ]?glcnac|"
                    r"galnac\(β1-4\)glcnac|"
                    r"n[,']?n'-diacetyl[- ]?lactosediamine|"
                    r"diacetyl[- ]?lactosediamine|"
                    r"β1,4[- ]?n-acetylgalactosaminyltransferase|"
                    r"b4galnac|"
                    r"wfa[- ]?binding"
                    r")"
                )
                if re.search(lacdiNAc_pattern, description):
                    return True
                else:
                    return False

            def judge_sialylation_gbp(description):
                sialylation_pattern = (
                    r"(?i)\b("
                    r"siglec|"
                    r"sialic[- ]?acid[- ]?binding.*lectin|"
                    r"sialic[- ]?acid[- ]?dependent|"
                    r"binds? (to )?([a-zA-Z0-9\-]* )?(alpha-)?(2,3|2,6|α2,3|α2,6|2-3|2-6)[- ]?linked sialic acid|"
                    r"binds? ([a-zA-Z0-9\-]* )?sialic[- ]?acid([ -]?derivatives)?|"
                    r"containing sialic[- ]?acid( engaged in a (2,3|2,6|2-3|2-6) linkage)?|"
                    r"sialic[- ]?acid recognition|"
                    r"sialylation|"
                    r"sialylated"
                    r")\b"
                )
                if re.search(sialylation_pattern, description):
                    return True
                else:
                    return False
                    
            def judge_neu5ac_sialylation_gbp(description):
                return False

            def judge_neu5gc_sialylation_gbp(description):
                return False
                
            def judge_neu5gc_sialylation_gbp(description):
                return False
                
            def judge_fucosylation_gbp(description):
                fucosylation_pattern = re.compile(
                    r"(fucosyltransferase|catalyzes[^.]{0,50}(transfer|addition)[^.]{0,20}fucos|"
                    r"transfer[sred]?[^.]{0,60}fucos|addition of fucose|fucosylates|fucosylation|"
                    r"GDP-?fucose|GDP-?β?-?L-?fucose|L-?fucose)",
                    re.I
                )
                N_glycan_pattern = re.compile(
                    r"("
                        r"N-?\s?linked|N-?\s?glycan|N-?\s?oligosaccharide|N-\s*or\s*O-linked|"
                        r"both\s+O-\s*and\s+N-linked|O-\s*and\s*N-linked|"
                        r"Asn[^.;,]{0,20}N[^P][ST]|N[^P][ST][^.;,]{0,20}Asn"
                        r"|"
                        r"(?:type\s*2\s*(?:lacto)?samine|LacNAc|polylactosamine|Galβ?1-4GlcNAc|"
                        r"type\s*II\s*(?:glycan|lactosamine)|neolactoside|terminal galactose residue)"
                        r"(?![^.]{0,120}\b(?:lactose|milk oligosaccharide)\b)"
                    r")",
                    re.I | re.S
                )
                N_exclusion = re.compile(r"\b(lactose|milk oligosaccharide)\b", re.I)
                if not fucosylation_pattern.search(description):
                    return False
                elif N_glycan_pattern.search(description) and not N_exclusion.search(description):
                    return True
                else:
                    return False

            def judge_fucosylation_gbp(description):
                fucosylation_pattern = (
                    r"(?i)\b("
                    r"binds?([^.;]*?)fucose|"
                    r"affinity for ([^.;]*?)fucose|"
                    r"fucosylated|"
                    r"fucosylation|"
                    r"sialyl[ -]?lewis[ -]?x|"
                    r"(alpha-)?1,6-linked fucose.*(N-acetylglucosamine|GlcNAc)|"
                    r"hydrolyz(es|ing).*(alpha-)?1,6-linked fucose|"
                    r"fucosidase"
                    r")\b"
                )
                exclude_pattern = r"(?i)(transporter|mediates? (the )?uptake of (glucose|galactose|mannose|xylose|fucose|dehydroascorbate|fructose)|transport of (glucose|galactose|mannose|xylose|fucose|dehydroascorbate|fructose))"
                exclude_pattern2 = r"(?i)(interconversion|epimerase|isomerase|mutarotase|mutase|metabolism|catabolism|anabolism)"
                if re.search(exclude_pattern, description) or re.search(exclude_pattern2, description):
                    return False
                elif re.search(fucosylation_pattern, description):
                    return True
                else:
                    return False

            def judge_core_fucosylation_gbp(description):
                return False
                
            def judge_antenna_fucosylation_gbp(description):
                return False

            glyco_feature_gbp_judge_dict = {
                'LacNAc': judge_lacnac_gbp,
                'GlcNAc': judge_glcnac_gbp,
                'Lewis': judge_lewis_gbp,
                'LacdiNAc': judge_lacdinac_gbp,
                'Sialylation': judge_sialylation_gbp,
                'Neu5Ac sialylation': judge_neu5ac_sialylation_gbp,
                'Neu5Gc sialylation': judge_neu5gc_sialylation_gbp,
                'Fucosylation': judge_fucosylation_gbp, 
                'Core fucosylation': judge_core_fucosylation_gbp, 
                'Antenna fucosylation': judge_antenna_fucosylation_gbp, 
            }
            
            # 2. 缓存Uniprot Function
            function_cache = {}

            def get_function_from_uniprot(accession):
                if accession in function_cache:
                    return function_cache[accession]
                else:
                    # 在线请求API
                    func_str = get_function_from_uniprot_api(accession)
                    function_cache[accession] = func_str
                    return func_str

            def get_function_from_uniprot_api(accession):
                url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
                try:
                    response = requests.get(url, timeout=10)
                    if response.status_code != 200:
                        return ""
                    data = response.json()
                    # 找到Function注释
                    comments = data.get("comments", [])
                    for c in comments:
                        if c.get("commentType") == "FUNCTION":
                            texts = c.get("texts", [])
                            if texts:
                                return texts[0].get("value", "")
                    return ""
                except Exception as e:
                    # print(f"Error fetching {accession}: {e}")
                    return ""
                
            def get_gene_name_from_uniprot_api(accession):
                """
                通过Uniprot API，根据蛋白accession获取gene name（如ST3GAL6）。
                找不到则返回空字符串。
                """
                url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
                try:
                    response = requests.get(url, timeout=10)
                    if response.status_code != 200:
                        return ""
                    data = response.json()
                    # gene信息在"genes"列表里，主gene名称一般在"geneName"
                    genes = data.get("genes", [])
                    for gene in genes:
                        gene_name_data = gene.get("geneName")
                        if gene_name_data and gene_name_data.get("value"):
                            return gene_name_data["value"]
                    return ""
                except Exception as e:
                    # print(f"Error fetching {accession}: {e}")
                    return ""
                
            network_dict = {}
            for reg_label, glycopeptide_df in [
                ('Up', upregulation_only_in_glycopeptide),
                ('Down', downregulation_only_in_glycopeptide)
            ]:
                network_dict[reg_label] = {}
                for feature in glyco_feature_judge_dict.keys():
                    glycan_rows = glycopeptide_df[
                        glycopeptide_df['GlycoFeature'].str.contains(feature)
                    ].copy()
                    if glycan_rows.empty:
                        continue
                    # 聚糖节点：你可以合并为一个中心节点，或每个聚糖为节点
                    glycan_center_fc = glycan_rows['fc_g'].mean()  # 举例，中心fc为均值
                    feature_average_level = glycan_rows[[*map(str, getattr(self.module_records['StrucGAP_Preprocess']['instance'], 'sample_group').index)]].copy()
                    feature_average_level = feature_average_level.mean()
                    labels = ['Abundances (Normalized): F1: 126, Control',
                                          'Abundances (Normalized): F1: 127N, Control',
                                          'Abundances (Normalized): F1: 127C, Control',
                                          'Abundances (Normalized): F1: 128N, Control',
                                          'Abundances (Normalized): F1: 128C, Control',
                                          'Abundances (Normalized): F1: 129N, Sample',
                                          'Abundances (Normalized): F1: 129C, Sample',
                                          'Abundances (Normalized): F1: 130N, Sample',
                                          'Abundances (Normalized): F1: 130C, Sample',
                                          'Abundances (Normalized): F1: 131, Sample']
                    filtered_labels = [label for ratio, label in zip(getattr(self.module_records['StrucGAP_GlycoNetwork']['instance'], 'abundance_ratio'), labels) if ratio != 0]
                    # B. 寻找相关糖基转移酶 
                    enzyme_nodes = []
                    for table_name in glyco_feature_to_enzyme_tables.get(feature, ['glycosyltransferases']):
                        # table_name = 'fucosyltransferase'
                        enzyme_df = getattr(self.module_records['StrucGAP_GlycoNetwork']['instance'], table_name).copy().reset_index()
                        for _, row in enzyme_df.iterrows():
                            accession = row['Accession']
                            function_str = row.get('Function', '') or get_function_from_uniprot(accession)
                            if glyco_feature_judge_dict[feature](function_str):
                                # 筛选显著性 accession='Q8VIB3'
                                # if row['pvalue_ttest'] < 0.05 and (row['fc'] > 1.5 or row['fc'] < 1/1.5):
                                if row['pvalue_ttest'] < 0.05:
                                    enzyme_expr = row[filtered_labels]
                                    # print(enzyme_expr)
                                    gene_name = get_gene_name_from_uniprot_api(accession)
                                    if not gene_name:  
                                        gene_name = accession
                                    rho, pval = pearsonr(enzyme_expr.values.tolist(), feature_average_level.values.tolist())
                                    enzyme_nodes.append({
                                        'id': gene_name,
                                        'type': 'enzyme',
                                        'name': row.get('Gene', accession),
                                        'fc': row['fc'],
                                        'pvalue': row['pvalue_ttest'],
                                        'edge_weight': rho,
                                        'edge_pvalue': pval,
                                    })
                    enzyme_nodes = [dict(t) for t in {frozenset(node.items()) for node in enzyme_nodes}]
                    # C. 寻找相关糖结合蛋白
                    gbp_nodes = []
                    gbp_df = getattr(self.module_records['StrucGAP_GlycoNetwork']['instance'], 'glycan_binding_protein').copy().reset_index()
                    gbp_pattern = glyco_feature_gbp_judge_dict.get(feature)
                    if gbp_pattern:
                        for _, row in gbp_df.iterrows():
                            accession = row['Accession']
                            function_str = row.get('Function', '') or get_function_from_uniprot(accession)
                            if glyco_feature_gbp_judge_dict[feature](function_str):
                                # if row['pvalue_ttest'] < 0.05 and (row['fc'] > 1.5 or row['fc'] < 1/1.5):
                                if row['pvalue_ttest'] < 0.05:
                                    gbp_expr = row[filtered_labels]
                                    gene_name = get_gene_name_from_uniprot_api(accession)
                                    if not gene_name:  
                                        gene_name = accession
                                    rho, pval = pearsonr(gbp_expr.values.tolist(), feature_average_level.values.tolist())
                                    gbp_nodes.append({
                                        'id': gene_name,
                                        'type': 'GBP',
                                        'name': row.get('Gene', accession),
                                        'fc': row['fc'],
                                        'pvalue': row['pvalue_ttest'],
                                        'edge_weight': rho,
                                        'edge_pvalue': pval,
                                    })
                    gbp_nodes = [dict(t) for t in {frozenset(node.items()) for node in gbp_nodes}]
                    # D. downstream pathway
                    pathway_nodes = []
                    module5 = StrucGAP_FunctionAnnotation(glycan_rows,data_manager=self)  
                    try:
                        module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,
                                    selected_terms=['KEGG']) 
                        pathway = module5.ora_no_background_both_proteins_result
                        pathway = pathway[pathway['P-value']<0.05]
                        pathway['Term'] = pathway['Term'].str.extract(r'^(.*) \(')
                        pathway['Overlap_num'] = pathway['Overlap'].astype(str).str.extract(r'^(\d+)/')
                        # 合成节点数据
                        pathway_nodes += [
                            {
                                'term': row['Term'],
                                'pvalue': row['P-value'],
                                'overlap': int(row['Overlap_num']) if pd.notnull(row['Overlap_num']) else None
                            }
                            for _, row in pathway.iterrows()
                        ]
                    except Exception:
                        pathway_nodes += []
                    network_dict[reg_label][feature] = {
                        'glycan_center': {
                            'id': feature + '_center',
                            'type': 'glycan',
                            'fc': glycan_center_fc,
                            # 可继续加统计
                        },
                        'enzymes': enzyme_nodes,
                        'gbps': gbp_nodes,
                        'pathway': pathway_nodes,
                    }
            self.StrucGAP_GlycoNetwork_key_information = network_dict.copy()
            # output the results as excel        
            with pd.ExcelWriter('./analysis_result/StrucGAP_GlycoNetwork_key_information.xlsx') as writer:
                for reg_label in network_dict:
                    for feature in network_dict[reg_label]:
                        net = network_dict[reg_label][feature]
                        nodes = []
                        edges = []
                        # 节点
                        nodes.append({
                            'id': net['glycan_center']['id'],
                            'name': feature,
                            'type': 'glycan',
                            'fc': net['glycan_center']['fc']
                        })
                        for node in net['enzymes']:
                            nodes.append({**node, 'type': 'enzyme'})
                        for node in net['gbps']:
                            nodes.append({**node, 'type': 'GBP'})
                        for node in net['pathway']:
                            nodes.append({
                                'id': node['term'],
                                'name': node['term'],
                                'type': 'pathway',
                                'fc': '',  # pathway没fc，用空字符串
                                'overlap': node.get('overlap', '')
                            })
                        # 边
                        for node in net['enzymes'] + net['gbps']:
                            edges.append({
                                'source': node['id'],
                                'target': net['glycan_center']['id'],
                                'weight': node.get('edge_weight', ''),
                                'pvalue': node.get('edge_pvalue', ''),
                                'color': '' # 可加edge颜色
                            })
                        for node in net['pathway']:
                            edges.append({
                                'source': net['glycan_center']['id'],
                                'target': node['term'],
                                'weight': node.get('overlap', ''),
                                'pvalue': node.get('pvalue', ''),
                                'color': ''
                            })
                        nodes_df = pd.DataFrame(nodes)
                        edges_df = pd.DataFrame(edges)
                        # Sheet名如 Up-Fucosylation
                        nodes_df.to_excel(writer, sheet_name=f"{reg_label.capitalize()}-{feature}-Nodes", index=False)
                        edges_df.to_excel(writer, sheet_name=f"{reg_label.capitalize()}-{feature}-Edges", index=False)
    
    def output_pickle(self):
        """
        Serializes and saves all picklable variables (excluding Flask-related objects) to a file.
    
        Returns:
            None
        """
        # 定义一个函数来检查对象是否可以被序列化
        def is_picklable(obj):
            try:
                pickle.dumps(obj)
            except (pickle.PicklingError, TypeError):
                return False
            return True
        # 过滤出可序列化且非Flask相关的变量
        with open('all_variables.pkl', 'wb') as f:
            variables_to_save = {}
            for k, v in globals().items():
                # 先过滤掉 Flask 相关的上下文代理对象和模块类型
                if isinstance(v, werkzeug.local.LocalProxy) or isinstance(v, types.ModuleType):
                    continue
                # 只保存可以序列化的对象
                if is_picklable(v):
                    variables_to_save[k] = v
            pickle.dump(variables_to_save, f)
            
    def read_pickle(self):
        """
        Loads and restores previously serialized variables from a pickle file.
    
        Returns:
            None
        """
        # 确保在加载前定义所需的自定义函数和类
        def is_picklable(obj):
            try:
                pickle.dumps(obj)
            except (pickle.PicklingError, TypeError):
                return False
            return True
        # 从保存的文件中加载变量
        with open('all_variables.pkl', 'rb') as f:
            loaded_variables = pickle.load(f)
        # 将变量重新加载到当前全局命名空间中
        globals().update(loaded_variables)
    
    
 