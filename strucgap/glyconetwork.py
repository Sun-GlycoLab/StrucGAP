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
## 多组学联合分析模块--39
class StrucGAP_GlycoNetwork:
    def __init__(self, gs_data, data_manager):
        self.abundance_ratio = gs_data.abundance_ratio
        self.glycopeptide_data = gs_data.data 
        self.sample_group = gs_data.sample_group
        self.gs_data = gs_data
        
        self.sample_size = None
        
        self.data_manager = data_manager    
        self.data_manager.register_module('StrucGAP_GlycoNetwork', self, {})
        self.data_manager.log_params('StrucGAP_GlycoNetwork', '', {})
        
    def normal_distribution_detect(self, data, sample_size):
        """
        An auxiliary function called by other functions to perform normal distribution detection.
        
        Parameters:
            data: the data to be tested.
            sample_size: sample size in each group (e.g. 5).
        
        Returns:
            self.normal_list
            self.lognormal_list
                
        Return type:
            dataframe
        
        """
        if sample_size == None:
            sample_size = input("Please enter the sample size (such as: 5 controls vs 5 samples, just enter 5): ")
            self.sample_size = int(sample_size) 
        self.sample_size = int(sample_size)
        # data = data1
        normal_list = {}
        lognormal_list = {}
        for i in list(data.index):
            list1 = data.loc[i][:self.sample_size]
            list2 = data.loc[i][self.sample_size:]
            if (kstest(list1, cdf = 'norm')[1] > 0.05) & (kstest(list2, cdf = 'norm')[1] > 0.05):
                normal_list[i] = 'Normal distribution'
            else:
                normal_list[i] = 'No normal distribution'
            list1 = np.log2(list1)
            list2 = np.log2(list2)
            if (kstest(list1, cdf = 'norm')[1] > 0.05) & (kstest(list2, cdf = 'norm')[1] > 0.05):
                lognormal_list[i] = 'Lognormal distribution'
            else:
                lognormal_list[i] = 'No lognormal distribution'
        self.normal_list = pd.DataFrame(list(normal_list.items()), columns=['Protein', 'Status'])
        self.lognormal_list = pd.DataFrame(list(lognormal_list.items()), columns=['Protein', 'Status'])
        return self
    
    def outliers_detect(self, data):
        """
        An auxiliary function called by other functions to perform outliers detection.
        
        Parameters:
            data: the data to be tested.
        
        Returns:
            no_outliers_data
                
        Return type:
            dataframe
        
        """
        # Tukey's method
        for i in list(data.index):
            list1 = data.loc[i][:self.sample_size]
            list2 = data.loc[i][self.sample_size:]
            #
            Q1 = list1.quantile(0.25)
            Q3 = list1.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            # outliers = data_series[(data_series < lower_bound) | (data_series > upper_bound)]
            list1[(list1 < lower_bound) | (list1 > upper_bound)] = np.nan
            #
            Q1 = list2.quantile(0.25)
            Q3 = list2.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            # outliers = data_series[(data_series < lower_bound) | (data_series > upper_bound)]
            list2[(list2 < lower_bound) | (list2 > upper_bound)] = np.nan
            #
            data.loc[i][:self.sample_size] = list1
            data.loc[i][self.sample_size:] = list2
            
        no_outliers_data = data
        return no_outliers_data
    
    def missing_values_imputation(self, data):
        """
        An auxiliary function called by other functions to perform missing values imputation.
        
        Parameters:
            data: the data to be tested.
        
        Returns:
            no_missing_value_data
                
        Return type:
            dataframe
        
        """
        protein_ids = data.index
        #
        control_data = data.iloc[:, :self.sample_size]  
        experiment_data = data.iloc[:, self.sample_size:]  # self.sample_size*2
        #
        knn_imputer = KNNImputer(n_neighbors=self.sample_size)
        control_filled = knn_imputer.fit_transform(control_data)
        experiment_filled = knn_imputer.fit_transform(experiment_data)
        #
        control_filled_df = pd.DataFrame(control_filled, columns=control_data.columns, index=protein_ids)
        experiment_filled_df = pd.DataFrame(experiment_filled, columns=experiment_data.columns, index=protein_ids)
        no_missing_value_data = pd.concat([control_filled_df, experiment_filled_df], axis=1)
        #
        return no_missing_value_data
    
    def cv_filter(self, data, threshold = None):
        """
        An auxiliary function called by other functions to perform cv filtering.
        
        Parameters:
            data: the data to be tested.
        
        Returns:
            cv_filter_data
                
        Return type:
            dataframe
        
        """
        if threshold == None:
            threshold = input("Please enter the cv threshold (such as: 0.1, 0.2, 0.3): ")
        threshold = float(threshold)
        self.threshold = threshold
        data_samplewise_normalized = copy.deepcopy(data)
        #
        control_data = data_samplewise_normalized.iloc[:, :self.sample_size]  
        experiment_data = data_samplewise_normalized.iloc[:, self.sample_size:self.sample_size*2]
        control_cv = control_data.std(axis=1) / control_data.mean(axis=1)
        experiment_cv = experiment_data.std(axis=1) / experiment_data.mean(axis=1) 
        data_samplewise_normalized['Control_CV'] = control_cv
        data_samplewise_normalized['Experiment_CV'] = experiment_cv
        cv_filter_data = data_samplewise_normalized[(data_samplewise_normalized['Control_CV']<threshold)&(data_samplewise_normalized['Experiment_CV']<threshold)]
        cv_filter_data = cv_filter_data.drop(columns=['Control_CV', 'Experiment_CV'])
        #
        return cv_filter_data
    
    def normalization_samplewise(self, data, method = None):
        """
        An auxiliary function called by other functions to perform sample-wise normalization.
        
        Parameters:
            data: the data to be tested.
            method: normalization method in ['robust normalization', 'quantile normalization', 'z-score', 'range normalization', 'pareto normalization', 'level normalization', 'no'].
        
        Returns:
            data_samplewise_normalized
                
        Return type:
            dataframe
        
        """
        if method == None:
            method = input("Please enter the normalization method (select from: 'robust normalization', 'quantile normalization', 'z-score', 'range normalization', 'pareto normalization', 'level normalization' or 'no' to skip) among sample wise: ")
            expected_options = ['robust normalization', 'quantile normalization', 'z-score', 
                                'range normalization', 'pareto normalization', 'level normalization', 'no']
            matches = get_close_matches(method, expected_options, n=1, cutoff=0.5)
            if matches:
                method = matches[0]
                print(f"Using '{method}' as the input.")
            else:
                print("No close match found. Using 'robust normalization' as the input.")
                method = 'robust normalization'
            
        self.normalization_samplewise_method = method
        if method in ['robust normalization', 'quantile normalization', 'z-score', 'range normalization', 'pareto normalization', 'level normalization', 'no']:
            if method == 'robust normalization':
                protein_ids = data.index
                scaler = RobustScaler()
                data_samplewise_normalized = scaler.fit_transform(data)
                data_samplewise_normalized = pd.DataFrame(data_samplewise_normalized, columns=data.columns, index=protein_ids)
                return data_samplewise_normalized
            
            if method == 'quantile normalization':
                data_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
                data_samplewise_normalized = data.rank(method='min').stack().astype(int).map(data_mean).unstack()
                return data_samplewise_normalized
            
            if method == 'z-score':
                data_samplewise_normalized = pd.DataFrame(preprocessing.scale(data,axis=0), index=data.index, columns=data.columns) 
                return data_samplewise_normalized
            
            if method == 'range normalization':
                data_samplewise_normalized = copy.deepcopy(data)
                for i in list(data_samplewise_normalized.columns):  
                    data_samplewise_normalized[i] = (data_samplewise_normalized[i]-data_samplewise_normalized[i].mean())/(data_samplewise_normalized[i].max()-data_samplewise_normalized[i].min())
                return data_samplewise_normalized
                
            if method == 'pareto normalization':
                data_samplewise_normalized = copy.deepcopy(data)
                for i in list(data_samplewise_normalized.index):  
                    data_samplewise_normalized[i] = (data_samplewise_normalized[i]-data_samplewise_normalized[i].mean())/pow(data_samplewise_normalized[i].std(),1/2)
                return data_samplewise_normalized
                
            if method == 'level normalization':
                data_samplewise_normalized = copy.deepcopy(data)
                for i in list(data_samplewise_normalized.index):  
                    data_samplewise_normalized[i] = (data_samplewise_normalized[i]-data_samplewise_normalized[i].mean())/data_samplewise_normalized[i].mean()
                return data_samplewise_normalized
                
            if method == 'no':
                data_samplewise_normalized = data
                return data_samplewise_normalized
            
    def normalization_featurewise(self, data, method=None):
        """
        An auxiliary function called by other functions to perform feature-wise normalization.
        
        Parameters:
            data: the data to be tested.
            method: normalization method in ['robust normalization', 'quantile normalization', 'z-score', 'range normalization', 'pareto normalization', 'level normalization', 'no'].
        
        Returns:
            data_featurewise_normalized
                
        Return type:
            dataframe
        
        """
        if method == None:
            method = input("Please enter the normalization method (select from: 'robust normalization', 'quantile normalization', 'z-score', 'range normalization', 'pareto normalization', 'level normalization' or 'no' to skip) among sample wise: ")
            expected_options = ['robust normalization', 'quantile normalization', 'z-score', 
                                'range normalization', 'pareto normalization', 'level normalization', 'no']
            matches = get_close_matches(method, expected_options, n=1, cutoff=0.5)
            if matches:
                method = matches[0]
                print(f"Using '{method}' as the input.")
            else:
                print("No close match found. Using 'robust normalization' as the input.")
                method = 'robust normalization'
                
        self.normalization_featurewise_method = method
                
        data = data.T
        if method in ['robust normalization', 'quantile normalization', 'z-score', 'range normalization', 'pareto normalization', 'level normalization', 'no']:
            if method == 'robust normalization':
                protein_ids = data.index
                scaler = RobustScaler()
                data_featurewise_normalized = scaler.fit_transform(data)
                data_featurewise_normalized = pd.DataFrame(data_featurewise_normalized, columns=data.columns, index=protein_ids)
                return data_featurewise_normalized.T
            
            if method == 'quantile normalization':
                data_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
                data_featurewise_normalized = data.rank(method='min').stack().astype(int).map(data_mean).unstack()
                return data_featurewise_normalized.T
            
            if method == 'z-score':
                data_featurewise_normalized = pd.DataFrame(preprocessing.scale(data,axis=0), index=data.index, columns=data.columns) 
                return data_featurewise_normalized.T
            
            if method == 'range normalization':
                data_featurewise_normalized = copy.deepcopy(data)
                for i in list(data_featurewise_normalized.columns):  
                    data_featurewise_normalized[i] = (data_featurewise_normalized[i]-data_featurewise_normalized[i].mean())/(data_featurewise_normalized[i].max()-data_featurewise_normalized[i].min())
                return data_featurewise_normalized.T
                
            if method == 'pareto normalization':
                data_featurewise_normalized = copy.deepcopy(data)
                for i in list(data_featurewise_normalized.index):  
                    data_featurewise_normalized[i] = (data_featurewise_normalized[i]-data_featurewise_normalized[i].mean())/pow(data_featurewise_normalized[i].std(),1/2)
                return data_featurewise_normalized.T
                
            if method == 'level normalization':
                data_featurewise_normalized = copy.deepcopy(data)
                for i in list(data_featurewise_normalized.index):  
                    data_featurewise_normalized[i] = (data_featurewise_normalized[i]-data_featurewise_normalized[i].mean())/data_featurewise_normalized[i].mean()
                return data_featurewise_normalized.T
                
            if method == 'no':
                data_featurewise_normalized = data
                return data_featurewise_normalized.T
    
    def geometric_mean_of_middle(self, data, n):
        """
        An auxiliary function called by other functions to calculates the geometric mean of the middle 'n' elements in a given dataset.
        
        Parameters:
            data (list or iterable): the dataset from which the middle 'n' elements will be selected.
            n (int): the number of elements from the center of the data to use for the calculation.
        
        Returns:
            The geometric mean of the selected middle 'n' elements.
            
        Return type:
            float
        
        """
        length = len(data)
        start = (length - n) // 2 if length % 2 != 0 else ((length - n) // 2) + 1
        end = start + n
        middle_n = data[start:end]
        product = functools.reduce(operator.mul, middle_n, 1)
        geom_mean = product ** (1 / n)
        return geom_mean
    
    def glycosylation_rate(self, data):
        """
        An auxiliary function called by other functions to calculates glycosylation rate.
        
        Parameters:
            data: the data to be used.
        
        Returns:
            fc_result
            
        Return type:
            float
        
        """
        # fc
        data_c = self.cv_filter_data.iloc[:, :self.sample_size]  
        data_s = self.cv_filter_data.iloc[:, self.sample_size:self.sample_size*2] 
        #
        fc_result = pd.DataFrame(index=self.cv_filter_data.index, columns=['fc','pvalue_mannwhitneyu','pvalue_ttest', 'ttest_applicable', 'pvalue_ttest_mannwhitneyu'])
        fc = []
        pvalue_mannwhitneyu = []
        pvalue_ttest = []
        ttest_applicable = []
        pvalue_ttest_mannwhitneyu = []
        #
        for i in self.cv_filter_data.index:
            list1 = list(data_c.loc[i])
            list2 = list(data_s.loc[i]) 
            fc.append(statistics.mean(list2)/statistics.mean(list1)) 
            # 计算P值
            if (sum(list1)/len(list1))/(sum(list2)/len(list2))==1:
                pvalue_mannwhitneyu.append(1)
                pvalue_ttest.append(1)
                ttest_applicable.append('')
                pvalue_ttest_mannwhitneyu.append(1)
            else:
                pvalue_mannwhitneyu.append(stats.mannwhitneyu(list1, list2,alternative='two-sided')[1])
                if (kstest(list1, cdf = 'norm')[1] > 0.05) & (kstest(list2, cdf = 'norm')[1] > 0.05):
                    if stats.levene(list1, list2)[1] > 0.05:
                        pvalue_ttest.append(stats.ttest_ind(list1, list2,equal_var=True)[1])
                        pvalue_ttest_mannwhitneyu.append(stats.ttest_ind(list1, list2,equal_var=True)[1])
                    else:
                        pvalue_ttest.append(stats.ttest_ind(list1, list2,equal_var=False)[1])  
                        pvalue_ttest_mannwhitneyu.append(stats.ttest_ind(list1, list2,equal_var=False)[1])  
                    ttest_applicable.append('yes')
                else:
                    pvalue_ttest.append(stats.ttest_ind(list1, list2, equal_var=False)[1])
                    pvalue_ttest_mannwhitneyu.append(stats.mannwhitneyu(list1, list2,alternative='two-sided')[1]) 
                    ttest_applicable.append('no')
        #
        fc_result['fc'] = fc     
        fc_result['pvalue_mannwhitneyu'] = pvalue_mannwhitneyu 
        fc_result['pvalue_ttest'] = pvalue_ttest
        fc_result['ttest_applicable'] = ttest_applicable
        fc_result['pvalue_ttest_mannwhitneyu'] = pvalue_ttest_mannwhitneyu
        fc_result = pd.concat([data, fc_result], axis=1) 
        #
        return fc_result
    
    def glycoprotein_site(self, data):
        """
        Auxiliary computation function based on glycosylation sites of glycoproteins.
        
        Parameters:
            data
        
        Returns:
            self.glycoprotein_site_count
            self.glycoprotein_site_value
            
        Return type:
            float
        
        """
        #
        temp_data = pd.DataFrame(data[['ProteinID', 'Glycosite_Position', 'fc_g', 'fc_p', 'p_g', 'p_p']])
        result = []
        for idx, row in temp_data.iterrows():
            peptide = row['ProteinID']
            glycosite = row['Glycosite_Position']
            filtered_data = temp_data[(temp_data['ProteinID'] == peptide) & (temp_data['Glycosite_Position'] == glycosite)]
            #
            g_up_count = filtered_data[(filtered_data['fc_g']>1)&(filtered_data['p_g']<0.05)].shape[0]
            g_down_count = filtered_data[(filtered_data['fc_g']<1)&(filtered_data['p_g']<0.05)].shape[0]
            if filtered_data[(filtered_data['fc_p']>1)&(filtered_data['p_p']<0.05)].shape[0] != 0:
                p_type = 'Up'
            elif filtered_data[(filtered_data['fc_p']<1)&(filtered_data['p_p']<0.05)].shape[0] != 0:
                p_type = 'Down'
            elif filtered_data[(filtered_data['p_p']>0.05)].shape[0] != 0:
                p_type = 'None significance'
            result.append([peptide, glycosite, g_up_count, g_down_count, p_type]) 
        glycoprotein_site_count = pd.DataFrame(result, columns=['ProteinID', 'Glycosite_Position', 'Up_glycopeptide_count', 'Down_glycopeptide_count', 'Protein_type'])
        self.glycoprotein_site_count = glycoprotein_site_count.drop_duplicates()
        #
        result = []
        for idx, row in temp_data.iterrows():
            peptide = row['ProteinID']
            glycosite = row['Glycosite_Position']
            filtered_data = temp_data[(temp_data['ProteinID'] == peptide) & (temp_data['Glycosite_Position'] == glycosite)]
            #
            g_up_value = filtered_data[(filtered_data['fc_g']>1)&(filtered_data['p_g']<0.05)]['fc_g'].mean()
            g_down_value = filtered_data[(filtered_data['fc_g']<1)&(filtered_data['p_g']<0.05)]['fc_g'].mean()
            p_value = filtered_data['fc_p'].mean()
            result.append([peptide, glycosite, g_up_value, g_down_value, p_value]) 
        glycoprotein_site_value = pd.DataFrame(result, columns=['ProteinID', 'Glycosite_Position', 'Up_glycopeptide_fc', 'Down_glycopeptide_fc', 'Protein_fc'])
        self.glycoprotein_site_value = glycoprotein_site_value.drop_duplicates()
        #
        return self
                
    def glycopeptide_site(self, data):
        """
        Auxiliary computation function based on glycopeptide sites of glycoproteins.
        
        Parameters:
            data
        
        Returns:
            self.glycopeptide_site_count
            self.glycopeptide_site_value
            
        Return type:
            float
        
        """
        #
        temp_data = pd.DataFrame(data[['PeptideSequence', 'Glycosite_Position', 'fc_g', 'fc_p', 'p_g', 'p_p']])
        result = []
        for idx, row in temp_data.iterrows():
            peptide = row['PeptideSequence']
            glycosite = row['Glycosite_Position']
            filtered_data = temp_data[(temp_data['PeptideSequence'] == peptide) & (temp_data['Glycosite_Position'] == glycosite)]
            #
            g_up_count = filtered_data[(filtered_data['fc_g']>1)&(filtered_data['p_g']<0.05)].shape[0]
            g_down_count = filtered_data[(filtered_data['fc_g']<1)&(filtered_data['p_g']<0.05)].shape[0]
            if filtered_data[(filtered_data['fc_p']>1)&(filtered_data['p_p']<0.05)].shape[0] != 0:
                p_type = 'Up'
            elif filtered_data[(filtered_data['fc_p']<1)&(filtered_data['p_p']<0.05)].shape[0] != 0:
                p_type = 'Down'
            elif filtered_data[(filtered_data['p_p']>0.05)].shape[0] != 0:
                p_type = 'None significance'
            result.append([peptide, glycosite, g_up_count, g_down_count, p_type]) 
        glycopeptide_site_count = pd.DataFrame(result, columns=['PeptideSequence', 'Glycosite_Position', 'Up_glycopeptide_count', 'Down_glycopeptide_count', 'Protein_type'])
        self.glycopeptide_site_count = glycopeptide_site_count.drop_duplicates()
        #
        result = []
        for idx, row in temp_data.iterrows():
            peptide = row['PeptideSequence']
            glycosite = row['Glycosite_Position']
            filtered_data = temp_data[(temp_data['PeptideSequence'] == peptide) & (temp_data['Glycosite_Position'] == glycosite)]
            #
            g_up_value = filtered_data[(filtered_data['fc_g']>1)&(filtered_data['p_g']<0.05)]['fc_g'].mean()
            g_down_value = filtered_data[(filtered_data['fc_g']<1)&(filtered_data['p_g']<0.05)]['fc_g'].mean()
            p_value = filtered_data['fc_p'].mean()
            result.append([peptide, glycosite, g_up_value, g_down_value, p_value]) 
        glycopeptide_site_value = pd.DataFrame(result, columns=['PeptideSequence', 'Glycosite_Position', 'Up_glycopeptide_fc', 'Down_glycopeptide_fc', 'Protein_fc'])
        self.glycopeptide_site_value = glycopeptide_site_value.drop_duplicates()
        #
        return self
    
    def transcriptomic(self):
        pass
    
    def proteomic(self, protein_data_dir, fc=1.5, pvalue=0.05, data_sheet_name=None, pvalue_type='pvalue_ttest_mannwhitneyu',
                  fdr = 'Medium', psm = 3, cv = 0.3, 
                  normalization_samplewise_method = 'robust normalization',
                  normalization_featurewise_method = 'robust normalization'):
        """
        Implements a comprehensive preprocessing pipeline for global proteomic datasets.
        
        Parameters:
            protein_data_dir: proteomic data file directory.
            data_sheet_name: data sheet name in proteomic data file.
            fc: FC threshold used for differential analysis.
            pvalue: P value used for differential analysis.
            pvalue_type: 'pvalue_ttest', 'pvalue_mannwhitneyu' or 'pvalue_ttest_mannwhitneyu'.
            fdr: Proteomics data fdr filtering threshold (select from: Low, Medium, or no).
            psm: Proteomics data psm filtering threshold (e.g. 3).
            cv: Proteomics data cv filtering threshold (e.g. 0.3).
            normalization_samplewise_method: Samplewise normalization methods.
            normalization_featurewise_method: Featurewise normalization methods.
        
        Returns:
            self.protein_raw_data.
            self.normal_distribution_result.
            self.no_outliers_data.
            self.no_missing_value_data.
            self.cv_filter_data.
            self.samplewise_normalized_data.
            self.samplewise_featurewise_normalized_data.
            self.proteomic_fc.
            self.glycopeptide_fc.
            self.pg_fc (combined global proteins fc and glycopeptides fc results).
            self.proteomic_protein_glycosite_count.
            self.proteomic_protein_glycosite_value.
            self.proteomic_glycopeptide_glycosite_count.
            self.proteomic_glycopeptide_glycosite_value.
            self.proteomic_protein_glycosite_same_direction (the glycopeptides from glycosylation sites on the same glycoprotein exhibit the same regulation direction).
            self.proteomic_protein_glycosite_different_direction (the glycopeptides from glycosylation sites on the same glycoprotein exhibit the different regulation direction).
            self.proteomic_glycopeptide_glycosite_same_direction (the glycopeptides from glycosylation sites on the same glycopeptide exhibit the same regulation direction).
            self.proteomic_glycopeptide_glycosite_different_direction (the glycopeptides from glycosylation sites on the same glycopeptide exhibit the different regulation direction).
            self.glycoprotein_site_count (the number of upregulated and downregulated glycopeptides from glycosylation sites on the same glycoprotein).
            self.glycoprotein_site_value (the FC of upregulated and downregulated glycopeptides from glycosylation sites on the same glycoprotein).
            self.glycopeptide_site_count (the number of upregulated and downregulated glycopeptides from glycosylation sites on the same glycopeptide).
            self.glycopeptide_site_value (the FC of upregulated and downregulated glycopeptides from glycosylation sites on the same glycopeptide).
            self.protein_up_glyco_up (upregulated glycopeptide on upregulated protein).
            self.protein_up_glyco_no.
            self.protein_up_glyco_down.
            self.protein_no_glyco_up.
            self.protein_no_glyco_no.
            self.protein_no_glyco_down.
            self.protein_down_glyco_up.
            self.protein_down_glyco_no.
            self.protein_down_glyco_down.
            self.protein_up_glyco_up_tree_structure (tree structure data of upregulated glycopeptides on upregulated proteins).
            self.protein_up_glyco_no_tree_structure.
            self.protein_up_glyco_down_tree_structure.
            self.protein_no_glyco_up_tree_structure.
            self.protein_no_glyco_no_tree_structure.
            self.protein_no_glyco_down_tree_structure.
            self.protein_down_glyco_up_tree_structure.
            self.protein_down_glyco_no_tree_structure.
            self.protein_down_glyco_down_tree_structure.
            
        Return type:
            float
        
        """
        protein_data = pd.read_excel(protein_data_dir, sheet_name=data_sheet_name or 0)
        #
        # protein_data['ProteinID'] = [x[2] for x in protein_data['FASTA Title Lines'].str.split('|')]
        # protein_data['ProteinID'] = [x[0] for x in protein_data['ProteinID'].str.split('_')]
        if 'Description' in list(protein_data.columns):
            protein_data['ProteinID'] = protein_data['Description'].str.extract(r'GN=([^\s]+)')[0].str.upper()
        #
        # protein_data = protein_data[~protein_data['ProteinID'].isnull()]
        #
        if 'Score Sequest HT: Sequest HT' in list(protein_data.columns):
            protein_data = protein_data[protein_data['Score Sequest HT: Sequest HT']!=0]
        #
        # fdr
        if 'Protein FDR Confidence: Combined' in list(protein_data.columns):
            if fdr is None:
                fdr = input('Please enter a level of fdr (select from: Low, Medium, or no) that is unacceptable to you: ')
            expected_options = ['Low', 'Medium', 'no']
            matches = get_close_matches(fdr, expected_options, n=1, cutoff=0.5)
            if matches:
                fdr = matches[0]
                print(f"Using '{fdr}' as the input.")
            else:
                print("No close match found. Using 'Low' as the input.")
                fdr = 'Low'
            
            if fdr in ['Low', 'Medium', 'no']:
                if fdr == 'Low':
                    protein_data = protein_data[protein_data['Protein FDR Confidence: Combined']!='Low']
                elif fdr == 'Medium':
                    protein_data = protein_data[protein_data['Protein FDR Confidence: Combined']=='High']
        
        # psm
        if '# PSMs' in list(protein_data.columns):
            if psm is None:
                psm = input('Please enter a level of psm (such as: 1, 2, 3...): ')
            protein_data = protein_data[protein_data['# PSMs'] >= int(psm)]
        
        # missing values
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
        filtered_labels = [label for ratio, label in zip(self.abundance_ratio, labels) if ratio != 0]
        protein_data = protein_data.dropna(subset=filtered_labels)  
        
        filtered_labels = ['Accession'] + filtered_labels
        protein_data = protein_data[filtered_labels]
        protein_data = protein_data.set_index('Accession', drop=True)
        
        self.protein_raw_data = protein_data.reset_index()
        # outlier
        self.normal_distribution_result = self.normal_distribution_detect(protein_data, sample_size = int((len(filtered_labels)-1)/2))
        print('This is normal distribution detection result, please check it and use outliers detect carefully: ', self.normal_distribution_result)
        self.no_outliers_data = self.outliers_detect(protein_data)
        
        # imputation
        self.no_missing_value_data = self.missing_values_imputation(self.no_outliers_data)
        
        # cv
        self.cv_filter_data = self.cv_filter(self.no_missing_value_data, threshold = cv)
        
        # normalization
        # samplewise
        self.samplewise_normalized_data = self.normalization_samplewise(self.cv_filter_data, method = normalization_samplewise_method)
        # featurewise
        self.samplewise_featurewise_normalized_data = self.normalization_featurewise(self.samplewise_normalized_data, method = normalization_featurewise_method)
        #
        ## analysis
        # glycosylation
        self.proteomic_fc = self.glycosylation_rate(self.samplewise_featurewise_normalized_data)
        # protein-glyco-fc
        if fc == None:
            fc = input('Please enter the fc value: ')
        if pvalue == None:
            pvalue = input('Please enter the pvalue: ')
        #
        proteomic_fc = pd.DataFrame(self.proteomic_fc[[pvalue_type,'fc']]).rename(columns={pvalue_type:'p_p', 'fc':'fc_p'})
        self.glycopeptide_fc = pd.DataFrame(self.gs_data.fc_result[[pvalue_type,'fc']]).rename(columns={pvalue_type:'p_g', 'fc':'fc_g'})
        glycopeptide_fc = self.glycopeptide_fc
        pg_fc = pd.merge(glycopeptide_fc, self.glycopeptide_data, left_index=True, right_index=True, how='left')
        pg_fc = pg_fc.set_index('ProteinID', drop=False)
        pg_fc = pd.merge(proteomic_fc, pg_fc, left_index=True, right_index=True, how='left')
        pg_fc = pg_fc[~pg_fc['fc_g'].isnull()]
        #
        def transform_to_tree_structure(df):
            new_data = []
            for protein, group in df.groupby(df.index):
                p_p = group['p_p'].iloc[0]
                fc_p = group['fc_p'].iloc[0]
                row = [protein, p_p, fc_p]
                for peptide, peptide_group in group.groupby('PeptideSequence'):
                    row.append(peptide)
                    for struct_code, struct_group in peptide_group.groupby('structure_coding'):
                        p_g = struct_group['p_g'].iloc[0]
                        fc_g = struct_group['fc_g'].iloc[0]
                        row.extend([struct_code, p_g, fc_g])
                new_data.append(row)
            new_df = pd.DataFrame(new_data)
            return new_df
        #
        self.pg_fc = pg_fc
        self.protein_up_glyco_up_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']>fc)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']>fc)])
        self.protein_up_glyco_no_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']>fc)&(pg_fc['p_g']>pvalue)])
        self.protein_up_glyco_down_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']>fc)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']<(1/fc))])
        self.protein_no_glyco_up_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']>pvalue)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']>fc)])
        self.protein_no_glyco_no_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']>pvalue)&(pg_fc['p_g']>pvalue)])
        self.protein_no_glyco_down_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']>pvalue)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']<(1/fc))])
        self.protein_down_glyco_up_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']<(1/fc))&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']>fc)])
        self.protein_down_glyco_no_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']<(1/fc))&(pg_fc['p_g']>pvalue)])
        self.protein_down_glyco_down_tree_structure = transform_to_tree_structure(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']<(1/fc))&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']<(1/fc))])
        #  data=module6.protein_up_glyco_up
        self.protein_up_glyco_up = pd.DataFrame(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']>fc)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']>fc)])
        self.protein_up_glyco_up['normalized_fc_g'] = self.protein_up_glyco_up['fc_g'] / self.protein_up_glyco_up['fc_p']
        self.protein_up_glyco_up.insert(4, 'normalized_fc_g', self.protein_up_glyco_up.pop('normalized_fc_g'))
        self.protein_up_glyco_no = pd.DataFrame(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']>fc)&(pg_fc['p_g']>pvalue)])
        self.protein_up_glyco_no['normalized_fc_g'] = self.protein_up_glyco_no['fc_g'] / self.protein_up_glyco_no['fc_p']
        self.protein_up_glyco_no.insert(4, 'normalized_fc_g', self.protein_up_glyco_no.pop('normalized_fc_g'))
        self.protein_up_glyco_down = pd.DataFrame(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']>fc)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']<(1/fc))])
        self.protein_up_glyco_down['normalized_fc_g'] = self.protein_up_glyco_down['fc_g'] / self.protein_up_glyco_down['fc_p']
        self.protein_up_glyco_down.insert(4, 'normalized_fc_g', self.protein_up_glyco_down.pop('normalized_fc_g'))
        self.protein_no_glyco_up = pd.DataFrame(pg_fc[(pg_fc['p_p']>pvalue)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']>fc)])
        self.protein_no_glyco_up['normalized_fc_g'] = self.protein_no_glyco_up['fc_g'] / self.protein_no_glyco_up['fc_p']
        self.protein_no_glyco_up.insert(4, 'normalized_fc_g', self.protein_no_glyco_up.pop('normalized_fc_g'))
        self.protein_no_glyco_no = pd.DataFrame(pg_fc[(pg_fc['p_p']>pvalue)&(pg_fc['p_g']>pvalue)])
        self.protein_no_glyco_no['normalized_fc_g'] = self.protein_no_glyco_no['fc_g'] / self.protein_no_glyco_no['fc_p']
        self.protein_no_glyco_no.insert(4, 'normalized_fc_g', self.protein_no_glyco_no.pop('normalized_fc_g'))
        self.protein_no_glyco_down = pd.DataFrame(pg_fc[(pg_fc['p_p']>pvalue)&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']<(1/fc))])
        self.protein_no_glyco_down['normalized_fc_g'] = self.protein_no_glyco_down['fc_g'] / self.protein_no_glyco_down['fc_p']
        self.protein_no_glyco_down.insert(4, 'normalized_fc_g', self.protein_no_glyco_down.pop('normalized_fc_g'))
        self.protein_down_glyco_up = pd.DataFrame(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']<(1/fc))&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']>fc)])
        self.protein_down_glyco_up['normalized_fc_g'] = self.protein_down_glyco_up['fc_g'] / self.protein_down_glyco_up['fc_p']
        self.protein_down_glyco_up.insert(4, 'normalized_fc_g', self.protein_down_glyco_up.pop('normalized_fc_g'))
        self.protein_down_glyco_no = pd.DataFrame(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']<(1/fc))&(pg_fc['p_g']>pvalue)])
        self.protein_down_glyco_no['normalized_fc_g'] = self.protein_down_glyco_no['fc_g'] / self.protein_down_glyco_no['fc_p']
        self.protein_down_glyco_no.insert(4, 'normalized_fc_g', self.protein_down_glyco_no.pop('normalized_fc_g'))
        self.protein_down_glyco_down = pd.DataFrame(pg_fc[(pg_fc['p_p']<pvalue)&(pg_fc['fc_p']<(1/fc))&(pg_fc['p_g']<pvalue)&(pg_fc['fc_g']<(1/fc))])
        self.protein_down_glyco_down['normalized_fc_g'] = self.protein_down_glyco_down['fc_g'] / self.protein_down_glyco_down['fc_p']
        self.protein_down_glyco_down.insert(4, 'normalized_fc_g', self.protein_down_glyco_down.pop('normalized_fc_g'))
        #
        # glycoprotein_site
        proteomic_fc = pd.DataFrame(self.proteomic_fc[self.proteomic_fc[pvalue_type]<0.05]['fc']).rename(columns={'fc':'fc_p'})
        glycopeptide_fc = pd.DataFrame(self.glycopeptide_fc[self.glycopeptide_fc['p_g']<0.05]['fc_g'])
        fc_gp = pd.merge(glycopeptide_fc, self.glycopeptide_data, left_index=True, right_index=True, how='left')
        fc_gp = fc_gp.set_index('ProteinID', drop=False)
        fc_gp = pd.merge(proteomic_fc, fc_gp, left_index=True, right_index=True, how='left')
        fc_gp = fc_gp[~fc_gp['fc_g'].isnull()]
        #
        self.fc_gp = fc_gp
        self.proteomic_protein_glycosite_count = self.glycoprotein_site(pg_fc).glycoprotein_site_count
        self.proteomic_protein_glycosite_value = self.glycoprotein_site(pg_fc).glycoprotein_site_value
        self.proteomic_glycopeptide_glycosite_count = self.glycopeptide_site(pg_fc).glycopeptide_site_count
        self.proteomic_glycopeptide_glycosite_value = self.glycopeptide_site(pg_fc).glycopeptide_site_value
        #
        proteomic_protein_glycosite_same_direction = []
        for i in set(self.proteomic_protein_glycosite_count['ProteinID']):
            temp_data = self.proteomic_protein_glycosite_count[self.proteomic_protein_glycosite_count['ProteinID']==i]
            temp_data = temp_data.replace(0, np.nan)
            if (temp_data[~temp_data['Up_glycopeptide_count'].isnull()].shape[0] == temp_data.shape[0]) & (temp_data[~temp_data['Down_glycopeptide_count'].isnull()].shape[0] == 0) & (len(temp_data['Protein_type'].unique())==1) & (temp_data['Protein_type'].unique()[0]=='Down'):
                proteomic_protein_glycosite_same_direction.append(i)
            if (temp_data[~temp_data['Down_glycopeptide_count'].isnull()].shape[0] == temp_data.shape[0]) & (temp_data[~temp_data['Up_glycopeptide_count'].isnull()].shape[0] == 0) & (len(temp_data['Protein_type'].unique())==1) & (temp_data['Protein_type'].unique()[0]=='Up'):
                proteomic_protein_glycosite_same_direction.append(i)
        self.proteomic_protein_glycosite_same_direction = fc_gp[fc_gp['ProteinID'].isin(list(proteomic_protein_glycosite_same_direction))].reset_index()
        proteomic_protein_glycosite_different_direction = []
        for i in set(self.proteomic_protein_glycosite_count['ProteinID']):
            temp_data = self.proteomic_protein_glycosite_count[self.proteomic_protein_glycosite_count['ProteinID']==i]
            temp_data = temp_data.replace(0, np.nan)
            if (temp_data[~temp_data['Up_glycopeptide_count'].isnull()].shape[0] != temp_data.shape[0]) & (temp_data[~temp_data['Down_glycopeptide_count'].isnull()].shape[0] != temp_data.shape[0]):
                proteomic_protein_glycosite_different_direction.append(i)
        self.proteomic_protein_glycosite_different_direction = fc_gp[fc_gp['ProteinID'].isin(list(proteomic_protein_glycosite_different_direction))]
        #
        proteomic_glycopeptide_glycosite_same_direction = []
        for i in set(self.proteomic_glycopeptide_glycosite_count['PeptideSequence']):
            temp_data = self.proteomic_glycopeptide_glycosite_count[self.proteomic_glycopeptide_glycosite_count['PeptideSequence']==i]
            temp_data = temp_data.replace(0, np.nan)
            if (temp_data[~temp_data['Up_glycopeptide_count'].isnull()].shape[0] == temp_data.shape[0]) & (temp_data[~temp_data['Down_glycopeptide_count'].isnull()].shape[0] == 0) & (len(temp_data['Protein_type'].unique())==1) & (temp_data['Protein_type'].unique()[0]=='Down'):
                proteomic_glycopeptide_glycosite_same_direction.append(i)
            if (temp_data[~temp_data['Down_glycopeptide_count'].isnull()].shape[0] == temp_data.shape[0]) & (temp_data[~temp_data['Up_glycopeptide_count'].isnull()].shape[0] == 0) & (len(temp_data['Protein_type'].unique())==1) & (temp_data['Protein_type'].unique()[0]=='Up'):
                proteomic_glycopeptide_glycosite_same_direction.append(i)
        self.proteomic_glycopeptide_glycosite_same_direction = fc_gp[fc_gp['PeptideSequence'].isin(list(proteomic_glycopeptide_glycosite_same_direction))]
        proteomic_glycopeptide_glycosite_different_direction = []
        for i in set(self.proteomic_glycopeptide_glycosite_count['PeptideSequence']):
            temp_data = self.proteomic_glycopeptide_glycosite_count[self.proteomic_glycopeptide_glycosite_count['PeptideSequence']==i]
            temp_data = temp_data.replace(0, np.nan)
            if (temp_data[~temp_data['Up_glycopeptide_count'].isnull()].shape[0] != temp_data.shape[0]) & (temp_data[~temp_data['Down_glycopeptide_count'].isnull()].shape[0] != temp_data.shape[0]):
                proteomic_glycopeptide_glycosite_different_direction.append(i)
        self.proteomic_glycopeptide_glycosite_different_direction = fc_gp[fc_gp['PeptideSequence'].isin(list(proteomic_glycopeptide_glycosite_different_direction))]
        #
        # 
        self.data_manager.log_params('StrucGAP_GlycoNetwork', 'proteomic', {'fdr':fdr,
                                                                 'psm':psm,
                                                                 'missing values':'drop any null',
                                                                 'outliers':'Tukey method',
                                                                 'missing values imputation':'KNNImputer',
                                                                 'cv filter threshold':self.threshold,
                                                                 'samplewise normalization':self.normalization_samplewise_method,
                                                                 'featurewise normalization':self.normalization_featurewise_method,
                                                                 })
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'protein raw data', self.protein_raw_data)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'features normal distribution list', self.normal_list)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'features lognormal distribution list', self.lognormal_list)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'no outliers data', self.no_outliers_data)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'no missing value data', self.no_missing_value_data)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'cv filter data', self.cv_filter_data)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'samplewise normalized data', self.samplewise_normalized_data)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'samplewise featurewise normalized data', self.samplewise_featurewise_normalized_data)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic fc', self.proteomic_fc)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein glycosite count', self.proteomic_protein_glycosite_count)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein glycosite value', self.proteomic_protein_glycosite_value)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic glycopeptide glycosite count', self.proteomic_glycopeptide_glycosite_count)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic glycopeptide glycosite value', self.proteomic_glycopeptide_glycosite_value)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein up glyco up', self.protein_up_glyco_up)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein up glyco no', self.protein_up_glyco_no)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein up glyco down', self.protein_up_glyco_down)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein no glyco up', self.protein_no_glyco_up)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein no glyco no', self.protein_no_glyco_no)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein no glyco down', self.protein_no_glyco_down)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein down glyco up', self.protein_down_glyco_up)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein down glyco no', self.protein_down_glyco_no)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein down glyco down', self.protein_down_glyco_down)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein up glyco up tree structure', self.protein_up_glyco_up_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein up glyco no tree structure', self.protein_up_glyco_no_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein up glyco down tree structure', self.protein_up_glyco_down_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein no glyco up tree structure', self.protein_no_glyco_up_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein no glyco no tree structure', self.protein_no_glyco_no_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein no glyco down tree structure', self.protein_no_glyco_down_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein down glyco up tree structure', self.protein_down_glyco_up_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein down glyco no tree structure', self.protein_down_glyco_no_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein down glyco down tree structure', self.protein_down_glyco_down_tree_structure)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein glycosite same direction', self.proteomic_protein_glycosite_same_direction)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic protein glycosite different direction', self.proteomic_protein_glycosite_different_direction)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic glycopeptide glycosite same direction', self.proteomic_glycopeptide_glycosite_same_direction)
        self.data_manager.log_output('StrucGAP_GlycoNetwork', 'proteomic glycopeptide glycosite different direction', self.proteomic_glycopeptide_glycosite_different_direction)
        #
        return self
    
    def glycosyltransferases(self, glycosyltransferases_data_dir=None, data_sheet_name=None):
        """
        Quantifies 161 glycosyltransferases (default) curated from previous literature based on global proteomics.
        
        Parameters:
            glycosyltransferases_data_dir: the file path of glycosyltransferases of interest.
            data_sheet_name: sheet name of uploaded data.
        
        Returns:
            self.glycosyltransferases
            
        Return type:
            dataframe
        
        """
        if glycosyltransferases_data_dir is not None:
            glycosyltransferases_list = pd.read_excel(glycosyltransferases_data_dir, sheet_name=data_sheet_name or 0)
        else:
            glycosyltransferases_list = pd.DataFrame()
            glycosyltransferases_list['Entry'] = ['Q5XPT3','Q8R1T4','Q8VCS3','Q920V1','Q9Z1M7','P23336','Q9EQC0','Q8BHA9','Q8VDB2','Q9EPL0',
                             'Q811B1','Q3U4G3','Q8BG19','Q8BRH0','Q56A06','Q3UV71','Q3TDQ1','P46978','O88829','P97325',
                             'Q76K27','Q64685','Q8K4T1','P70126','Q64692','Q64689','O35696','Q64687','Q9JM95','Q9QYJ1',
                             'Q9R2B6','Q9WUV2','P70277','Q9QZ39','Q91Y74','Q11204','P54751','Q8VIB3','Q8VDX6','O09009',
                             'Q61144','Q8BGQ4','Q8R2R1','Q8BW41','Q91X88','Q8BTP0','Q7TPN3','Q8C2R7','Q9JJQ0','Q64323',
                             'Q8BYB9','P12815','Q8CGY8','Q8VHI3','Q91ZW2','Q765H6','Q8R4G6','Q9D306','Q812F8','Q812G0',
                             'Q10470','Q921V5','P27808','O09008','O09010','Q810K9','Q3UHH8','Q6NVG7','Q8K297','Q9D4M9',
                             'Q8K1B9','Q7TT15','Q9JJ61','Q9D2N8','Q8BVG5','Q8CF93','Q8BGT9','Q921L8','Q6P9S7','Q3V3K7',
                             'E9Q649','Q5JCT0','P97402','Q09324','Q80VA0','Q8C7U7','Q8C102','O08832','P70419','Q6PB93',
                             'O08912','O88819','Q9WTS2','Q11131','Q11127','Q9JL27','O09160','Q8BHC9','Q5F2L2','Q8R507',
                             'Q8CG64','Q9WVL6','Q9ES89','Q9JKV7','P70428','P97464','Q8BYW9','A2AJQ3','Q71B07','P0CW70',
                             'A6X919','Q5DTK1','Q6IQX7','Q6ZQ11','Q8C1F4','Q8BJQ9','O88693','Q9JJ06','P38649','Q9WVK5',
                             'Q9JMK0','Q9JJ04','Q91YY2','Q9Z2Y2','P15535','Q766D5','Q6L8S8','Q09199','Q09200','Q8BWP8',
                             'Q91Z92','Q9JI67','Q9Z0F0','O54905','O54904','Q8VI16','Q8R3I9','Q8K0J2','Q3USF0','Q8BGY6',
                             'Q1RLK6','Q5JCS9','Q9Z222','Q8BHT6','Q8BG28','P58158','P59270','Q9CW73','Q8VDI9','Q6P8H8',
                             'Q3TAE8','Q9DB25','Q8K2A8','Q9DBE8','Q921Q3','Q9D081','Q9D8C3','Q3TZM9','Q3UGP8','Q14BT6',
                             'Q67BJ4']
        
        self.glycosyltransferases = self.proteomic_fc[self.proteomic_fc.index.isin(glycosyltransferases_list['Entry'])]
        
    def glycosidases(self, glycosidases_data_dir=None, data_sheet_name=None):
        """
        Profiles 63 glycosidases annotated in UniProt based on global proteomics.
        
        Parameters:
            glycosidases_data_dir: the file path of glycosidases of interest.
            data_sheet_name: sheet name of uploaded data.
        
        Returns:
            self.glycosidases
            
        Return type:
            dataframe
        
        """
        if glycosidases_data_dir is not None:
            glycosidases_list = pd.read_excel(glycosidases_data_dir, sheet_name=data_sheet_name or 0)
        else:
            glycosidases_list = pd.DataFrame()
            glycosidases_list['Entry'] = ['A2AJ15','O35082','O35632','O35744','P00687','P00688','P10852','P12265','P17439','P18826',
                                          'P20060','P23780','P29416','P39098','P45700','P48441','P48794','P51569','P54818','P70699',
                                          'P82343','Q2HXL6','Q3U4H6','Q61362','Q62010','Q69ZF3','Q69ZQ1','Q6YGZ1','Q80UM7','Q812F3',
                                          'Q8BHN3','Q8K1F9','Q8K2I4','Q8VEI3','Q91WV7','Q91XA9','Q91ZJ9','Q925U4','Q99N32','Q9D7Q1',
                                          'Q9EQQ9','Q05A56','Q8BJT9','Q91Z98','Q99KR8','Q99LJ1','Q9D6Y9','Q9JLT2','Q9QWR8','B2RY83',
                                          'Q7TSH2','Q80YT5','Q8BP56','Q8BVW0','Q8BWJ3','Q8BX80','Q8R242','Q8VC60','Q922Q9','A2RSQ1',
                                          'Q3UPY5','Q6NXH2','Q6P1J0']
            
        self.glycosidases = self.proteomic_fc[self.proteomic_fc.index.isin(glycosidases_list['Entry'])]
    
    def sialyltransferases(self, sialyltransferases_data_dir=None, data_sheet_name=None):
        """
        Covers 21 sialyltransferases from UniProt based on global proteomics.
        
        Parameters:
            sialyltransferases_data_dir: the file path of sialyltransferases of interest.
            data_sheet_name: sheet name of uploaded data.
        
        Returns:
            self.sialyltransferases
            
        Return type:
            dataframe
        
        """
        if sialyltransferases_data_dir is not None:
            sialyltransferases_list = pd.read_excel(sialyltransferases_data_dir, sheet_name=data_sheet_name or 0)
        else:
            sialyltransferases_list = pd.DataFrame()
            sialyltransferases_list['Entry'] = ['O35696','O88829','P54751','P70126','P97325','Q11204','Q60994','Q64685','Q64687','Q64689',
                                                'Q64692','Q76K27','Q8K4T1','Q91Y74','Q9JM95','Q9QYJ1','Q9QZ39','Q9R2B6','Q9WUV2','P70277',
                                                'Q8VIB3']
            
        self.sialyltransferases = self.proteomic_fc[self.proteomic_fc.index.isin(sialyltransferases_list['Entry'])]
        
    def fucosyltransferase(self, fucosyltransferase_data_dir=None, data_sheet_name=None):
        """
        Tracks 51 fucosyltransferases from UniProt based on global proteomics.
        
        Parameters:
            fucosyltransferase_data_dir: the file path of fucosyltransferase of interest.
            data_sheet_name: sheet name of uploaded data.
        
        Returns:
            self.fucosyltransferase
            
        Return type:
            dataframe
        
        """
        if fucosyltransferase_data_dir is not None:
            fucosyltransferase_list = pd.read_excel(fucosyltransferase_data_dir, sheet_name=data_sheet_name or 0)
        else:
            fucosyltransferase_list = pd.DataFrame()
            fucosyltransferase_list['Entry'] = ['O09160','O88819','P97353','Q11127','Q11131','Q5F2L2','Q8VHI3','Q91ZW2','Q9JL27','Q9WTS2',
                                                'A0A1B0GRD2','A6H6C9','Q14AE3','Q91VF0','Q920W3','B2RPT3','B2RV73','E2D0W5','Q32MG3','Q3UYN7',
                                                'Q544B8','Q8BHC9','Q91V20','Q91VB5','Q920V7','Q920V8','Q920V9','Q920W0','Q920W1','Q920W5',
                                                'Q9JL28','A0A1W2P6H4','A0A1W2P7Z1','A2AMC3','B2RRR9','E9PZ15','E9Q686','Q3SWS0','Q3UHF9','Q3UXG7',
                                                'Q3V1L7','Q3ZB27','Q496T8','Q5DU44','Q8CDC9','A0A1W2P844','A2AJ24','B2X2D7','B2X2D8','E9Q154',
                                                'E9Q7A1']
            
        self.fucosyltransferase = self.proteomic_fc[self.proteomic_fc.index.isin(fucosyltransferase_list['Entry'])]
        
    def glycan_binding_protein(self, glycan_binding_protein_data_dir=None, data_sheet_name=None):
        """
        Annotates 276 glycan-binding proteins curated from the GlyCosmos database based on global proteomics.
        
        Parameters:
            glycan_binding_protein_data_dir: the file path of glycan binding protein of interest.
            data_sheet_name: sheet name of uploaded data.
        
        Returns:
            self.glycan_binding_protein
            
        Return type:
            dataframe
        
        """
        if glycan_binding_protein_data_dir is not None:
            glycan_binding_protein_list = pd.read_excel(glycan_binding_protein_data_dir, sheet_name=data_sheet_name or 0)
        else:
            glycan_binding_protein_list = pd.DataFrame()
            glycan_binding_protein_list['Entry'] = ['Q80TR1', 'Q80TS3', 'Q61282', 'Q9WU60', 'Q6A051', 'P21855', 'P35329', 'Q80W49', 'Q61878', 'Q61361',
                                                     'Q8BWY2', 'Q9JL99', 'O88200', 'Q504P2', 'Q149M0', 'Q8VCP9', 'Q91V08', 'Q80XD9', 'Q8VI21', 'Q9D676',
                                                     'Q8C1T8', 'Q9WVF9', 'P0C7M9', 'Q9EPW4', 'Q9QZ15', 'Q9Z2H6', 'Q9R0Q8', 'P70194', 'Q8BNX1', 'Q8VBX4',
                                                     'Q9R007', 'Q9JKF4', 'Q6QLQ4', 'Q8BRU4', 'Q64449', 'P35564', 'P14211', 'Q9DCG2', 'Q8BI06', 'Q5FWI3',
                                                     'O35744', 'Q9CXM0', 'Q8K4Q8', 'O89103', 'P37217', 'Q91V98', 'Q80UW2', 'O70165', 'Q684R7', 'P16045',
                                                     'P97400', 'Q9CQW5', 'P16110', 'Q8K419', 'O54891', 'O54974', 'Q9JL15', 'O08573', 'Q7TPX9', 'Q8VED9',
                                                     'Q9D1U0', 'P98154', 'O88310', 'Q60660', 'Q64329', 'Q60651', 'Q60652', 'Q60653', 'Q60654', 'Q60682',
                                                     'Q0ZUP1', 'P27811', 'P27812', 'Q99JB4', 'P27814', 'Q8VD98', 'Q8CJC7', 'O88713', 'Q3UM83', 'B2KG20',
                                                     'Q5DT36', 'Q8C351', 'P20693', 'Q60767', 'P20917', 'Q63994', 'Q80VA0', 'O54707', 'P55066', 'O70340',
                                                     'O54709', 'Q9EQ09', 'O08852', 'Q7TN88', 'Q2EG98', 'O08912', 'Q6P9S7', 'Q921L8', 'Q8BGT9', 'Q8CF93',
                                                     'Q8BVG5', 'Q9D2N8', 'Q9JJ61', 'Q7TT15', 'Q8K1B9', 'Q6PB93', 'P70419', 'O08832', 'Q8C102', 'Q8C7U7',
                                                     'P15501', 'P58659', 'Q9D8T0', 'Q9D309', 'Q91VU0', 'P97805', 'Q91X88', 'Q8K2C7', 'Q9JL95', 'P35242',
                                                     'Q62028', 'P12246', 'Q80ZE3', 'Q62230', 'P20937', 'P43025', 'Q62059', 'Q9DBH5', 'Q8K0C5', 'Q9D6Y9',
                                                     'Q9DCD0', 'Q8JZZ7', 'Q921V5', 'P27046', 'Q91W89', 'Q8BRK9', 'P47857', 'P12265', 'P20060', 'Q6GQT9',
                                                     'D3Z7H8', 'Q91ZX1', 'Q8VDP6', 'P12960', 'Q9R0H2', 'O54782', 'Q9EP72', 'Q9QZN4', 'Q8K157', 'Q8CFX1',
                                                     'Q91X44', 'P97324', 'Q00612', 'P47856', 'Q8CI94', 'Q9ET01', 'Q9WUB3', 'Q64314', 'P04351', 'P04104',
                                                     'Q9WUA5', 'P70699', 'Q6ZQI3', 'Q9JKF6', 'Q8BHN3', 'Q8BVW0', 'P47968', 'Q9D7G0', 'Q9CS42', 'Q8VEE0',
                                                     'P24807', 'Q99LJ1', 'Q93092', 'P40142', 'Q91YN5', 'P10761', 'Q5YIR8', 'P34927', 'P24721', 'P49300',
                                                     'Q7TSQ1', 'Q91ZW7', 'Q00690', 'O70497', 'Q91VD1', 'Q80ZA0', 'P18337', 'P43137', 'Q08731', 'Q61830',
                                                     'Q01102', 'Q9D0F3', 'Q8VCD3', 'P50404', 'O09037', 'P35230', 'O09049', 'Q9D8G5', 'Q91Y57', 'Q920G3',
                                                     'P59481', 'Q99LC4', 'Q60675', 'Q8CJ91', 'P48759', 'P41317', 'Q63961', 'Q9R0N0', 'Q9Z1E4', 'P17710',
                                                     'Q91ZJ5', 'Q9R062', 'P32037', 'Q3TRM8', 'Q91W97', 'Q9JIF3', 'O08528', 'Q8VCB3', 'P52792', 'P39039',
                                                     'P06802', 'Q9R1E6', 'Q9JM99', 'P29788', 'Q3V188', 'Q3SXB8', 'Q922Q9', 'P58022', 'Q8BRJ4', 'Q9JIG4',
                                                     'Q8C7E7', 'Q7TMB3', 'Q9CW07', 'Q8C767', 'Q99MR9', 'Q8C0L9', 'Q91ZW8', 'P23578', 'Q8R2K1', 'Q9CXB8', 
                                                     'Q9QXD6', 'P12382', 'P53657', 'P06745', 'Q9WUA3', 'Q9JK38', 'Q9CQ60', 'Q9WV38', 'Q91Y97', 'Q07113',
                                                     'O09159', 'Q8CF98', 'P18572', 'Q3UMW8', 'Q91ZW9', 'O70152', 'Q8K2I4', 'Q9R0E2', 'Q9CQ04', 'Q6W3F0',
                                                     'Q3U0K8', 'Q60716', 'Q91YE3', 'Q8BG58', 'Q91YE2', 'Q91UZ4', 'Q9R0E1', 'Q3V1T4', 'Q9R0B9', 'Q8CG70',
                                                     'Q9D136', 'Q8CG71', 'P97467', 'O35386', 'Q60715', 'Q64237']
        self.glycan_binding_protein = self.proteomic_fc[self.proteomic_fc.index.isin(glycan_binding_protein_list['Entry'])]
    
    def phosphorylation(self):
        pass
    
    def convert_accession_to_gene(self, df, protein_col, species=10090):
        """
        An auxiliary function called by other functions to map protein to gene.
        
        Parameters:
            df: the data that needs to be processed.
            protein_col: The column in the data containing protein accession.
        
        Returns:
            Column named 'gene_id'.
            
        Return type:
            series
        
        """
        base_url = "https://version-12-0.string-db.org/api/json/get_string_ids"
        def fetch_gene_id(accession):
            params = {
                "identifiers": accession,
                "species": species,
                "limit": 1,
                "caller_identity": "my_app"  
            }
            try:
                response = requests.get(base_url, params=params, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    if data:
                        preferred_name = data[0].get("preferredName")
                        if preferred_name and "_" in preferred_name:
                            return preferred_name.split("_")[0]
                        return preferred_name
                    else:
                        return None
                else:
                    print(f"请求 {accession} 返回状态码：{response.status_code}")
                    return None
            except Exception as e:
                print(f"请求 {accession} 时发生错误：{e}")
                return None
        df['gene_id'] = df[protein_col].apply(fetch_gene_id)
        return df
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        
        with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_GlycoNetwork.xlsx'), engine='xlsxwriter') as writer:
            self.pg_fc.to_excel(writer, sheet_name='protein_glycan_fc')
            self.protein_raw_data.to_excel(writer, sheet_name='protein_raw_data')
            self.normal_list.to_excel(writer, sheet_name='normal_list')
            self.lognormal_list.to_excel(writer, sheet_name='lognormal_list')
            self.no_outliers_data.to_excel(writer, sheet_name='no_outliers_data')
            self.no_missing_value_data.to_excel(writer, sheet_name='no_missing_value_data'[:31])
            self.cv_filter_data.to_excel(writer, sheet_name='cv_filter_data')
            self.samplewise_normalized_data.to_excel(writer, sheet_name='samplewise_normalized_data'[:31])
            self.samplewise_featurewise_normalized_data.to_excel(writer, sheet_name='samplewise_featurewise_normalized_data'[:31])
            self.proteomic_fc.to_excel(writer, sheet_name='proteomic_fc')
            self.proteomic_protein_glycosite_count.to_excel(writer, sheet_name='protein_glycosite_count'[:31])
            self.proteomic_protein_glycosite_value.to_excel(writer, sheet_name='protein_glycosite_value'[:31])
            self.proteomic_glycopeptide_glycosite_count.to_excel(writer, sheet_name='glycopeptide_glycosite_count'[:31])
            self.proteomic_glycopeptide_glycosite_value.to_excel(writer, sheet_name='glycopeptide_glycosite_value'[:31])
            self.proteomic_protein_glycosite_same_direction.to_excel(writer, sheet_name='protein_glycosite_same_direction'[:31])
            self.proteomic_protein_glycosite_different_direction.to_excel(writer, sheet_name='protein_glycosite_different_direction'[:31])
            self.proteomic_glycopeptide_glycosite_same_direction.to_excel(writer, sheet_name='glycopeptide_glycosite_same_direction'[:31])
            self.proteomic_glycopeptide_glycosite_different_direction.to_excel(writer, sheet_name='glycopeptide_glycosite_different_direction'[:31])
            
            self.glycoprotein_site_count.to_excel(writer, sheet_name='glycoprotein_site_count')
            self.glycoprotein_site_value.to_excel(writer, sheet_name='glycoprotein_site_value')
            
            self.glycopeptide_site_count.to_excel(writer, sheet_name='glycopeptide_site_count')
            self.glycopeptide_site_value.to_excel(writer, sheet_name='glycopeptide_site_value')
            
            self.protein_up_glyco_up.to_excel(writer, sheet_name='protein_up_glyco_up'[:31])
            self.protein_up_glyco_no.to_excel(writer, sheet_name='protein_up_glyco_no'[:31])
            self.protein_up_glyco_down.to_excel(writer, sheet_name='protein_up_glyco_down'[:31])
            self.protein_no_glyco_up.to_excel(writer, sheet_name='protein_no_glyco_up'[:31])
            self.protein_no_glyco_no.to_excel(writer, sheet_name='protein_no_glyco_no'[:31])
            self.protein_no_glyco_down.to_excel(writer, sheet_name='protein_no_glyco_down'[:31])
            self.protein_down_glyco_up.to_excel(writer, sheet_name='protein_down_glyco_up'[:31])
            self.protein_down_glyco_no.to_excel(writer, sheet_name='protein_down_glyco_no'[:31])
            self.protein_down_glyco_down.to_excel(writer, sheet_name='protein_down_glyco_down'[:31])
            
            self.protein_up_glyco_up_tree_structure.to_excel(writer, sheet_name='protein_up_glyco_up_tree_structure'[:31])
            self.protein_up_glyco_no_tree_structure.to_excel(writer, sheet_name='protein_up_glyco_no_tree_structure'[:31])
            self.protein_up_glyco_down_tree_structure.to_excel(writer, sheet_name='protein_up_glyco_down_tree_structure'[:31])
            self.protein_no_glyco_up_tree_structure.to_excel(writer, sheet_name='protein_no_glyco_up_tree_structure'[:31])
            self.protein_no_glyco_no_tree_structure.to_excel(writer, sheet_name='protein_no_glyco_no_tree_structure'[:31])
            self.protein_no_glyco_down_tree_structure.to_excel(writer, sheet_name='protein_no_glyco_down_tree_structure'[:31])
            self.protein_down_glyco_up_tree_structure.to_excel(writer, sheet_name='protein_down_glyco_up_tree_structure'[:31])
            self.protein_down_glyco_no_tree_structure.to_excel(writer, sheet_name='protein_down_glyco_no_tree_structure'[:31])
            self.protein_down_glyco_down_tree_structure.to_excel(writer, sheet_name='protein_down_glyco_down_tree_structure'[:31])
            
            self.glycosyltransferases.to_excel(writer, sheet_name='glycosyltransferases')
            self.glycosidases.to_excel(writer, sheet_name='glycosidases')
            self.sialyltransferases.to_excel(writer, sheet_name='sialyltransferases')
            self.fucosyltransferase.to_excel(writer, sheet_name='fucosyltransferase')
            self.glycan_binding_protein.to_excel(writer, sheet_name='glycan_binding_protein')
    
    

