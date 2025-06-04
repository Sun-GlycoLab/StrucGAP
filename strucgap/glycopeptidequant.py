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
## 糖肽定量分析模块--86
class StrucGAP_GlycoPeptideQuant:
    def __init__(self, gs_data, data_manager, data_type = 'psm_filtered'):
        self.gs_data = gs_data
        self.sample_group = self.gs_data.sample_group
        if hasattr(self.gs_data, 'abundance_ratio'):
            self.abundance_ratio = self.gs_data.abundance_ratio
        #
        if data_type == 'psm_filtered':
            self.data = self.gs_data.data_psm_filtered
        elif data_type == 'cv_filtered':
            self.data = self.gs_data.data_cv_filtered
        elif data_type == 'outliers_filtered':
            self.data = self.gs_data.data_outliers_filtered
        elif data_type == 'data':
            self.data = self.gs_data.data
        else:
            raise ValueError("Invalid data_type provided.")
        
        if str(self.sample_group.index[0]) in self.data.columns:
            #
            self.data_quant = self.data[['PeptideSequence+structure_coding+ProteinID',*map(str, self.sample_group.index)]]
            self.data_quant = self.data_quant.loc[:, ~self.data_quant.T.duplicated(keep='first')].drop_duplicates()
            #
            self.data_quant = self.data_quant.dropna(how='any')
        else:
            self.data_quant = pd.DataFrame()
        
        self.data_manager = data_manager    
        self.data_manager.register_module('StrucGAP_GlycoPeptideQuant', self, {})
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'input_data', {'data_type': data_type})
        
    def statistics(self):
        """
        Provides a comprehensive overview of glycoproteomic identification results.
        
        Parameters:
            None.
        
        Returns:
            self.number_of_identified_spectrum
            self.number_of_identified_peptide
            self.type_of_identified_peptide
            self.number_of_identified_glycopeptide
            self.type_of_identified_glycopeptide
            self.number_of_identified_glycoprotein
            self.type_of_identified_glycoprotein
            
        Return type:
            dataframe
        
        """
        #
        self.number_of_identified_spectrum = self.data['MS2Scan'].count()
        #
        self.number_of_identified_peptide = self.data['PeptideSequence'].value_counts().count()
        #
        self.type_of_identified_peptide = pd.DataFrame(self.data['PeptideSequence'].value_counts().reset_index())
        self.type_of_identified_peptide.columns = ['peptide', 'count']
        #
        self.number_of_identified_glycopeptide = self.data[['PeptideSequence', 'structure_coding']].dropna().drop_duplicates().count()[0]
        #
        self.type_of_identified_glycopeptide = self.data.groupby(['PeptideSequence', 'structure_coding']).size().reset_index(name='Count')
        #
        self.number_of_identified_glycoprotein = self.data[['ProteinID', 'structure_coding']].dropna().drop_duplicates().count()[0]
        #
        self.type_of_identified_glycoprotein = self.data.groupby(['ProteinID', 'structure_coding']).size().reset_index(name='Count')
        #
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'statistics', {})
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'number_of_identified_spectrum', self.number_of_identified_spectrum)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'number_of_identified_peptide', self.number_of_identified_peptide)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'type_of_identified_peptide', self.type_of_identified_peptide)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'number_of_identified_glycopeptide', self.number_of_identified_glycopeptide)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'type_of_identified_glycopeptide', self.type_of_identified_glycopeptide)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'number_of_identified_glycoprotein', self.number_of_identified_glycoprotein)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'type_of_identified_glycoprotein', self.type_of_identified_glycoprotein)
        
        return self
    
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

    def statistics_index(self):
        """
        Computes multiple indices for differential analysis,including:
            FC with significance tests (P values from t-tests and Mann-Whitney U tests), 
            area under the ROC curve (AUC) with P values, 
            feature importance scores from pca, random forest and XGBoost,
            Other indices include ANOVA F-scores and chi-squared scores with corresponding P values.
        
        Parameters:
            None.
        
        Returns:
            self.fc_result.
            self.roc_result.
            self.ml_result (random forest and XGBoost).
            self.pca_result.
            self.anova_result.
            self.chi2_result.
            
        Return type:
            dataframe
        
        """
        # fc
        data_c = self.data_quant[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[0]].index)]]
        data_s = self.data_quant[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[1]].index)]]
        fc_result = pd.DataFrame(index=self.data_quant.index,columns=['fc','pvalue_mannwhitneyu','pvalue_ttest', 'ttest_applicable', 'pvalue_ttest_mannwhitneyu'])
        fc = []
        pvalue_mannwhitneyu = []
        pvalue_ttest = []
        ttest_applicable = []
        pvalue_ttest_mannwhitneyu = []
        
        for i in self.data_quant.index:
            list1 = list(data_c.loc[i])
            list2 = list(data_s.loc[i]) 
            fc_all = []
            for s in list2:
                for c in list1:
                    fc_all.append(s/c)
            fc_all = [x for x in fc_all if not np.isnan(x)]
            fc_all.sort()
            fc.append(self.geometric_mean_of_middle(fc_all,int((self.data_quant.shape[1])/2)))
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
        
        fc_result['fc'] = fc     
        fc_result['pvalue_mannwhitneyu'] = pvalue_mannwhitneyu 
        fc_result['pvalue_ttest'] = pvalue_ttest
        fc_result['ttest_applicable'] = ttest_applicable
        fc_result['pvalue_ttest_mannwhitneyu'] = pvalue_ttest_mannwhitneyu
        self.fc_result = fc_result   
        # roc
        auclist = []
        pvaluelist = []
        roc_result = pd.DataFrame(index=self.data_quant.index,columns=['auc','auc_pvalue'])
        for i in self.data_quant.index:
            control = self.data_quant[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[0]].index)]].loc[i].values
            sample = self.data_quant[[*map(str, self.sample_group[self.sample_group['group']==self.sample_group['group'].unique()[1]].index)]].loc[i].values
            X = np.concatenate([control, sample])
            y = np.concatenate([np.zeros(len(control)), np.ones(len(sample))])
            fpr, tpr, thresholds = roc_curve(y, X)
            roc_auc = auc(fpr, tpr)
            if roc_auc < 0.5:
                roc_auc = 1 - roc_auc
                fpr, tpr = 1 - fpr, 1 - tpr
            auclist.append(roc_auc)
            n1 = len(control)
            n2 = len(sample)
            SE1 = np.sqrt((roc_auc * (1 - roc_auc) + (n1 - 1) * ((roc_auc / (2-roc_auc)) - roc_auc**2) + (n2 - 1) * ((2 * roc_auc**2 / (1+roc_auc)) - roc_auc**2)) / (n1 * n2))
            z1 = (roc_auc - 0.5) / SE1
            p_value1 = 2 * (1 - stats.norm.cdf(abs(z1)))
            pvaluelist.append(p_value1)

        roc_result['auc'] = auclist
        roc_result['auc_pvalue'] = pvaluelist
        self.roc_result = roc_result
        # ml
        # rf xgbclassifier
        feature = self.data_quant[[*map(str, self.sample_group.index)]].T
        label = pd.factorize(self.sample_group['group'])[0]
        importances_rf = pd.DataFrame(index = self.data_quant.index)
        importances_xgb = pd.DataFrame(index = self.data_quant.index)
        random_integers = [random.randint(0, 100) for _ in range(10)]
        for i in random_integers:
            rf = RandomForestClassifier(n_estimators=500, random_state=i, n_jobs=-1,
                                       max_depth=30,min_samples_split=5,min_samples_leaf=1,
                                       max_features='sqrt')
            rf.fit(feature, label)
            importances_rf = pd.concat([importances_rf,pd.DataFrame(rf.feature_importances_,index = self.data_quant.index,columns=[r'rf_score_%s'%i])],axis=1)
            xgb = XGBClassifier(n_estimators=500,seed=i,learning_rate=0.01)
            xgb.fit(np.array(feature), label)
            importances_xgb = pd.concat([importances_xgb,pd.DataFrame(xgb.feature_importances_,index = self.data_quant.index,columns=[r'xgb_score_%s'%i])],axis=1)
        
        importances_rf['randomforest_features_importance_means'] = importances_rf.mean(axis=1)
        importances_xgb['xgbclassifier_features_importance_means'] = importances_xgb.mean(axis=1)
        self.ml_result = pd.DataFrame(pd.concat([importances_rf['randomforest_features_importance_means'], importances_xgb['xgbclassifier_features_importance_means']],axis=1))
        self.ml_result = self.ml_result.sort_values(by='randomforest_features_importance_means', ascending=False)
        # pca features importance
        scaler = StandardScaler()
        feature_scaled = scaler.fit_transform(feature.values)
        pca = PCA()
        feature_pca = pca.fit_transform(feature_scaled)
        loading_matrix = np.abs(pca.components_[:2])
        explained_variance_ratio = pca.explained_variance_ratio_[:2]
        self.pca_result = pd.DataFrame(index = self.data_quant.index, columns = ['pca_features_importance'], data = (loading_matrix.T * explained_variance_ratio).sum(axis=1))
        self.pca_result = self.pca_result.sort_values(by='pca_features_importance', ascending=False)
        # anova
        fval = f_classif(feature, label)
        print('Please apply this function with caution, as your data is most likely not normally distributed!')
        self.anova_result = pd.DataFrame({'f score': fval[0],'anova_pvalue': fval[1]}, index=self.data_quant.index)
        # chi2
        chi_scores = chi2(feature, label)
        self.chi2_result = pd.DataFrame({'chi2 score': chi_scores[0],'chi2_pvalue': chi_scores[1]}, index=self.data_quant.index)
        #
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'statistics_index', {})
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'fc_result', self.fc_result)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'roc_result', self.roc_result)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'ml_result', self.ml_result)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'pca_result', self.pca_result)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'anova_result', self.anova_result)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'chi2_result', self.chi2_result)
        
        return self
    
    def identify_core_structure(self, data): 
        """An auxiliary function called by other functions to identify core structure."""
        temp_data = data[(data['GlycanComposition']!='N2H2')&(data['GlycanComposition']!='N2H2F1')]
        temp_data = temp_data[['structure_coding','Corefucose', 'Bisection']]
        core_structure_list = []
        for i, row in temp_data.iterrows():
            if row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD1'
            elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD1dcbB5'
            elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD2dD1'
            elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD2dD1dcbB5'
            else:
                continue
            core_structure_list.append([i, structure])
        core_structure = pd.DataFrame(core_structure_list, columns=['PeptideSequence+structure_coding+ProteinID', 'Core_structure']).set_index('PeptideSequence+structure_coding+ProteinID', drop=True)
    
        return core_structure
    
    def differential_analysis(self, index = 'fc', pvalue = 0.05, fc = 1.5, pvalue_type='pvalue_ttest'):
        """
        Identifies significantly regulated glycopeptides under user-defined thresholds.
        
        Parameters:
            index: statistic indices used for differential analysis ('fc', 'roc', 'ml', 'pca', 'anova', 'chi2').
            pvalue: P value used for differential analysis.
            fc: FC used for differential analysis.
            pvalue_type: 'pvalue_ttest', 'pvalue_mannwhitneyu' or 'pvalue_ttest_mannwhitneyu' used for FC-based differential analysis.
        
        Returns:
            self.differential_analysis_data.
            self.up_data (upregulated IGPs).
            self.down_data (downregulated IGPs).
            self.differential_analysis_overview.
            self.differential_analysis_structure_coding (the number of different glycan structures in upregulated and downregulated glycopeptides).
            self.differential_analysis_structure_coding_no_oligomannose. (the number of different glycan structures in upregulated and downregulated glycopeptides except oligo mannose).
            self.differential_analysis_core_structure (the number of different core structures in upregulated and downregulated glycopeptides).
            self.differential_analysis_glycan_type (the number of different glycan types in upregulated and downregulated glycopeptides).
            self.differential_analysis_branches_structure (the number of different branch structures in upregulated and downregulated glycopeptides).
            self.differential_analysis_branches_count (the number of different branch counts in upregulated and downregulated glycopeptides).
            self.differential_analysis_sialicacid_count (the number of different sialic acid counts in upregulated and downregulated glycopeptides).
            self.differential_analysis_fucose_count (the number of different fucose counts in upregulated and downregulated glycopeptides).
            self.differential_analysis_sialicacid_structure (the number of different sialic acid structures in upregulated and downregulated glycopeptides).
            self.differential_analysis_fucose_structure (the number of different fucose structures in upregulated and downregulated glycopeptides).
            self.differential_analysis_lacdinac (the number of different LacdiNAc in upregulated and downregulated glycopeptides).
            self.differential_analysis_fucosylated_type (the number of different fucosylated types in upregulated and downregulated glycopeptides).
            self.differential_analysis_acgc (the number of different sialylated types in upregulated and downregulated glycopeptides).
            
        Return type:
            dataframe
        
        """
        index = input("Please enter 1 or more statistic indices used for differential analysis ('fc', 'roc', 'ml', 'pca', 'anova', 'chi2'): ")
        index = [x.strip() for x in index.split(',')]
        print(f"Using '{index}' as the input.")
        index_list = {
            'fc': self.fc_result,
            'roc': self.roc_result,
            'ml': self.ml_result,
            'pca': self.pca_result,
            'anova': self.anova_result,
            'chi2': self.chi2_result
        }
        for i in index:
            if i not in index_list:
                raise ValueError(f"Invalid index '{i}'. Available options are {list(index_list.keys())}")
        
        auxiliary_filter_index = input("Please enter other auxiliary filtering index ( 'roc', 'ml', 'pca', 'anova', 'chi2') for further filtering after fc filtering, the primary filter index has been set to 'fc' by default: ")
        auxiliary_filter_index = [x.strip() for x in auxiliary_filter_index.split(',')]
        print(f"Using '{auxiliary_filter_index}' as the input.")
        auxiliary_filter_index_list = {
            'roc': "auc_pvalue<0.05",
            'ml': "randomforest_features_importance_means!=0",
            'pca': "chi2_pvalue<0.05",
            'anova': "anova_pvalue<0.05",
            'chi2': "chi2_pvalue<0.05",
            'no':'null'
            }

        result = self.data  # data1 = test1.ml_result
        for i in index:
            result = pd.merge(result, index_list[i], left_index=True, right_index=True, how='left')
        self.differential_analysis_data = result
        #
        user_pvalue = input(f"Please enter the pvalue used for differential analysis (default pvalue was {pvalue}): ")
        if user_pvalue:
            try:
                pvalue = float(user_pvalue)
            except ValueError:
                print("Invaild input, the pvalue was set to 0.05")
                pvalue = 0.05
        user_fc = input(f"Please enter the fc used for differential analysis (default fc was {fc}): ")
        if user_fc:
            try:
                fc = float(user_fc)
            except ValueError:
                print("Invaild input, the fc was set to 1.5")
                fc = 1.5
        self.differential_analysis_data = self.differential_analysis_data[self.differential_analysis_data[pvalue_type]<pvalue]
        if 'no' not in auxiliary_filter_index:
            for i in auxiliary_filter_index:
                self.differential_analysis_data = self.differential_analysis_data.query(auxiliary_filter_index_list[i])
        #
        self.up_data = self.differential_analysis_data[self.differential_analysis_data['fc'] > fc]
        self.down_data = self.differential_analysis_data[self.differential_analysis_data['fc'] < (1/fc)]
        #
        up_data_overview = pd.DataFrame()
        up_data_overview['item'] = ['Glycan', 'Glycanpeptide', 'Glycanprotein', 'Glycosite']
        temp_data1 = self.up_data[['ProteinID','structure_coding']]
        temp_data1['ProteinID'] = temp_data1['ProteinID'].str.split(';')
        temp_data1 = temp_data1.explode(['ProteinID'])
        temp_data2 = pd.DataFrame(self.up_data[['ProteinID', 'Glycosite_Position']])
        temp_data2 = temp_data2[~temp_data2['Glycosite_Position'].isnull()]
        temp_data2['ProteinID'] = temp_data2['ProteinID'].str.split(';')
        temp_data2['Glycosite_Position'] = temp_data2['Glycosite_Position'].str.split(';')
        temp_data2 = temp_data2.explode(['ProteinID', 'Glycosite_Position'])
        up_data_overview['Up_data_item_count'] = [self.up_data['structure_coding'].nunique(),
                                                  (self.up_data['PeptideSequence'] + self.up_data['structure_coding']).nunique(),
                                                  (temp_data1['ProteinID'] + temp_data1['structure_coding']).nunique(),
                                                  (temp_data2['ProteinID'] + temp_data2['Glycosite_Position']).nunique()
                                                  ]
        down_data_overview = pd.DataFrame()
        down_data_overview['item'] = ['Glycan', 'Glycanpeptide', 'Glycanprotein', 'Glycosite']
        temp_data1 = self.down_data[['ProteinID','structure_coding']]
        temp_data1['ProteinID'] = temp_data1['ProteinID'].str.split(';')
        temp_data1 = temp_data1.explode(['ProteinID'])
        temp_data2 = pd.DataFrame(self.down_data[['ProteinID', 'Glycosite_Position']])
        temp_data2 = temp_data2[~temp_data2['Glycosite_Position'].isnull()]
        temp_data2['ProteinID'] = temp_data2['ProteinID'].str.split(';')
        temp_data2['Glycosite_Position'] = temp_data2['Glycosite_Position'].str.split(';')
        temp_data2 = temp_data2.explode(['ProteinID', 'Glycosite_Position'])
        down_data_overview['Down_data_item_count'] = [self.down_data['structure_coding'].nunique(),
                                                  (self.down_data['PeptideSequence'] + self.down_data['structure_coding']).nunique(),
                                                  (temp_data1['ProteinID'] + temp_data1['structure_coding']).nunique(),
                                                  (temp_data2['ProteinID'] + temp_data2['Glycosite_Position']).nunique()
                                                  ]
        self.differential_analysis_overview = pd.concat([up_data_overview, down_data_overview], axis=1)
        #
        up_data_structure_coding = pd.DataFrame(self.up_data['structure_coding'].value_counts())
        down_data_structure_coding = pd.DataFrame(self.down_data['structure_coding'].value_counts())
        self.differential_analysis_structure_coding = pd.concat([up_data_structure_coding, down_data_structure_coding], axis=1).reset_index()
        self.differential_analysis_structure_coding.columns = ['structure_coding', 'Up_data_structure_coding', 'Down_data_structure_coding']
        df1 = self.differential_analysis_structure_coding.copy()
        back = pd.DataFrame(self.data['structure_coding'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_structure_coding']   / df1['structure_coding'].map(back)
        df1['Down_ratio'] = df1['Down_data_structure_coding'] / df1['structure_coding'].map(back)
        self.differential_analysis_structure_coding = df1
        #
        up_data_structure_coding = pd.DataFrame(self.up_data[self.up_data['Glycan_type'] != 'Oligo mannose']['structure_coding'].value_counts())
        down_data_structure_coding = pd.DataFrame(self.down_data[self.down_data['Glycan_type'] != 'Oligo mannose']['structure_coding'].value_counts())
        self.differential_analysis_structure_coding_no_oligomannose = pd.concat([up_data_structure_coding, down_data_structure_coding], axis=1).reset_index()
        self.differential_analysis_structure_coding_no_oligomannose.columns = ['structure_coding_no_oligomannose', 'Up_data_structure_coding', 'Down_data_structure_coding']
        df1 = self.differential_analysis_structure_coding_no_oligomannose.copy()
        back = pd.DataFrame(self.data[self.data['Glycan_type'] != 'Oligo mannose']['structure_coding'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_structure_coding']   / df1['structure_coding_no_oligomannose'].map(back)
        df1['Down_ratio'] = df1['Down_data_structure_coding'] / df1['structure_coding_no_oligomannose'].map(back)
        self.differential_analysis_structure_coding_no_oligomannose = df1
        #
        up_data_core_structure = self.identify_core_structure(self.up_data)['Core_structure'].value_counts()
        down_data_core_structure = self.identify_core_structure(self.down_data)['Core_structure'].value_counts()
        self.differential_analysis_core_structure = pd.concat([up_data_core_structure, down_data_core_structure], axis=1).reset_index()
        self.differential_analysis_core_structure.columns = ['Core_structure', 'Up_data_Core_structure_count', 'Down_data_Core_structure_count']
        df1 = self.differential_analysis_core_structure.copy()
        back = pd.DataFrame(self.identify_core_structure(self.data)['Core_structure'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_Core_structure_count']   / df1['Core_structure'].map(back)
        df1['Down_ratio'] = df1['Down_data_Core_structure_count'] / df1['Core_structure'].map(back)
        self.differential_analysis_core_structure = df1
        #
        up_data_glycan_type = pd.DataFrame(self.up_data['Glycan_type'].value_counts())
        down_data_glycan_type = pd.DataFrame(self.down_data['Glycan_type'].value_counts())
        self.differential_analysis_glycan_type = pd.concat([up_data_glycan_type, down_data_glycan_type], axis=1).reset_index()
        self.differential_analysis_glycan_type.columns = ['Glycan_type', 'Up_data_Glycan_type_count', 'Down_data_Glycan_type_count']
        df1 = self.differential_analysis_glycan_type.copy()
        back = pd.DataFrame(self.data['Glycan_type'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_Glycan_type_count']   / df1['Glycan_type'].map(back)
        df1['Down_ratio'] = df1['Down_data_Glycan_type_count'] / df1['Glycan_type'].map(back)
        self.differential_analysis_glycan_type = df1
        #
        up_data_branches_structure = pd.DataFrame(self.up_data['Branches'].value_counts().reset_index())
        up_data_branches_structure.columns = ['Branches', 'Up_data_Branches_count']
        count_dict = {}
        for index, row in up_data_branches_structure.iterrows():
            branch_list = literal_eval(row['Branches'])
            count = row['Up_data_Branches_count']
            for branch in branch_list:
                if branch:
                    if branch in count_dict:
                        count_dict[branch] += count
                    else:
                        count_dict[branch] = count
        up_data_branches_structure = pd.DataFrame(list(count_dict.items()), columns=['Branches', 'Up_data_Branches_count'])
        up_data_branches_structure = up_data_branches_structure.set_index('Branches',drop=True)
        down_data_branches_structure = pd.DataFrame(self.down_data['Branches'].value_counts().reset_index())
        down_data_branches_structure.columns = ['Branches', 'Down_data_Branches_count']
        count_dict = {}
        for index, row in down_data_branches_structure.iterrows():
            branch_list = literal_eval(row['Branches'])
            count = row['Down_data_Branches_count']
            for branch in branch_list:
                if branch:
                    if branch in count_dict:
                        count_dict[branch] += count
                    else:
                        count_dict[branch] = count
        down_data_branches_structure = pd.DataFrame(list(count_dict.items()), columns=['Branches', 'Down_data_Branches_count'])
        down_data_branches_structure = down_data_branches_structure.set_index('Branches',drop=True)
        self.differential_analysis_branches_structure = pd.concat([up_data_branches_structure, down_data_branches_structure], axis=1).reset_index()
        df1 = self.differential_analysis_branches_structure.copy()
        back = pd.DataFrame(self.data['Branches'].value_counts().reset_index())
        back.columns = ['Branches', 'Branches_count']
        count_dict = {}
        for index, row in back.iterrows():
            branch_list = literal_eval(row['Branches'])
            count = row['Branches_count']
            for branch in branch_list:
                if branch:
                    if branch in count_dict:
                        count_dict[branch] += count
                    else:
                        count_dict[branch] = count
        back = pd.DataFrame(list(count_dict.items()), columns=['index', 'Branches'])
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_Branches_count'] / df1['Branches'].map(back)
        df1['Down_ratio'] = df1['Down_data_Branches_count'] / df1['Branches'].map(back)
        self.differential_analysis_branches_structure = df1
        #
        up_data_branches_count = pd.DataFrame(self.up_data['BranchNumber'].value_counts())
        down_data_branches_count = pd.DataFrame(self.down_data['BranchNumber'].value_counts())
        self.differential_analysis_branches_count = pd.concat([up_data_branches_count, down_data_branches_count], axis=1).reset_index()
        self.differential_analysis_branches_count.columns = ['BranchNumber', 'Up_data_BranchNumber_count', 'Down_data_BranchNumber_count']
        df1 = self.differential_analysis_branches_count.copy()
        back = pd.DataFrame(self.data['BranchNumber'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_BranchNumber_count'] / df1['BranchNumber'].map(back)
        df1['Down_ratio'] = df1['Down_data_BranchNumber_count'] / df1['BranchNumber'].map(back)
        self.differential_analysis_branches_count = df1
        #
        up_data_sialicacid_count = pd.DataFrame(self.up_data['structure_coding'].str.count('3').value_counts())
        down_data_sialicacid_count = pd.DataFrame(self.down_data['structure_coding'].str.count('3').value_counts())
        self.differential_analysis_sialicacid_count = pd.concat([up_data_sialicacid_count, down_data_sialicacid_count], axis=1).reset_index()
        self.differential_analysis_sialicacid_count.columns = ['Sialicacid_count', 'Up_data_Number', 'Down_data_Number']
        df1 = self.differential_analysis_sialicacid_count.copy()
        back = pd.DataFrame(self.data['structure_coding'].str.count('3').value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_Number'] / df1['Sialicacid_count'].map(back)
        df1['Down_ratio'] = df1['Down_data_Number'] / df1['Sialicacid_count'].map(back)
        self.differential_analysis_sialicacid_count = df1
        #
        up_data_fucose_count = pd.DataFrame(self.up_data['structure_coding'].str.count('5').value_counts())
        down_data_fucose_count = pd.DataFrame(self.down_data['structure_coding'].str.count('5').value_counts())
        self.differential_analysis_fucose_count = pd.concat([up_data_fucose_count, down_data_fucose_count], axis=1).reset_index()
        self.differential_analysis_fucose_count.columns = ['Fucose_count', 'Up_data_Number', 'Down_data_Number']
        df1 = self.differential_analysis_fucose_count.copy()
        back = pd.DataFrame(self.data['structure_coding'].str.count('5').value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_Number'] / df1['Fucose_count'].map(back)
        df1['Down_ratio'] = df1['Down_data_Number'] / df1['Fucose_count'].map(back)
        self.differential_analysis_fucose_count = df1
        #
        up_data_sialicacid_structure = pd.DataFrame(self.up_data[self.up_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts())
        down_data_sialicacid_structure = pd.DataFrame(self.down_data[self.down_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts())
        self.differential_analysis_sialicacid_structure = pd.concat([up_data_sialicacid_structure, down_data_sialicacid_structure], axis=1).reset_index()
        self.differential_analysis_sialicacid_structure.columns = ['structure_coding', 'Up_data_sialicacid_structure_count', 'Down_data_sialicacid_structure_count']
        df1 = self.differential_analysis_sialicacid_structure.copy()
        back = pd.DataFrame(self.data[self.data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_sialicacid_structure_count'] / df1['structure_coding'].map(back)
        df1['Down_ratio'] = df1['Down_data_sialicacid_structure_count'] / df1['structure_coding'].map(back)
        self.differential_analysis_sialicacid_structure = df1
        #
        up_data_fucose_structure = pd.DataFrame(self.up_data[self.up_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts())
        down_data_fucose_structure = pd.DataFrame(self.down_data[self.down_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts())
        self.differential_analysis_fucose_structure = pd.concat([up_data_fucose_structure, down_data_fucose_structure], axis=1).reset_index()
        self.differential_analysis_fucose_structure.columns = ['structure_coding', 'Up_data_fucose_structure_count', 'Down_data_fucose_structure_count']
        df1 = self.differential_analysis_fucose_structure.copy()
        back = pd.DataFrame(self.data[self.data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().reset_index())
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_fucose_structure_count'] / df1['structure_coding'].map(back)
        df1['Down_ratio'] = df1['Down_data_fucose_structure_count'] / df1['structure_coding'].map(back)
        self.differential_analysis_fucose_structure = df1
        #
        up_data_lacdinac = pd.DataFrame(self.up_data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts())
        down_data_lacdinac = pd.DataFrame(self.down_data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts())
        self.differential_analysis_lacdinac = pd.concat([up_data_lacdinac, down_data_lacdinac], axis=1).reset_index()
        self.differential_analysis_lacdinac.columns = ['lacdinac', 'Up_data_lacdinac', 'Down_data_lacdinac']
        self.differential_analysis_lacdinac = self.differential_analysis_lacdinac[self.differential_analysis_lacdinac['lacdinac']!=' ']
        df1 = self.differential_analysis_lacdinac.copy()
        back = pd.DataFrame(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().reset_index())
        back = back[back['lacdinac']!=' ']
        back = back.set_index('lacdinac')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_lacdinac'] / df1['lacdinac'].map(back)
        df1['Down_ratio'] = df1['Down_data_lacdinac'] / df1['lacdinac'].map(back)
        self.differential_analysis_lacdinac = df1
        #
        up_data_fucosylated_type = pd.DataFrame(self.up_data['fucosylated type'].value_counts())
        down_data_fucosylated_type = pd.DataFrame(self.down_data['fucosylated type'].value_counts())
        self.differential_analysis_fucosylated_type = pd.concat([up_data_fucosylated_type, down_data_fucosylated_type], axis=1).reset_index()
        self.differential_analysis_fucosylated_type.columns = ['fucosylated type', 'Up_data_fucosylated_type', 'Down_data_fucosylated_type']
        self.differential_analysis_fucosylated_type = self.differential_analysis_fucosylated_type[self.differential_analysis_fucosylated_type['fucosylated type']!=' ']
        df1 = self.differential_analysis_fucosylated_type.copy()
        back = pd.DataFrame(self.data['fucosylated type'].value_counts().reset_index())
        back = back[back['index']!=' ']
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_fucosylated_type'] / df1['fucosylated type'].map(back)
        df1['Down_ratio'] = df1['Down_data_fucosylated_type'] / df1['fucosylated type'].map(back)
        self.differential_analysis_fucosylated_type = df1
        #
        up_data_acgc = pd.DataFrame(self.up_data['Ac/Gc'].value_counts())
        down_data_acgc = pd.DataFrame(self.down_data['Ac/Gc'].value_counts())
        self.differential_analysis_acgc = pd.concat([up_data_acgc, down_data_acgc], axis=1).reset_index()
        self.differential_analysis_acgc.columns = ['Ac/Gc', 'Up_data_Ac/Gc', 'Down_data_Ac/Gc']
        self.differential_analysis_acgc = self.differential_analysis_acgc[self.differential_analysis_acgc['Ac/Gc']!=' ']
        df1 = self.differential_analysis_acgc.copy()
        back = pd.DataFrame(self.data['Ac/Gc'].value_counts().reset_index())
        back = back[back['index']!=' ']
        back = back.set_index('index')[back.columns[1]].to_dict()
        df1['Up_ratio'] = df1['Up_data_Ac/Gc'] / df1['Ac/Gc'].map(back)
        df1['Down_ratio'] = df1['Down_data_Ac/Gc'] / df1['Ac/Gc'].map(back)
        self.differential_analysis_acgc = df1
        #
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'differential_analysis', {'index':index, 'pvalue':pvalue, 'fc':fc})
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_overview)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_core_structure)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_branches_count)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_sialicacid_structure)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_fucose_structure)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'differential_analysis', self.differential_analysis_acgc)
        return self
    
    def identify_glycan_composition(self, data):  
        """An auxiliary function called by other functions to identify glycan composition."""
        temp_data = data[['GlycanComposition', 'structure_coding']] 
        glycancomposition = {}
        for i in temp_data['GlycanComposition'].unique(): 
            glycancomposition[i] = temp_data[temp_data['GlycanComposition'] == i]['structure_coding'].value_counts().max()
        GlycanComposition_rank = pd.DataFrame(list(glycancomposition.items()), columns=['GlycanComposition', 'GlycanComposition_count'])
        GlycanComposition_rank = GlycanComposition_rank.sort_values(by='GlycanComposition_count', ascending=False)     
        
        return GlycanComposition_rank
    
    def threshold_variation_analysis(self, statistic_index=None, pvalue_type='pvalue_ttest'):
        """
        Explores how glycan substructural regulation varies across a range of index thresholds.
        
        Parameters:
            statistic_index: the statistic index ('fc', 'roc', 'ml', 'pca', 'anova', 'chi2') would like to use for threshold variation analysis.
            pvalue_type: 'pvalue_ttest', 'pvalue_mannwhitneyu' or 'pvalue_ttest_mannwhitneyu'.
        
        Returns:
            self.result_branches_structure
            self.result_branches_structure_up
            self.result_branches_structure_down
            self.result_branches_structure_split
            self.result_branches_structure_ratio
            self.result_core_structure
            self.result_core_structure_up
            self.result_core_structure_down
            self.result_core_structure_split
            self.result_core_structure_ratio
            self.result_glycan_type
            self.result_glycan_type_up
            self.result_glycan_type_down
            self.result_glycan_type_split
            self.result_glycan_type_ratio
            self.result_branches_count
            self.result_branches_count_up
            self.result_branches_count_down
            self.result_branches_count_split
            self.result_branches_count_ratio
            self.result_glycan_composition
            self.result_glycan_composition_up
            self.result_glycan_composition_down
            self.result_glycan_composition_split
            self.result_glycan_composition_ratio
            self.result_lacdinac
            self.result_lacdinac_up
            self.result_lacdinac_down
            self.result_lacdinac_split
            self.result_lacdinac_ratio
            self.result_fucosylated_type
            self.result_fucosylated_type_up
            self.result_fucosylated_type_down
            self.result_fucosylated_type_split
            self.result_fucosylated_type_ratio
            self.result_acgc
            self.result_acgc_up
            self.result_acgc_down
            self.result_acgc_split
            self.result_acgc_ratio
                
        Return type:
            dataframe
        
        """
        if statistic_index == None:
            statistic_index = input("Please enter the statistic index ('fc', 'roc', 'ml', 'pca', 'anova', 'chi2') you would like to use for threshold variation analysis: ")
            expected_options = ['fc', 'roc', 'ml', 'pca', 'anova', 'chi2']
            matches = get_close_matches(statistic_index, expected_options, n=1, cutoff=0.5)
            if matches:
                statistic_index = matches[0]
                print(f"Using '{statistic_index}' as the input.")
            else:
                print("No close match found. Using 'fc' as the input.")
                statistic_index = 'fc'
                
        if pvalue_type is None:
            pvalue_type = input("Please enter the pvalue type ('pvalue_ttest_mannwhitneyu', 'pvalue_ttest') you would like to use for threshold variation analysis: ")
            
        if not statistic_index:
            raise ValueError("Invalid statistic index provided, the statistic index was set to 'fc'. ") 
            statistic_index = 'fc'
        else:
            statistic_index = statistic_index
            
        index_list = {
            'fc': self.fc_result,
            'roc': self.roc_result,
            'ml': self.ml_result,
            'pca': self.pca_result,
            'anova': self.anova_result,
            'chi2': self.chi2_result
        }
        #
        result = pd.merge(self.data, index_list[statistic_index], left_index=True, right_index=True, how='left')
        #
        if statistic_index == 'fc':
            result = result[result[pvalue_type]<0.05]
            
            unique_branches = set()
            for branches_str in result['Branches'].value_counts().index:
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            result_branches_structure = pd.DataFrame(index=unique_branches_list)
            
            result_core_structure = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            result_lacdinac = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            result_fucosylated_type = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            result_acgc = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)
            
            result_branches_structure_up = pd.DataFrame(index=unique_branches_list)
            result_core_structure_up = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type_up = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count_up = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition_up = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            result_lacdinac_up = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            result_fucosylated_type_up = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            result_acgc_up = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)
            
            # result_branches_structure_up_ratio = pd.DataFrame(index=unique_branches_list)
            # result_core_structure_up_ratio = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            # result_glycan_type_up_ratio = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            # result_branches_count_up_ratio = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            # result_glycan_composition_up_ratio = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            # result_lacdinac_up_ratio = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            # result_fucosylated_type_up_ratio = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            # result_acgc_up_ratio = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)

            result_branches_structure_down = pd.DataFrame(index=unique_branches_list)
            result_core_structure_down = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type_down = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count_down = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition_down = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            result_lacdinac_down = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            result_fucosylated_type_down = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            result_acgc_down = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)
            
            # result_branches_structure_down_ratio = pd.DataFrame(index=unique_branches_list)
            # result_core_structure_down_ratio = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            # result_glycan_type_down_ratio = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            # result_branches_count_down_ratio = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            # result_glycan_composition_down_ratio = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            # result_lacdinac_down_ratio = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            # result_fucosylated_type_down_ratio = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            # result_acgc_down_ratio = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)

            result_branches_structure_split = pd.DataFrame(index=unique_branches_list)
            result_core_structure_split = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type_split = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count_split = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition_split = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            result_lacdinac_split = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            result_fucosylated_type_split = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            result_acgc_split = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)
            
            result_branches_structure_ratio = pd.DataFrame(index=unique_branches_list)
            result_core_structure_ratio = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type_ratio = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count_ratio = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition_ratio = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            result_lacdinac_ratio = pd.DataFrame(index=self.data[self.data['lacdinac'] != ' ']['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts().index)
            result_fucosylated_type_ratio = pd.DataFrame(index=self.data[self.data['fucosylated type'] != ' ']['fucosylated type'].value_counts().index)
            result_acgc_ratio = pd.DataFrame(index=self.data[self.data['Ac/Gc'] != ' ']['Ac/Gc'].value_counts().index)

            for i in [1.2, 1.5, 2, 2.5, 3]:
                unique_branches = set()
                for branches_str in result['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                
                temp_data = pd.DataFrame(result[(result['fc']>i)|(result['fc']<(1/i))]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure = pd.concat([result_branches_structure, resultlist], axis=1)
                result_core_structure = pd.concat([result_core_structure, self.identify_core_structure(result[(result['fc']>i)|(result['fc']<(1/i))])['Core_structure'].value_counts()], axis=1)
                result_glycan_type = pd.concat([result_glycan_type, result[(result['fc']>i)|(result['fc']<(1/i))]['Glycan_type'].value_counts()], axis=1)
                result_branches_count = pd.concat([result_branches_count, result[(result['fc']>i)|(result['fc']<(1/i))]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition = pd.concat([result_glycan_composition, result[(result['fc']>i)|(result['fc']<(1/i))]['GlycanComposition'].value_counts()], axis=1)
                result_lacdinac = pd.concat([result_lacdinac, result[(result['fc']>i)|(result['fc']<(1/i))]['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts()], axis=1)
                result_fucosylated_type = pd.concat([result_fucosylated_type, result[(result['fc']>i)|(result['fc']<(1/i))]['fucosylated type'].value_counts()], axis=1)
                result_acgc = pd.concat([result_acgc, result[(result['fc']>i)|(result['fc']<(1/i))]['Ac/Gc'].value_counts()], axis=1)
                
                temp_data = pd.DataFrame(result[result['fc']>i]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure_up = pd.concat([result_branches_structure_up, resultlist], axis=1)
                result_core_structure_up = pd.concat([result_core_structure_up, self.identify_core_structure(result[result['fc']>i])['Core_structure'].value_counts()], axis=1)
                result_glycan_type_up = pd.concat([result_glycan_type_up, result[result['fc']>i]['Glycan_type'].value_counts()], axis=1)
                result_branches_count_up = pd.concat([result_branches_count_up, result[result['fc']>i]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition_up = pd.concat([result_glycan_composition_up, result[result['fc']>i]['GlycanComposition'].value_counts()], axis=1)
                result_lacdinac_up = pd.concat([result_lacdinac_up, result[result['fc']>i]['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts()], axis=1)
                result_fucosylated_type_up = pd.concat([result_fucosylated_type_up, result[result['fc']>i]['fucosylated type'].value_counts()], axis=1)
                result_acgc_up = pd.concat([result_acgc_up, result[result['fc']>i]['Ac/Gc'].value_counts()], axis=1)
                
                temp_data = pd.DataFrame(result[result['fc']<(1/i)]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure_down = pd.concat([result_branches_structure_down, resultlist], axis=1)
                result_core_structure_down = pd.concat([result_core_structure_down, self.identify_core_structure(result[result['fc']<(1/i)])['Core_structure'].value_counts()], axis=1)
                result_glycan_type_down = pd.concat([result_glycan_type_down, result[result['fc']<(1/i)]['Glycan_type'].value_counts()], axis=1)
                result_branches_count_down = pd.concat([result_branches_count_down, result[result['fc']<(1/i)]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition_down = pd.concat([result_glycan_composition_down, result[result['fc']<(1/i)]['GlycanComposition'].value_counts()], axis=1)
                result_lacdinac_down = pd.concat([result_lacdinac_down, result[result['fc']<(1/i)]['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac').value_counts()], axis=1)
                result_fucosylated_type_down = pd.concat([result_fucosylated_type_down, result[result['fc']<(1/i)]['fucosylated type'].value_counts()], axis=1)
                result_acgc_down = pd.concat([result_acgc_down, result[result['fc']<(1/i)]['Ac/Gc'].value_counts()], axis=1)
                
                result_branches_structure_split = pd.concat([result_branches_structure_split, result_branches_structure_up.iloc[:,-1], result_branches_structure_down.iloc[:,-1]], axis=1)
                result_core_structure_split = pd.concat([result_core_structure_split, result_core_structure_up.iloc[:,-1], result_core_structure_down.iloc[:,-1]], axis=1)
                result_glycan_type_split = pd.concat([result_glycan_type_split, result_glycan_type_up.iloc[:,-1], result_glycan_type_down.iloc[:,-1]], axis=1)
                result_branches_count_split = pd.concat([result_branches_count_split, result_branches_count_up.iloc[:,-1], result_branches_count_down.iloc[:,-1]], axis=1)
                result_glycan_composition_split = pd.concat([result_glycan_composition_split, result_glycan_composition_up.iloc[:,-1], result_glycan_composition_down.iloc[:,-1]], axis=1)
                result_lacdinac_split = pd.concat([result_lacdinac_split, result_lacdinac_up.iloc[:,-1], result_lacdinac_down.iloc[:,-1]], axis=1)
                result_fucosylated_type_split = pd.concat([result_fucosylated_type_split, result_fucosylated_type_up.iloc[:,-1], result_fucosylated_type_down.iloc[:,-1]], axis=1)
                result_acgc_split = pd.concat([result_acgc_split, result_acgc_up.iloc[:,-1], result_acgc_down.iloc[:,-1]], axis=1)
                
                result_branches_structure_ratio = pd.concat([result_branches_structure_ratio, result_branches_structure_up.iloc[:,-1] / result_branches_structure_down.iloc[:,-1]], axis=1)
                result_core_structure_ratio = pd.concat([result_core_structure_ratio, result_core_structure_up.iloc[:,-1] / result_core_structure_down.iloc[:,-1]], axis=1)
                result_glycan_type_ratio = pd.concat([result_glycan_type_ratio, result_glycan_type_up.iloc[:,-1] / result_glycan_type_down.iloc[:,-1]], axis=1)
                result_branches_count_ratio = pd.concat([result_branches_count_ratio, result_branches_count_up.iloc[:,-1] / result_branches_count_down.iloc[:,-1]], axis=1)
                result_glycan_composition_ratio = pd.concat([result_glycan_composition_ratio, result_glycan_composition_up.iloc[:,-1] / result_glycan_composition_down.iloc[:,-1]], axis=1)
                result_lacdinac_ratio = pd.concat([result_lacdinac_ratio, result_lacdinac_up.iloc[:,-1] / result_lacdinac_down.iloc[:,-1]], axis=1)
                result_fucosylated_type_ratio = pd.concat([result_fucosylated_type_ratio, result_fucosylated_type_up.iloc[:,-1] / result_fucosylated_type_down.iloc[:,-1]], axis=1)
                result_acgc_ratio = pd.concat([result_acgc_ratio, result_acgc_up.iloc[:,-1] / result_acgc_down.iloc[:,-1]], axis=1)
                
                def process_data(data):
                    # data = data.replace([0, np.inf, np.nan], 1)
                    # data = np.log2(data)
                    return data
                
                self.result_branches_structure_ratio = process_data(result_branches_structure_ratio)
                self.result_core_structure_ratio = process_data(result_core_structure_ratio)
                self.result_glycan_type_ratio = process_data(result_glycan_type_ratio)
                self.result_branches_count_ratio = process_data(result_branches_count_ratio)
                self.result_glycan_composition_ratio = process_data(result_glycan_composition_ratio)
                self.result_lacdinac_ratio = process_data(result_lacdinac_ratio)
                self.result_fucosylated_type_ratio = process_data(result_fucosylated_type_ratio)
                self.result_acgc_ratio = process_data(result_acgc_ratio)

            self.result_branches_structure = result_branches_structure.reset_index()
            self.result_branches_structure.columns = ['Branches', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_core_structure = result_core_structure.reset_index()
            self.result_core_structure.columns = ['Core structure', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_glycan_type = result_glycan_type.reset_index()
            self.result_glycan_type.columns = ['Glycan type', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_branches_count = result_branches_count.reset_index()
            self.result_branches_count.columns = ['Branches count', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_glycan_composition = result_glycan_composition.reset_index()
            self.result_glycan_composition.columns = ['Glycan composition', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_lacdinac = result_lacdinac.reset_index()
            self.result_lacdinac.columns = ['Lacdinac', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_lacdinac = self.result_lacdinac[self.result_lacdinac['Lacdinac']!=' ']
            self.result_fucosylated_type = result_fucosylated_type.reset_index()
            self.result_fucosylated_type.columns = ['Fucosylated type', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_fucosylated_type = self.result_fucosylated_type[self.result_fucosylated_type['Fucosylated type']!=' ']
            self.result_acgc = result_acgc.reset_index()
            self.result_acgc.columns = ['Ac/Gc', 'fc>1.2|fc<1/1.2', 'fc>1.5|fc<1/1.5', 'fc>2|fc<1/2', 'fc>2.5|fc<1/2.5', 'fc>3|fc<1/3']
            self.result_acgc = self.result_acgc[self.result_acgc['Ac/Gc']!=' ']
            
            self.result_branches_structure_up = result_branches_structure_up.reset_index()
            self.result_branches_structure_up.columns = ['Branches', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_core_structure_up = result_core_structure_up.reset_index()
            self.result_core_structure_up.columns = ['Core structure', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_glycan_type_up = result_glycan_type_up.reset_index()
            self.result_glycan_type_up.columns = ['Glycan type', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_branches_count_up = result_branches_count_up.reset_index()
            self.result_branches_count_up.columns = ['Branches count', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_glycan_composition_up = result_glycan_composition_up.reset_index()
            self.result_glycan_composition_up.columns = ['Glycan composition', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_lacdinac_up = result_lacdinac_up.reset_index()
            self.result_lacdinac_up.columns = ['Lacdinac', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_lacdinac_up = self.result_lacdinac_up[self.result_lacdinac_up['Lacdinac']!=' ']
            self.result_fucosylated_type_up = result_fucosylated_type_up.reset_index()
            self.result_fucosylated_type_up.columns = ['Fucosylated type', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_fucosylated_type_up = self.result_fucosylated_type_up[self.result_fucosylated_type_up['Fucosylated type']!=' ']
            self.result_acgc_up = result_acgc_up.reset_index()
            self.result_acgc_up.columns = ['Ac/Gc', 'fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_acgc_up = self.result_acgc_up[self.result_acgc_up['Ac/Gc']!=' ']
            
            self.result_branches_structure_up_ratio = self.result_branches_structure_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_branches_structure_up_ratio[columns_to_normalize] = self.result_branches_structure_up_ratio[columns_to_normalize].div(self.result_branches_structure_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_core_structure_up_ratio = self.result_core_structure_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_core_structure_up_ratio[columns_to_normalize] = self.result_core_structure_up_ratio[columns_to_normalize].div(self.result_core_structure_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_glycan_type_up_ratio = self.result_glycan_type_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_glycan_type_up_ratio[columns_to_normalize] = self.result_glycan_type_up_ratio[columns_to_normalize].div(self.result_glycan_type_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_branches_count_up_ratio = self.result_branches_count_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_branches_count_up_ratio[columns_to_normalize] = self.result_branches_count_up_ratio[columns_to_normalize].div(self.result_branches_count_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_glycan_composition_up_ratio = self.result_glycan_composition_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_glycan_composition_up_ratio[columns_to_normalize] = self.result_glycan_composition_up_ratio[columns_to_normalize].div(self.result_glycan_composition_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_lacdinac_up_ratio = self.result_lacdinac_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_lacdinac_up_ratio[columns_to_normalize] = self.result_lacdinac_up_ratio[columns_to_normalize].div(self.result_lacdinac_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_fucosylated_type_up_ratio = self.result_fucosylated_type_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_fucosylated_type_up_ratio[columns_to_normalize] = self.result_fucosylated_type_up_ratio[columns_to_normalize].div(self.result_fucosylated_type_up_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_acgc_up_ratio = self.result_acgc_up.copy()
            columns_to_normalize = ['fc>1.2', 'fc>1.5', 'fc>2', 'fc>2.5', 'fc>3']
            self.result_acgc_up_ratio[columns_to_normalize] = self.result_acgc_up_ratio[columns_to_normalize].div(self.result_acgc_up_ratio[columns_to_normalize].sum(axis=0), axis=1)

            self.result_branches_structure_down = result_branches_structure_down.reset_index()
            self.result_branches_structure_down.columns = ['Branches', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_core_structure_down = result_core_structure_down.reset_index()
            self.result_core_structure_down.columns = ['Core structure', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_glycan_type_down = result_glycan_type_down.reset_index()
            self.result_glycan_type_down.columns = ['Glycan type', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_branches_count_down = result_branches_count_down.reset_index()
            self.result_branches_count_down.columns = ['Branches count', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_glycan_composition_down = result_glycan_composition_down.reset_index()
            self.result_glycan_composition_down.columns = ['Glycan composition', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_lacdinac_down = result_lacdinac_down.reset_index()
            self.result_lacdinac_down.columns = ['Lacdinac', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_lacdinac_down = self.result_lacdinac_down[self.result_lacdinac_down['Lacdinac']!=' ']
            self.result_fucosylated_type_down = result_fucosylated_type_down.reset_index()
            self.result_fucosylated_type_down.columns = ['Fucosylated type', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_fucosylated_type_down = self.result_fucosylated_type_down[self.result_fucosylated_type_down['Fucosylated type']!=' ']
            self.result_acgc_down = result_acgc_down.reset_index()
            self.result_acgc_down.columns = ['Ac/Gc', 'fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_acgc_down = self.result_acgc_down[self.result_acgc_down['Ac/Gc']!=' ']
            
            self.result_branches_structure_down_ratio = self.result_branches_structure_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_branches_structure_down_ratio[columns_to_normalize] = self.result_branches_structure_down_ratio[columns_to_normalize].div(self.result_branches_structure_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_core_structure_down_ratio = self.result_core_structure_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_core_structure_down_ratio[columns_to_normalize] = self.result_core_structure_down_ratio[columns_to_normalize].div(self.result_core_structure_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_glycan_type_down_ratio = self.result_glycan_type_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_glycan_type_down_ratio[columns_to_normalize] = self.result_glycan_type_down_ratio[columns_to_normalize].div(self.result_glycan_type_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_branches_count_down_ratio = self.result_branches_count_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_branches_count_down_ratio[columns_to_normalize] = self.result_branches_count_down_ratio[columns_to_normalize].div(self.result_branches_count_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_glycan_composition_down_ratio = self.result_glycan_composition_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_glycan_composition_down_ratio[columns_to_normalize] = self.result_glycan_composition_down_ratio[columns_to_normalize].div(self.result_glycan_composition_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_lacdinac_down_ratio = self.result_lacdinac_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_lacdinac_down_ratio[columns_to_normalize] = self.result_lacdinac_down_ratio[columns_to_normalize].div(self.result_lacdinac_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_fucosylated_type_down_ratio = self.result_fucosylated_type_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_fucosylated_type_down_ratio[columns_to_normalize] = self.result_fucosylated_type_down_ratio[columns_to_normalize].div(self.result_fucosylated_type_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            self.result_acgc_down_ratio = self.result_acgc_down.copy()
            columns_to_normalize = ['fc<1/1.2', 'fc<1/1.5', 'fc<1/2', 'fc<1/2.5', 'fc<1/3']
            self.result_acgc_down_ratio[columns_to_normalize] = self.result_acgc_down_ratio[columns_to_normalize].div(self.result_acgc_down_ratio[columns_to_normalize].sum(axis=0), axis=1)
            
            self.result_branches_structure_split = result_branches_structure_split.reset_index()
            self.result_branches_structure_split.columns = ['Branches', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_core_structure_split = result_core_structure_split.reset_index()
            self.result_core_structure_split.columns = ['Core structure', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_glycan_type_split = result_glycan_type_split.reset_index()
            self.result_glycan_type_split.columns = ['Glycan type', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_branches_count_split = result_branches_count_split.reset_index()
            self.result_branches_count_split.columns = ['Branches count', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_glycan_composition_split = result_glycan_composition_split.reset_index()
            self.result_glycan_composition_split.columns = ['Glycan composition', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_lacdinac_split = result_lacdinac_split.reset_index()
            self.result_lacdinac_split.columns = ['Lacdinac', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_lacdinac_split = self.result_lacdinac_split[self.result_lacdinac_split['Lacdinac']!=' ']
            self.result_fucosylated_type_split = result_fucosylated_type_split.reset_index()
            self.result_fucosylated_type_split.columns = ['Fucosylated type', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_fucosylated_type_split = self.result_fucosylated_type_split[self.result_fucosylated_type_split['Fucosylated type']!=' ']
            self.result_acgc_split = result_acgc_split.reset_index()
            self.result_acgc_split.columns = ['Ac/Gc', 'fc>1.2', 'fc<1/1.2', 'fc>1.5', 'fc<1/1.5', 'fc>2', 'fc<1/2', 'fc>2.5', 'fc<1/2.5', 'fc>3', 'fc<1/3']
            self.result_acgc_split = self.result_acgc_split[self.result_acgc_split['Ac/Gc']!=' ']
            
            self.result_branches_structure_ratio = self.result_branches_structure_ratio.reset_index()
            self.result_branches_structure_ratio.columns = ['Branches', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_core_structure_ratio = self.result_core_structure_ratio.reset_index()
            self.result_core_structure_ratio.columns = ['Core structure', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_glycan_type_ratio = self.result_glycan_type_ratio.reset_index()
            self.result_glycan_type_ratio.columns = ['Glycan type', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_branches_count_ratio = self.result_branches_count_ratio.reset_index()
            self.result_branches_count_ratio.columns = ['Branches count', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_glycan_composition_ratio = self.result_glycan_composition_ratio.reset_index()
            self.result_glycan_composition_ratio.columns = ['Glycan composition', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_lacdinac_ratio = self.result_lacdinac_ratio.reset_index()
            self.result_lacdinac_ratio.columns = ['Lacdinac', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_lacdinac_ratio = self.result_lacdinac_ratio[self.result_lacdinac_ratio['Lacdinac']!=' ']
            self.result_fucosylated_type_ratio = self.result_fucosylated_type_ratio.reset_index()
            self.result_fucosylated_type_ratio.columns = ['Fucosylated type', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_fucosylated_type_ratio = self.result_fucosylated_type_ratio[self.result_fucosylated_type_ratio['Fucosylated type']!=' ']
            self.result_acgc_ratio = self.result_acgc_ratio.reset_index()
            self.result_acgc_ratio.columns = ['Ac/Gc', 'fc>1.2 / fc<1/1.2', 'fc>1.5 / fc<1/1.5', 'fc>2 / fc<1/2', 'fc>2.5 / fc<1/2.5', 'fc>3 / fc<1/3']
            self.result_acgc_ratio = self.result_acgc_ratio[self.result_acgc_ratio['Ac/Gc']!=' ']

        elif statistic_index == 'roc':
            result = result[result['auc_pvalue']<0.05]
            
            unique_branches = set()
            for branches_str in result['Branches'].value_counts().index:
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            result_branches_structure = pd.DataFrame(index=unique_branches_list)
            
            result_core_structure = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            
            for i in [0.6, 0.7, 0.8, 0.9, 1]:
                
                unique_branches = set()
                for branches_str in result['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                temp_data = pd.DataFrame(result[result['auc']>=i]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure = pd.concat([result_branches_structure, resultlist], axis=1)
                
                result_core_structure = pd.concat([result_core_structure, self.identify_core_structure(result[result['auc']>=i])['Core_structure'].value_counts()], axis=1)
                result_glycan_type = pd.concat([result_glycan_type, result[result['auc']>=i]['Glycan_type'].value_counts()], axis=1)
                result_branches_count = pd.concat([result_branches_count, result[result['auc']>=i]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition = pd.concat([result_glycan_composition, result[result['auc']>=i]['GlycanComposition'].value_counts()], axis=1)
                
            self.result_branches_structure = result_branches_structure.reset_index()
            self.result_branches_structure.columns = ['Branches', 'auc>=0.6', 'auc>=0.7', 'auc>=0.8', 'auc>=0.9', 'auc>=1']
            self.result_core_structure = result_core_structure.reset_index()
            self.result_core_structure.columns = ['Core structure', 'auc>=0.6', 'auc>=0.7', 'auc>=0.8', 'auc>=0.9', 'auc>=1']
            self.result_glycan_type = result_glycan_type.reset_index()
            self.result_glycan_type.columns = ['Glycan type', 'auc>=0.6', 'auc>=0.7', 'auc>=0.8', 'auc>=0.9', 'auc>=1']
            self.result_branches_count = result_branches_count.reset_index()
            self.result_branches_count.columns = ['Branches count', 'auc>=0.6', 'auc>=0.7', 'auc>=0.8', 'auc>=0.9', 'auc>=1']
            self.result_glycan_composition = result_glycan_composition.reset_index()
            self.result_glycan_composition.columns = ['GlycanComposition', 'auc>=0.6', 'auc>=0.7', 'auc>=0.8', 'auc>=0.9', 'auc>=1']

        elif statistic_index == 'ml':
            
            unique_branches = set()
            for branches_str in result['Branches'].value_counts().index:
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            result_branches_structure = pd.DataFrame(index=unique_branches_list)
            
            quantiles = np.quantile(sorted(result['randomforest_features_importance_means'].dropna()), np.linspace(0, 1, 5))
            column_names = ['Branches'] + [f'rf_features_importance_score>={q:.5f}' for q in quantiles]
            result_core_structure = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            
            for i in quantiles:
                
                unique_branches = set()
                for branches_str in result['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                temp_data = pd.DataFrame(result[result['randomforest_features_importance_means']>=i]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure = pd.concat([result_branches_structure, resultlist], axis=1)

                result_core_structure = pd.concat([result_core_structure, self.identify_core_structure(result[result['randomforest_features_importance_means']>=i])['Core_structure'].value_counts()], axis=1)
                result_glycan_type = pd.concat([result_glycan_type, result[result['randomforest_features_importance_means']>=i]['Glycan_type'].value_counts()], axis=1)
                result_branches_count = pd.concat([result_branches_count, result[result['randomforest_features_importance_means']>=i]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition = pd.concat([result_glycan_composition, result[result['randomforest_features_importance_means']>=i]['GlycanComposition'].value_counts()], axis=1)
                
            self.result_branches_structure = result_branches_structure.reset_index()
            self.result_branches_structure.columns = column_names
            column_names = ['Core structure'] + [f'rf_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_core_structure = result_core_structure.reset_index()
            self.result_core_structure.columns = column_names
            column_names = ['Glycan type'] + [f'rf_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_glycan_type = result_glycan_type.reset_index()
            self.result_glycan_type.columns = column_names
            column_names = ['Branches count'] + [f'rf_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_branches_count = result_branches_count.reset_index()
            self.result_branches_count.columns = column_names
            column_names = ['GlycanComposition'] + [f'rf_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_glycan_composition = result_glycan_composition.reset_index()
            self.result_glycan_composition.columns = column_names
        
        elif statistic_index == 'pca':
            
            unique_branches = set()
            for branches_str in result['Branches'].value_counts().index:
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            result_branches_structure = pd.DataFrame(index=unique_branches_list)
            
            quantiles = np.quantile(sorted(result['pca_features_importance'].dropna()), np.linspace(0, 1, 5))
            column_names = ['Branches'] + [f'pca_features_importance_score>={q:.5f}' for q in quantiles]
            result_core_structure = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            
            for i in quantiles:
                
                unique_branches = set()
                for branches_str in result['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                temp_data = pd.DataFrame(result[result['pca_features_importance']>=i]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure = pd.concat([result_branches_structure, resultlist], axis=1)
                
                result_core_structure = pd.concat([result_core_structure, self.identify_core_structure(result[result['pca_features_importance']>=i])['Core_structure'].value_counts()], axis=1)
                result_glycan_type = pd.concat([result_glycan_type, result[result['pca_features_importance']>=i]['Glycan_type'].value_counts()], axis=1)
                result_branches_count = pd.concat([result_branches_count, result[result['pca_features_importance']>=i]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition = pd.concat([result_glycan_composition, result[result['pca_features_importance']>=i]['GlycanComposition'].value_counts()], axis=1)
                
            self.result_branches_structure = result_branches_structure.reset_index()
            self.result_branches_structure.columns = column_names
            column_names = ['Core structure'] + [f'pca_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_core_structure = result_core_structure.reset_index()
            self.result_core_structure.columns = column_names
            column_names = ['Glycan type'] + [f'pca_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_glycan_type = result_glycan_type.reset_index()
            self.result_glycan_type.columns = column_names
            column_names = ['Branches count'] + [f'pca_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_branches_count = result_branches_count.reset_index()
            self.result_branches_count.columns = column_names
            column_names = ['GlycanComposition'] + [f'pca_features_importance_score>={q:.5f}' for q in quantiles]
            self.result_glycan_composition = result_glycan_composition.reset_index()
            self.result_glycan_composition.columns = column_names
            
        elif statistic_index == 'anova':
            result = result[result['anova_pvalue']<0.05]
            
            unique_branches = set()
            for branches_str in result['Branches'].value_counts().index:
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            result_branches_structure = pd.DataFrame(index=unique_branches_list)
            
            quantiles = np.quantile(sorted(result['f score'].dropna()), np.linspace(0, 1, 5))
            column_names = ['Branches'] + [f'f score>={q:.2f}' for q in quantiles]
            result_core_structure = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            
            for i in quantiles:
                
                unique_branches = set()
                for branches_str in result['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                temp_data = pd.DataFrame(result[result['f score']>=i]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure = pd.concat([result_branches_structure, resultlist], axis=1)
                
                result_core_structure = pd.concat([result_core_structure, self.identify_core_structure(result[result['f score']>=i])['Core_structure'].value_counts()], axis=1)
                result_glycan_type = pd.concat([result_glycan_type, result[result['f score']>=i]['Glycan_type'].value_counts()], axis=1)
                result_branches_count = pd.concat([result_branches_count, result[result['f score']>=i]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition = pd.concat([result_glycan_composition, result[result['f score']>=i]['GlycanComposition'].value_counts()], axis=1)
                
            self.result_branches_structure = result_branches_structure.reset_index()
            self.result_branches_structure.columns = column_names
            column_names = ['Core structure'] + [f'f score>={q:.2f}' for q in quantiles]
            self.result_core_structure = result_core_structure.reset_index()
            self.result_core_structure.columns = column_names
            column_names = ['Glycan type'] + [f'f score>={q:.2f}' for q in quantiles]
            self.result_glycan_type = result_glycan_type.reset_index()
            self.result_glycan_type.columns = column_names
            column_names = ['Branches count'] + [f'f score>={q:.2f}' for q in quantiles]
            self.result_branches_count = result_branches_count.reset_index()
            self.result_branches_count.columns = column_names
            column_names = ['GlycanComposition'] + [f'f score>={q:.2f}' for q in quantiles]
            self.result_glycan_composition = result_glycan_composition.reset_index()
            self.result_glycan_composition.columns = column_names
            
        elif statistic_index == 'chi2':
            result = result[result['chi2_pvalue']<0.05]
            
            unique_branches = set()
            for branches_str in result['Branches'].value_counts().index:
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            result_branches_structure = pd.DataFrame(index=unique_branches_list)
            
            quantiles = np.quantile(sorted(result['chi2 score'].dropna()), np.linspace(0, 1, 5))
            column_names = ['Branches'] + [f'chi2 score>={q:.2f}' for q in quantiles]
            result_core_structure = pd.DataFrame(index = self.identify_core_structure(self.data)['Core_structure'].value_counts().index)
            result_glycan_type = pd.DataFrame(index = self.data['Glycan_type'].value_counts().index)
            result_branches_count = pd.DataFrame(index = self.data['BranchNumber'].value_counts().index)
            result_glycan_composition = pd.DataFrame(index = self.data['GlycanComposition'].value_counts().index)
            
            for i in quantiles:
                
                unique_branches = set()
                for branches_str in result['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                temp_data = pd.DataFrame(result[result['chi2 score']>=i]['Branches'])
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                resultlist = pd.DataFrame(index=resultlist.keys(),data=resultlist.values())
                result_branches_structure = pd.concat([result_branches_structure, resultlist], axis=1)

                result_core_structure = pd.concat([result_core_structure, self.identify_core_structure(result[result['chi2 score']>=i])['Core_structure'].value_counts()], axis=1)
                result_glycan_type = pd.concat([result_glycan_type, result[result['chi2 score']>=i]['Glycan_type'].value_counts()], axis=1)
                result_branches_count = pd.concat([result_branches_count, result[result['chi2 score']>=i]['BranchNumber'].value_counts()], axis=1)
                result_glycan_composition = pd.concat([result_glycan_composition, result[result['chi2 score']>=i]['GlycanComposition'].value_counts()], axis=1)
                
            self.result_branches_structure = result_branches_structure.reset_index()
            self.result_branches_structure.columns = column_names
            column_names = ['Core structure'] + [f'chi2 score>={q:.2f}' for q in quantiles]
            self.result_core_structure = result_core_structure.reset_index()
            self.result_core_structure.columns = column_names
            column_names = ['Glycan type'] + [f'chi2 score>={q:.2f}' for q in quantiles]
            self.result_glycan_type = result_glycan_type.reset_index()
            self.result_glycan_type.columns = column_names
            column_names = ['Branches count'] + [f'chi2 score>={q:.2f}' for q in quantiles]
            self.result_branches_count = result_branches_count.reset_index()
            self.result_branches_count.columns = column_names
            column_names = ['GlycanComposition'] + [f'chi2 score>={q:.2f}' for q in quantiles]
            self.result_glycan_composition = result_glycan_composition.reset_index()
            self.result_glycan_composition.columns = column_names
        #
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'threshold_variation_analysis', {'statistic_index':statistic_index})
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_core_structure', self.result_core_structure)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_branches_structure', self.result_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_glycan_type', self.result_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_branches_count', self.result_branches_count)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_glycan_composition', self.result_glycan_composition)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_lacdinac', self.result_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_fucosylated_type', self.result_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_acgc', self.result_acgc)

        return self
    
    def glycopeptide_glycosite_glycan_variation(self, fc = 1.5):
        """
        Quantifies the number and distribution of upregulated and downregulated glycan structures associated with each glycosylation site at glycopeptide level.
        
        Parameters:
            fc: FC threshold used for analysis.
        
        Returns:
            self.result_glycopeptide_glycosite_glycan_variation
                
        Return type:
            dataframe
        
        """
        #
        if fc == None:
            fc = input(f"Please enter the fc used for glycopeptide glycosite variation analysis (default fc was {fc}): ")

        temp_data = pd.merge(self.data[['PeptideSequence', 'Glycosite_Position']], self.fc_result, left_index=True, right_index=True, how='left')
        temp_data = temp_data[~temp_data['fc'].isnull()]
        result = []
        for idx, row in temp_data.iterrows():
            peptide = row['PeptideSequence']
            glycosite = row['Glycosite_Position']
            filtered_data = temp_data[(temp_data['PeptideSequence'] == peptide) & (temp_data['Glycosite_Position'] == glycosite)]
            up_num = filtered_data[filtered_data['fc'] > fc].shape[0] 
            down_num = filtered_data[filtered_data['fc'] < 1/fc].shape[0] 
            if not (1/fc < row['fc'] <fc):
                result.append([idx, peptide, glycosite, up_num, down_num]) 
        self.result_glycopeptide_glycosite_glycan_variation = pd.DataFrame(result, columns=['Index', 'Glycopeptide', 'Glycosite_Position', 'Up count', 'Down count'])
        self.result_glycopeptide_glycosite_glycan_variation = self.result_glycopeptide_glycosite_glycan_variation.drop_duplicates()
        
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'glycopeptide_glycosite_glycan_variation', {})
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_glycopeptide_glycosite_glycan_variation', self.result_glycopeptide_glycosite_glycan_variation)
        
        return self
    
    def glycoprotein_glycosite_glycan_variation(self, fc = 1.5):
        """
        Quantifies the number and distribution of upregulated and downregulated glycan structures associated with each glycosylation site at glycoprotein level.
        
        Parameters:
            fc: FC threshold used for analysis.
        
        Returns:
            self.result_glycoprotein_glycosite_glycan_variation
                
        Return type:
            dataframe
        
        """
        #
        if fc == None:
            fc = input(f"Please enter the fc used for glycoprotein glycosite variation analysis (default fc was {fc}): ")
 
        temp_data = pd.merge(self.data[['ProteinID', 'Glycosite_Position']], self.fc_result, left_index=True, right_index=True, how='left')
        temp_data = temp_data[~temp_data['fc'].isnull()]
        result = []
        for idx, row in temp_data.iterrows():
            peptide = row['ProteinID']
            glycosite = row['Glycosite_Position']
            filtered_data = temp_data[(temp_data['ProteinID'] == peptide) & (temp_data['Glycosite_Position'] == glycosite)]
            up_num = filtered_data[filtered_data['fc'] > fc].shape[0] 
            down_num = filtered_data[filtered_data['fc'] < 1/fc].shape[0] 
            if not (1/fc < row['fc'] <fc):
                result.append([idx, peptide, glycosite, up_num, down_num]) 
        self.result_glycoprotein_glycosite_glycan_variation = pd.DataFrame(result, columns=['Index', 'Glycoprotein', 'Glycosite_Position', 'Up count', 'Down count'])
        self.result_glycoprotein_glycosite_glycan_variation = self.result_glycoprotein_glycosite_glycan_variation.drop_duplicates()
        
        self.data_manager.log_params('StrucGAP_GlycoPeptideQuant', 'glycoprotein_glycosite_glycan_variation', {})
        self.data_manager.log_output('StrucGAP_GlycoPeptideQuant', 'result_glycoprotein_glycosite_glycan_variation', self.result_glycoprotein_glycosite_glycan_variation)

        return self
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        
        with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_GlycoPeptideQuant.xlsx'), engine='xlsxwriter') as writer:
            self.data_quant.reset_index(drop=True).to_excel(writer, sheet_name='data_quant')
            
            self.differential_analysis_data.reset_index(drop=True).to_excel(writer, sheet_name='differential_analysis_data'[:31])
            self.up_data.reset_index(drop=True).to_excel(writer, sheet_name='up_data'[:31])
            self.down_data.reset_index(drop=True).to_excel(writer, sheet_name='down_data'[:31])
            
            pd.DataFrame([self.number_of_identified_spectrum]).to_excel(writer, sheet_name='number_of_identified_spectrum'[:31])
            pd.DataFrame([self.number_of_identified_peptide]).to_excel(writer, sheet_name='number_of_identified_peptide'[:31])
            self.type_of_identified_peptide.to_excel(writer, sheet_name='type_of_identified_peptide'[:31])
            pd.DataFrame([self.number_of_identified_glycopeptide]).to_excel(writer, sheet_name='number_of_identified_glycopeptide'[:31])
            self.type_of_identified_glycopeptide.to_excel(writer, sheet_name='type_of_identified_glycopeptide'[:31])
            pd.DataFrame([self.number_of_identified_glycoprotein]).to_excel(writer, sheet_name='number_of_identified_glycoprotein'[:31])
            self.type_of_identified_glycoprotein.to_excel(writer, sheet_name='type_of_identified_glycoprotein'[:31])
            
            self.fc_result.to_excel(writer, sheet_name='fc_result')
            self.roc_result.to_excel(writer, sheet_name='roc_result')
            self.ml_result.to_excel(writer, sheet_name='ml_result')
            self.pca_result.to_excel(writer, sheet_name='pca_result')
            self.anova_result.to_excel(writer, sheet_name='anova_result')
            self.chi2_result.to_excel(writer, sheet_name='chi2_result')
            
            self.differential_analysis_overview.to_excel(writer, sheet_name='da_overview'[:31])
            self.differential_analysis_structure_coding.to_excel(writer, sheet_name='da_glycan'[:31])
            self.differential_analysis_structure_coding_no_oligomannose.to_excel(writer, sheet_name='da_glycan_no_oligomannose'[:31])
            self.differential_analysis_core_structure.to_excel(writer, sheet_name='da_core_structure'[:31])
            self.differential_analysis_glycan_type.to_excel(writer, sheet_name='da_glycan_type'[:31])
            self.differential_analysis_branches_structure.to_excel(writer, sheet_name='da_branch_structure'[:31])
            self.differential_analysis_branches_count.to_excel(writer, sheet_name='da_branch_count'[:31])
            self.differential_analysis_sialicacid_count.to_excel(writer, sheet_name='da_sialicacid_count'[:31])
            self.differential_analysis_fucose_count.to_excel(writer, sheet_name='da_fucose_count'[:31])
            self.differential_analysis_sialicacid_structure.to_excel(writer, sheet_name='da_sialicacid_structure'[:31])
            self.differential_analysis_fucose_structure.to_excel(writer, sheet_name='da_fucose_structure'[:31])
            self.differential_analysis_lacdinac.to_excel(writer, sheet_name='da_lacdinac'[:31])
            self.differential_analysis_fucosylated_type.to_excel(writer, sheet_name='da_fucosylated_type'[:31])
            self.differential_analysis_acgc.to_excel(writer, sheet_name='da_acgc'[:31])
            
            self.result_core_structure.to_excel(writer, sheet_name='core_structure'[:31])
            self.result_core_structure_split.to_excel(writer, sheet_name='core_structure_split'[:31])
            self.result_core_structure_ratio.to_excel(writer, sheet_name='core_structure_ratio'[:31])
            self.result_core_structure_up.to_excel(writer, sheet_name='core_structure_up'[:31])
            self.result_core_structure_up_ratio.to_excel(writer, sheet_name='core_structure_up_ratio'[:31])
            self.result_core_structure_down.to_excel(writer, sheet_name='core_structure_down'[:31])
            self.result_core_structure_down_ratio.to_excel(writer, sheet_name='core_structure_down_ratio'[:31])
            
            self.result_branches_structure.to_excel(writer, sheet_name='branches_structure'[:31])
            self.result_branches_structure_split.to_excel(writer, sheet_name='branches_structure_split'[:31])
            self.result_branches_structure_ratio.to_excel(writer, sheet_name='branches_structure_ratio'[:31])
            self.result_branches_structure_up.to_excel(writer, sheet_name='branches_structure_up'[:31])
            self.result_branches_structure_up_ratio.to_excel(writer, sheet_name='branches_structure_up_ratio'[:31])
            self.result_branches_structure_down.to_excel(writer, sheet_name='branches_structure_down'[:31])
            self.result_branches_structure_down_ratio.to_excel(writer, sheet_name='branches_structure_down_ratio'[:31])
            
            self.result_glycan_type.to_excel(writer, sheet_name='glycan_type'[:31])
            self.result_glycan_type_split.to_excel(writer, sheet_name='glycan_type_split'[:31])
            self.result_glycan_type_ratio.to_excel(writer, sheet_name='glycan_type_ratio'[:31])
            self.result_glycan_type_up.to_excel(writer, sheet_name='glycan_type_up'[:31])
            self.result_glycan_type_up_ratio.to_excel(writer, sheet_name='glycan_type_up_ratio'[:31])
            self.result_glycan_type_down.to_excel(writer, sheet_name='glycan_type_down'[:31])
            self.result_glycan_type_down_ratio.to_excel(writer, sheet_name='glycan_type_down_ratio'[:31])
            
            self.result_branches_count.to_excel(writer, sheet_name='branches_count'[:31])
            self.result_branches_count_split.to_excel(writer, sheet_name='branches_count_split'[:31])
            self.result_branches_count_ratio.to_excel(writer, sheet_name='branches_count_ratio'[:31])
            self.result_branches_count_up.to_excel(writer, sheet_name='branches_count_up'[:31])
            self.result_branches_count_up_ratio.to_excel(writer, sheet_name='branches_count_up_ratio'[:31])
            self.result_branches_count_down.to_excel(writer, sheet_name='branches_count_down'[:31])
            self.result_branches_count_down_ratio.to_excel(writer, sheet_name='branches_count_down_ratio'[:31])
            
            self.result_glycan_composition.to_excel(writer, sheet_name='glycan_composition'[:31])
            self.result_glycan_composition_split.to_excel(writer, sheet_name='glycan_composition_split'[:31])
            self.result_glycan_composition_ratio.to_excel(writer, sheet_name='glycan_composition_ratio'[:31])
            self.result_glycan_composition_up.to_excel(writer, sheet_name='glycan_composition_up'[:31])
            self.result_glycan_composition_up_ratio.to_excel(writer, sheet_name='glycan_composition_up_ratio'[:31])
            self.result_glycan_composition_down.to_excel(writer, sheet_name='glycan_composition_down'[:31])
            self.result_glycan_composition_down_ratio.to_excel(writer, sheet_name='glycan_composition_down_ratio'[:31])
            
            self.result_lacdinac.to_excel(writer, sheet_name='lacdinac'[:31])
            self.result_lacdinac_split.to_excel(writer, sheet_name='lacdinac_split'[:31])
            self.result_lacdinac_ratio.to_excel(writer, sheet_name='lacdinac_ratio'[:31])
            self.result_lacdinac_up.to_excel(writer, sheet_name='lacdinac_up'[:31])
            self.result_lacdinac_up_ratio.to_excel(writer, sheet_name='lacdinac_up_ratio'[:31])
            self.result_lacdinac_down.to_excel(writer, sheet_name='lacdinac_down'[:31])
            self.result_lacdinac_down_ratio.to_excel(writer, sheet_name='lacdinac_down_ratio'[:31])
            
            self.result_fucosylated_type.to_excel(writer, sheet_name='fucosylated_type'[:31])
            self.result_fucosylated_type_split.to_excel(writer, sheet_name='fucosylated_type_split'[:31])
            self.result_fucosylated_type_ratio.to_excel(writer, sheet_name='fucosylated_type_ratio'[:31])
            self.result_fucosylated_type_up.to_excel(writer, sheet_name='fucosylated_type_up'[:31])
            self.result_fucosylated_type_up_ratio.to_excel(writer, sheet_name='fucosylated_type_up_ratio'[:31])
            self.result_fucosylated_type_down.to_excel(writer, sheet_name='fucosylated_type_down'[:31])
            self.result_fucosylated_type_down_ratio.to_excel(writer, sheet_name='fucosylated_type_down_ratio'[:31])
            
            self.result_acgc.to_excel(writer, sheet_name='acgc'[:31])
            self.result_acgc_split.to_excel(writer, sheet_name='acgc_split'[:31])
            self.result_acgc_ratio.to_excel(writer, sheet_name='acgc_ratio'[:31])
            self.result_acgc_up.to_excel(writer, sheet_name='acgc_up'[:31])
            self.result_acgc_up_ratio.to_excel(writer, sheet_name='acgc_up_ratio'[:31])
            self.result_acgc_down.to_excel(writer, sheet_name='acgc_down'[:31])
            self.result_acgc_down_ratio.to_excel(writer, sheet_name='acgc_down_ratio'[:31])
            
            self.result_glycopeptide_glycosite_glycan_variation.to_excel(writer, sheet_name='glycopeptide_glycosite_glycan_variation'[:31])
            self.result_glycoprotein_glycosite_glycan_variation.to_excel(writer, sheet_name='glycoprotein_glycosite_glycan_variation'[:31])


 
