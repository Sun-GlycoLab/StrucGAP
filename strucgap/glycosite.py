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
## 位点特征分析模块--16
class StrucGAP_GlycoSite:
    def __init__(self, gs_data, data_manager, data_type = 'psm_filtered'):
        self.gs_data = gs_data
        self.sample_group = self.gs_data.sample_group
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
            
        self.data_manager = data_manager    
        self.data_manager.register_module('StrucGAP_GlycoSite', self, {})
        self.data_manager.log_params('StrucGAP_GlycoSite', 'input_data', {'data_type': data_type})
    
    def map_protein_to_gene(self, data, protein_column):
        """
        An auxiliary function called by other functions to map protein to gene.
        
        Parameters:
            data: the data that needs to be processed.
            protein_column: The column in the data containing protein accession.
        
        Returns:
            Column named 'gene_name'.
            
        Return type:
            series
        
        """
        # 将 ProteinID 和 GeneName 列拆分并映射为字典
        protein_to_gene = {}
        for idx, row in self.gs_data.data.iterrows():
            protein_ids = row['ProteinID'].split(';')
            gene_names = row['GeneName'].split(';')
            
            # 构建 ProteinID 到 GeneName 的映射
            for protein_id, gene_name in zip(protein_ids, gene_names):
                protein_to_gene[protein_id] = gene_name
        
        # 定义映射函数，将 protein_column 中的 ProteinID 替换为相应的 GeneName
        def map_genes(protein_ids):
            proteins = protein_ids.split(';')
            genes = [protein_to_gene.get(protein, "") for protein in proteins]
            return ';'.join(gene for gene in genes if gene)  # 忽略映射不到的 ProteinID
    
        # 将映射的 GeneName 添加到新的一列中，命名为 `gene_name`
        data['gene_name'] = data[protein_column].apply(map_genes)
        
    def glycoprotein_site(self):
        """
        Summarizes glycosylation at the protein level.
        
        Parameters:
            None.
        
        Returns:
            self.glycoprotein_glycosite_count
            self.glycoprotein_glycan_count
            self.glycoprotein_glycan_type
            self.glycoprotein_glycosite
            
        Return type:
            dataframe
        
        """
        #
        result = {}
        temp_data = pd.DataFrame(self.data[['ProteinID', 'Glycosite_Position']])
        temp_data = temp_data[~temp_data['Glycosite_Position'].isnull()]
        temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
        type_counts = pd.DataFrame(temp_data['ProteinID'].value_counts())
        for i in type_counts.index:
            filtered_data = temp_data[temp_data['ProteinID'] == i]['Glycosite_Position']
            item = filtered_data.unique()
            result[i] = len(item)
        self.glycoprotein_glycosite_count = pd.DataFrame(list(result.items()), columns=['glycoprotein', 'glycosite_count'])
        self.glycoprotein_glycosite_count = self.glycoprotein_glycosite_count.sort_values(by='glycosite_count', ascending=False)
        self.map_protein_to_gene(self.glycoprotein_glycosite_count, 'glycoprotein')
        #
        self.glycoprotein_glycan_count = pd.DataFrame(self.data['ProteinID'].value_counts().reset_index())
        self.glycoprotein_glycan_count.columns = ['glycoprotein', 'glycan_count']
        self.map_protein_to_gene(self.glycoprotein_glycan_count, 'glycoprotein')
        #
        if 'structure_coding' in self.data.columns:
            result = {}
            temp_data = pd.DataFrame(self.data[['ProteinID', 'structure_coding']])
            temp_data = temp_data[~temp_data['structure_coding'].isnull()]
            temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
            temp_data = temp_data.explode('ProteinID')
            type_counts = pd.DataFrame(temp_data['ProteinID'].value_counts())
            for i in type_counts.index:
                filtered_data = temp_data[temp_data['ProteinID'] == i]['structure_coding']
                item = filtered_data.unique()
                result[i] = len(item)
            self.glycoprotein_glycan_type = pd.DataFrame(list(result.items()), columns=['glycoprotein', 'glycan_type'])
            self.glycoprotein_glycan_type = self.glycoprotein_glycan_type.sort_values(by='glycan_type', ascending=False)
            self.map_protein_to_gene(self.glycoprotein_glycan_type, 'glycoprotein')
        else:
            self.glycoprotein_glycan_type = pd.DataFrame()
        #
        temp_data = pd.DataFrame(self.data[['ProteinID', 'Glycosite_Position']])
        temp_data = temp_data[~temp_data['Glycosite_Position'].isnull()]
        temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
        temp_data = temp_data.drop_duplicates()
        temp_data.columns = ['glycoprotein', 'glycosite']
        self.glycoprotein_glycosite = temp_data.reset_index(drop=True)
        self.map_protein_to_gene(self.glycoprotein_glycosite, 'glycoprotein')
        #
        self.data_manager.log_params('StrucGAP_GlycoSite', 'glycoprotein_site', {})
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycoprotein_glycosite_count', self.glycoprotein_glycosite_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycoprotein_glycan_count', self.glycoprotein_glycan_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycoprotein_glycan_type', self.glycoprotein_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycoprotein_glycosite', self.glycoprotein_glycosite)

        return self
    
    def glycopeptide_site(self):
        """
        Summarizes glycosylation at the glycopeptide level.
        
        Parameters:
            None.
        
        Returns:
            self.glycopeptide_glycosite_count
            self.glycopeptide_glycan_count
            self.glycopeptide_glycan_type
            self.glycopeptide_glycosite
            
        Return type:
            dataframe
        
        """
        #
        result = {}
        temp_data = pd.DataFrame(self.data[['PeptideSequence', 'Glycosite_Position']])
        temp_data = temp_data[~temp_data['Glycosite_Position'].isnull()]
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode('Glycosite_Position')
        type_counts = pd.DataFrame(temp_data['PeptideSequence'].value_counts())
        for i in type_counts.index:
            filtered_data = temp_data[temp_data['PeptideSequence'] == i]['Glycosite_Position']
            item = filtered_data.unique()
            result[i] = len(item)
        self.glycopeptide_glycosite_count = pd.DataFrame(list(result.items()), columns=['glycopeptide', 'glycosite_count'])
        self.glycopeptide_glycosite_count = self.glycopeptide_glycosite_count.sort_values(by='glycosite_count', ascending=False)
        #
        self.glycopeptide_glycan_count = pd.DataFrame(self.data['PeptideSequence'].value_counts().reset_index())
        self.glycopeptide_glycan_count.columns = ['glycopeptide', 'glycan_count']
        #
        if 'structure_coding' in self.data.columns:
            result = {}
            temp_data = pd.DataFrame(self.data[['PeptideSequence', 'structure_coding']])
            temp_data = temp_data[~temp_data['structure_coding'].isnull()]
            temp_data['PeptideSequence'] = temp_data['PeptideSequence'].str.split(';')
            temp_data = temp_data.explode('PeptideSequence')
            type_counts = pd.DataFrame(temp_data['PeptideSequence'].value_counts())
            for i in type_counts.index:
                filtered_data = temp_data[temp_data['PeptideSequence'] == i]['structure_coding']
                item = filtered_data.unique()
                result[i] = len(item)
            self.glycopeptide_glycan_type = pd.DataFrame(list(result.items()), columns=['glycopeptide', 'glycan_type'])
            self.glycopeptide_glycan_type = self.glycopeptide_glycan_type.sort_values(by='glycan_type', ascending=False)
        else:
            self.glycopeptide_glycan_type = pd.DataFrame()
        #
        temp_data = pd.DataFrame(self.data[['PeptideSequence', 'Glycosite_Position']])
        temp_data = temp_data[~temp_data['Glycosite_Position'].isnull()]
        # temp_data['PeptideSequence'] = temp_data['PeptideSequence'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['Glycosite_Position'])
        temp_data = temp_data.drop_duplicates()
        temp_data.columns = ['glycopeptide', 'glycosite']
        self.glycopeptide_glycosite = temp_data.reset_index(drop=True)
        #
        self.data_manager.log_params('StrucGAP_GlycoSite', 'glycopeptide_site', {})
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycopeptide_glycosite_count', self.glycopeptide_glycosite_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycopeptide_glycosite_count', self.glycopeptide_glycan_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycopeptide_glycosite_count', self.glycopeptide_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'glycopeptide_glycosite', self.glycopeptide_glycosite)

        return self
    
    def specific_site(self):
        """
        Conducts in-depth analysis of individual glycosylation sites.
        
        Parameters:
            None.
        
        Returns:
            self.protein_glycosite_glycan_count
            self.protein_glycosite_glycan_type
            self.protein_glycosite_glycan_composition_count
            self.protein_glycosite_isoforms_count
            self.peptide_glycosite_glycan_count
            self.peptide_glycosite_glycan_type
            self.peptide_glycosite_glycan_composition_count
            self.peptide_glycosite_isoforms_count
            
        Return type:
            dataframe
        
        """
        #
        temp_data = pd.DataFrame(self.data[['ProteinID', 'Glycosite_Position', 'structure_coding']])
        temp_data = temp_data[~temp_data['structure_coding'].isnull()]
        temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
        temp_data = temp_data.dropna(subset=['ProteinID', 'Glycosite_Position', 'structure_coding'])
        result_df = (
            temp_data.groupby(['ProteinID', 'Glycosite_Position'])['structure_coding']
            .nunique()
            .reset_index(name='glycan_count')
        )
        result = {
            (row['ProteinID'], row['Glycosite_Position']): row['glycan_count']
            for _, row in result_df.iterrows()
        }
        self.protein_glycosite_glycan_count = pd.DataFrame(list(result.items()), columns=['Protein_Glycosite', 'glycan_type_count'])
        self.protein_glycosite_glycan_count[['ProteinID', 'Glycosite_Position']] = pd.DataFrame(self.protein_glycosite_glycan_count['Protein_Glycosite'].tolist(), index=self.protein_glycosite_glycan_count.index)
        self.protein_glycosite_glycan_count = self.protein_glycosite_glycan_count.drop(columns='Protein_Glycosite') 
        self.protein_glycosite_glycan_count = self.protein_glycosite_glycan_count[['ProteinID', 'Glycosite_Position', 'glycan_type_count']]
        self.protein_glycosite_glycan_count = self.protein_glycosite_glycan_count.sort_values(by='glycan_type_count', ascending=False)
        self.map_protein_to_gene(self.protein_glycosite_glycan_count, 'ProteinID')
        #
        temp_data = pd.DataFrame(self.data[['ProteinID', 'Glycosite_Position', 'structure_coding']])
        temp_data = temp_data[~temp_data['structure_coding'].isnull()]
        temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
        temp_data = temp_data.dropna(subset=['ProteinID', 'Glycosite_Position', 'structure_coding'])
        result_df = (
            temp_data.groupby(['ProteinID', 'Glycosite_Position', 'structure_coding'])
            .size()
            .reset_index(name='count')
        )
        result = result_df.values.tolist()
        self.protein_glycosite_glycan_type = pd.DataFrame(result, columns=['ProteinID', 'Glycosite_Position', 'glycan_type', 'count'])
        self.protein_glycosite_glycan_type = self.protein_glycosite_glycan_type.sort_values(by='count', ascending=False)
        self.protein_glycosite_glycan_type = self.protein_glycosite_glycan_type.drop_duplicates()
        self.map_protein_to_gene(self.protein_glycosite_glycan_type, 'ProteinID')
        #
        temp_data = pd.DataFrame(self.data[['ProteinID', 'Glycosite_Position', 'GlycanComposition']])
        temp_data = temp_data[~temp_data['GlycanComposition'].isnull()]
        temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
        temp_data = temp_data.dropna(subset=['ProteinID', 'Glycosite_Position', 'GlycanComposition'])
        result_df = (
            temp_data.groupby(['ProteinID', 'Glycosite_Position', 'GlycanComposition'])
            .size()
            .reset_index(name='count')
        )
        result = result_df.values.tolist()
        self.protein_glycosite_glycan_composition_count = pd.DataFrame(result, columns=['ProteinID', 'Glycosite_Position', 'GlycanComposition', 'count'])
        self.protein_glycosite_glycan_composition_count = self.protein_glycosite_glycan_composition_count.sort_values(by='count', ascending=False)
        self.protein_glycosite_glycan_composition_count = self.protein_glycosite_glycan_composition_count.drop_duplicates()
        self.map_protein_to_gene(self.protein_glycosite_glycan_composition_count, 'ProteinID')
        #
        temp_data = pd.DataFrame(self.data[['ProteinID', 'Glycosite_Position', 'GlycanComposition', 'structure_coding']])
        temp_data = temp_data[~temp_data['structure_coding'].isnull()]
        temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
        temp_data = temp_data.dropna(subset=['ProteinID', 'Glycosite_Position', 'GlycanComposition', 'structure_coding'])
        result_df = (
            temp_data.groupby(['ProteinID', 'Glycosite_Position', 'GlycanComposition', 'structure_coding'])
            .size()
            .reset_index(name='count')
        )
        result = result_df.values.tolist()
        self.protein_glycosite_isoforms_count = pd.DataFrame(result, columns=['ProteinID', 'Glycosite_Position', 'GlycanComposition', 'isoforms', 'count'])
        self.protein_glycosite_isoforms_count = self.protein_glycosite_isoforms_count.sort_values(by='count', ascending=False)
        self.protein_glycosite_isoforms_count = self.protein_glycosite_isoforms_count.drop_duplicates()
        self.map_protein_to_gene(self.protein_glycosite_isoforms_count, 'ProteinID')
        
        #
        temp_data = pd.DataFrame(self.data[['PeptideSequence', 'Glycosite_Position', 'structure_coding']])
        temp_data = temp_data[~temp_data['structure_coding'].isnull()]
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode('Glycosite_Position')
        temp_data = temp_data.dropna(subset=['PeptideSequence', 'Glycosite_Position', 'structure_coding'])
        result_df = (
            temp_data.groupby(['PeptideSequence', 'Glycosite_Position'])['structure_coding']
            .nunique()
            .reset_index(name='glycan_count')
        )
        result = {
            (row['PeptideSequence'], row['Glycosite_Position']): row['glycan_count']
            for _, row in result_df.iterrows()
        }
        self.peptide_glycosite_glycan_count = pd.DataFrame(list(result.items()), columns=['Peptide_Glycosite', 'glycan_type_count'])
        self.peptide_glycosite_glycan_count[['PeptideSequence', 'Glycosite_Position']] = pd.DataFrame(self.peptide_glycosite_glycan_count['Peptide_Glycosite'].tolist(), index=self.peptide_glycosite_glycan_count.index)
        self.peptide_glycosite_glycan_count = self.peptide_glycosite_glycan_count.drop(columns='Peptide_Glycosite') 
        self.peptide_glycosite_glycan_count = self.peptide_glycosite_glycan_count[['PeptideSequence', 'Glycosite_Position', 'glycan_type_count']]
        self.peptide_glycosite_glycan_count = self.peptide_glycosite_glycan_count.sort_values(by='glycan_type_count', ascending=False)
        #
        temp_data = pd.DataFrame(self.data[['PeptideSequence', 'Glycosite_Position', 'structure_coding']])
        temp_data = temp_data[~temp_data['structure_coding'].isnull()]
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode('Glycosite_Position')
        temp_data = temp_data.dropna(subset=['PeptideSequence', 'Glycosite_Position', 'structure_coding'])
        result_df = (
            temp_data.groupby(['PeptideSequence', 'Glycosite_Position', 'structure_coding'])
            .size()
            .reset_index(name='count')
        )
        result = result_df.values.tolist()
        self.peptide_glycosite_glycan_type = pd.DataFrame(result, columns=['PeptideSequence', 'Glycosite_Position', 'glycan_type', 'count'])
        self.peptide_glycosite_glycan_type = self.peptide_glycosite_glycan_type.sort_values(by='count', ascending=False)
        self.peptide_glycosite_glycan_type = self.peptide_glycosite_glycan_type.drop_duplicates()
        #
        temp_data = pd.DataFrame(self.data[['PeptideSequence', 'Glycosite_Position', 'GlycanComposition']])
        temp_data = temp_data[~temp_data['GlycanComposition'].isnull()]
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode('Glycosite_Position')
        temp_data = temp_data.dropna(subset=['PeptideSequence', 'Glycosite_Position', 'GlycanComposition'])
        result_df = (
            temp_data.groupby(['PeptideSequence', 'Glycosite_Position', 'GlycanComposition'])
            .size()
            .reset_index(name='count')
        )
        result = result_df.values.tolist()
        self.peptide_glycosite_glycan_composition_count = pd.DataFrame(result, columns=['PeptideSequence', 'Glycosite_Position', 'GlycanComposition', 'count'])
        self.peptide_glycosite_glycan_composition_count = self.peptide_glycosite_glycan_composition_count.sort_values(by='count', ascending=False)
        self.peptide_glycosite_glycan_composition_count = self.peptide_glycosite_glycan_composition_count.drop_duplicates()
        #
        temp_data = pd.DataFrame(self.data[['PeptideSequence', 'Glycosite_Position', 'GlycanComposition', 'structure_coding']])
        temp_data = temp_data[~temp_data['structure_coding'].isnull()]
        temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
        temp_data = temp_data.explode('Glycosite_Position')
        temp_data = temp_data.dropna(subset=['PeptideSequence', 'Glycosite_Position', 'GlycanComposition', 'structure_coding'])
        result_df = (
            temp_data.groupby(['PeptideSequence', 'Glycosite_Position', 'GlycanComposition', 'structure_coding'])
            .size()
            .reset_index(name='count')
        )
        result = result_df.values.tolist()
        self.peptide_glycosite_isoforms_count = pd.DataFrame(result, columns=['PeptideSequence', 'Glycosite_Position', 'GlycanComposition', 'isoforms', 'count'])
        self.peptide_glycosite_isoforms_count = self.peptide_glycosite_isoforms_count.sort_values(by='count', ascending=False)
        self.peptide_glycosite_isoforms_count = self.peptide_glycosite_isoforms_count.drop_duplicates()
        #
        self.data_manager.log_params('StrucGAP_GlycoSite', 'specific_site', {})
        self.data_manager.log_output('StrucGAP_GlycoSite', 'protein_glycosite_glycan_count', self.protein_glycosite_glycan_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'protein_glycosite_glycan_type', self.protein_glycosite_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'protein_glycosite_glycan_composition_count', self.protein_glycosite_glycan_composition_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'protein_glycosite_isoforms_count', self.protein_glycosite_isoforms_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'peptide_glycosite_glycan_count', self.peptide_glycosite_glycan_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'peptide_glycosite_glycan_type', self.peptide_glycosite_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'peptide_glycosite_glycan_composition_count', self.peptide_glycosite_glycan_composition_count)
        self.data_manager.log_output('StrucGAP_GlycoSite', 'peptide_glycosite_isoforms_count', self.peptide_glycosite_isoforms_count)

        return self
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_GlycoSite.xlsx'), engine='xlsxwriter') as writer:
            self.glycoprotein_glycosite_count.to_excel(writer, sheet_name='glycoprotein_glycosite_count'[:31])
            self.glycoprotein_glycan_count.to_excel(writer, sheet_name='glycoprotein_glycan_count'[:31])
            self.glycoprotein_glycan_type.to_excel(writer, sheet_name='glycoprotein_glycan_type'[:31])
            self.glycoprotein_glycosite.to_excel(writer, sheet_name='glycoprotein_glycosite'[:31])
            
            self.glycopeptide_glycosite_count.to_excel(writer, sheet_name='glycopeptide_glycosite_count'[:31])
            self.glycopeptide_glycan_count.to_excel(writer, sheet_name='glycopeptide_glycan_count'[:31])
            self.glycopeptide_glycan_type.to_excel(writer, sheet_name='glycopeptide_glycan_type'[:31])
            self.glycopeptide_glycosite.to_excel(writer, sheet_name='glycopeptide_glycosite'[:31])
            
            if 'structure_coding' in self.data.columns:
                self.protein_glycosite_glycan_count.to_excel(writer, sheet_name='protein_glycosite_glycan_count'[:31])
                self.protein_glycosite_glycan_type.to_excel(writer, sheet_name='protein_glycosite_glycan_type'[:31])
                self.protein_glycosite_glycan_composition_count.to_excel(writer, sheet_name='protein_glycosite_glycan_composition_count'[:31])
                self.protein_glycosite_isoforms_count.to_excel(writer, sheet_name='protein_glycosite_isoforms_count'[:31])
                
                self.peptide_glycosite_glycan_count.to_excel(writer, sheet_name='peptide_glycosite_glycan_count'[:31])
                self.peptide_glycosite_glycan_type.to_excel(writer, sheet_name='peptide_glycosite_glycan_type'[:31])
                self.peptide_glycosite_glycan_composition_count.to_excel(writer, sheet_name='peptide_glycosite_glycan_composition_count'[:31])
                self.peptide_glycosite_isoforms_count.to_excel(writer, sheet_name='peptide_glycosite_isoforms_count'[:31])



