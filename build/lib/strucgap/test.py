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

os.chdir('D:\\doctor\\analysisys\\StrucGAP')

##
from strucgap.preprocess import StrucGAP_Preprocess
from strucgap.glycanstructure import StrucGAP_GlycanStructure
from strucgap.glycosite import StrucGAP_GlycoSite
from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
from strucgap.functionannotation import StrucGAP_FunctionAnnotation
from strucgap.glyconetwork import StrucGAP_GlycoNetwork
from strucgap.datavisualization import StrucGAP_DataVisualization
from strucgap.insighttracker import StrucGAP_InsightTracker


data_manager = StrucGAP_InsightTracker()
# 
module1 = StrucGAP_Preprocess(data_dir="D:\\doctor\\zzd\\20250627\\Rat_aging_review_thymus_psm.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'D:\\doctor\\analysisys\\data\\sample_group.xlsx',
                      branch_list_dir = "D:\\doctor\\wyq\\branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager)
module1.data_cleaning(data_type='tmt')
# 1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464
module1.fdr(feature_type='no')
module1.outliers(abundance_ratio=[1.240003449, 0, 1.344387558, 0, 1.576533442, 0, 1, 0, 1.956346409, 1.517000766])
module1.cv(threshold = 'no')
module1.psm(psm_number = 'no')
module1.output() 

# # 
# module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
# module2.statistics(remove_oligo_mannose = False) 
# module2.structure_statistics()
# module2.lacdinac()
# module2.cor()
# module2.isoforms()
# module2.output()

# #
# module3 = StrucGAP_GlycoSite(module1, data_manager=data_manager)
# module3.glycoprotein_site()
# module3.glycopeptide_site()
# module3.specific_site()
# module3.output()

#
module4 = StrucGAP_GlycoPeptideQuant(module1, data_type = 'psm_filtered', data_manager=data_manager)
module4.statistics()
module4.statistics_index()
module4.differential_analysis(pvalue_type='pvalue_ttest')
module4.threshold_variation_analysis(pvalue_type='pvalue_ttest',statistic_index='fc')
module4.glycopeptide_glycosite_glycan_variation()
module4.glycoprotein_glycosite_glycan_variation()
# module4.output()

#
module6 = StrucGAP_GlycoNetwork(module4, data_manager=data_manager)
module6.proteomic(protein_data_dir="D:\\doctor\\zzd\\20250716\\Protein_all database.xlsx",
                  data_sheet_name = 'Sheet1', cv = 'no')
module6.phosphorylation(phospho_data_dir="D:\doctor\zzd\_Rat_Phospho_mixThymus_TMT6c_.xlsx",
                        data_sheet_name='PeptideGroups', cv = 'no')
module6.glycosyltransferases(glycosyltransferases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", 
                             data_sheet_name="glycosyltransferases")
module6.glycosidases(glycosidases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", 
                     data_sheet_name='glycosidases')
module6.sialyltransferases()
module6.fucosyltransferase()
module6.glycan_binding_protein()
module6.output()

