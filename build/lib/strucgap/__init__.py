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

# os.chdir('D:\\doctor\\StrucGAP\\strucgap')

##
from .preprocess import StrucGAP_Preprocess
from .glycanstructure import StrucGAP_GlycanStructure
from .glycosite import StrucGAP_GlycoSite
from .glycopeptidequant import StrucGAP_GlycoPeptideQuant
from .functionannotation import StrucGAP_FunctionAnnotation
from .glyconetwork import StrucGAP_GlycoNetwork
from .datavisualization import StrucGAP_DataVisualization
from .insighttracker import StrucGAP_InsightTracker

