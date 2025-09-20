Get Started
===========

Dependencies
------------------------

This package relies on the following standard Python libraries and third-party packages:

.. code-block:: python

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

Note: We recommend that users verify the availability of all required dependencies before using this package to avoid runtime errors.


Get Started
--------------------

Here is a basic example of how to use StrucGAP. We strongly recommend using the output results from StrucGP as the input for StrucGAP. Although we also provide data processing pipelines for search results from MSFragger-Glyco, pGlyco3, and Glyco-Decipher, the information contained in the outputs of these three search engines is very limited. As a result, only the StrucGAP_GlycanStructure or StrucGAP_GlycoSite modules can be used for subsequent analyses.

.. code-block:: python

    from strucgap.preprocess import StrucGAP_Preprocess
    from strucgap.glycanstructure import StrucGAP_GlycanStructure
    from strucgap.glycosite import StrucGAP_GlycoSite
    from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
    from strucgap.functionannotation import StrucGAP_FunctionAnnotation
    from strucgap.glyconetwork import StrucGAP_GlycoNetwork
    from strucgap.datavisualization import StrucGAP_DataVisualization
    from strucgap.insighttracker import StrucGAP_InsightTracker
    import os
    # Initialization
    data_manager = StrucGAP_InsightTracker()
    # Setting the result storage path (folder)
    os.chdir('tests/')
    # Read it if you've already done the analysis
    data_manager.read_pickle()
    # Import data
    module1 = StrucGAP_Preprocess(data_dir="tests/mouse uterus.xlsx",
                      data_sheet_name = '1 PSM',
                      sample_group_data_dir = 'tests/sample_group.xlsx',
                      branch_list_dir = "tests/branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager)
    module1.data_cleaning(data_type='tmt')
    module1.cv_raw(threshold='no')
    module1.fdr(feature_type='no')
    module1.outliers(abundance_ratio=[1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464],
                 samplewise_normalization = False)
    module1.cv(threshold = 'no')
    module1.psm(psm_number = 'no')
    # Using glytoucan = True and biosynthetic_pathways = True is a very time-consuming task, due to the limitations of the GlyTouCan and KEGG APIs. Please be patient when enabling these two annotations. If you prefer faster execution, set both options to False.
    module1.annotation(glytoucan = True, glytoucan_structure = True, glytoucan_wurcs_file = "tests/glycosmos_glycans_wurcs.csv", biosynthetic_pathways = True, glycobiology_filter = True)
    module1.output() 
    # ... Other analysis

Each module can be instantiated and run independently depending on your workflow.
