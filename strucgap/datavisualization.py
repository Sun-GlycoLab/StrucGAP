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
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from packaging import version
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'
## 数据可视化和报告生成模块
class StrucGAP_DataVisualization:
    """
    Parameters:
        data_manager: Data manager instance, such as 'data_manager' if data_manager = StrucGAP_InsightTracker().
    
    """
    def __init__(self, data_manager):
        self.data_manager = data_manager    
        self.data_manager.register_module('StrucGAP_DataVisualization', self, {})
        self.data_manager.log_params('StrucGAP_DataVisualization', '', {}) 
        
        # 初始化其他原有数据存储结构
        self.analysis_records = []  # 假设已有分析记录存储
        # 新增图表组合管理结构
        self.figure_collections: Dict[str, List[dict]] = {}  # {figure_name: [figure_meta1, ...]}
        self.current_figure = None  # 当前操作的figure名称
        
    def trim_white_border(self, image):
        """
        An auxiliary function called by other functions to trim the white border from an image.

        Parameters:
            image: The image to be processed, which may contain a white border.
    
        Returns:
            A cropped image with the white border removed.
    
        Return type:
            PIL Image object
            
        """
        if image.mode != 'RGB':
            image = image.convert('RGB')
        
        bg = Image.new(image.mode, image.size, image.getpixel((0,0)))
        diff = ImageChops.difference(image, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
        if bbox:
            return image.crop(bbox)
    
    def select_palette(self, data, child_column):
        """
        An interactive auxiliary function to select a qualitative color palette based on the number of unique values in a specified column.
    
        Parameters:
            data: The dataset containing categorical values.
            child_column: The column name whose unique value count determines the palette size.
    
        Returns:
            The selected qualitative color palette object, which includes color definitions.
    
        Return type:
            qualitative color palette (e.g., from colorbrewer or cartocolors)
            
        """
        # 获取 brewer_qualitative 和 carto_qualitative 的调色板列表
        brewer_palettes = [attr for attr in dir(brewer_qualitative) if not attr.startswith("__")]
        carto_palettes = [attr for attr in dir(carto_qualitative) if not attr.startswith("__")]
    
        # 合并两个调色板列表
        all_palettes = brewer_palettes + carto_palettes
    
        # 输出所有可用的调色板
        print("Available palettes (from both brewer and cartocolors):")
        print(all_palettes)
        
        # 动态根据数据生成提示信息
        print(f'You should select palette that ends with {len(data[child_column].unique())}')
        
        # 用户输入调色板名称
        palette_input = input("Please input the palette name (e.g., 'Set3_9'): ")
    
        # 根据输入判断并导入相应调色板
        if palette_input in brewer_palettes:
            palette = getattr(brewer_qualitative, palette_input)
            print(f'You selected {palette_input} from colorbrewer.qualitative with colors: {palette.colors}')
        elif palette_input in carto_palettes:
            palette = getattr(carto_qualitative, palette_input)
            print(f'You selected {palette_input} from cartocolors.qualitative with colors: {palette.colors}')
        else:
            raise ValueError(f"Palette {palette_input} not found in available palettes.")
        
        return palette
    
    def radar(self, *data_names, columns, colors = None, screen_column = None, screen_values = None,
              shape = 'circle', 
              text_color = 'black',
              text_font_weight = 'bold',
              text_font_size = 12,
              text_split = 10,
              splitline_color = '#001F78',
              splitline_type = 'dashed',
              splitline_width = 3,
              splitline_opacity = 0.5,
              axisline_width = 2,
              axisline_opacity = 0.3,
              axisline_color = 'grey',
              background_opacity = 0.2,
              background_color = 'white',
              symbol = 'none',
              area_opacity = 0.5,
              if_unique = True,
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              subfolder='plot',
              filename = None,
              figure_description = 'Radar plot',
              ):
        """
        Generate a radar chart visualization from one or multiple DataFrames.
    
        This method computes values per column (either unique counts or raw counts, 
        optionally filtered by conditions) and visualizes them on a radar plot. 
        The output figure is saved in multiple formats (HTML, SVG, PDF, PNG), 
        and plotting parameters are logged using `data_manager`.
    
        Args:
            *data_names: One or more input DataFrames or variable names (str). 
                Strings are `eval`-ed into DataFrames; objects can be passed directly.
            columns (list of str): Columns to be used as radar axes.
            colors (list of str, optional): Colors for each dataset. Defaults to a 
                predefined color palette.
            screen_column (str, optional): Column name for filtering data. If provided, 
                `if_unique` is automatically set to False.
            screen_values (list, optional): Values of `screen_column` to include in 
                the radar plot.
            shape (str, optional): Radar shape, "circle" or "polygon". Default is "circle".
            text_color (str, optional): Color of axis labels. Default is "black".
            text_font_weight (str, optional): Font weight for axis labels. 
                One of {"normal", "bold", "bolder", "lighter"}. Default is "bold".
            text_font_size (int, optional): Font size for axis labels. Default is 12.
            text_split (int, optional): Maximum character length per line for axis labels. 
                Labels are wrapped if > 0. Default is 10.
            splitline_color (str, optional): Color of grid split lines. Default "#001F78".
            splitline_type (str, optional): Grid line style, one of {"solid", "dashed", "dotted"}. 
                Default is "dashed".
            splitline_width (int, optional): Width of split lines. Default is 3.
            splitline_opacity (float, optional): Opacity of split lines [0–1]. Default is 0.5.
            axisline_width (int, optional): Width of axis lines. Default is 2.
            axisline_opacity (float, optional): Opacity of axis lines [0–1]. Default is 0.3.
            axisline_color (str, optional): Color of axis lines. Default is "grey".
            background_opacity (float, optional): Opacity of background area [0–1]. Default is 0.2.
            background_color (str, optional): Color of background area. Default is "white".
            symbol (str, optional): Marker symbol for radar lines. 
                Options: {"circle","rect","roundRect","triangle","diamond","pin","arrow","none"}. 
                Default is "none".
            area_opacity (float, optional): Opacity of filled radar area [0–1]. Default is 0.5.
            if_unique (bool, optional): If True, radar values are based on counts of unique 
                values per column; otherwise use total counts. Default is True.
            legend (list of str or str, optional): Legend labels. If None, they are inferred 
                from `data_names` or `screen_values`.
            legend_font_size (int, optional): Font size for legend text. Default is 12.
            plot_title (str, optional): Title of the radar plot. Default is None.
            subfolder (str, optional): Subfolder under `./plot/` for saving output. Default is "plot".
            filename (str, optional): Prefix for saved file names. A timestamp is appended automatically.
            figure_description (str, optional): Description string recorded with the output. 
                Default is "Radar plot".
    
        Returns:
            dict: A dictionary with:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.
    
        Side Effects:
            - Saves radar chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output files via `data_manager`.
    
        Raises:
            ValueError: If invalid parameter values are provided or filtering fails.
        """
        if isinstance(data_names, str):
            data_names = [data_names]
            
        # dataframes = [
        #     eval(name) if isinstance(name, str) else name
        #     for name in data_names
        # ]
        # dataframe_names = [name.split('.')[-1] for name in data_names]
        
        dataframes = []
        dataframe_names = []
        for idx, name in enumerate(data_names):
            if isinstance(name, str):
                df = eval(name)
                dataframes.append(df)
                dataframe_names.append(name.split('.')[-1])
            else:  # 直接是 DataFrame
                dataframes.append(name)
                dataframe_names.append(f"df{idx}")
        
        if screen_column is not None:
            if_unique = False
        
        #
        radar_data = []
        schema = []
        #
        if shape not in ['circle', 'polygon']:
            shape = 'circle'
        
        if text_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            text_font_weight = 'bold' 
        
        if splitline_type not in ['solid', 'dashed', 'dotted']:
            splitline_type = 'dashed'
            
        splitline_width = int(splitline_width)  
        
        splitline_opacity = float(splitline_opacity)
        if (splitline_opacity < 0) or (splitline_opacity > 1):
            splitline_opacity = float(0.5)
            
        axisline_width = int(axisline_width)  
        
        axisline_opacity = float(axisline_opacity)
        if (axisline_opacity < 0) or (axisline_opacity > 1):
            axisline_opacity = float(0.5)
        
        background_opacity = float(background_opacity)
        if (background_opacity < 0) or (background_opacity > 1):
            background_opacity = float(0.5)
        
        if symbol not in ['circle','rect','roundRect','triangle','diamond','pin','arrow','none']:
            symbol = 'none'
        
        area_opacity = float(area_opacity)
        if (area_opacity < 0) or (area_opacity > 1):
            area_opacity = float(0.5)
            
        if if_unique not in [True, False]:
            if_unique = True
        #
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]

        max_values = {}
        for column in columns:
            if if_unique:
                max_values[column] = max([df[column].dropna().unique().shape[0] * 1.1 for df in dataframes])
            else:
                if screen_column is None:
                    max_values[column] = max([df[column].dropna().shape[0] * 1.1 for df in dataframes])
                else:
                    max_values[column] = [max(df[df[screen_column].isin(screen_values)][column]) * 1.1 for df in dataframes][0]

        if text_split != 0:
            wrapped_labels = ['\n'.join([label[i:i+text_split] for i in range(0, len(label), text_split)]) 
                      for label in columns]
        
        for column in wrapped_labels:
            schema.append(opts.RadarIndicatorItem(name=column, max_=max(max_values.values()))) # max_ = max_values[column]

        for df in dataframes:
            if if_unique:
                row_counts = [df[col].dropna().unique().shape[0] for col in columns]
                radar_data.append([row_counts])
            else:
                if screen_column is None:
                    row_counts = [df[col].dropna().shape[0] for col in columns]
                    radar_data.append([row_counts])
                else:
                    for value in screen_values:
                        radar_data.append([df[df[screen_column]==value][columns].values[0].tolist()])

        radar = (
            Radar(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add_schema(
                schema = schema,
                shape = shape,  
                radius = "60%",
                textstyle_opts = opts.TextStyleOpts(color = text_color, 
                                                    font_family = 'Arial',
                                                    font_weight = text_font_weight,
                                                    font_size = text_font_size),
                # 分割线：同心圆
                splitline_opt = opts.SplitLineOpts(                            
                    is_show = True,
                    linestyle_opts = opts.LineStyleOpts(
                        is_show = True,
                        width = splitline_width,
                        opacity = splitline_opacity,
                        curve = 1,
                        type_ = splitline_type, 
                        color = splitline_color
                    )
                ),
                # 坐标轴轴线
                axisline_opt = opts.AxisLineOpts(                              
                    is_show = True,
                    linestyle_opts = opts.LineStyleOpts(
                        width = axisline_width,
                        opacity = axisline_opacity,
                        color = axisline_color
                        )
                    ),
                # 坐标轴刻度
                axistick_opt = opts.AxisTickOpts(                              
                    is_show = False
                    ),
                # 坐标轴刻度标签
                axislabel_opt = opts.LabelOpts(                                
                    is_show = False
                    ),
                # 被分割的区域
                splitarea_opt = opts.SplitAreaOpts(
                    is_show = True, 
                    areastyle_opts = opts.AreaStyleOpts(
                        opacity = background_opacity,
                        color = {
                            'type': 'radial',
                            'x': 0.5,
                            'y': 0.5,
                            'r': 0.5,
                            'colorStops': [
                                {'offset': 0, 'color': background_color},  # 0%处的颜色
                                {'offset': 1, 'color': background_color}  # 100%处的颜色
                            ],
                            'global': False  
                        }
                    )
                )
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
        )
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        if legend is None:
            if screen_column is None:
                legend = dataframe_names
            if screen_column is not None:
                legend = screen_values
        
        for idx, data in enumerate(radar_data):
            radar.add(
                series_name = legend[idx],
                data = data,
                linestyle_opts = opts.LineStyleOpts(color = colors[idx], width = 3),  
                symbol = symbol, # 'circle','rect','roundRect','triangle','diamond','pin','arrow','none'
                color_by = 'series',
                color = colors[idx],
                areastyle_opts = opts.AreaStyleOpts(opacity = area_opacity)
            )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_radar_{timestamp}.html")
        radar.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_radar_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_radar_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_radar_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'radar', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'radar', f"{'_'.join(dataframe_names)}_radar.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def polar1(self, *data_names, columns, number_column=None, colors=None,
                  splitline_color = ['#001F78'],
                  splitline_type = 'dashed',
                  splitline_width = 3,
                  splitline_opacity = 0.5,
                  axisline_width = 2,
                  axisline_opacity = 0.3,
                  axisline_color = 'grey',
                  background_opacity = 0.2,
                  background_color = 'white',
                  symbol = 'none',
                  symbol_size = 0,
                  area_opacity = 0.5,
                  if_unique = True,
                  radiusaxis_label_show = False,
                  radiusaxis_label_color = 'black',
                  radiusaxis_label_font_weight = 'normal',
                  radiusaxis_label_font_size = 12,
                  angleaxis_label_show = True,
                  angleaxis_label_color = 'black',
                  angleaxis_label_font_weight = 'bold',
                  angleaxis_label_font_size = 12,
                  legend = None,
                  legend_font_size = 12,
                  plot_title = None,
                  subfolder='plot',
                  filename = None,
                  figure_description = 'Polar plot',
                  ):
        """
        Generate a polar plot visualization from one or more DataFrames.
    
        This method creates a polar coordinate chart using the specified `columns`
        (categorical variables) or `number_column` (numeric values). It supports 
        both multiple dataset comparison and single dataset visualization. The 
        figure is exported in multiple formats (HTML, SVG, PDF, PNG), and plotting 
        parameters are logged using `data_manager`.
    
        Behavior:
            - Multiple DataFrames:
                Each DataFrame is summarized by counting values per column 
                (unique counts if `if_unique=True`, otherwise raw counts). 
                Each dataset is plotted as a separate series in polar coordinates.
            - Single DataFrame + `number_column`:
                Plots the numeric values in `number_column` against the categories 
                in `columns[0]`.
    
        Args:
            *data_names: Input DataFrames or variable names (str). Strings are 
                evaluated with `eval()`, or DataFrames can be passed directly.
            columns (list of str): Categorical columns to use for polar axes. 
                If length is 1, `number_column` must be provided.
            number_column (str, optional): Column containing numeric values when 
                `columns` has only one entry. Default is None.
            colors (list of str, optional): Colors for series. Defaults to a preset palette.
            splitline_color (list of str, optional): Colors for radial split lines. Default ["#001F78"].
            splitline_type (str, optional): Style of split lines. One of {"solid", "dashed", "dotted"}. 
                Default "dashed".
            splitline_width (int, optional): Width of split lines. Default 3.
            splitline_opacity (float, optional): Opacity of split lines [0–1]. Default 0.5.
            axisline_width (int, optional): Width of axis lines. Default 2.
            axisline_opacity (float, optional): Opacity of axis lines [0–1]. Default 0.3.
            axisline_color (str, optional): Color of axis lines. Default "grey".
            background_opacity (float, optional): Opacity of background fill [0–1]. Default 0.2.
            background_color (str, optional): Background color. Default "white".
            symbol (str, optional): Marker symbol for data points. 
                Options: {"circle","rect","roundRect","triangle","diamond","pin","arrow","none"}. 
                Default "none".
            symbol_size (int, optional): Size of marker symbols. Default 0.
            area_opacity (float, optional): Opacity of filled area [0–1]. Default 0.5.
            if_unique (bool, optional): Whether to count unique values (`True`) 
                or total counts (`False`). Default True.
            radiusaxis_label_show (bool, optional): Show/hide radius axis labels. Default False.
            radiusaxis_label_color (str, optional): Color of radius axis labels. Default "black".
            radiusaxis_label_font_weight (str, optional): Font weight for radius axis labels. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            radiusaxis_label_font_size (int, optional): Font size for radius axis labels. Default 12.
            angleaxis_label_show (bool, optional): Show/hide angle axis labels. Default True.
            angleaxis_label_color (str, optional): Color of angle axis labels. Default "black".
            angleaxis_label_font_weight (str, optional): Font weight for angle axis labels. 
                One of {"normal","bold","bolder","lighter"}. Default "bold".
            angleaxis_label_font_size (int, optional): Font size for angle axis labels. Default 12.
            legend (list of str or str, optional): Legend labels. If None, inferred from input names. Default None.
            legend_font_size (int, optional): Font size for legend labels. Default 12.
            plot_title (str, optional): Title of the polar plot. Default None.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. Default "Polar plot".
    
        Returns:
            dict: 
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.
    
        Side Effects:
            - Saves polar chart as HTML, SVG, PDF, and PNG under `./plot/<subfolder>/`.
            - Logs parameters and outputs using `data_manager`.
    
        Raises:
            ValueError: If parameter values are invalid or required inputs are missing.
        """
        # 读取传入的dataframe
        if isinstance(data_names, str):
            data_names = [data_names]
        dataframes = [
            eval(name) if isinstance(name, str) else name
            for name in data_names
        ]
        dataframe_names = [name.split('.')[-1] for name in data_names]

        polar_data = []
        schema = []

        # 处理颜色和其他参数
        # splitline_color.insert(1, 'white')
        
        if splitline_type not in ['solid', 'dashed', 'dotted']:
            splitline_type = 'dashed'
            
        splitline_width = int(splitline_width)  
        splitline_opacity = float(splitline_opacity)
        if (splitline_opacity < 0) or (splitline_opacity > 1):
            splitline_opacity = float(0.5)
            
        axisline_width = int(axisline_width)  
        axisline_opacity = float(axisline_opacity)
        if (axisline_opacity < 0) or (axisline_opacity > 1):
            axisline_opacity = float(0.5)

        background_opacity = float(background_opacity)
        if (background_opacity < 0) or (background_opacity > 1):
            background_opacity = float(0.5)

        if symbol not in ['circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none']:
            symbol = 'none'

        area_opacity = float(area_opacity)
        if (area_opacity < 0) or (area_opacity > 1):
            area_opacity = float(0.5)

        if if_unique not in [True, False]:
            if_unique = True
            
        if radiusaxis_label_show not in [True, False]:
            radiusaxis_label_show = False
            
        if radiusaxis_label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            radiusaxis_label_font_weight = 'normal'   
        
        if angleaxis_label_show not in [True, False]:
            angleaxis_label_show = True
            
        if angleaxis_label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            angleaxis_label_font_weight = 'bold' 

        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                        "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(dataframes)]
        
        max_values = {}
        if len(columns) != 1:
            for column in columns:
                if if_unique:
                    max_values[column] = max([df[column].dropna().unique().shape[0] * 1.1 for df in dataframes])
                else:
                    max_values[column] = max([df[column].dropna().shape[0] * 1.1 for df in dataframes])
        else:
            max_values[columns[0]] = max(list(dataframes[0][number_column])) * 1.1

        # 处理 schema 并收集数据
        if len(dataframes) != 1:
            for df in dataframes:
                if if_unique:
                    row_counts = [df[col].dropna().unique().shape[0] for col in columns]
                else:
                    row_counts = [df[col].dropna().shape[0] for col in columns]
                polar_data.append([row_counts])
        else:
            polar_data = list(dataframes[0][number_column])
            
        if len(columns) == 1:
            columns = list(dataframes[0][columns[0]])

        # 创建极坐标图
        polar = (
            Polar(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add_schema(
                radiusaxis_opts=opts.RadiusAxisOpts(type_='category', 
                                                    data = columns,
                                                    name_location = 'middle',
                                                    splitline_opts = opts.SplitLineOpts(                            
                                                                    is_show = True,
                                                                    linestyle_opts = opts.LineStyleOpts(
                                                                        is_show = True,
                                                                        width = splitline_width,
                                                                        opacity = splitline_opacity,
                                                                        curve = 1,
                                                                        type_ = splitline_type, 
                                                                        color = splitline_color
                                                                    )
                                                                ),
                                                    splitarea_opts = opts.SplitAreaOpts(
                                                                    is_show = True, 
                                                                    areastyle_opts = opts.AreaStyleOpts(
                                                                        opacity = background_opacity,
                                                                        color = {
                                                                            'type': 'radial',
                                                                            'x': 0.5,
                                                                            'y': 0.5,
                                                                            'r': 0.5,
                                                                            'colorStops': [
                                                                                {'offset': 0, 'color': background_color},  # 0%处的颜色
                                                                                {'offset': 1, 'color': background_color}  # 100%处的颜色
                                                                            ],
                                                                            'global': False  
                                                                        }
                                                                    )
                                                                ),
                                                    axisline_opts = opts.AxisLineOpts(                              
                                                                    is_show = True,
                                                                    linestyle_opts = opts.LineStyleOpts(
                                                                        width = axisline_width,
                                                                        opacity = axisline_opacity,
                                                                        color = axisline_color
                                                                        )
                                                                    ),
                                                    axislabel_opts = opts.LabelOpts(                                
                                                                    is_show = radiusaxis_label_show,
                                                                    color = radiusaxis_label_color,
                                                                    font_family = 'Arial',
                                                                    font_weight = radiusaxis_label_font_weight,
                                                                    interval = 0,
                                                                    font_size = radiusaxis_label_font_size,
                                                                    ),
                                                    ),
                angleaxis_opts=opts.AngleAxisOpts(is_clockwise=True, 
                                                    max_=max(max_values.values()),
                                                    splitline_opts = opts.SplitLineOpts(                            
                                                                    is_show = True,
                                                                    linestyle_opts = opts.LineStyleOpts(
                                                                                    width = axisline_width,
                                                                                    opacity = axisline_opacity,
                                                                                    color = axisline_color
                                                                                    )
                                                                    ),
                                                    axisline_opts = opts.AxisLineOpts(                              
                                                                    is_show = False
                                                                    ),
                                                    axislabel_opts = opts.LabelOpts(                                
                                                                    is_show = angleaxis_label_show,
                                                                    color = angleaxis_label_color,
                                                                    font_family = 'Arial',
                                                                    font_weight = angleaxis_label_font_weight,
                                                                    interval = 0,
                                                                    font_size = angleaxis_label_font_size,
                                                                    ),
                                                    
                                                    
                                                    ),
                
            )
            .set_global_opts(
            title_opts=opts.TitleOpts(
                title = plot_title,
                pos_left="center",  # 标题居中
                pos_top="0%"       # 标题位置靠上
            ),
            legend_opts=opts.LegendOpts(
                pos_right="center",  # 图例靠右
                pos_top="top",
                textstyle_opts = opts.TextStyleOpts(
                                                font_family = 'Arial',
                                                font_weight = 'normal',
                                                font_size = legend_font_size),
            ),
            )
        )
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        if legend is None:
            if len(dataframes) != 1:
                legend = dataframe_names
            if len(dataframes) == 1:
                legend = [dataframe_names[0]]
                
        if len(dataframes) != 1:
            # 为每个数据集添加线条
            for idx, data in enumerate(polar_data):
                polar.add(
                    series_name = legend[idx],
                    data = data[0],
                    type_='bar',
                    # linestyle_opts = opts.LineStyleOpts(color=colors[idx], width=3),  
                    symbol = symbol,
                    symbol_size = symbol_size,
                    # color = colors[idx],
                    label_opts = opts.LabelOpts(                                
                        is_show = False
                        ),
                    areastyle_opts = opts.AreaStyleOpts(opacity = area_opacity)
                )
            
        if len(dataframes) == 1:
            polar.add(
                    series_name = legend[0],
                    data = polar_data,
                    type_='bar',
                    # linestyle_opts = opts.LineStyleOpts(color=colors[idx], width=3),  
                    symbol = symbol,
                    symbol_size = symbol_size,
                    # color = colors[idx],
                    label_opts = opts.LabelOpts(                                
                        is_show = False
                        ),
                    areastyle_opts = opts.AreaStyleOpts(opacity = area_opacity)
                )
            
        polar.set_colors(colors)
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        file_name = os.path.join(output_dir, f"{filename}_polar1_{timestamp}.html")
        polar.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_polar1_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_polar1_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_polar1_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'polar1', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'polar1',  f"{'_'.join(dataframe_names)}_polar1.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def polar2(self, *data_names, columns, colors=None,
              splitline_color = ['#001F78'],
              splitline_type = 'dashed',
              splitline_width = 3,
              splitline_opacity = 0.5,
              axisline_width = 2,
              axisline_opacity = 0.3,
              axisline_color = 'grey',
              background_opacity = 0.2,
              background_color = 'white',
              symbol = 'none',
              symbol_size = 0,
              area_opacity = 0.5,
              if_unique = True,
              radiusaxis_label_show = True,
              radiusaxis_label_color = 'black',
              radiusaxis_label_font_weight = 'normal',
              radiusaxis_label_font_size = 12,
              angleaxis_label_show = True,
              angleaxis_label_color = 'black',
              angleaxis_label_font_weight = 'bold',
              angleaxis_label_font_size = 12,
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              subfolder='plot',
              filename = None,
              figure_description = 'Polar plot',
              ):
        """
        Generate a stacked polar plot comparing multiple DataFrames across specified columns.

        This method creates a polar coordinate chart where:
        - The angle axis represents different input DataFrames.
        - The radius axis represents the selected `columns`.
        - Each column is plotted as a separate stacked bar series, showing either
            unique value counts or total counts depending on `if_unique`.

        The figure is exported in multiple formats (HTML, SVG, PDF, PNG), and plotting
        parameters are logged using `data_manager`.

        Args:
            *data_names: Input DataFrames or variable names (str). Strings are 
                evaluated with `eval()`, or DataFrames can be passed directly.
            columns (list of str): Column names to compare across datasets.
            colors (list of str, optional): Series colors. Defaults to a preset palette.
            splitline_color (list of str, optional): Colors of radial split lines. Default ["#001F78"].
            splitline_type (str, optional): Split line style, one of {"solid","dashed","dotted"}. Default "dashed".
            splitline_width (int, optional): Width of split lines. Default 3.
            splitline_opacity (float, optional): Opacity of split lines [0–1]. Default 0.5.
            axisline_width (int, optional): Width of axis lines. Default 2.
            axisline_opacity (float, optional): Opacity of axis lines [0–1]. Default 0.3.
            axisline_color (str, optional): Axis line color. Default "grey".
            background_opacity (float, optional): Background opacity [0–1]. Default 0.2.
            background_color (str, optional): Background color. Default "white".
            symbol (str, optional): Marker symbol for data points. 
                Options: {"circle","rect","roundRect","triangle","diamond","pin","arrow","none"}. 
                Default "none".
            symbol_size (int, optional): Marker symbol size. Default 0.
            area_opacity (float, optional): Opacity of filled areas [0–1]. Default 0.5.
            if_unique (bool, optional): Whether to count unique values (`True`) 
                or total counts (`False`). Default True.
            radiusaxis_label_show (bool, optional): Show/hide labels on radius axis. Default True.
            radiusaxis_label_color (str, optional): Radius axis label color. Default "black".
            radiusaxis_label_font_weight (str, optional): Font weight of radius labels. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            radiusaxis_label_font_size (int, optional): Radius axis label font size. Default 12.
            angleaxis_label_show (bool, optional): Show/hide labels on angle axis. Default True.
            angleaxis_label_color (str, optional): Angle axis label color. Default "black".
            angleaxis_label_font_weight (str, optional): Font weight of angle labels. Default "bold".
            angleaxis_label_font_size (int, optional): Angle axis label font size. Default 12.
            legend (list of str or str, optional): Legend labels. Defaults to `columns` if None.
            legend_font_size (int, optional): Legend font size. Default 12.
            plot_title (str, optional): Title of the polar plot. Default None.
            subfolder (str, optional): Subfolder under `./plot/` for outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp appended automatically.
            figure_description (str, optional): Description string logged with output. Default "Polar plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves polar chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Raises:
            ValueError: If parameter values are invalid or required inputs are missing.
        """
        # 读取传入的dataframe
        if isinstance(data_names, str):
            data_names = [data_names]
        dataframes = [
            eval(name) if isinstance(name, str) else name
            for name in data_names
        ]
        dataframe_names = [name.split('.')[-1] for name in data_names]

        polar_data = []
        schema = []

        # 处理颜色和其他参数
        # splitline_color.insert(1, 'white')
        
        if splitline_type not in ['solid', 'dashed', 'dotted']:
            splitline_type = 'dashed'
            
        splitline_width = int(splitline_width)  
        splitline_opacity = float(splitline_opacity)
        if (splitline_opacity < 0) or (splitline_opacity > 1):
            splitline_opacity = float(0.5)
            
        axisline_width = int(axisline_width)  
        axisline_opacity = float(axisline_opacity)
        if (axisline_opacity < 0) or (axisline_opacity > 1):
            axisline_opacity = float(0.5)

        background_opacity = float(background_opacity)
        if (background_opacity < 0) or (background_opacity > 1):
            background_opacity = float(0.5)

        if symbol not in ['circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none']:
            symbol = 'none'

        area_opacity = float(area_opacity)
        if (area_opacity < 0) or (area_opacity > 1):
            area_opacity = float(0.5)

        if if_unique not in [True, False]:
            if_unique = True
            
        if radiusaxis_label_show not in [True, False]:
            radiusaxis_label_show = False
            
        if radiusaxis_label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            radiusaxis_label_font_weight = 'normal'   
        
        if angleaxis_label_show not in [True, False]:
            angleaxis_label_show = True
            
        if angleaxis_label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            angleaxis_label_font_weight = 'bold' 

        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(columns)]

        max_values = {}
        for column in columns:
            if if_unique:
                max_values[column] = max([df[column].dropna().unique().shape[0] * 1.1 for df in dataframes])
            else:
                max_values[column] = max([df[column].dropna().shape[0] * 1.1 for df in dataframes])

        # 处理 schema 并收集数据
        for col in columns:
            if if_unique:
                row_counts = [df[col].dropna().unique().shape[0] for df in dataframes]
            else:
                row_counts = [df[col].dropna().shape[0] for df in dataframes]
            polar_data.append([row_counts])    

        # 创建极坐标图
        polar = (
            Polar(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add_schema(
                angleaxis_opts=opts.AngleAxisOpts(type_='category', 
                                                  data = dataframe_names,
                                                  is_clockwise=True, 
                                                  splitline_opts = opts.SplitLineOpts(                            
                                                                    is_show = True,
                                                                    linestyle_opts = opts.LineStyleOpts(
                                                                                    width = axisline_width,
                                                                                    opacity = axisline_opacity,
                                                                                    color = axisline_color
                                                                                    )
                                                                    ),
                                                  axisline_opts = opts.AxisLineOpts(                              
                                                                    is_show = False
                                                                    ),
                                                  axislabel_opts = opts.LabelOpts(                                
                                                                    is_show = angleaxis_label_show,
                                                                    color = angleaxis_label_color,
                                                                    font_family = 'Arial',
                                                                    font_weight = angleaxis_label_font_weight,
                                                                    interval = 0,
                                                                    font_size = angleaxis_label_font_size,
                                                                    ),
                                                    ),
                radiusaxis_opts=opts.RadiusAxisOpts(
                                                  splitline_opts = opts.SplitLineOpts(                            
                                                                    is_show = True,
                                                                    linestyle_opts = opts.LineStyleOpts(
                                                                                    width = splitline_width,
                                                                                    opacity = splitline_opacity,
                                                                                    color = splitline_color,
                                                                                    type_ = splitline_type,
                                                                                    )
                                                                    ),
                                                  axisline_opts = opts.AxisLineOpts(                              
                                                                    is_show = False
                                                                    ),
                                                  axislabel_opts = opts.LabelOpts(                                
                                                                    is_show = radiusaxis_label_show,
                                                                    color = radiusaxis_label_color,
                                                                    font_family = 'Arial',
                                                                    font_weight = radiusaxis_label_font_weight,
                                                                    interval = 1,
                                                                    font_size = radiusaxis_label_font_size,
                                                                    ),
                                                  splitarea_opts = opts.SplitAreaOpts(
                                                                    is_show = True, 
                                                                    areastyle_opts = opts.AreaStyleOpts(
                                                                        opacity = background_opacity,
                                                                        color = {
                                                                            'type': 'radial',
                                                                            'x': 0.5,
                                                                            'y': 0.5,
                                                                            'r': 0.5,
                                                                            'colorStops': [
                                                                                {'offset': 0, 'color': background_color},  # 0%处的颜色
                                                                                {'offset': 1, 'color': background_color}  # 100%处的颜色
                                                                            ],
                                                                            'global': False  
                                                                        }
                                                                    )
                                                                ),
                                                  ),
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
        )
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        if legend is None:
            legend = columns
        # 为每个数据集添加线条
        for idx, data in enumerate(polar_data):
            polar.add(
                series_name = legend[idx],
                data = data[0],
                type_='bar',
                stack = 'stack0',
                # linestyle_opts=opts.LineStyleOpts(color=colors[idx], width=3),  
                symbol = symbol,
                symbol_size = symbol_size,
                # color=colors[idx]
                label_opts = opts.LabelOpts(                                
                    is_show = False
                    ),
                areastyle_opts = opts.AreaStyleOpts(opacity = area_opacity)
            )
            
        polar.set_colors(colors)

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_polar2_{timestamp}.html")
        polar.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_polar2_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_polar2_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_polar2_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'polar2', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'polar2',  f"{'_'.join(dataframe_names)}_polar2.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def funnel(self, data, item_column, number_column, colors=None,
               top = 10,
               orient = 'vertical',
               sort = 'ascending',
               gap = 1,
               label_font_size = 12,
               label_font_weight = 'normal',
               plot_title = None,
               legend_font_size = 12,
               subfolder='plot',
               filename = None,
               figure_description = 'Funnel plot',
               ):
        """
        Generate a funnel chart visualization from a DataFrame.

        The funnel chart displays values from `number_column` grouped by categories 
        in `item_column`. Items are truncated to the top N (`top`) and can be sorted 
        and oriented as specified. The chart is exported in multiple formats 
        (HTML, SVG, PDF, PNG), and plotting parameters are logged using `data_manager`.

        Args:
            data (pd.DataFrame): Input dataset.
            item_column (str): Column containing item/category names.
            number_column (str): Column containing numeric values to plot.
            colors (list of str, optional): Custom colors for funnel layers. 
                Defaults to a preset palette.
            top (int, optional): Number of top items to display. Default is 10.
            orient (str, optional): Orientation of the funnel. One of {"vertical", "horizontal"}. 
                Default "vertical".
            sort (str, optional): Sorting order of items. One of {"ascending", "descending", "none"}. 
                Default "ascending".
            gap (int, optional): Gap size between funnel layers. Default 1.
            label_font_size (int, optional): Font size of labels inside the funnel. Default 12.
            label_font_weight (str, optional): Font weight of labels. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            plot_title (str, optional): Title of the funnel plot. Default None.
            legend_font_size (int, optional): Font size of legend text. Default 12.
            subfolder (str, optional): Subfolder under `./plot/` for saving output files. Default "plot".
            filename (str, optional): Base filename prefix. A timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. Default "Funnel plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves funnel chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file using `data_manager`.

        Raises:
            ValueError: If parameter values are invalid.
        """
        data = data.iloc[:top,:]
        max_value = max(data[number_column])
        min_value = 0
        
        if top > data.shape[0]:
            top = data.shape[0]
        
        if orient not in ['vertical', 'horizontal']:
            orient = 'vertical'
            
        if sort not in ['ascending', 'descending', 'none']:
            sort = 'ascending'
        
        if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            label_font_weight = 'normal'
        
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:top]
        
        funnel = (
            Funnel(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add(
                item_column,
                [list(z) for z in zip(data[item_column], data[number_column])],
                orient = orient,
                sort_ = sort,
                gap = gap,
                label_opts=opts.LabelOpts(position="inside",
                                          font_family = 'Arial',
                                          font_weight = label_font_weight,
                                          font_size = label_font_size,
                                          ),
                min_=min_value,  
                max_=max_value, 
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
        )
        
        funnel.set_colors(colors)

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_funnel_{timestamp}.html")
        funnel.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_funnel_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_funnel_{timestamp}.pdf")
        drawing = svg2rlg(svg_file) 
        renderPDF.drawToFile(drawing, pdf_file) 
        
        png_file = os.path.join(output_dir, f"{filename}_funnel_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'funnel', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'funnel',  f"{item_column}_{number_column}_funnel.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def parallel(self, *data_names, columns, colors = None, 
                 axis_width = 2,
                 axis_color = 'black',
                 axis_label_color = 'black',
                 axis_label_font_weight = 'normal',
                 axis_label_font_size = 12,
                 line_width = 3,
                 is_smooth = False,
                 if_unique = True,
                 legend = None,
                 legend_font_size = 12,
                 plot_title = None,
                 subfolder='plot',
                 filename = None,
                 figure_description = 'Parallel plot',
                 ):
        """
        Generate a parallel coordinates plot comparing multiple DataFrames.

        Each column in `columns` is treated as a dimension in the parallel plot. 
        For each input DataFrame, values are computed as either unique counts 
        (`if_unique=True`) or total counts (`if_unique=False`). The datasets are 
        drawn as separate polylines across the parallel axes, allowing comparison 
        across multiple categorical dimensions.

        The plot is exported in multiple formats (HTML, SVG, PDF, PNG), and plotting 
        parameters are logged using `data_manager`.

        Args:
            *data_names: Input DataFrames or variable names (str). Strings are 
                evaluated with `eval()`, or DataFrames can be passed directly.
            columns (list of str): Column names to use as parallel dimensions.
            colors (list of str, optional): Line colors for each dataset. 
                Defaults to a preset palette.
            axis_width (int, optional): Width of axis lines. Default 2.
            axis_color (str, optional): Color of axis lines. Default "black".
            axis_label_color (str, optional): Color of axis labels. Default "black".
            axis_label_font_weight (str, optional): Font weight of axis labels. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            axis_label_font_size (int, optional): Font size of axis labels. Default 12.
            line_width (int, optional): Line width for polylines. Default 3.
            is_smooth (bool, optional): Whether to smooth polylines. Default False.
            if_unique (bool, optional): Whether to count unique values (`True`) 
                or total counts (`False`). Default True.
            legend (list of str or str, optional): Legend labels. If None, inferred 
                from dataset names. Default None.
            legend_font_size (int, optional): Font size of legend text. Default 12.
            plot_title (str, optional): Title of the parallel plot. Default None.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. 
                Default "plot".
            filename (str, optional): Base filename prefix. Timestamp appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Parallel plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves parallel plot as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file using `data_manager`.

        Raises:
            ValueError: If parameter values are invalid.
        """
        dataframes = [
                eval(name) if isinstance(name, str) else name
                for name in data_names
            ]
        dataframe_names = [name.split('.')[-1] for name in data_names]
        
        parallel_data = []
        schema = []
        
        if is_smooth not in [True, False]:
            is_smooth = False
            
        if axis_label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            axis_label_font_weight = 'normal'  
            
        #
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(dataframes)]

        max_values = {}
        for column in columns:
            if if_unique:
                max_values[column] = max([df[column].dropna().unique().shape[0] * 1.1 for df in dataframes])
            else:
                max_values[column] = max([df[column].dropna().shape[0] * 1.1 for df in dataframes])

        # 处理 schema 并收集数据
        for df in dataframes:
            if if_unique:
                row_counts = [df[col].dropna().unique().shape[0] for col in columns]
            else:
                row_counts = [df[col].dropna().shape[0] for col in columns]
            parallel_data.append([row_counts])  
        
        #
        parallel_axis = [
            opts.ParallelAxisOpts(
                dim = index, 
                name = col,
                axisline_opts = opts.AxisLineOpts(  
                    is_show = True,
                    linestyle_opts = opts.LineStyleOpts(
                        width = axis_width,      
                        color = axis_color, 
                    )
                ),
                axislabel_opts = opts.LabelOpts(
                    color = axis_label_color,
                    font_family = 'Arial',
                    font_weight = axis_label_font_weight,
                    font_size = axis_label_font_size,
                    )
            ) 
            for index, col in enumerate(columns)
        ]
        
        #
        parallel = (
            Parallel(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add_schema(
                schema = parallel_axis,  
                parallel_opts = opts.ParallelOpts(),
                
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
        )        
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        if legend is None:
            legend = dataframe_names
        
        for idx, data in enumerate(parallel_data):
            parallel.add(
                series_name = legend[idx],
                data = data,
                linestyle_opts=opts.LineStyleOpts(color = colors[idx], 
                                                  width = line_width,
                                                  ),  
                is_smooth = is_smooth,
            )

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_parallel_{timestamp}.html")
        parallel.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_parallel_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_parallel_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_parallel_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'parallel', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'parallel',  f"{'_'.join(dataframe_names)}_parallel.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def pie(self, data, item_column, number_column, colors = None,
             top = 10,
             end = None,
             radius = ["10%", "70%"],
             rosetype = 'radius',
             label_show = True,
             label_color = 'black',
             label_font_weight = 'normal',
             label_font_size = 12,
             label_text_width = 10,
             plot_title = None,
             legend_show = True,
             legend_font_size = 12,
             subfolder='plot',
             filename = None,
             figure_description = 'Pie chart',
             ):
        """
        Generate a pie chart (or rose chart) from a DataFrame.

        This method plots categories from `item_column` against numeric values 
        from `number_column`. Users can limit the number of displayed items 
        (via `top` or `end`), choose pie or rose style, customize label/legend 
        appearance, and adjust chart radius. The chart is exported in multiple 
        formats (HTML, SVG, PDF, PNG), and plotting parameters are logged using 
        `data_manager`.

        Args:
            data (pd.DataFrame): Input dataset.
            item_column (str): Column containing category/item labels.
            number_column (str): Column containing numeric values.
            colors (list of str, optional): Colors for pie slices. 
                Defaults to a preset palette.
            top (int, optional): Display the first N items. Default 10.
            end (int, optional): Display the last N items. If both `top` and `end` 
                are None, all items are used.
            radius (list of str, optional): Inner and outer radius for the pie. 
                Default ["10%", "70%"].
            rosetype (str or None, optional): Rose chart mode. 
                One of {"radius", "area", None}. Default "radius".
            label_show (bool, optional): Whether to show labels on slices. Default True.
            label_color (str, optional): Label text color. Default "black".
            label_font_weight (str, optional): Label font weight. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            label_font_size (int, optional): Label font size. Default 12.
            label_text_width (int, optional): Max characters per line in label text. Default 10.
            plot_title (str, optional): Title of the pie chart. Default None.
            legend_show (bool, optional): Whether to display legend. Default True.
            legend_font_size (int, optional): Font size of legend text. Default 12.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Pie chart".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves pie chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Raises:
            ValueError: If invalid parameter values are provided.
        """
        if top is not None and end is not None:
            data = data
        if top is not None and end is None:
            data = data.iloc[:top,:]
        if top is None and end is not None:
            data = data.iloc[-end:,:]
        if top is None and end is None:
            data = data
            
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:top]

        if rosetype not in ['radius', 'area', None]:
            rosetype = 'radius'
            
        if label_show not in [True, False]:
            label_show = True
            
        if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            label_font_weight = 'normal'  

        pie = (
            Pie(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add("", [list(z) for z in zip(data[item_column], data[number_column])],
                 radius = radius,
                 rosetype = rosetype,
                 label_opts=opts.LabelOpts(is_show = label_show,
                                           color = label_color, 
                                           font_family = 'Arial',
                                           font_weight = label_font_weight,
                                           font_size = label_font_size,
                                           formatter="{b}:\n {c} ({d}%)"
                                           ),
                 )
            .set_colors(colors)
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    is_show = legend_show,
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
        )

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_pie_{timestamp}.html")
        pie.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_pie_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_pie_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file) 
        
        png_file = os.path.join(output_dir, f"{filename}_pie_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'pie', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'pie',f"{item_column}_{number_column}_pie.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def nested_pie(self, data, item_column, number_column, value_counts_column = None, colors=None, 
                   top=15,
                   split=3, 
                   inner_data_first=True, 
                   inner_radius=["0%", "30%"], 
                   inner_rosetype = 'radius',
                   outer_radius=["40%", "70%"],
                   outer_rosetype = 'radius',
                   inner_label_show=True, 
                   outer_label_show=True,
                   label_color="black", 
                   label_font_weight="normal", 
                   label_font_size=12,
                   plot_title = None,
                   legend_font_size = 12,
                   subfolder='plot',
                   filename = None,
                   figure_description = 'Nested pie chart',
                   ):
        """
        Generate a nested (two-level) pie/rose chart from a DataFrame.

        The chart consists of an inner ring and an outer ring, which together 
        represent a split of the input data. Data can be provided directly 
        (`item_column` + `number_column`) or generated by counting values 
        in `value_counts_column`. The top N items are selected, then divided 
        into two groups based on `split` and displayed in concentric layers. 

        Args:
            data (pd.DataFrame): Input dataset.
            item_column (str): Column containing category labels.
            number_column (str): Column containing numeric values.
            value_counts_column (str, optional): If provided, use value counts of 
                this column as data instead of `item_column` and `number_column`. 
                Default None.
            colors (list of str, optional): Custom colors for slices. Defaults to 
                a preset palette sized to `top`.
            top (int, optional): Number of top items to include. Default 15.
            split (int, optional): Index to split items into inner and outer rings. 
                Default 3.
            inner_data_first (bool, optional): Whether the first split goes to the 
                inner ring (`True`) or outer ring (`False`). Default True.
            inner_radius (list of str, optional): Radius range for inner ring. 
                Default ["0%", "30%"].
            inner_rosetype (str, optional): Rose chart type for inner ring. 
                One of {"radius","area"}. Default "radius".
            outer_radius (list of str, optional): Radius range for outer ring. 
                Default ["40%", "70%"].
            outer_rosetype (str, optional): Rose chart type for outer ring. 
                One of {"radius","area"}. Default "radius".
            inner_label_show (bool, optional): Whether to show labels on inner ring. Default True.
            outer_label_show (bool, optional): Whether to show labels on outer ring. Default True.
            label_color (str, optional): Label text color. Default "black".
            label_font_weight (str, optional): Label font weight. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            label_font_size (int, optional): Label font size. Default 12.
            plot_title (str, optional): Title of the chart. Default None.
            legend_font_size (int, optional): Font size for legend text. Default 12.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. 
                Default "plot".
            filename (str, optional): Base filename prefix. Timestamp appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Nested pie chart".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves nested pie chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Raises:
            ValueError: If invalid parameter values are provided.
        """
        if value_counts_column == None:
            data = data.iloc[:top,:]
        else:
            data = data[value_counts_column].value_counts().reset_index()
            data['index'] = data['index'].astype(str)
        
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:top]
            
        if split > data.shape[0]:
            split = data.shape[0] // 2

        if inner_rosetype not in ['radius', 'area']:
            inner_rosetype = 'radius'
            
        if outer_rosetype not in ['radius', 'area']:
            outer_rosetype = 'radius'
            
        if inner_label_show not in [True, False]:
            inner_label_show = True
            
        if outer_label_show not in [True, False]:
            outer_label_show = True
            
        if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            label_font_weight = 'normal'  
        
        # 数据分割
        data_pairs = [list(z) for z in zip(data.iloc[:,0], data.iloc[:,1])]
        inner_data = data_pairs[:split]
        outer_data = data_pairs[split:]

        # 用户选择 inner 和 outer 的数据
        if not inner_data_first:
            inner_data, outer_data = outer_data, inner_data

        # 设置颜色
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:top]

        # 确保颜色数量和数据数量一致
        inner_color = colors[:len(inner_data)]
        outer_color = colors[len(inner_data):]

        # 绘制嵌套饼图
        pie = (
            Pie(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            # 添加 inner 数据
            .add(
                series_name="Inner Data",
                data_pair=inner_data,
                radius=inner_radius,
                rosetype = inner_rosetype,
                label_opts=opts.LabelOpts(
                    is_show=inner_label_show,
                    color=label_color,
                    font_weight=label_font_weight,
                    font_size=label_font_size,
                    formatter="{b}: {c}",  # 这里使用pyecharts的内置百分比显示
                ),
            )
            # 添加 outer 数据
            .add(
                series_name="Outer Data",
                data_pair=outer_data,
                radius=outer_radius,
                rosetype = outer_rosetype,
                label_opts=opts.LabelOpts(
                    is_show=outer_label_show,
                    color=label_color,
                    font_weight=label_font_weight,
                    font_size=label_font_size,
                    formatter="{b}: {c}",  # 这里也同样使用pyecharts的内置百分比显示
                ),
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
            
        )
        
        pie.set_colors(colors)
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        file_name = os.path.join(output_dir, f"{filename}_nested_pie_{timestamp}.html")
        pie.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_nested_pie_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_nested_pie_{timestamp}.pdf")
        drawing = svg2rlg(svg_file) 
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_nested_pie_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'nested_pie', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'nested_pie',f"{item_column}_{number_column}_nested_pie.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def sankey_test(self, data, node1_column, node2_column,
               top = 3,
               subfolder='plot',
               plot_title = None,
               legend_font_size = 12,
               filename = None,
               figure_description = 'Sankey plot',
               ):
        """
        Generate a Sankey diagram for term–gene relationships.

        The function expects an input DataFrame with at least two columns:
        - "Term": category/annotation names (e.g., enriched pathways). 
            Names may include counts in parentheses (e.g., "Pathway (10)").
        - "Genes": semicolon-separated list of genes associated with each term.

        The top N terms are selected (`top`), then each gene is linked to its 
        corresponding term to construct a Sankey diagram. The chart is exported 
        in multiple formats (HTML, SVG, PDF, PNG), and plotting parameters are 
        logged using `data_manager`.

        Args:
            data (pd.DataFrame): Input dataset containing at least "Term" and "Genes".
            node1_column (str): Name of the first node column (kept for logging).
            node2_column (str): Name of the second node column (kept for logging).
            top (int, optional): Number of top terms to visualize. Default 3.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. 
                Default "plot".
            plot_title (str, optional): Title of the Sankey diagram. Default None.
            legend_font_size (int, optional): Font size of legend text. Default 12.
            filename (str, optional): Base filename prefix. A timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Sankey plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves Sankey diagram as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Input `node1_column` and `node2_column` are not used in plotting logic, 
            but are included in output logs for traceability.
            - Terms are truncated at `(` if counts are included in their names.
            - Each gene–term pair is drawn as a link with value=1.
        """
        data = data[:top]
        
        nodes = []
        for term in data['Term']:
            term_name = term.split(' (')[0]
            nodes.append({"name": term_name})
        for genes in data['Genes']:
            gene_list = genes.split(';')
            for gene in gene_list:
                nodes.append({"name": gene})
        nodes = [dict(t) for t in {tuple(d.items()) for d in nodes}]
        
        links = []
        for index, row in data.iterrows():
            term = row['Term'].split(' (')[0]  
            genes = row['Genes'].split(';')  
            for gene in genes:
                link = {"source": gene, "target": term, "value": 1}
                links.append(link)
        
        c = (
            Sankey(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add(
                "sankey",
                nodes,
                links,
                linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source"),
                label_opts=opts.LabelOpts(position="right"),
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                )
        )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_sankey_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_sankey_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        # 将 SVG 转换为 PDF
        pdf_file = os.path.join(output_dir, f"{filename}_sankey_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  # 读取 SVG 文件
        renderPDF.drawToFile(drawing, pdf_file)  # 保存为 PDF 文件
        
        png_file = os.path.join(output_dir, f"{filename}_sankey_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'sankey', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'sankey',f"{node1_column}_{node2_column}_sankey.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def sunburst_mulit(self, data, root_column, columns, colors = None,
                       radius=["30%", "90%"],
                       label_show=True,
                       label_color='black',
                       label_font_weight='normal',
                       label_font_size=12,
                       hide_root_label=True, 
                       plot_title = None,
                       subfolder='plot',
                       filename = None,
                       figure_description = 'Sunburst plot',
                       ):
        """
        Generate a multi-level sunburst plot from hierarchical data.

        This method constructs a sunburst chart where:
        - The `root_column` defines the root grouping (outermost category).
        - The subsequent `columns` define hierarchical values that expand 
            inward to outward.
        - Each row in the DataFrame becomes a hierarchical path through 
            the specified columns.

        The chart is exported in multiple formats (HTML, SVG, PDF, PNG), 
        and plotting parameters are logged using `data_manager`.

        Args:
            data (pd.DataFrame): Input dataset.
            root_column (str): Column used as the root node (first grouping).
            columns (list of str): Columns representing hierarchical levels 
                to be expanded in the sunburst plot.
            colors (list of str, optional): Colors for root categories. 
                Defaults to a preset palette sized by unique values of `root_column`.
            radius (list of str, optional): Inner and outer radius of the sunburst. 
                Default ["30%", "90%"].
            label_show (bool, optional): Whether to display labels on nodes. Default True.
            label_color (str, optional): Label text color. Default "black".
            label_font_weight (str, optional): Label font weight. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            label_font_size (int, optional): Label font size. Default 12.
            hide_root_label (bool, optional): Whether to hide root node labels. Default True.
            plot_title (str, optional): Title of the sunburst plot. Default None.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Sunburst plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves sunburst plot as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Each path from `root_column` through `columns` is converted into 
            hierarchical `SunburstItem` nodes.
            - If `hide_root_label=True`, root node names are hidden but their 
            children remain visible.
            - Values in non-root columns are expected to be numeric and 
            are converted into int before plotting.
        """
        data = data[columns]
        
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(data[root_column].unique())]

        def create_sunburst_structure(df):
            sunburst_data = []

            # 获取 root_column 中的唯一值作为图例
            unique_core_structures = df[root_column].unique()

            for core_structure in unique_core_structures:
                subset = df[df[root_column] == core_structure]
                children = []
                for _, row in subset.iterrows():
                    values = row[1:].dropna().astype(int).tolist()

                    def create_children(values):
                        if not values:
                            return []
                        return [
                            opts.SunburstItem(
                                name=str(values[0]),
                                value=values[0],
                                children=create_children(values[1:])
                            )
                        ]

                    children.extend(create_children(values))

                # 将每个 core_structure 作为一个 series
                sunburst_data.append(
                    opts.SunburstItem(
                        name=core_structure if not hide_root_label else '',  # 根据参数控制是否显示根节点名称
                        children=children
                    )
                )

            return sunburst_data

        # 构建数据结构
        sunburst_data = create_sunburst_structure(data)

        # 创建 Sunburst 图
        sunburst = (
            Sunburst(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add(
                series_name = list(data[root_column].unique()),  # 使用 root_column 作为 series_name
                data_pair=sunburst_data,
                radius=radius,
                label_opts=opts.LabelOpts(
                    is_show=label_show,
                    color=label_color,
                    font_family='Arial',
                    font_weight=label_font_weight,
                    font_size=label_font_size,
                    formatter="{b}",  # 使用 {b} 来显示节点的名称
                    distance=10,
                ),
            )
            .set_colors(colors)
            # 设置全局选项，包括图例
            .set_global_opts(
                tooltip_opts=opts.TooltipOpts(formatter="{b}: {c}"),
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
            )
        )

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_sunburst_mulit_{timestamp}.html")
        sunburst.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_sunburst_mulit_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_sunburst_mulit_{timestamp}.pdf")
        drawing = svg2rlg(svg_file) 
        renderPDF.drawToFile(drawing, pdf_file) 
        
        png_file = os.path.join(output_dir, f"{filename}_sunburst_mulit_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'sunburst_mulit', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'sunburst_mulit',f"{root_column}_sunburst_mulit.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def sunburst(self, data, root_column, child_column, child_column_value, colors = None,
                       radius=["30%", "90%"],
                       label_show=True,
                       label_color='black',
                       label_font_weight='normal',
                       label_font_size=12,
                       legend_show=True,
                       hide_root_label=True,
                       plot_title = None,
                       subfolder='plot',
                       filename = None,
                       figure_description = 'Sunburst plot',
                       ):
        """
        Generate a two-level sunburst plot from a DataFrame.

        The chart shows hierarchical relationships between `root_column` 
        (outer categories) and `child_column` (subcategories) with their 
        corresponding numeric values (`child_column_value`). Each root node 
        is colored from `colors`, while children receive distinct colors from 
        a generated palette.

        Args:
            data (pd.DataFrame): Input dataset.
            root_column (str): Column defining root-level categories.
            child_column (str): Column defining child-level categories under each root.
            child_column_value (str): Column containing numeric values for child nodes.
            colors (list of str, optional): Colors for root nodes. Defaults to a preset palette.
            radius (list of str, optional): Inner and outer radius of the sunburst. 
                Default ["30%", "90%"].
            label_show (bool, optional): Whether to display labels. Default True.
            label_color (str, optional): Label text color. Default "black".
            label_font_weight (str, optional): Label font weight. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            label_font_size (int, optional): Label font size. Default 12.
            legend_show (bool, optional): Whether to display legend. Default True.
            hide_root_label (bool, optional): Whether to hide root node labels. Default True.
            plot_title (str, optional): Title of the sunburst plot. Default None.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Sunburst plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves sunburst plot as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Root nodes are assigned colors from `colors`.
            - Child nodes (`child_column`) are assigned distinct colors from a palette 
            generated by `self.select_palette`.
            - Tooltip displays both name and value (`{b}: {c}`).
        """
        palette = self.select_palette(data, child_column)

        # 生成Core_structure的颜色
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(data[root_column].unique())]
 
        # 使用 palettable 为 Glycan_type 分配颜色
        unique_glycan_types = data[child_column].unique()
        glycan_colors_palette = palette.colors  
        glycan_colors = {
            glycan_type: f"rgb({c[0]}, {c[1]}, {c[2]})" 
            for glycan_type, c in zip(unique_glycan_types, glycan_colors_palette[:len(unique_glycan_types)])
        }
        
        def create_sunburst_data(df):
            sunburst_data = []
            
            # 获取每一个Core_structure的唯一值
            unique_core_structures = df[root_column].unique()
            
            # 遍历每一个Core_structure
            for idx, core_structure in enumerate(unique_core_structures):
                # 获取该Core_structure对应的所有行数据
                subset = df[df[root_column] == core_structure]
                
                # 创建children
                children = []
                for _, row in subset.iterrows():
                    glycan_type = row[child_column]
                    count_value = row[child_column_value]
                    glycan_color = glycan_colors.get(glycan_type, "#000000")  # 如果没有找到颜色，默认为黑色
                    
                    # 为每个Glycan_type创建字典
                    children.append({
                        "name": glycan_type,
                        "value": count_value,
                        "itemStyle": {"color": glycan_color}
                    })
                
                # 为每个Core_structure分配 colors 中的颜色
                sunburst_data.append({
                    "name": core_structure if not hide_root_label else '',
                    "itemStyle": {"color": colors[idx]},  # 分配 colors 中的颜色
                    "children": children
                })
            
            return sunburst_data

        # 构建数据结构
        sunburst_data = create_sunburst_data(data)

        # 创建 Sunburst 图
        sunburst = (
            Sunburst(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add(
                series_name = list(data[root_column].unique()),  # 使用 root_column 作为 series_name
                data_pair=sunburst_data,
                radius=radius,
                label_opts=opts.LabelOpts(
                    is_show=label_show,
                    color=label_color,
                    font_family='Arial',
                    font_weight=label_font_weight,
                    font_size=label_font_size,
                    formatter="{b}",  # 使用 {b} 来显示节点的名称
                    distance=10,
                ),
            )
            .set_colors(colors)  # 设置全局颜色
            # 设置全局选项，包括图例
            .set_global_opts(
                legend_opts=opts.LegendOpts(
                    is_show=legend_show,  
                    orient="vertical",  
                    pos_left="left",  
                    pos_top="middle",  
                    legend_icon="circle",  
                    textstyle_opts=opts.TextStyleOpts(color="black"),  
                ),
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                tooltip_opts=opts.TooltipOpts(formatter="{b}: {c}")  # 设置工具提示
            )
        )

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_sunburst_{timestamp}.html")
        sunburst.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_sunburst_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_sunburst_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_sunburst_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'sunburst', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'sunburst',f"{root_column}_sunburst.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def bar(self, data, x_column=None, y_column=None, y_column_value=None, colors = None,
            top = 10,
            end = None,
                  if_stack = True,
                  transform_ratio = False,
                  label_show = False,
                  label_color = 'black',
                  label_font_weight = 'normal',
                  label_font_size = 12,
                  bar_width = "10%",
                  category_gap = "40%",
                  gap='80%',
                  xaxis_label_show = True,
                  xaxis_label_color = 'black',
                  xaxis_label_font_weight = 'normal',
                  xaxis_label_font_size = 14,
                  xaxis_label_rotate = -15,
                  xaxis_label_margin = 15,
                  xaxis_label_text_split = 10,
                  xaxis_line_show = True,
                  xaxis_line_width = 2,
                  xaxis_line_color = 'black',
                  xaxis_tick_show = True,
                  xaxis_tick_length = 8,
                  xaxis_tick_width = 2,
                  xaxis_tick_color = 'black',
                  xaxis_splitline_show = True,
                  xaxis_splitline_width = 0.5,
                  xaxis_splitline_color = 'grey',
                  top_xaxis_line_show = True,
                  y_max = None,
                  yaxis_label_show = True,
                  yaxis_label_color = 'black',
                  yaxis_label_font_weight = 'normal',
                  yaxis_label_font_size = 14,
                  yaxis_label_margin = 15,
                  yaxis_line_show = True,
                  yaxis_line_width = 2,
                  yaxis_line_color = 'black',
                  yaxis_tick_show = True,
                  yaxis_tick_length = 8,
                  yaxis_tick_width = 2,
                  yaxis_tick_color = 'black',
                  yaxis_splitline_show = True,
                  yaxis_splitline_width = 0.5,
                  yaxis_splitline_color = 'grey',
                  right_yaxis_line_show = True,
                  legend = None,
                  legend_font_size = 12,
                  plot_title = None,
                  xaxis_title = None,
                  xaxis_title_gap = 25,
                  yaxis_title = None,
                  yaxis_title_gap = 40,
                  subfolder='plot',
                  filename = None,
                  figure_description = 'Bar chart',
                  ):
            """
            Generate a bar chart (stacked or simple) from a DataFrame.

            This method supports:
            - Standard bar charts (single category vs value).
            - Stacked bar charts (x_column vs y_column groups, stacked by category).
            - Ratio transformation (`transform_ratio=True`) to normalize values 
                by total (e.g., percentage bars).

            It provides extensive customization for axes, labels, legends, colors, 
            and layout. The chart is exported in multiple formats (HTML, SVG, PDF, PNG), 
            and plotting parameters are logged using `data_manager`.

            Args:
                data (pd.DataFrame): Input dataset.
                x_column (str, optional): Column for x-axis categories. If None, uses `y_column`.
                y_column (str, optional): Column for grouping/legend categories.
                y_column_value (str, optional): Column containing numeric values for bar heights.
                colors (list of str, optional): Colors for bars or groups. Defaults to preset palette.
                top (int, optional): Number of top rows to display. Default 10.
                end (int, optional): Number of last rows to display. Default None.
                if_stack (bool, optional): Whether to use stacked bar chart. Default True.
                transform_ratio (bool, optional): Whether to normalize values as ratios. Default False.
                label_show (bool, optional): Whether to display labels on bars. Default False.
                label_color (str, optional): Label text color. Default "black".
                label_font_weight (str, optional): Label font weight. 
                    One of {"normal","bold","bolder","lighter"}. Default "normal".
                label_font_size (int, optional): Label font size. Default 12.
                bar_width (str or int, optional): Width of bars. Default "10%".
                category_gap (str, optional): Gap between categories. Default "40%".
                gap (str, optional): Gap between bars. Default "80%".
                xaxis_* , yaxis_* (various): Options for axis labels, ticks, lines, 
                    and split lines (see code for full details).
                y_max (float, optional): Maximum value for y-axis. Default None.
                right_yaxis_line_show (bool, optional): Whether to display right-side y-axis line. Default True.
                legend (list of str or str, optional): Custom legend labels. Defaults to `y_column` unique values.
                legend_font_size (int, optional): Legend font size. Default 12.
                plot_title (str, optional): Title of the chart. Default None.
                xaxis_title (str, optional): Title for x-axis. Default None.
                xaxis_title_gap (int, optional): Gap between x-axis title and axis. Default 25.
                yaxis_title (str, optional): Title for y-axis. Default None.
                yaxis_title_gap (int, optional): Gap between y-axis title and axis. Default 40.
                subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
                filename (str, optional): Base filename prefix. Timestamp is appended automatically.
                figure_description (str, optional): Description string logged with output. 
                    Default "Bar chart".

            Returns:
                dict:
                    - "file_path": Path to the saved PNG file.
                    - "legend": Figure description string.

            Side Effects:
                - Saves bar chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
                - Logs plotting parameters and output file via `data_manager`.

            Notes:
                - If both `top` and `end` are provided, no truncation is applied.
                - If `transform_ratio=True`, values are normalized either within 
                each x-category or across y-categories depending on stacking.
            """
            if top is not None and end is not None:
                data = data
            if top is not None and end is None:
                data = data.iloc[:top,:]
            if top is None and end is not None:
                data = data.iloc[-end:,:]
            if top is None and end is None:
                data = data
        
            if y_column_value is not None:
                data = data.sort_values(by=y_column_value, ascending=False) 
        
            if colors is None:
                colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                          "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(data[y_column].unique())]
            
            if if_stack in [True, False]:
                if if_stack:
                    stack = 'stack1'
                    bar_width = None
                else:
                    stack = False
                    bar_width = bar_width
                    
            if transform_ratio in [True, False]:
                if transform_ratio is not False:
                    data[y_column_value] = data[y_column_value]/data[y_column_value].sum()
                else:
                    pass
                    
            if label_show not in [True, False]:
                label_show = False
            
            if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
                label_font_weight = 'normal' 
            
            # 用于处理堆积柱状图的数据生成函数
            def generate_data_for_bar(df, core_column=None, glycan_column=None, value_column=None, transform_ratio=False):
                if core_column:
                    if transform_ratio == False:
                        # 堆积柱状图
                        unique_core_structures = df[core_column].unique()
                        unique_glycan_types = df[glycan_column].unique()
                        all_lists = []
                        for glycan_type in unique_glycan_types:
                            glycan_list = []
                            for core_structure in unique_core_structures:
                                subset = df[(df[core_column] == core_structure) & (df[glycan_column] == glycan_type)]
                                value = subset[value_column].values[0] if not subset.empty else 0
                                total = df[df[core_column] == core_structure][value_column].sum()
                                glycan_list.append({
                                    "value": value,
                                    "percent": value / total if total > 0 else 0
                                })
                            all_lists.append(glycan_list)
                        return all_lists
                    else:
                        # 堆积柱状图
                        unique_core_structures = df[core_column].unique()
                        unique_glycan_types = df[glycan_column].unique()
                        all_lists = []
                        for glycan_type in unique_glycan_types:
                            glycan_list = []
                            for core_structure in unique_core_structures:
                                subset = df[(df[core_column] == core_structure) & (df[glycan_column] == glycan_type)]
                                value = subset[value_column].values[0] if not subset.empty else 0
                                total = df[df[glycan_column] == glycan_type][value_column].sum()
                                glycan_list.append({
                                    "value": value / total if total > 0 else 0,
                                    "percent": value / total if total > 0 else 0
                                })
                            all_lists.append(glycan_list)
                        return all_lists
                else:
                    if transform_ratio == False:
                        # 普通柱状图
                        glycan_list = [{"value": value, "percent": 0} for value in df[value_column]]
                        return [glycan_list]
                    else:
                        # ratio
                        df[value_column] = df[value_column]/df[value_column].sum()
                        glycan_list = [{"value": value, "percent": 0} for value in df[value_column]]
                        return [glycan_list]
            
            all_lists = generate_data_for_bar(data, x_column, y_column, y_column_value, transform_ratio=transform_ratio)
    
            c = Bar(init_opts=opts.InitOpts(
                    renderer=RenderType.SVG,bg_color='#fff', 
                ))
            # 检查是否为普通柱状图还是堆积柱状图
            if x_column:
                c.add_xaxis(list(data[x_column].unique()))
            else:
                c.add_xaxis(list(data[y_column]))
            # 动态添加 y 轴数据
            glycan_types = data[y_column].unique()
            
            if legend is not None:
                if isinstance(legend, str):
                    legend = [legend]
            
            if legend is None:
                legend = glycan_types
            
            for idx, glycan_list in enumerate(all_lists):
                c.add_yaxis(str(legend[idx]), 
                            glycan_list, 
                            stack = stack, 
                            bar_width = bar_width,
                            category_gap = category_gap,
                            gap = gap,
                            )
            c.set_series_opts(
                label_opts=opts.LabelOpts(
                    is_show=label_show,
                    color=label_color,
                    font_family='Arial',
                    font_weight=label_font_weight,
                    font_size=label_font_size,
                    position="right",
                    formatter=JsCode(
                        "function(x){return Number(x.data.percent * 100).toFixed(2) + '%';}"
                    ),
                )
            )
            
            if x_column is not None:
                xaxis_list = list(data[x_column].unique())
            else:
                xaxis_list = list(data[y_column])
                
            c.extend_axis(xaxis_list, 
                xaxis=opts.AxisOpts(
                    type_="category",
                    position='top',
                    axisline_opts=opts.AxisLineOpts(
                        is_show = top_xaxis_line_show,
                        is_on_zero = False, 
                        linestyle_opts = opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                    ),
                    axislabel_opts=opts.LabelOpts(
                         is_show=False),
                    axistick_opts = opts.AxisTickOpts(
                        is_show = False,
                        ),
                    splitline_opts = opts.SplitLineOpts(is_show=False),
                )
            )
            
            c.extend_axis(yaxis=opts.AxisOpts(position="right",
                                              axisline_opts=opts.AxisLineOpts(
                                                    is_show = right_yaxis_line_show,
                                                    linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                                                ),
                                              axislabel_opts = opts.LabelOpts(
                                                    is_show=False,
                                                    ),
                                              splitline_opts = opts.SplitLineOpts(is_show=False),
                                              )
                          )
            
            c.set_colors(colors)
            
            if xaxis_label_text_split != 0:
                if x_column is not None:
                    wrapped_labels = ['\n'.join([str(label)[i:i+xaxis_label_text_split] for i in range(0, len(str(label)), xaxis_label_text_split)]) 
                              for label in list(data[x_column].unique())]
                    c.add_xaxis(wrapped_labels)
                else:
                    wrapped_labels = ['\n'.join([str(label)[i:i+xaxis_label_text_split] for i in range(0, len(str(label)), xaxis_label_text_split)]) 
                              for label in list(data[y_column].unique())]
                    c.add_xaxis(wrapped_labels)
            
            c.set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center", 
                    pos_top="0%"       
                ),
                legend_opts=opts.LegendOpts(
                    pos_right="center",  # 图例靠右
                    pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
                ),
                xaxis_opts=opts.AxisOpts(
                    name = xaxis_title,
                    name_location='center',
                    name_gap = xaxis_title_gap,
                    name_textstyle_opts = opts.TextStyleOpts(color = xaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = xaxis_label_font_weight,
                                                    font_size = xaxis_label_font_size),
                    position='bottom',
                    axislabel_opts = opts.LabelOpts(
                        is_show = xaxis_label_show,
                        font_size = xaxis_label_font_size,  
                        color = xaxis_label_color,  
                        font_family = 'Arial',
                        font_weight = xaxis_label_font_weight,
                        rotate = xaxis_label_rotate,
                        margin = xaxis_label_margin,
                    ),
                    axisline_opts=opts.AxisLineOpts(
                        is_show = xaxis_line_show,
                        linestyle_opts=opts.LineStyleOpts(
                                                          width = xaxis_line_width, 
                                                          color = xaxis_line_color) 
                    ),
                    axistick_opts = opts.AxisTickOpts(
                        is_show = xaxis_tick_show,
                        length = xaxis_tick_length,
                        linestyle_opts=opts.LineStyleOpts(
                                                          width = xaxis_tick_width , 
                                                          color = xaxis_tick_color)
                        ),
                    splitline_opts = opts.SplitLineOpts(
                        is_show = xaxis_splitline_show,
                        linestyle_opts=opts.LineStyleOpts(
                                                          width = xaxis_splitline_width, 
                                                          color = xaxis_splitline_color)
                        )
                ),
                yaxis_opts=opts.AxisOpts(
                    name = yaxis_title,
                    name_location='center',
                    name_gap = yaxis_title_gap,
                    name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                    position='left',
                    max_ = y_max,
                    axislabel_opts=opts.LabelOpts(
                        is_show = yaxis_label_show,
                        font_size = yaxis_label_font_size,  
                        color = yaxis_label_color,  
                        font_family = 'Arial',
                        font_weight = yaxis_label_font_weight,
                        margin = yaxis_label_margin,
                    ),
                    axisline_opts=opts.AxisLineOpts(
                        is_show = yaxis_line_show,
                        linestyle_opts = opts.LineStyleOpts(
                                                          width = yaxis_line_width, 
                                                          color = yaxis_line_color)  
                    ),
                    axistick_opts = opts.AxisTickOpts(
                        is_show = yaxis_tick_show,
                        length =  yaxis_tick_length,
                        linestyle_opts=opts.LineStyleOpts(
                                                          width = yaxis_tick_width, 
                                                          color = yaxis_tick_color)
                        ),
                    splitline_opts = opts.SplitLineOpts(
                        is_show = yaxis_splitline_show,
                        linestyle_opts=opts.LineStyleOpts(
                                                          width = yaxis_splitline_width, 
                                                          color = yaxis_splitline_color)
                        )
                ),
    
            )
            
            grid = Grid(init_opts = opts.InitOpts(width="1600px",height="800px",
                                                  renderer=RenderType.SVG,bg_color='#fff', ))
            grid.add(c, grid_opts = opts.GridOpts(pos_left='30%',pos_bottom='30%'))
    
            output_dir = os.path.join('./plot', subfolder)
            os.makedirs(output_dir, exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            file_name = os.path.join(output_dir, f"{filename}_bar_{timestamp}.html")
            grid.render(file_name)
            svg_file = os.path.join(output_dir, f"{filename}_bar_{timestamp}.svg")
            make_snapshot(snapshot, file_name, svg_file)
            
            pdf_file = os.path.join(output_dir, f"{filename}_bar_{timestamp}.pdf")
            drawing = svg2rlg(svg_file)  
            renderPDF.drawToFile(drawing, pdf_file)  
            
            png_file = os.path.join(output_dir, f"{filename}_bar_{timestamp}.png")
            pdf_document = fitz.open(pdf_file)
            page = pdf_document.load_page(0)
            pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
            pix.save(png_file)
            
            png_file = png_file
            image = Image.open(png_file)
            trimmed_image = self.trim_white_border(image)
            trimmed_image.save(png_file)
    
            # 自动获取所有参数
            params = locals()
            # 去掉不需要记录的局部变量 'self'
            params.pop('self')
            # 使用data_manager记录这些参数
            self.data_manager.log_params('StrucGAP_DataVisualization', 'bar', params)
            self.data_manager.log_output('StrucGAP_DataVisualization', 'bar',f"{x_column}_{y_column}_bar.pdf")
    
            return {'file_path': png_file, 'legend': figure_description}
    
    def butterfly_plot(self,data1,data2,item_column=None,count_column=None,colors=None,
                   top = 10,
                   order_by = 'descending', 
                   bar_width = 0.8,
                   label_font_size = 10,
                   xaxis_title = None,
                   xaxis_title_font_size = 20,
                   plot_title = None,
                   plot_title_font_size = 20,
                   legend = None,
                   legend_fontsize = 12,
                   legend_loc = 'best',
                   subfolder='plot',
                   filename = None,
                   figure_description = 'Butterfly plot',
                   ):
        """
        Generate a butterfly plot (mirrored horizontal bar chart) to compare two datasets.

        This chart is typically used to compare the distribution of two groups 
        (e.g., male vs female, case vs control) across the same categories.

        Args:
            data1 (pd.DataFrame): First dataset containing category and value columns.
            data2 (pd.DataFrame): Second dataset with the same structure as `data1`.
            item_column (str): Column containing category labels.
            count_column (str): Column containing numeric values for bar lengths.
            colors (list of str, optional): Two colors for left and right bars. 
                Must be length 2. Default None.
            top (int, optional): Number of top categories to display. Default 10.
            order_by (str, optional): Sorting order. One of {"descending","ascending"}. 
                Default "descending".
            bar_width (float, optional): Thickness of horizontal bars. Default 0.8.
            label_font_size (int, optional): Font size of category labels. Default 10.
            xaxis_title (str, optional): Title for the x-axis. Default None.
            xaxis_title_font_size (int, optional): Font size of x-axis title. Default 20.
            plot_title (str, optional): Title of the plot. Default None.
            plot_title_font_size (int, optional): Font size of plot title. Default 20.
            legend (list of str or str, optional): Labels for the two datasets in the legend. 
                If string, converted to list. Default None.
            legend_fontsize (int, optional): Font size of legend text. Default 12.
            legend_loc (str, optional): Location of the legend (matplotlib style). Default "best".
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Butterfly plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves butterfly plot as PDF and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Left bars represent `data1` (plotted as negative values).
            - Right bars represent `data2`.
            - Category labels are drawn outside of bars for readability.
            - X-axis tick labels are shown as absolute values.
        """
        if order_by not in ['descending', 'ascending']:
            print('The data was ranked by descending!')
        if top is not None:
            if order_by is 'descending':
                data1 = data1[[item_column, count_column]].sort_values(by=count_column,ascending=False).iloc[:top,:]
                data2 = data2[[item_column, count_column]].sort_values(by=count_column,ascending=False).iloc[:top,:]
            elif order_by is 'ascending':
                data1 = data1[[item_column, count_column]].sort_values(by=count_column,ascending=True).iloc[:top,:]
                data2 = data2[[item_column, count_column]].sort_values(by=count_column,ascending=True).iloc[:top,:]
        else:
            if order_by is 'descending':
                data1 = data1[[item_column, count_column]].sort_values(by=count_column,ascending=False)
                data2 = data2[[item_column, count_column]].sort_values(by=count_column,ascending=False)
            elif order_by is 'ascending':
                data1 = data1[[item_column, count_column]].sort_values(by=count_column,ascending=True)
                data2 = data2[[item_column, count_column]].sort_values(by=count_column,ascending=True)
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        # 提取数据
        left_values = data1[count_column].values
        right_values = data2[count_column].values
        left_labels = data1[item_column].values
        right_labels = data2[item_column].values
        y_positions = np.arange(len(left_values))
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(12, 8))  # 加大图形尺寸以适应标签
        
        # 绘制左右两侧的条形
        ax.barh(y_positions, -left_values, height=bar_width, color=colors[0], align='center')
        ax.barh(y_positions, right_values, height=bar_width, color=colors[1], align='center')
        
        # 添加左侧标签
        for i, (value, label) in enumerate(zip(left_values, left_labels)):
            ax.text(-value - 50, i, label, 
                    ha='right', va='center',
                    fontsize=label_font_size)
        
        # 添加右侧标签
        for i, (value, label) in enumerate(zip(right_values, right_labels)):
            ax.text(value + 50, i, label,
                    ha='left', va='center',
                    fontsize=label_font_size)
        
        plt.tick_params(width=2,length=5,color='k')
        ax = plt.gca()
        ax.spines['bottom'].set_linewidth('2')
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_linewidth('0')
        ax.spines['top'].set_color('k')
        ax.spines['left'].set_linewidth('0')
        ax.spines['left'].set_color('k')
        ax.spines['right'].set_linewidth('0')
        ax.spines['right'].set_color('k')
        
        # 设置x轴标签和标题
        ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        plt.xticks(fontproperties = 'Arial', size = xaxis_title_font_size)
        plt.xticks()
        # 调整x轴显示
        max_value = max(max(left_values), max(right_values))
        margin = 300  # 边距
        ax.set_xlim(-(max_value + margin), max_value + margin)
        target_steps = 3  
        raw_step = (max_value + margin) / target_steps
        magnitude = 10 ** np.floor(np.log10(raw_step))
        tick_step = np.ceil(raw_step / magnitude) * magnitude
        max_tick = ((max_value + margin) // tick_step + 1) * tick_step
        xticks = np.arange(-max_tick, max_tick + tick_step, tick_step)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(abs(int(x))) for x in xticks])
        ax.set_yticks([])
        
        plt.legend(legend, loc=legend_loc, fontsize=legend_fontsize)
    
        # 调整布局
        plt.tight_layout()
        # plt.show()
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_butterfly_plot_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_butterfly_plot_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'butterfly_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'butterfly_plot',f"{filename}_butterfly_plot_{timestamp}.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_butterfly_plot_{timestamp}.png"), 
                'legend': figure_description}
    
    def multi_bar(self, data, x_column, y_column=None, y_column_value=None, colors=None, 
                  label_show=False,
                  label_color='black',
                  label_font_weight='normal',
                  label_font_size=12,
                  bar_width="10%",
                  category_gap="40%",
                  xaxis_label_show=False,
                  xaxis_label_color='black',
                  xaxis_label_font_weight='normal',
                  xaxis_label_font_size=5,
                  xaxis_label_rotate=-15,
                  xaxis_label_margin=15,
                  xaxis_label_text_split=10,
                  xaxis_line_show=True,
                  xaxis_line_width=2,
                  xaxis_line_color='black',
                  xaxis_tick_show=True,
                  xaxis_tick_length=8,
                  xaxis_tick_width=2,
                  xaxis_tick_color='black',
                  xaxis_splitline_show=True,
                  xaxis_splitline_width=0.5,
                  xaxis_splitline_color='grey',
                  yaxis_label_show=True,
                  yaxis_label_color='black',
                  yaxis_label_font_weight='normal',
                  yaxis_label_font_size=10,
                  yaxis_label_margin=15,
                  yaxis_line_show=True,
                  yaxis_line_width=2,
                  yaxis_line_color='black',
                  yaxis_tick_show=True,
                  yaxis_tick_length=8,
                  yaxis_tick_width=2,
                  yaxis_tick_color='black',
                  yaxis_splitline_show=True,
                  yaxis_splitline_width=0.5,
                  yaxis_splitline_color='grey',
                  right_yaxis_line_show=True,
                  plot_title = None,
                  subfolder='plot',
                  filename = None,
                  figure_description = 'Multi bar chart',
                  ):
        """
        Generate a multi-panel bar chart where each row of the DataFrame is drawn as a separate bar chart.

        This method is useful for visualizing multiple entities (rows) across the same set 
        of categories (columns). Each row becomes its own horizontal subplot stacked vertically. 
        The x-axis contains the non-`x_column` columns of the DataFrame, while the y-values 
        are taken from each row. Each subplot has its row label shown on the y-axis.

        Args:
            data (pd.DataFrame): Input dataset.
            x_column (str): Column to be treated as row identifier (e.g., sample, branch).
            y_column (str, optional): If provided, filters to a subset of data by column. Default None.
            y_column_value (str, optional): Not directly used in this method, reserved for compatibility.
            colors (list of str, optional): Colors for bars. Defaults to a preset palette.

            label_show (bool, optional): Whether to display labels on bars. Default False.
            label_color (str, optional): Bar label color. Default "black".
            label_font_weight (str, optional): Bar label font weight. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            label_font_size (int, optional): Bar label font size. Default 12.
            bar_width (str, optional): Width of bars. Default "10%".
            category_gap (str, optional): Gap between bar categories. Default "40%".

            xaxis_* , yaxis_* (various): Options for axis labels, ticks, lines, and split lines 
                (see function definition for full parameter list).
            plot_title (str, optional): Title of the whole figure. Default None.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Multi bar chart".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves the multi-bar chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Each row in the DataFrame becomes a separate subplot stacked vertically.
            - X-axis = all columns except `x_column`. Y-axis values = row values.
            - The y-axis label of each subplot is set to the value in `x_column`.
            - Missing values are replaced with 0 before plotting.
        """
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]
            
        if y_column is not None:
            data = data[y_column]
     
        # 过滤掉所有数值为 NaN 的行
        data_filtered = data.dropna(how='all', subset=data.columns.difference([x_column]))
        data_filtered = data_filtered.replace(np.nan, 0)
        bar_charts = []

        # 遍历每一行，生成单独的柱状图
        for index, row in data_filtered.iterrows():
            # 收集 y 轴的数值，跳过 NaN，并转换为 Python 列表
            y_values = row.drop(labels=[x_column]).dropna().tolist()
            branch_label = row[x_column]

            # X轴标签也转换为 Python 列表
            xaxis_labels = list(data.columns.drop(x_column))

            # 每个柱状图都是单独的，不使用堆积
            c = Bar(init_opts=opts.InitOpts(width="800px", height="200px"))
            c.add_xaxis(xaxis_labels)  # X 轴是除了 `BranchNumber` 的列
            c.add_yaxis(branch_label, y_values, bar_width=bar_width)

            # 设置样式
            c.set_series_opts(
                label_opts=opts.LabelOpts(
                    is_show=False,
                )
            )
            c.set_global_opts(
                xaxis_opts=opts.AxisOpts(
                    axislabel_opts=opts.LabelOpts(
                        is_show=xaxis_label_show,
                        font_size=xaxis_label_font_size,  
                        color=xaxis_label_color,  
                        font_family='Arial',
                        font_weight=xaxis_label_font_weight,
                        rotate=xaxis_label_rotate,
                        margin=xaxis_label_margin,
                    ),
                    axisline_opts=opts.AxisLineOpts(
                        is_show=xaxis_line_show,
                        linestyle_opts=opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                    ),
                    axistick_opts=opts.AxisTickOpts(
                        is_show=xaxis_tick_show,
                        length=xaxis_tick_length,
                        linestyle_opts=opts.LineStyleOpts(width=xaxis_tick_width, color=xaxis_tick_color)
                    ),
                    splitline_opts=opts.SplitLineOpts(
                        is_show=xaxis_splitline_show,
                        linestyle_opts=opts.LineStyleOpts(width=xaxis_splitline_width, color=xaxis_splitline_color)
                    )
                ),
                yaxis_opts=opts.AxisOpts(
                    name = data_filtered['Branches'][index],
                    name_location = 'middle',
                    name_gap = 15,
                    name_rotate = 90,
                    name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                    axislabel_opts=opts.LabelOpts(
                        is_show=yaxis_label_show,
                        font_size=yaxis_label_font_size,  
                        color=yaxis_label_color,  
                        font_family='Arial',
                        font_weight=yaxis_label_font_weight,
                        margin=yaxis_label_margin,
                    ),
                    axisline_opts=opts.AxisLineOpts(
                        is_show=yaxis_line_show,
                        linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                    ),
                    axistick_opts=opts.AxisTickOpts(
                        is_show=yaxis_tick_show,
                        length=yaxis_tick_length,
                        linestyle_opts=opts.LineStyleOpts(width=yaxis_tick_width, color=yaxis_tick_color)
                    ),
                    splitline_opts=opts.SplitLineOpts(
                        is_show=yaxis_splitline_show,
                        linestyle_opts=opts.LineStyleOpts(width=yaxis_splitline_width, color=yaxis_splitline_color)
                    )
                ),
                legend_opts=opts.LegendOpts(
                    is_show=False,  
                ),
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
            )
            c.set_colors(colors[0])
            bar_charts.append(c)
        
        # 使用 Grid 来添加多个柱状图，取消堆积效果
        grid = Grid(init_opts=opts.InitOpts(width="800px", 
                                            height=f'{data_filtered.shape[0]*300}px',renderer=RenderType.SVG,bg_color='#fff' ))  # 设置整体宽度和高度 "4000px"
        
        # 循环加入每个图表，并调整每个图表的位置
        top_offset = 10
        for chart in bar_charts:
            grid.add(chart, grid_opts=opts.GridOpts(pos_top=f"{top_offset}", height="200px"))  # 控制每个图表的高度和位置
            top_offset += 230  # 每个图表向下移动一定距离

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_multi_bar_{timestamp}.html")
        grid.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_multi_bar_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_multi_bar_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file) 
        
        png_file = os.path.join(output_dir, f"{filename}_multi_bar_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        Image.MAX_IMAGE_PIXELS = None
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'multi_bar', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'multi_bar',f"{x_column}_multi_bar.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def bar_up_down(self, data, x_column, up_column, down_column, colors = None, 
              if_stack = True,
              label_show = False,
              label_color = 'black',
              label_font_weight = 'normal',
              label_font_size = 12,
              bar_width = "40%",
              category_gap = "40%",
              xaxis_label_show = True,
              xaxis_label_color = 'black',
              xaxis_label_font_weight = 'normal',
              xaxis_label_font_size = 14,
              xaxis_label_rotate = -15,
              xaxis_label_margin = 15,
              xaxis_label_text_split = 0,
              xaxis_line_show = True,
              xaxis_line_width = 2,
              xaxis_line_color = 'black',
              xaxis_tick_show = True,
              xaxis_tick_length = 8,
              xaxis_tick_width = 2,
              xaxis_tick_color = 'black',
              xaxis_splitline_show = True,
              xaxis_splitline_width = 0.5,
              xaxis_splitline_color = 'grey',
              top_xaxis_line_show = True,
              yaxis_label_show = True,
              yaxis_label_color = 'black',
              yaxis_label_font_weight = 'normal',
              yaxis_label_font_size = 14,
              yaxis_label_margin = 15,
              yaxis_line_show = True,
              yaxis_line_width = 2,
              yaxis_line_color = 'black',
              yaxis_tick_show = True,
              yaxis_tick_length = 8,
              yaxis_tick_width = 2,
              yaxis_tick_color = 'black',
              yaxis_splitline_show = True,
              yaxis_splitline_width = 0.5,
              yaxis_splitline_color = 'grey',
              right_yaxis_line_show = True,
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              xaxis_title = None,
              xaxis_title_gap = 25,
              yaxis_title = None,
              yaxis_title_gap = 40,
              subfolder='plot',
              filename = None,
              figure_description = 'Up down bar chart',
              ):
        """
        Generate an up-down bar chart with symmetric positive and negative bars.

        This chart is designed to compare two related quantities (e.g., 
        up-regulated vs down-regulated, gains vs losses). The "up" values 
        are drawn as positive bars, while the "down" values are inverted 
        (multiplied by -1) to appear on the opposite side.

        Args:
            data (pd.DataFrame): Input dataset containing x-axis categories and 
                two numeric columns (`up_column`, `down_column`).
            x_column (str): Column containing category labels for the x-axis.
            up_column (str): Column containing "up" (positive) values.
            down_column (str): Column containing "down" (negative) values.
            colors (list of str, optional): Colors for up and down bars. 
                Must be length 2. Defaults to preset palette.
            if_stack (bool, optional): Whether to stack bars. Default True.
            label_show (bool, optional): Whether to show labels on bars. Default False.
            label_color (str, optional): Label text color. Default "black".
            label_font_weight (str, optional): Label font weight. 
                One of {"normal","bold","bolder","lighter"}. Default "normal".
            label_font_size (int, optional): Font size of bar labels. Default 12.
            bar_width (str or float, optional): Width of bars. Default "40%".
            category_gap (str, optional): Gap between bar categories. Default "40%".

            xaxis_* , yaxis_* (various): Options for axis labels, ticks, lines, 
                and split lines (see function definition for full parameter list).

            legend (list of str or str, optional): Custom legend labels for up/down. 
                Defaults to `[up_column, down_column]`.
            legend_font_size (int, optional): Font size of legend text. Default 12.
            plot_title (str, optional): Title of the chart. Default None.
            xaxis_title (str, optional): Title for x-axis. Default None.
            xaxis_title_gap (int, optional): Gap between x-axis title and axis. Default 25.
            yaxis_title (str, optional): Title for y-axis. Default None.
            yaxis_title_gap (int, optional): Gap between y-axis title and axis. Default 40.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description string logged with output. 
                Default "Up down bar chart".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Handles duplicate column names by merging automatically.
            - "Down" values are negated for symmetric visualization.
            - Useful for visualizing changes (e.g., up vs down regulation).
        """
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:2]
            
        if if_stack in [True, False]:
            if if_stack:
                stack = 'stack1'
                # bar_width = None
            else:
                stack = False
                # bar_width = bar_width

        def handle_duplicate_columns(data, x_column, up_column, down_column):
            # 检测重复列并自动重命名
            cols = pd.Series(data.columns)
            duplicated_columns = cols[cols.duplicated()].unique()
            
            # 如果发现重复列，则进行重命名
            if len(duplicated_columns) > 0:
                for dup_col in duplicated_columns:
                    indices = cols[cols == dup_col].index.tolist()
                    cols.loc[indices] = [f"{dup_col}_{i}" if i != 0 else dup_col for i in range(len(indices))]
                data.columns = cols
        
            # 动态合并重复列（比如 x_column 和 x_column_1）
            merged_column_name = f"{x_column}_merged"  # 动态生成合并列的名字
            if f"{x_column}_1" in data.columns:
                data[merged_column_name] = data[x_column].combine_first(data[f"{x_column}_1"])
            else:
                # 如果没有重复的列，则直接使用原始列
                data[merged_column_name] = data[x_column]
        
            # 保留合并后的列以及上下数据列，并重命名
            data_cleaned = data[[merged_column_name, up_column, down_column]]
            data_cleaned.columns = [x_column, up_column, down_column]  # 重命名为原始列名
            
            return data_cleaned
        # 生成 all_lists
        data = handle_duplicate_columns(data, x_column, up_column, down_column)
        
        def generate_all_lists(df, core_column, up_column, down_column):
            # 获取唯一的 Core_structure
            unique_core_structures = df[core_column].tolist()
            
            # 初始化大列表
            all_lists = []
            
            # 创建每个 Glycan_type 对应的 list
            up_list = []
            down_list = []
            
            for core_structure in unique_core_structures:
                # 获取对应的 up 和 down 值
                raw_up = df.loc[df[core_column] == core_structure, up_column].iloc[0]
                raw_down = df.loc[df[core_column] == core_structure, down_column].iloc[0]
                def normalize_number(x):
                    if pd.isnull(x):
                        return 0
                    if isinstance(x, (int, np.integer)):
                        return int(x)
                    if isinstance(x, (float, np.floating)):
                        if x.is_integer():
                            return int(x)
                        return x
                    return x
                
                up_value = normalize_number(raw_up)
                down_value = normalize_number(raw_down)                  
                # 填充 up_list 和 down_list
                up_list.append({"value": up_value})
                down_list.append({"value": -down_value})  # down 的数值取负
            
            # 将两个 list 加入到 all_lists 中
            all_lists.append(up_list)
            all_lists.append(down_list)
            
            return all_lists

        all_lists = generate_all_lists(data, x_column, up_column, down_column)
        
        c = Bar(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
        c.add_xaxis(list(data[x_column]))  
        up_down = [up_column, down_column]
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        if legend is None:
            legend = up_down
        
        # 动态添加 y 轴数据
        for idx, glycan_list in enumerate(all_lists):
            c.add_yaxis(legend[idx], 
                        glycan_list, 
                        stack = stack,
                        bar_width = bar_width,
                        category_gap = category_gap,
                        )
        c.set_series_opts(
            label_opts=opts.LabelOpts(
                is_show=label_show,
                color=label_color,
                font_family='Arial',
                font_weight=label_font_weight,
                font_size=label_font_size,
                position="right",
                formatter=JsCode(
                    "function(x){return Number(x.data.percent * 100).toFixed(2) + '%';}"
                ),
            )
        )
        
        c.extend_axis(list(data[x_column]), 
            xaxis=opts.AxisOpts(
                type_="category",
                position='top',
                axisline_opts=opts.AxisLineOpts(
                    is_show = top_xaxis_line_show,
                    is_on_zero = False, 
                    linestyle_opts = opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                ),
                axislabel_opts=opts.LabelOpts(
                     is_show=False),
                axistick_opts = opts.AxisTickOpts(
                    is_show = False,
                    ),
                splitline_opts = opts.SplitLineOpts(is_show=False),
            )
        )
        
        c.extend_axis(list(data[x_column]), 
            xaxis=opts.AxisOpts(
                type_="category",
                position='bottom',
                axisline_opts=opts.AxisLineOpts(
                    is_show = top_xaxis_line_show,
                    is_on_zero = False, 
                    linestyle_opts = opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                ),
                axislabel_opts=opts.LabelOpts(
                     is_show=False),
                axistick_opts = opts.AxisTickOpts(
                    is_show = False,
                    ),
                splitline_opts = opts.SplitLineOpts(is_show=False),
            )
        )
        
        c.extend_axis(yaxis=opts.AxisOpts(position="right",
                                          axisline_opts=opts.AxisLineOpts(
                                                is_show = right_yaxis_line_show,
                                                linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                                            ),
                                          axislabel_opts = opts.LabelOpts(
                                                is_show=False,
                                                ),
                                          splitline_opts = opts.SplitLineOpts(is_show=False),
                                          )
                      )
        
        c.set_colors(colors)
        
        if xaxis_label_text_split != 0:
            wrapped_labels = ['\n'.join([label[i:i+xaxis_label_text_split] for i in range(0, len(str(label)), xaxis_label_text_split)]) 
                      for label in list(data[x_column].unique())]
            c.add_xaxis(wrapped_labels)
        
        c.set_global_opts(
            title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
            ),
            legend_opts=opts.LegendOpts(
                pos_right="center",  # 图例靠右
                pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
            ),
            xaxis_opts=opts.AxisOpts(
                name = xaxis_title,
                name_location='center',
                name_gap = xaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = xaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = xaxis_label_font_weight,
                                                    font_size = xaxis_label_font_size),
                axislabel_opts = opts.LabelOpts(
                    is_show = xaxis_label_show,
                    font_size = xaxis_label_font_size,  
                    color = xaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = xaxis_label_font_weight,
                    rotate = xaxis_label_rotate,
                    margin = xaxis_label_margin,
                    
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = xaxis_line_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_line_width, 
                                                      color = xaxis_line_color) 
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = xaxis_tick_show,
                    length = xaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_tick_width , 
                                                      color = xaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = xaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_splitline_width, 
                                                      color = xaxis_splitline_color)
                    )
            ),
            yaxis_opts=opts.AxisOpts(
                name = yaxis_title,
                name_location='center',
                name_gap = yaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                axislabel_opts=opts.LabelOpts(
                    is_show = yaxis_label_show,
                    font_size = yaxis_label_font_size,  
                    color = yaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = yaxis_label_font_weight,
                    margin = yaxis_label_margin,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = yaxis_line_show,
                    linestyle_opts = opts.LineStyleOpts(
                                                      width = yaxis_line_width, 
                                                      color = yaxis_line_color)  
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = yaxis_tick_show,
                    length =  yaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_tick_width, 
                                                      color = yaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = yaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_splitline_width, 
                                                      color = yaxis_splitline_color)
                    )
            )
        )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_bar_up_down_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_bar_up_down_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_bar_up_down_{timestamp}.pdf")
        drawing = svg2rlg(svg_file) 
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_bar_up_down_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'bar_up_down', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'bar_up_down',f"{x_column}_bar_up_down.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def bar_up_down_ratio(self, feature, colors,
                      screen_feature = [],
                      axis_line_width = 3,
                      axis_line_color = 'black',
                      plot_title = None,
                      plot_title_font_size = 20,
                      xaxis_font_size = 20,
                      yaxis_font_size = 20,
                      legend_fontsize = 15,
                      subfolder='plot',
                      filename = None,
                      figure_description = 'Volcano plot'):
        """
        Generate an up/down ratio bar chart with dual y-axes.

        This plot combines:
            - Bar charts showing the proportion of up-regulated (positive) and 
            down-regulated (negative) IGPs for each glycan feature across multiple 
            fold-change (FC) thresholds.
            - Line plots showing the up/down ratio at each threshold.

        Args:
            feature (str): Glycan feature to visualize. Must be one of:
                ['core_structure','branches_structure','glycan_type',
                'branches_count','glycan_composition','lacdinac',
                'fucosylated_type','acgc'].
            colors (list of str): List of hex color codes for glycan types. If empty, 
                random colors will be generated.
            screen_feature (list, optional): Subset of glycan features to display. 
                If empty, all features are included.
            axis_line_width (int, optional): Width of axis lines. Default 3.
            axis_line_color (str, optional): Color of axis lines. Default "black".
            plot_title (str, optional): Title of the chart. Default None.
            plot_title_font_size (int, optional): Font size of the title. Default 20.
            xaxis_font_size (int, optional): Font size for x-axis labels. Default 20.
            yaxis_font_size (int, optional): Font size for y-axis labels. Default 20.
            legend_fontsize (int, optional): Font size for legend text. Default 15.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. 
                Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended 
                automatically.
            figure_description (str, optional): Description logged with output. 
                Default "Volcano plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves chart as PDF and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Up bars are plotted with full opacity; down bars are plotted with 
            reduced opacity.
            - A secondary y-axis shows the up/down ratio.
            - For `feature='core_structure'`, shorthand names (Core-I, Core-II, …) 
            are substituted automatically.
        """
        print("You can only choose feature from ['core_structure','branches_structure','glycan_type','branches_count','glycan_composition','lacdinac','fucosylated_type','acgc']")
    
        data_up = getattr(self.data_manager.module_records['StrucGAP_GlycoPeptideQuant']['instance'], f"result_{feature}_up_ratio") 
        data_down = getattr(self.data_manager.module_records['StrucGAP_GlycoPeptideQuant']['instance'], f"result_{feature}_down_ratio")     
        data_down.columns = data_up.columns
        data_ratio = getattr(self.data_manager.module_records['StrucGAP_GlycoPeptideQuant']['instance'], f"result_{feature}_ratio")     
        data_ratio.columns = data_up.columns
        
        if feature == 'core_structure':
            data_up = data_up.replace('A2B2C1D1dD1','Core-I') \
                     .replace('A2B2C1D1dD1dcbB5','Core-II') \
                     .replace('A2B2C1D1dD2dD1','Core-III') \
                     .replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
            data_down = data_down.replace('A2B2C1D1dD1','Core-I') \
                     .replace('A2B2C1D1dD1dcbB5','Core-II') \
                     .replace('A2B2C1D1dD2dD1','Core-III') \
                     .replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
            data_ratio = data_ratio.replace('A2B2C1D1dD1','Core-I') \
                     .replace('A2B2C1D1dD1dcbB5','Core-II') \
                     .replace('A2B2C1D1dD2dD1','Core-III') \
                     .replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
        
        # 阈值列表，顺序与 DataFrame 列一致
        thresholds = [1.2, 1.5, 2, 2.5, 3]
        # x 轴坐标位置，每个阈值一个位置
        x_positions = np.arange(len(thresholds))
        
        if len(screen_feature) == 0:
            glycan_types = list(data_up[data_up.columns[0]])   
            n_types = len(glycan_types)
        else:
            glycan_types = screen_feature   
            n_types = len(glycan_types)
        # 每个聚糖类型的柱子宽度
        total_group_width = 0.8
        bar_width = total_group_width / n_types
        offsets = (np.arange(n_types) - (n_types - 1) / 2) * bar_width
        
        colorslist = {}
        if not colors:
            colors = []
            for _ in range(len(glycan_types)):
                color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
                colors.append(color)
        for i, j in zip(glycan_types, colors):
            colorslist[i] = j
            
        fig, ax1 = plt.subplots(figsize=(10, 6))
        # 副 y 轴用于绘制上调/下调比值
        ax2 = ax1.twinx()
        
        # 逐个绘制每种聚糖类型的数据
        for i, glycan in enumerate(glycan_types):
            # 提取每个 glycan 类型在各阈值下的数据
            up_vals = data_up.loc[data_up[data_up.columns[0]] == glycan].iloc[0, 1:].values.astype(float)
            down_vals = data_down.loc[data_down[data_up.columns[0]] == glycan].iloc[0, 1:].values.astype(float)
            ratio_vals = data_ratio.loc[data_ratio[data_up.columns[0]] == glycan].iloc[0, 1:].values.astype(float)
            
            # 计算每个 glycan 类型在各阈值下的横向位置
            # 对于每个阈值，在相同位置上绘制上下对齐的柱子
            group_center = x_positions + offsets[i]  # 这里用 i 来确定每个 glycan 的偏移位置
            # 绘制上调柱状图（保持正值，显示在 x 轴上方）
            ax1.bar(group_center, up_vals, width=bar_width, color=colorslist[glycan], alpha=1,
                    label=f'{glycan} up')
            # 绘制下调柱状图，将数据取负值显示在 x 轴下方
            ax1.bar(group_center, -down_vals, width=bar_width, color=colorslist[glycan], alpha=0.3,
                    label=f'{glycan} down')
            
            # 在副 y 轴上绘制上调/下调比值折线图（保持正值），位置取组中心
            ax2.plot(group_center, ratio_vals, marker='o', linestyle='-',
                     color=colorslist[glycan], label=f'{glycan} ratio')
        
        ax1.tick_params(width=axis_line_width,length=5,color=axis_line_color) 
        ax2.tick_params(width=axis_line_width,length=5,color=axis_line_color)
        ax = plt.gca()
        ax.spines['bottom'].set_linewidth(axis_line_width)
        ax.spines['bottom'].set_color(axis_line_color)
        ax.spines['top'].set_linewidth(axis_line_width)
        ax.spines['top'].set_color(axis_line_color)
        ax.spines['left'].set_linewidth(axis_line_width)
        ax.spines['left'].set_color(axis_line_color)
        ax.spines['right'].set_linewidth(axis_line_width)
        ax.spines['right'].set_color(axis_line_color)
        
        # 设置 x 轴刻度与标签
        ax1.set_xticks(x_positions)
        ax1.set_xticklabels([str(t) for t in thresholds],fontproperties = 'Arial', size = xaxis_font_size)
        ax1.set_xlabel('Threshold',fontproperties = 'Arial', size = xaxis_font_size)
        ax1.set_ylabel('Proportion of regulated IGPs\n(per glycan type at each FC)',fontproperties = 'Arial', size = yaxis_font_size)
        ax2.set_ylabel('Up / Down Ratio',fontproperties = 'Arial', size = yaxis_font_size)
        ax1.tick_params(axis='both', labelsize = xaxis_font_size)
        ax2.tick_params(axis='both', labelsize = xaxis_font_size)
        # 添加一条水平线标识 x 轴
        ax1.axhline(0, color='black', linewidth=axis_line_width*0.5)
        
        # 合并图例，并确保显示所有的图例项
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        all_handles = handles1 + handles2
        all_labels  = labels1  + labels2
        max_per_col = 14
        n = len(all_labels)
        n_cols = (n + max_per_col - 1) // max_per_col

        fig.legend(
            all_handles,
            all_labels,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            fontsize=legend_fontsize,
            ncol=n_cols,
            frameon=False
        )
        plt.tight_layout()
        # plt.show()
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_bar_up_down_ratio_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_bar_up_down_ratio_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'bar_up_down_ratio', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'bar_up_down_ratio',f"{feature}_bar_up_down_ratio.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_bar_up_down_ratio_{timestamp}.png"), 
                'legend': figure_description}
    
    def bar_multi_columns(self, data, y_column, x_columns, colors = None, 
              if_stack = True,
              label_show = False,
              label_color = 'black',
              label_font_weight = 'normal',
              label_font_size = 12,
              bar_width = "10%",
              category_gap = "40%",
              xaxis_label_show = True,
              xaxis_label_color = 'black',
              xaxis_label_font_weight = 'normal',
              xaxis_label_font_size = 9,
              xaxis_label_rotate = -15,
              xaxis_label_margin = 15,
              xaxis_label_text_split = 0,
              xaxis_line_show = True,
              xaxis_line_width = 2,
              xaxis_line_color = 'black',
              xaxis_tick_show = True,
              xaxis_tick_length = 8,
              xaxis_tick_width = 2,
              xaxis_tick_color = 'black',
              xaxis_splitline_show = True,
              xaxis_splitline_width = 0.5,
              xaxis_splitline_color = 'grey',
              top_xaxis_line_show = True,
              y_max = None,
              yaxis_label_show = True,
              yaxis_label_color = 'black',
              yaxis_label_font_weight = 'normal',
              yaxis_label_font_size = 14,
              yaxis_label_margin = 15,
              yaxis_line_show = True,
              yaxis_line_width = 2,
              yaxis_line_color = 'black',
              yaxis_tick_show = True,
              yaxis_tick_length = 8,
              yaxis_tick_width = 2,
              yaxis_tick_color = 'black',
              yaxis_splitline_show = True,
              yaxis_splitline_width = 0.5,
              yaxis_splitline_color = 'grey',
              right_yaxis_line_show = True,
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              xaxis_title = None,
              xaxis_title_gap = 25,
              yaxis_title = None,
              yaxis_title_gap = 40,
              subfolder='plot',
              filename = None,
              figure_description = 'Multi columns bar chart',
              ):
        """
        Generate a grouped or stacked bar chart from multiple numeric columns.

        Each category defined by `y_column` is represented by a group of bars 
        corresponding to values from `x_columns`. Bars can be drawn side-by-side 
        or stacked depending on `if_stack`.

        Args:
            data (pd.DataFrame): Input DataFrame.
            y_column (str): Column name representing categories (e.g., groups).
            x_columns (list of str): Column names containing numeric values.
            colors (list of str, optional): Colors for each category. Defaults to a preset palette.
            if_stack (bool, optional): Whether to stack bars. Default True.
            label_show (bool, optional): Whether to show bar labels. Default False.
            label_color (str, optional): Label color. Default "black".
            label_font_weight (str, optional): Label font weight. Default "normal".
            label_font_size (int, optional): Label font size. Default 12.
            bar_width (str or int, optional): Width of bars. Default "10%".
            category_gap (str, optional): Gap between categories. Default "40%".
            xaxis_label_show (bool, optional): Show/hide x-axis labels. Default True.
            xaxis_label_color (str, optional): X-axis label color. Default "black".
            xaxis_label_font_weight (str, optional): X-axis label weight. Default "normal".
            xaxis_label_font_size (int, optional): X-axis label font size. Default 9.
            xaxis_label_rotate (int, optional): Rotation of x-axis labels. Default -15.
            xaxis_label_margin (int, optional): Margin of x-axis labels. Default 15.
            xaxis_label_text_split (int, optional): Split long labels into multiple lines. 
                Default 0 (no split).
            xaxis_line_show (bool, optional): Show x-axis line. Default True.
            xaxis_line_width (int, optional): X-axis line width. Default 2.
            xaxis_line_color (str, optional): X-axis line color. Default "black".
            xaxis_tick_show (bool, optional): Show x-axis ticks. Default True.
            xaxis_tick_length (int, optional): Length of ticks. Default 8.
            xaxis_tick_width (int, optional): Tick line width. Default 2.
            xaxis_tick_color (str, optional): Tick color. Default "black".
            xaxis_splitline_show (bool, optional): Show grid splitlines. Default True.
            xaxis_splitline_width (float, optional): Splitline width. Default 0.5.
            xaxis_splitline_color (str, optional): Splitline color. Default "grey".
            top_xaxis_line_show (bool, optional): Show top x-axis line. Default True.
            y_max (int, optional): Maximum value of y-axis. Default None.
            yaxis_label_show (bool, optional): Show y-axis labels. Default True.
            yaxis_label_color (str, optional): Y-axis label color. Default "black".
            yaxis_label_font_weight (str, optional): Y-axis label weight. Default "normal".
            yaxis_label_font_size (int, optional): Y-axis label font size. Default 14.
            yaxis_label_margin (int, optional): Margin for y-axis labels. Default 15.
            yaxis_line_show (bool, optional): Show y-axis line. Default True.
            yaxis_line_width (int, optional): Y-axis line width. Default 2.
            yaxis_line_color (str, optional): Y-axis line color. Default "black".
            yaxis_tick_show (bool, optional): Show y-axis ticks. Default True.
            yaxis_tick_length (int, optional): Length of y-axis ticks. Default 8.
            yaxis_tick_width (int, optional): Width of y-axis ticks. Default 2.
            yaxis_tick_color (str, optional): Y-axis tick color. Default "black".
            yaxis_splitline_show (bool, optional): Show horizontal splitlines. Default True.
            yaxis_splitline_width (float, optional): Splitline width. Default 0.5.
            yaxis_splitline_color (str, optional): Splitline color. Default "grey".
            right_yaxis_line_show (bool, optional): Show right y-axis line. Default True.
            legend (list of str, optional): Custom legend labels. Default None (use y_column values).
            legend_font_size (int, optional): Legend font size. Default 12.
            plot_title (str, optional): Title of the chart. Default None.
            xaxis_title (str, optional): Label for the x-axis. Default None.
            xaxis_title_gap (int, optional): Gap for x-axis title. Default 25.
            yaxis_title (str, optional): Label for the y-axis. Default None.
            yaxis_title_gap (int, optional): Gap for y-axis title. Default 40.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Description logged with output. 
                Default "Multi columns bar chart".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves chart as HTML, SVG, PDF, and PNG in `./plot/<subfolder>/`.
            - Logs plotting parameters and output file via `data_manager`.

        Notes:
            - Each `y_column` value produces one series across `x_columns`.
            - `if_stack=True` stacks bars of different series at each x position.
            - Labels can be auto-split across multiple lines with `xaxis_label_text_split`.
        """
        data = data.replace(np.nan, 0)
        
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(data[y_column].unique())]
        
        if if_stack in [True, False]:
            if if_stack:
                stack = 'stack1'
                bar_width = None
            else:
                stack = False
                bar_width = bar_width
                
        if label_show not in [True, False]:
            label_show = False
        
        if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            label_font_weight = 'normal' 
        
        def generate_all_lists(data, y_column, x_columns):
            # 获取唯一的 Core_structure
            unique_core_structures = data[y_column].tolist()
            
            # 初始化大列表
            all_lists = []

            for core_structure in unique_core_structures:
                # 获取对应的 up 和 down 值
                up_value = list(data[data[y_column] == core_structure].values[0][1:])  
                up_list = []
                for i in up_value:
                    # 填充 up_list 和 down_list
                    up_list.append({"value": i})
                # 将两个 list 加入到 all_lists 中
                all_lists.append(up_list)
            
            return all_lists

        all_lists = generate_all_lists(data, y_column, x_columns)
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        if legend is None:
            legend = list(data[y_column])
        
        c = Bar(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
        c.add_xaxis(list(x_columns))  
        # 动态添加 y 轴数据
        for idx, glycan_list in enumerate(all_lists):
            c.add_yaxis(legend[idx], 
                        glycan_list, 
                        stack = stack,
                        bar_width = bar_width,
                        category_gap = category_gap,
                        )
        c.set_series_opts(
            label_opts=opts.LabelOpts(
                is_show=label_show,
                color=label_color,
                font_family='Arial',
                font_weight=label_font_weight,
                font_size=label_font_size,
                position="right",
                formatter=JsCode(
                    "function(x){return Number(x.data.percent * 100).toFixed(2) + '%';}"
                ),
            )
        )
        
        c.extend_axis(list(x_columns), 
            xaxis=opts.AxisOpts(
                type_="category",
                position='top',
                axisline_opts=opts.AxisLineOpts(
                    is_show = top_xaxis_line_show,
                    is_on_zero = False, 
                    linestyle_opts = opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                ),
                axislabel_opts=opts.LabelOpts(
                      is_show=False),
                axistick_opts = opts.AxisTickOpts(
                    is_show = False,
                    ),
                splitline_opts = opts.SplitLineOpts(is_show=False),
            )
        )
        
        c.extend_axis(yaxis=opts.AxisOpts(position="right",
                                          axisline_opts=opts.AxisLineOpts(
                                                is_show = right_yaxis_line_show,
                                                linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                                            ),
                                          axislabel_opts = opts.LabelOpts(
                                                is_show=False,
                                                ),
                                          splitline_opts = opts.SplitLineOpts(is_show=False),
                                          )
                      )

        c.set_colors(colors)
        
        if xaxis_label_text_split != 0:
            wrapped_labels = ['\n'.join([label[i:i+xaxis_label_text_split] for i in range(0, len(label), xaxis_label_text_split)]) 
                      for label in list(x_columns)]
            c.add_xaxis(wrapped_labels)
        
        c.set_global_opts(
            title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
            ),
            legend_opts=opts.LegendOpts(
                pos_right="center",  # 图例靠右
                pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
            ),
            xaxis_opts=opts.AxisOpts(
                name = xaxis_title,
                name_location='center',
                name_gap = xaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = xaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = xaxis_label_font_weight,
                                                    font_size = xaxis_label_font_size),
                position='bottom',
                axislabel_opts = opts.LabelOpts(
                    is_show = xaxis_label_show,
                    font_size = xaxis_label_font_size,  
                    color = xaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = xaxis_label_font_weight,
                    rotate = xaxis_label_rotate,
                    margin = xaxis_label_margin,
                    interval = 0,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = xaxis_line_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_line_width, 
                                                      color = xaxis_line_color) 
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = xaxis_tick_show,
                    length = xaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_tick_width , 
                                                      color = xaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = xaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_splitline_width, 
                                                      color = xaxis_splitline_color)
                    )
            ),
            yaxis_opts=opts.AxisOpts(
                name = yaxis_title,
                name_location='center',
                name_gap = yaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                position='left',
                max_ = y_max,
                axislabel_opts=opts.LabelOpts(
                    is_show = yaxis_label_show,
                    font_size = yaxis_label_font_size,  
                    color = yaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = yaxis_label_font_weight,
                    margin = yaxis_label_margin,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = yaxis_line_show,
                    linestyle_opts = opts.LineStyleOpts(
                                                      width = yaxis_line_width, 
                                                      color = yaxis_line_color)  
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = yaxis_tick_show,
                    length =  yaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_tick_width, 
                                                      color = yaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = yaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_splitline_width, 
                                                      color = yaxis_splitline_color)
                    )
            ),

        )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
       
        file_name = os.path.join(output_dir, f"{filename}_bar_multi_columns_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_bar_multi_columns_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_bar_multi_columns_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)
        
        png_file = os.path.join(output_dir, f"{filename}_bar_multi_columns_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'bar_multi_columns', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'bar_multi_columns',f"{y_column}_bar_multi_columns.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def frequency_bar(self, data, columns, colors=None, 
                       log2_transformation = True,
                       bins=50, 
                       density=True, 
                       alpha=0.7, 
                       edge_color='white', 
                       ref_line_style='--', 
                       ref_line_width=1.5, 
                       ref_lines_value = 1.5, 
                       xaxis_label_color='black', 
                       xaxis_label_font_weight='normal', 
                       xaxis_label_font_size=9, 
                       yaxis_label_color='black', 
                       yaxis_label_font_weight='normal', 
                       yaxis_label_font_size=14, 
                       axis_line_width=2, 
                       axis_tick_length=5, 
                       axis_tick_width=2, 
                       add_x_grid=False,  
                       add_y_grid=True,  
                       x_grid_line_color='grey', 
                       x_grid_line_style='-', 
                       x_grid_line_width=0.5,
                       y_grid_line_color='grey', 
                       y_grid_line_style='-', 
                       y_grid_line_width=0.5,
                       xaxis_title = None,
                       xaxis_title_font_size = 20,
                       yaxis_title = None,
                       yaxis_title_font_size = 20,
                       plot_title = None,
                       plot_title_font_size = 20,
                       legend = None,
                       legend_fontsize = 12,
                       subfolder='plot',
                       filename = None,
                       figure_description = 'Frequency bar chart',
                       ):
        """
        Plot frequency distribution histograms for selected columns.

        Each selected column in `data` is plotted as a histogram, optionally log2-
        transformed. Reference vertical lines (e.g., ±log2(1.5)) can be added to 
        indicate fold-change thresholds. Multiple columns are overlaid for 
        comparison.

        Args:
            data (pd.DataFrame): Input DataFrame.
            columns (list of str): Column names to plot.
            colors (list of str, optional): List of colors for each column. Defaults 
                to preset palette.
            log2_transformation (bool, optional): Apply log2 transform before plotting. 
                Default True.
            bins (int, optional): Number of histogram bins. Default 50.
            density (bool, optional): Normalize histogram to density. Default True.
            alpha (float, optional): Bar transparency. Default 0.7.
            edge_color (str, optional): Bar edge color. Default "white".
            ref_line_style (str, optional): Line style for reference lines. Default "--".
            ref_line_width (float, optional): Reference line width. Default 1.5.
            ref_lines_value (float, optional): Fold-change value for reference lines 
                (drawn at log2(ref_lines_value) and -log2(ref_lines_value)). Default 1.5.
            xaxis_label_color, yaxis_label_color (str, optional): Axis label colors. 
                Default "black".
            xaxis_label_font_weight, yaxis_label_font_weight (str, optional): Axis label 
                font weights. Default "normal".
            xaxis_label_font_size, yaxis_label_font_size (int, optional): Axis label font 
                sizes. Defaults: 9 (x), 14 (y).
            axis_line_width (int, optional): Axis line thickness. Default 2.
            axis_tick_length (int, optional): Length of axis ticks. Default 5.
            axis_tick_width (int, optional): Width of axis ticks. Default 2.
            add_x_grid (bool, optional): Show vertical grid lines. Default False.
            add_y_grid (bool, optional): Show horizontal grid lines. Default True.
            x_grid_line_color, y_grid_line_color (str, optional): Grid line colors. Default "grey".
            x_grid_line_style, y_grid_line_style (str, optional): Grid line styles. Default "-".
            x_grid_line_width, y_grid_line_width (float, optional): Grid line widths. Default 0.5.
            xaxis_title, yaxis_title, plot_title (str, optional): Axis and plot titles. 
                Default None.
            xaxis_title_font_size, yaxis_title_font_size, plot_title_font_size (int, optional): 
                Title font sizes. Defaults: 20.
            legend (list of str, optional): Legend labels for each column. Default None.
            legend_fontsize (int, optional): Legend font size. Default 12.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs. 
                Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str, optional): Figure description. Default "Frequency bar chart".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves histogram to PDF and PNG in `./plot/<subfolder>/`.
            - Logs parameters and output file via `data_manager`.

        Notes:
            - If `log2_transformation=True`, data is transformed before plotting.
            - Reference lines are symmetric about 0 on the log2 scale.
            - Multiple histograms are overlaid for comparison.
        """
        data = data[columns]
        if log2_transformation == True:
            data = np.log2(data)
            
        xmax = math.ceil(data.abs().max().max())

        # 创建绘图区域
        plt.figure(figsize=(8, 6))

        # 确认颜色列表长度匹配列数
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]

        # 绘制每列的直方图
        for idx, col in enumerate(columns):
            plt.hist(data[col], bins=bins, alpha=alpha, label=col, color=colors[idx], 
                     edgecolor=edge_color, density=density)

        # 添加参考线
        if ref_lines_value is not None:
            plt.axvline(x=np.log2(ref_lines_value), color='black', linestyle=ref_line_style, linewidth=ref_line_width)
            plt.axvline(x=np.log2(1/ref_lines_value), color='black', linestyle=ref_line_style, linewidth=ref_line_width)

        # 设置图例和标签
        plt.xticks(fontproperties = 'Arial', size = xaxis_label_font_size)
        plt.yticks(fontproperties = 'Arial', size = yaxis_label_font_size)
        plt.xlabel('', fontsize=xaxis_label_font_size, color=xaxis_label_color, weight=xaxis_label_font_weight)
        plt.ylabel('', fontsize=yaxis_label_font_size, color=yaxis_label_color, weight=yaxis_label_font_weight)
        plt.xlim([-xmax, xmax])
        
        ax = plt.gca()
        ax.spines['top'].set_linewidth(axis_line_width)
        ax.spines['right'].set_linewidth(axis_line_width)
        ax.spines['left'].set_linewidth(axis_line_width)
        ax.spines['bottom'].set_linewidth(axis_line_width)

        ax.tick_params(axis='both', which='major', length=axis_tick_length, width=axis_tick_width)
        ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
        ax.set_ylabel(yaxis_title,fontproperties = 'Arial', size = yaxis_title_font_size) 
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)

        plt.legend(legend, loc='best', fontsize=legend_fontsize)
        # 添加x轴网格线
        if add_x_grid:
            ax.xaxis.grid(True, color=x_grid_line_color, linestyle=x_grid_line_style, linewidth=x_grid_line_width)

        # 添加y轴网格线
        if add_y_grid:
            ax.yaxis.grid(True, color=y_grid_line_color, linestyle=y_grid_line_style, linewidth=y_grid_line_width)

        # 设置输出路径
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        # 生成时间戳
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 保存图像
        pdf_name = f"{filename}_frequency_bar_{timestamp}.pdf"
        plt.savefig(os.path.join(output_dir, pdf_name), dpi=900, bbox_inches='tight')
        png_file = f"{filename}_frequency_bar_{timestamp}.png"
        plt.savefig(os.path.join(output_dir, png_file), dpi=900, bbox_inches='tight')
        plt.close('all')
        
        image = Image.open(os.path.join(output_dir, png_file))
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(os.path.join(output_dir, png_file))

        # 自动记录参数和输出
        params = locals()
        params.pop('self')
        self.data_manager.log_params('StrucGAP_DataVisualization', 'histogram_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'histogram_plot', pdf_name)

        return {'file_path': os.path.join(output_dir, png_file), 'legend': figure_description}
    
    def boxplot(self, data, item_column, item_name, group1_columns, group2_columns, 
              p_data, p_column, colors = None, 
              label_show = False,
              label_color = 'black',
              label_font_weight = 'normal',
              label_font_size = 12,
              bar_width = "10%",
              category_gap = "40%",
              xaxis_label_show = True,
              xaxis_label_color = 'black',
              xaxis_label_font_weight = 'normal',
              xaxis_label_font_size = 14,
              xaxis_label_rotate = -15,
              xaxis_label_margin = 15,
              xaxis_line_show = True,
              xaxis_line_width = 2,
              xaxis_line_color = 'black',
              xaxis_tick_show = True,
              xaxis_tick_length = 8,
              xaxis_tick_width = 2,
              xaxis_tick_color = 'black',
              xaxis_splitline_show = True,
              xaxis_splitline_width = 2,
              xaxis_splitline_color = 'black',
              yaxis_label_show = True,
              yaxis_label_color = 'black',
              yaxis_label_font_weight = 'normal',
              yaxis_label_font_size = 14,
              yaxis_label_margin = 15,
              yaxis_line_show = True,
              yaxis_line_width = 2,
              yaxis_line_color = 'black',
              yaxis_tick_show = True,
              yaxis_tick_length = 8,
              yaxis_tick_width = 2,
              yaxis_tick_color = 'black',
              yaxis_splitline_show = True,
              yaxis_splitline_width = 0.5,
              yaxis_splitline_color = 'grey',
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              xaxis_title = None,
              xaxis_title_gap = 25,
              yaxis_title = None,
              yaxis_title_gap = 40,
              subfolder='plot',
              filename = None,
              figure_description = 'Box plot',
              ):
        """
        Generate side-by-side boxplots for two groups across multiple items.

        For each `item_name` in `item_column`, this function extracts values from 
        `group1_columns` and `group2_columns` and plots paired boxplots. Statistical 
        significance values (`p_column`) from `p_data` are displayed on the top axis.

        Args:
            data (pd.DataFrame): Input DataFrame containing numeric values for groups.
            item_column (str): Column name identifying item categories (e.g., glycan types).
            item_name (list of str): Items to include in the boxplot.
            group1_columns (list of str): Column names for group 1 values.
            group2_columns (list of str): Column names for group 2 values.
            p_data (pd.DataFrame): DataFrame containing p-values per item.
            p_column (str): Column name in `p_data` that holds p-values.
            colors (list of str, optional): Colors for group1 and group2 boxes. 
                Defaults to preset palette.
            label_show (bool, optional): Show data labels inside boxes. Default False.
            label_color, label_font_weight, label_font_size: Style options for labels.
            bar_width (str, optional): Width of box elements. Default "10%".
            category_gap (str, optional): Gap between categories. Default "40%".
            xaxis_label_show, yaxis_label_show (bool, optional): Show/hide axis labels.
            xaxis_label_color, yaxis_label_color (str, optional): Axis label colors.
            xaxis_label_font_weight, yaxis_label_font_weight (str, optional): Font weight.
            xaxis_label_font_size, yaxis_label_font_size (int, optional): Font sizes.
            ... (many axis customization options)
            legend (list of str, optional): Legend labels for the two groups. 
                Default ["Group1", "Group2"].
            legend_font_size (int, optional): Legend font size. Default 12.
            plot_title (str, optional): Title of the plot.
            xaxis_title, yaxis_title (str, optional): Axis titles.
            subfolder (str, optional): Subfolder under `./plot/` for saving outputs.
                Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended.
            figure_description (str, optional): Figure description. Default "Box plot".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Saves plot in HTML, SVG, PDF, and PNG formats under `./plot/<subfolder>/`.
            - Crops whitespace from PNG output.
            - Logs parameters and output path with `data_manager`.

        Notes:
            - p-values from `p_data[p_column]` are aligned to each `item_name` and 
            displayed along the top axis.
            - Group1 and Group2 are plotted side-by-side with separate colors.
            - Long item labels are wrapped every 20 characters for readability.
        """
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]

        v1 = []
        v2 = []
        for i in item_name:
            v1.append(list(data[data[item_column] == i][group1_columns].values[0]))
            v2.append(list(data[data[item_column] == i][group2_columns].values[0]))
        
        max_value = round(max(list(chain(*v1, *v2)))*1.1, 1)
        min_value = (math.floor(min(list(chain(*v1, *v2))) * 10) / 10)*0.9
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        if legend is None:
            legend = ['Group1', 'Group2']

        p_list = []
        for i in item_name:
            p_list.append(round(p_data[p_data.index == i][p_column].values[0], 5))
 
        c = Boxplot(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
        c.add_xaxis(item_name)
        c.add_yaxis(legend[0], c.prepare_data(v1), yaxis_index = 0,
                    itemstyle_opts=opts.ItemStyleOpts(
                        color=colors[0],
                        border_color = 'black',
                        border_width = 2,
                        )
            )
        c.add_yaxis(legend[0], c.prepare_data(v2), yaxis_index = 1,
                    itemstyle_opts=opts.ItemStyleOpts(
                        color=colors[1],
                        border_color = 'black',
                        border_width = 2,
                        )
            )
        
        c.extend_axis(p_list, 
                      xaxis=opts.AxisOpts(
                            position='top',
                            axisline_opts=opts.AxisLineOpts(
                                is_show=True,
                                linestyle_opts=opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                            ),
                            axislabel_opts = opts.LabelOpts(
                                is_show = xaxis_label_show,
                                font_size = xaxis_label_font_size,  
                                color = xaxis_label_color,  
                                font_family = 'Arial',
                                font_weight = xaxis_label_font_weight,
                                # rotate = xaxis_label_rotate,
                                margin = xaxis_label_margin,
                            ),
                            splitline_opts = opts.SplitLineOpts(is_show=False),
                          ),
                      )
        
        c.extend_axis(yaxis=opts.AxisOpts(position="right",
                                          min_ = min_value,
                                          max_ = max_value,
                                          axisline_opts=opts.AxisLineOpts(
                                                is_show=True,
                                                linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                                            ),
                                          axislabel_opts = opts.LabelOpts(
                                                is_show=False,
                                                ),
                                          splitline_opts = opts.SplitLineOpts(is_show=False),
                                          )
                      )
        
        c.set_colors(colors)
        wrapped_labels = ['\n'.join([label[i:i+20] for i in range(0, len(label), 20)]) 
                  for label in item_name]
        c.add_xaxis(wrapped_labels)

        c.set_global_opts(
            xaxis_opts=opts.AxisOpts(
                name = xaxis_title,
                name_location='center',
                name_gap = xaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = xaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = xaxis_label_font_weight,
                                                    font_size = xaxis_label_font_size),
                position='bottom',
                axislabel_opts = opts.LabelOpts(
                    is_show = xaxis_label_show,
                    font_size = xaxis_label_font_size,  
                    color = xaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = xaxis_label_font_weight,
                    # rotate = xaxis_label_rotate,
                    margin = xaxis_label_margin,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = xaxis_line_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_line_width, 
                                                      color = xaxis_line_color) 
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = xaxis_tick_show,
                    length = xaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_tick_width , 
                                                      color = xaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = xaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_splitline_width, 
                                                      color = xaxis_splitline_color)
                    )
            ),
            yaxis_opts=opts.AxisOpts(
                name = yaxis_title,
                name_location='center',
                name_gap = yaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                position='left',
                min_ = min_value,
                max_ = max_value,
                axislabel_opts=opts.LabelOpts(
                    is_show = yaxis_label_show,
                    font_size = yaxis_label_font_size,  
                    color = yaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = yaxis_label_font_weight,
                    margin = yaxis_label_margin,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = yaxis_line_show,
                    linestyle_opts = opts.LineStyleOpts(
                                                      width = yaxis_line_width, 
                                                      color = yaxis_line_color)  
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = yaxis_tick_show,
                    length =  yaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_tick_width, 
                                                      color = yaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = yaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_splitline_width, 
                                                      color = yaxis_splitline_color)
                    )
            ),
            title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
            ),
            legend_opts=opts.LegendOpts(
                pos_right="center",  # 图例靠右
                pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
            ),
        )

        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 根据 item_name 的长度和 p_column 生成文件名
        item_count = len(item_name)
        file_name = os.path.join(output_dir, f"{filename}_boxplot_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_boxplot_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_boxplot_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)  
        
        png_file = os.path.join(output_dir, f"{filename}_boxplot_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'boxplot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'boxplot',f"{item_column}_item_{p_column}_boxplot.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def heatmap1(self, data, columns, colors = None, 
                cluster = 'both',
                cluster_method = 'weighted',
                log = True,
                minvalue = -1,
                maxvalue = 1,
                xaxis_label_show = True,
                xaxis_label_color = 'black',
                xaxis_label_font_weight = 'normal',
                xaxis_label_font_size = 14,
                yaxis_label_show = False,
                yaxis_label_color = 'black',
                yaxis_label_font_weight = 'normal',
                yaxis_label_font_size = 14,
                label_show = False,
                label_color = 'black',
                label_font_weight ='normal',
                label_font_size = 12,
                plot_title = None,
                subfolder='plot',
                filename = None,
                figure_description = 'Heatmap',
                ):
        """
        Generate a clustered heatmap with optional log-transformation and custom styling.

        Args:
            data (pd.DataFrame): Input DataFrame containing numeric values.
            columns (list of str): Column names to include in the heatmap.
            colors (list of str, optional): Color gradient for the heatmap. 
                Defaults to a red-blue diverging palette.
            cluster (str, optional): Clustering mode: 'row', 'col', 'both', or 'none'. Default 'both'.
            cluster_method (str, optional): Linkage method for hierarchical clustering. 
                One of ['single','complete','average','weighted','centroid','median','ward'].
                Default 'weighted'.
            log (bool, optional): Apply log2 transformation to data. Default True.
            minvalue (float, optional): Minimum value for color scale. Default -1.
            maxvalue (float, optional): Maximum value for color scale. Default 1.
            xaxis_label_show, yaxis_label_show (bool): Show/hide axis labels.
            xaxis_label_color, yaxis_label_color (str): Label font colors.
            xaxis_label_font_weight, yaxis_label_font_weight (str): Label font weights.
            xaxis_label_font_size, yaxis_label_font_size (int): Label font sizes.
            label_show (bool, optional): Whether to display heatmap values inside cells.
            label_color, label_font_weight, label_font_size: Styling for value labels.
            plot_title (str, optional): Title for the heatmap.
            subfolder (str, optional): Subfolder under `./plot/` to save outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended.
            figure_description (str, optional): Figure description string. Default "Heatmap".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Applies hierarchical clustering on rows/columns if requested.
            - Saves heatmap in HTML, SVG, PDF, and PNG formats under `./plot/<subfolder>/`.
            - Crops whitespace from PNG output.
            - Logs parameters and output path with `data_manager`.

        Notes:
            - Values are rounded to two decimals when displayed inside cells.
            - When `cluster='both'`, both rows and columns are reordered by hierarchical clustering.
            - Default color palette is a diverging scheme suitable for log fold changes.
        """
        df = data[columns]
        
        if colors == None:
            colors = ["#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"]
        
        if cluster not in ['row', 'col', 'both', 'none']:
            cluster = 'both'
            
        if cluster_method not in ['single', 'complete', 'average', 'weighted', 
                           'centroid', 'median', 'ward']:
            cluster_method = 'ward'
        
        if log not in [True, False]:
            log = True
        
        if log == True:
            df = df.apply(pd.to_numeric, errors='coerce')
            scaled_data = np.log2(df)
        elif log == False:
            scaled_data = df
        
        if cluster == 'row':
            # 对行进行层次聚类
            row_linkage = sch.linkage(df, method=cluster_method)
            row_dendro = sch.dendrogram(row_linkage, no_plot=True)
            row_order = row_dendro['leaves']
            # 只对行进行重新排列
            clustered_data = scaled_data.iloc[row_order, :]
        elif cluster == 'col':
            # 对列进行层次聚类
            col_linkage = sch.linkage(df.T, method=cluster_method)
            col_dendro = sch.dendrogram(col_linkage, no_plot=True)
            col_order = col_dendro['leaves']
            # 只对列进行重新排列
            clustered_data = scaled_data.T.iloc[col_order, :].T
        elif cluster == 'both':
            # 对行进行层次聚类
            row_linkage = sch.linkage(df, method=cluster_method)
            row_dendro = sch.dendrogram(row_linkage, no_plot=True)
            row_order = row_dendro['leaves']
            # 对列进行层次聚类
            col_linkage = sch.linkage(df.T, method=cluster_method)
            col_dendro = sch.dendrogram(col_linkage, no_plot=True)
            col_order = col_dendro['leaves']
            # 根据行和列聚类结果重新排列数据
            clustered_data = scaled_data.iloc[row_order, :].T.iloc[col_order, :].T
        elif cluster == 'none':
            # 不进行聚类，保留原始顺序
            clustered_data = scaled_data
        
        # minvalue = clustered_data.min().min()
        # maxvalue = clustered_data.max().max()
        
        # 处理为 pyecharts 的输入格式
        heatmap_data = [[j, i, round(clustered_data.values[i, j], 2)] for j in range(clustered_data.shape[1]) for i in range(clustered_data.shape[0])]
        
        # 使用pyecharts绘制热图
        heatmap = HeatMap(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
        heatmap.add_xaxis(list(clustered_data.columns))
        heatmap.add_yaxis("",
                          list(clustered_data.index), 
                          heatmap_data,)
        
        # 设置全局选项和显示标签，只保留两位小数
        heatmap.set_global_opts(
            visualmap_opts=opts.VisualMapOpts(min_ = minvalue,
                                              max_ = maxvalue,
                                              orient = "horizontal", 
                                              pos_left = "center",
                                              range_color = colors,
                                              ),
            xaxis_opts=opts.AxisOpts(
                axislabel_opts = opts.LabelOpts(
                    is_show = xaxis_label_show,
                    font_size = xaxis_label_font_size,  
                    color = xaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = xaxis_label_font_weight,
                    interval = 0,
                ),
                axisline_opts=opts.AxisLineOpts(is_show = False),
                axistick_opts = opts.AxisTickOpts(is_show = True),
            ),
            yaxis_opts=opts.AxisOpts(
                type_="category",
                axislabel_opts=opts.LabelOpts(
                    is_show = yaxis_label_show,
                    font_size = yaxis_label_font_size,  
                    color = yaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = yaxis_label_font_weight,
                ),
                axisline_opts=opts.AxisLineOpts(is_show = False),
                axistick_opts = opts.AxisTickOpts(is_show = False),
            ),
            title_opts=opts.TitleOpts(
                title = plot_title,
                pos_left="center",  # 标题居中
                pos_top="0%"       # 标题位置靠上
            ),
        )
        
        heatmap.set_series_opts(
            label_opts=opts.LabelOpts(
                is_show=label_show,
                color=label_color,
                font_family='Arial',
                font_weight=label_font_weight,
                font_size=label_font_size,
                formatter=JsCode("function(params) { return Number(params.value[2]).toFixed(2); }") 
            )
        )
                
        # heatmap.render("double_clustered_heatmap.html")
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        file_name = os.path.join(output_dir, f"{filename}_heatmap1_{timestamp}.html")
        heatmap.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_heatmap1_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        pdf_file = os.path.join(output_dir, f"{filename}_heatmap1_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  
        renderPDF.drawToFile(drawing, pdf_file)
        
        png_file = os.path.join(output_dir, f"{filename}_heatmap1_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'heatmap1', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'heatmap1',f"{columns}_heatmap1.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def heatmap2(self, data, columns, colors=None, 
                 filter_data = None,
                 filter_columns = None,
                 filter_values = None,
                 cluster='both', 
                 cluster_method='ward', 
                 log=False, 
                 text_annotation = False,
                 text_size = 10,
                 text_color = None,
                 minvalue=-1, 
                 centervalue = 0,
                 maxvalue=1,
                 xaxis_label_show=True, 
                 xaxis_label_color='black', 
                 xaxis_label_font_weight='normal', 
                 xaxis_label_font_size=14, 
                 xaxis_label_rotate = 90,
                 yaxis_label_show=False, 
                 yaxis_label_color='black',
                 yaxis_label_font_weight='normal', 
                 yaxis_label_font_size=14,
                 clusterline_width = 2,
                 splitline_width = 2,
                 splitline_color = 'white',
                 xaxis_title = None,
                 xaxis_title_font_size = 20,
                 yaxis_title = None,
                 yaxis_title_font_size = 20,
                 plot_title = None,
                 plot_title_font_size = 20,
                 subfolder='plot',
                 filename = None,
                 z_score=None,
                 figsize=(10, 8),  # 新增的参数：图形尺寸，默认为 10x8 英寸
                 figure_description = 'Heatmap',
                 ):
        """
        Generate a clustered heatmap using seaborn.clustermap with flexible filtering, 
        clustering, annotation, and styling options.

        Args:
            data (pd.DataFrame): Input DataFrame containing numeric values.
            columns (list of str): Columns to include in the heatmap.
            colors (str or list, optional): Colormap name (string from `plt.colormaps()`) 
                or a list of colors. Default 'coolwarm'.
            filter_data (pd.DataFrame, optional): External DataFrame for filtering rows.
            filter_columns (list, optional): Column(s) in `filter_data` to filter by.
            filter_values (list, optional): Threshold(s) for filtering.
                - If column == "fc", keeps rows with fold-change > val or < 1/val.
                - Else, keeps rows with values < val.
            cluster (str, optional): Clustering mode: 'row', 'col', 'both', or 'none'. Default 'both'.
            cluster_method (str, optional): Linkage method for hierarchical clustering.
                One of ['single','complete','average','weighted','centroid','median','ward'].
            log (bool, optional): Apply log2 transformation to data. Default False.
            text_annotation (bool, optional): Whether to display values inside cells. Default False.
            text_size (int, optional): Font size for annotations.
            text_color (str, optional): Font color for annotations.
            minvalue (float): Minimum value for colormap.
            centervalue (float): Center value for diverging colormaps. Default 0.
            maxvalue (float): Maximum value for colormap.
            xaxis_label_show, yaxis_label_show (bool): Whether to display tick labels.
            xaxis_label_color, yaxis_label_color (str): Axis label colors.
            xaxis_label_font_weight, yaxis_label_font_weight (str): Font weights.
            xaxis_label_font_size, yaxis_label_font_size (int): Font sizes.
            xaxis_label_rotate (int): Rotation angle for x-axis tick labels. Default 90.
            clusterline_width (float): Line width of dendrograms. Default 2.
            splitline_width (float): Line width of cell borders. Default 2.
            splitline_color (str): Border color for heatmap cells. Default 'white'.
            xaxis_title (str, optional): X-axis label.
            yaxis_title (str, optional): Y-axis label.
            xaxis_title_font_size, yaxis_title_font_size (int): Axis title font sizes.
            plot_title (str, optional): Plot title.
            plot_title_font_size (int): Title font size.
            subfolder (str): Subfolder under `./plot/` to save outputs. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended.
            z_score (int, optional): Standardize rows or columns before clustering (seaborn option).
            figsize (tuple, optional): Figure size in inches. Default (10, 8).
            figure_description (str): Figure description string. Default "Heatmap".

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description string.

        Side Effects:
            - Optionally filters rows based on `filter_data`.
            - Performs hierarchical clustering on rows/columns if requested.
            - Saves heatmap in PDF and PNG formats under `./plot/<subfolder>/`.
            - Crops whitespace from PNG output.
            - Logs parameters and output path with `data_manager`.

        Notes:
            - If `text_annotation=True`, integers are shown as whole numbers,
            floats as 2 decimal places.
            - Unlike `heatmap1` (pyecharts), this function uses seaborn and
            produces static plots suitable for publications.
        """
        print('You can choose colors in following range: ', plt.colormaps())
        
        # 取选定的列
        df = data[columns].apply(pd.to_numeric, errors='coerce')
        
        if filter_data is not None:
            # 如果filter_columns或filter_values不是列表，则转换为列表
            if not isinstance(filter_columns, list):
                filter_columns = [filter_columns]
            if not isinstance(filter_values, list):
                filter_values = [filter_values]
            # 初始条件设为True
            condition = pd.Series(True, index=filter_data.index)
            for col, val in zip(filter_columns, filter_values):
                if col == 'fc':
                    condition &= ((filter_data[col] > val) | (filter_data[col] < 1/val))
                else:
                    condition &= (filter_data[col] < val)
            df = df.loc[df.index.isin(filter_data[condition].index)]

        if text_annotation == False:
            # text_annotation = None
            annot = None
            annot_kws = None
        if text_annotation == True:
            def format_value(val):
                if isinstance(val, (int, np.integer)) or (isinstance(val, float) and val.is_integer()):
                    return f"{int(val)}"
                else:
                    return f"{val:.2f}"
            annot = np.vectorize(format_value)(data)
            annot_kws = {'size': text_size, 'weight': 'normal', 'color': text_color, 'fontfamily': None} if text_annotation else None
        
        # 如果颜色未设置，提供默认颜色
        if colors is None:
            colors = 'coolwarm'
    
        # 转换为数值类型并取对数（如果需要）
        if log:
            df = df.apply(pd.to_numeric, errors='coerce')
            df = np.log2(df)
            
        df = df.replace(np.nan, 0)
    
        # 根据聚类选项进行聚类处理
        if cluster == 'row':
            # 仅对行进行层次聚类
            row_linkage = linkage(df, method=cluster_method)
            plt.figure(figsize=figsize)  # 设置图形尺寸
            g = sns.clustermap(df, row_linkage=row_linkage, col_cluster=False,
                               cmap=colors, vmin=minvalue, vmax=maxvalue, center=centervalue,
                           xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,
                           linewidths=splitline_width, linecolor=splitline_color, z_score=z_score,
                           annot=annot, annot_kws=annot_kws, fmt='')
        elif cluster == 'col':
            # 仅对列进行层次聚类
            col_linkage = linkage(df.T, method=cluster_method)
            plt.figure(figsize=figsize)  # 设置图形尺寸
            g = sns.clustermap(df, col_linkage=col_linkage, row_cluster=False,
                               cmap=colors, vmin=minvalue, vmax=maxvalue, center=centervalue,
                           xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,
                           linewidths=splitline_width, linecolor=splitline_color, z_score=z_score,
                           annot=annot, annot_kws=annot_kws, fmt='')
        elif cluster == 'both':
            # 对行和列都进行层次聚类
            row_linkage = linkage(df, method=cluster_method)
            col_linkage = linkage(df.T, method=cluster_method)
            plt.figure(figsize=figsize)  # 设置图形尺寸
            g = sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage, cmap=colors, 
                           vmin=minvalue, vmax=maxvalue, xticklabels=xaxis_label_show, center=centervalue,
                           yticklabels=yaxis_label_show,
                           linewidths=splitline_width, linecolor=splitline_color, z_score=z_score,
                           annot=annot, annot_kws=annot_kws, fmt='')
        else:
            # 不进行聚类
            plt.figure(figsize=figsize)  # 设置图形尺寸
            g = sns.clustermap(df,row_cluster=False,col_cluster=False,cmap=colors, vmin=minvalue, vmax=maxvalue, center=centervalue,
                        xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,
                        linewidths=splitline_width, linecolor=splitline_color,
                        annot=annot, annot_kws=annot_kws, fmt='')
            
        g.ax_heatmap.set_xlabel(
            xaxis_title,
            fontsize=xaxis_title_font_size,
            fontweight=xaxis_label_font_weight,
            fontfamily='Arial'
        )
        g.ax_heatmap.set_ylabel(
            yaxis_title,
            fontsize=yaxis_title_font_size,
            fontweight=yaxis_label_font_weight,
            fontfamily='Arial'
        )
        
        # 2) 设置刻度（tick）标签的字体大小、颜色、旋转
        g.ax_heatmap.tick_params(
            axis='x',
            labelrotation=xaxis_label_rotate,
            labelsize=xaxis_label_font_size,
            labelcolor=xaxis_label_color
        )
        g.ax_heatmap.tick_params(
            axis='y',
            labelrotation=0,
            labelsize=yaxis_label_font_size,
            labelcolor=yaxis_label_color
        )
        
        # （可选）如果你还想单独微调每个 tick 对象：
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(),
            rotation=xaxis_label_rotate,
            fontsize=xaxis_label_font_size,
            color=xaxis_label_color,
            fontweight=xaxis_label_font_weight
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(),
            fontsize=yaxis_label_font_size,
            color=yaxis_label_color,
            fontweight=yaxis_label_font_weight
        )
    
        # # 设置x轴标签样式
        # plt.xticks(rotation=90, fontsize=xaxis_label_font_size, color=xaxis_label_color, 
        #            weight=xaxis_label_font_weight)
        # # 设置y轴标签样式
        # plt.yticks(fontsize=yaxis_label_font_size, color=yaxis_label_color, 
        #            weight=yaxis_label_font_weight)
        
        ax = plt.gca()       
        # ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
        # ax.set_ylabel(yaxis_title,fontproperties = 'Arial', size = yaxis_title_font_size) 
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        
        if cluster in ['col', 'row', 'both']:
            for a in g.ax_row_dendrogram.collections:
                a.set_linewidth(clusterline_width)  
            for a in g.ax_col_dendrogram.collections:
                a.set_linewidth(clusterline_width) 
    
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_heatmap2_{timestamp}.pdf"), dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_heatmap2_{timestamp}.png")
        plt.savefig(png_file, dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
    
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'heatmap2', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'heatmap2', f"{columns}_heatmap2.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_heatmap2_{timestamp}.png"), 
                'legend': figure_description}
    
    def correlation_heatmap(self, data_list, colors = 'coolwarm',
                            minvalue = -1, 
                            centervalue = 0,
                             maxvalue = 1,
                             xaxis_label_show=True, 
                             xaxis_label_color='black', 
                             xaxis_label_font_weight='normal', 
                             xaxis_label_font_size=14, 
                             yaxis_label_show=False, 
                             yaxis_label_color='black',
                             yaxis_label_font_weight='normal', 
                             yaxis_label_font_size=14,
                             xaxis_title = None,
                             xaxis_title_font_size = 20,
                             yaxis_title = None,
                             yaxis_title_font_size = 20,
                             plot_title = None,
                             plot_title_font_size = 20,
                             subfolder='plot',
                             filename = None,
                             figure_description = 'Heatmap',
                             ):
        """
        Generate a circle-style correlation heatmap with Spearman correlation coefficients 
        and p-values, using multiple input DataFrames.

        Workflow:
            1. All input DataFrames are aligned to have the same column names.
            2. DataFrames are concatenated (row-wise) and transposed so that 
            each variable is in a column and each row corresponds to samples.
            3. Pairwise Spearman correlations and p-values are calculated.
            4. Heatmap is drawn:
            - Upper triangle: left blank (transparent squares).
            - Lower triangle: circles whose radius reflects correlation strength, 
                color reflects correlation sign/magnitude.
            - Each circle is annotated with correlation coefficient (r) and p-value.

        Args:
            data_list (list of pd.DataFrame): List of DataFrames to combine and analyze.
                All DataFrames must have the same number of columns.
            colors (str or Colormap): Colormap for circle coloring. Default 'coolwarm'.
            minvalue (float): Minimum value for color scale. Default -1.
            centervalue (float): Center point for diverging colormaps. Default 0.
            maxvalue (float): Maximum value for color scale. Default 1.
            xaxis_label_show, yaxis_label_show (bool): Whether to display axis tick labels.
            xaxis_label_color, yaxis_label_color (str): Color of axis labels.
            xaxis_label_font_weight, yaxis_label_font_weight (str): Font weight for axis labels.
            xaxis_label_font_size, yaxis_label_font_size (int): Font size for axis labels.
            xaxis_title, yaxis_title (str, optional): Titles for x and y axes.
            xaxis_title_font_size, yaxis_title_font_size (int): Axis title font sizes.
            plot_title (str, optional): Overall plot title.
            plot_title_font_size (int): Title font size.
            subfolder (str): Subfolder under './plot/' to save results. Default "plot".
            filename (str, optional): Base filename prefix. Timestamp is appended automatically.
            figure_description (str): Short description of the figure for logging.

        Returns:
            dict:
                - "file_path": Path to saved PNG file.
                - "legend": Figure description.

        Side Effects:
            - Saves correlation heatmap as PDF and PNG in ./plot/<subfolder>/.
            - Crops whitespace from PNG output.
            - Logs parameters and outputs via `data_manager`.

        Notes:
            - Correlation values range between [-1, 1].
            - Circle radius ∝ sqrt(|r|), so stronger correlations produce larger circles.
            - Each circle shows both r (above) and p (below) for interpretability.
            - For consistent input, all DataFrames in `data_list` are forced to share 
            identical column names before merging.
        """
        # 1. 统一重设每个 DataFrame 的列名称
        new_cols = [f"col{i}" for i in range(data_list[0].shape[1])]
        data_list_updated = []
        for df in data_list:
            df = df.copy()
            df.columns = new_cols
            data_list_updated.append(df)
        
        # 2. 合并数据（忽略原始行索引），并转置
        df_combined = pd.concat(data_list_updated, axis=0)
        df_combined = df_combined.T  # 转置后：每列为一个变量，每行为 10 个样本
    
        # 3. 计算相关性矩阵和 p 值矩阵（使用 Spearman 相关系数）
        labels = df_combined.columns
        corr_matrix = pd.DataFrame(np.zeros((len(labels), len(labels))), index=labels, columns=labels)
        pval_matrix = pd.DataFrame(np.zeros((len(labels), len(labels))), index=labels, columns=labels)
        
        for i in range(len(labels)):
            for j in range(len(labels)):
                r, p = spearmanr(df_combined.iloc[:, i], df_combined.iloc[:, j])
                corr_matrix.iloc[i, j] = r
                pval_matrix.iloc[i, j] = p
        print(f'You can choose minvalue and maxvalue from {corr_matrix.min().min()} to {corr_matrix.max().max()}.')
        # 4. 绘制圆圈相关性热图
        # 仅绘制上三角部分的背景
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
        fig, ax = plt.subplots(figsize=(6, 6))
        
        # 使用 sns.heatmap 绘制背景正方形，设置透明填充，仅显示黑色边框和颜色条
        hm = sns.heatmap(corr_matrix, mask=mask, cmap=colors, square=True, linewidths=0, 
                         linecolor='black', cbar=True, vmin=minvalue, center=centervalue, vmax=maxvalue,
                         annot=False, ax=ax,
                         xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,)
        
        # 将背景正方形填充设为透明（仅保留边框）
        for collection in hm.collections:
            collection.set_alpha(0)
        
        # 提取 heatmap 返回的 norm 对象用于颜色映射
        norm = hm.collections[0].norm
        cmap_used = plt.get_cmap(colors)
        max_radius = 0.4  # 最大圆圈半径
        
        # 绘制下三角区域的圆圈
        for i in range(len(corr_matrix)):
            for j in range(len(corr_matrix)):
                if j > i:
                    continue  # 仅绘制下三角区域（含对角线）
                r = corr_matrix.iloc[i, j]
                p = pval_matrix.iloc[i, j]
                radius = max_radius * np.sqrt(abs(r))
                color = cmap_used(norm(r))
                circle = plt.Circle((j + 0.5, i + 0.5), radius=radius, color=color)
                ax.add_artist(circle)
                # 标注相关性系数（在上面）
                ax.text(j + 0.5, i + 0.5 + 0.12, f"r={r:.2f}", ha="center", va="center", 
                        fontsize=13)
                # 标注 p 值（在下面）
                ax.text(j + 0.5, i + 0.5 - 0.12, f"p={p:.2f}", ha="center", va="center", 
                        fontsize=13)
        
        # 设置坐标轴刻度与标签
        ax.set_xticks(np.arange(len(corr_matrix)) + 0.5)
        ax.set_yticks(np.arange(len(corr_matrix)) + 0.5)
        ax.set_xticklabels(corr_matrix.columns, rotation=90, fontsize=xaxis_label_font_size, color=xaxis_label_color, 
                       weight=xaxis_label_font_weight)
        ax.set_yticklabels(corr_matrix.index,fontsize=yaxis_label_font_size, color=yaxis_label_color, 
                       weight=yaxis_label_font_weight)
        ax = plt.gca()       
        ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
        ax.set_ylabel(yaxis_title,fontproperties = 'Arial', size = yaxis_title_font_size) 
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        
        plt.tight_layout()
        # plt.show()
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_correlation_heatmap_{timestamp}.pdf"), dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_correlation_heatmap_{timestamp}.png")
        plt.savefig(png_file, dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
    
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'correlation_heatmap', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'correlation_heatmap', f"{filename}_correlation_heatmap.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_correlation_heatmap_{timestamp}.png"), 
                'legend': figure_description}
    
    def heatmap_multi_data(self, *data_names, columns, statistical_methods, colors = None, 
                count_top = 10,
                cluster = 'none',
                cluster_method = 'weighted',
                minvalue = -1,
                maxvalue = 1,
                xaxis_label_show = True,
                xaxis_label_color = 'black',
                xaxis_label_font_weight = 'normal',
                xaxis_label_font_size = 14,
                yaxis_label_show = True,
                yaxis_label_color = 'black',
                yaxis_label_font_weight = 'normal',
                yaxis_label_font_size = 14,
                clusterline_width = 2,
                splitline_width = 2,
                splitline_color = 'white',
                xaxis_title = None,
                xaxis_title_font_size = 20,
                yaxis_title = None,
                yaxis_title_font_size = 20,
                plot_title = None,
                plot_title_font_size = 20,
                annotation_font_size = 20,
                subfolder='plot',
                filename = None,
                figure_description = 'Multi data heatmap',
                ):
        """
        Create heatmaps to compare multiple datasets, either by top frequency counts
        or by user-specified statistical methods.

        Modes:
            1. **Frequency mode** (when `statistical_methods` does not contain "both"/"unique"):
            - For each dataset, compute `value_counts()` of the specified column(s).
            - Take top `count_top` values (default=10).
            - Render side-by-side heatmaps with a shared colorbar.

            2. **Statistical mode** (when `statistical_methods` contains "both"/"unique"):
            - For each dataset and column:
                * "both": count total non-empty entries.
                * "unique": count number of unique non-empty values.
                * If column contains ';', it is exploded before counting.
                * If column is a combination (e.g., "col1+col2"), values are concatenated
                    across multiple columns before counting.
            - Combine results into a single DataFrame.
            - Render a heatmap or clustered heatmap (`cluster` = 'row', 'col', 'both', 'none').

        Args:
            *data_names: str or DataFrame
                Either variable names as strings (evaluated with `eval`) or DataFrame objects directly.
            columns (list[str]): Columns to analyze (single or combined with '+').
            statistical_methods (list[str]): Statistical methods applied to each column.
                Options: "both", "unique". Must align with `columns` in length.
            colors (str or Colormap): Colormap for the heatmap. Default "coolwarm".
            count_top (int): Number of top values to include in frequency mode. Default 10.
            cluster (str): Clustering option: 'row', 'col', 'both', 'none'. Default 'none'.
            cluster_method (str): Linkage method for hierarchical clustering. Default 'weighted'.
            minvalue, maxvalue (float): Color scale limits.
            xaxis_label_show, yaxis_label_show (bool): Whether to show axis labels.
            xaxis_label_color, yaxis_label_color (str): Label colors.
            xaxis_label_font_weight, yaxis_label_font_weight (str): Label font weights.
            xaxis_label_font_size, yaxis_label_font_size (int): Label font sizes.
            clusterline_width (float): Line width for dendrograms in clustered heatmaps.
            splitline_width (float): Width of grid lines in the heatmap.
            splitline_color (str): Color of grid lines in the heatmap.
            xaxis_title, yaxis_title (str): Titles for axes.
            xaxis_title_font_size, yaxis_title_font_size (int): Font sizes for axis titles.
            plot_title (str): Overall plot title.
            plot_title_font_size (int): Title font size.
            annotation_font_size (int): Font size for annotated values in non-cluster heatmaps.
            subfolder (str): Output subfolder under `./plot/`. Default "plot".
            filename (str): Base name for saved files. Timestamp appended automatically.
            figure_description (str): Short description of the figure for logging.

        Returns:
            dict:
                - "file_path": Path to saved PNG file.
                - "legend": Figure description.

        Side Effects:
            - Saves heatmap(s) as PDF and PNG in `./plot/<subfolder>/`.
            - Crops whitespace from PNG output.
            - Logs parameters and outputs via `data_manager`.

        Notes:
            - If `columns` contains a combination (e.g., "col1+col2"), values across columns
            are concatenated before counting.
            - Exploded ';'-separated entries are treated as separate values.
            - In frequency mode, heatmaps are rendered side-by-side with a shared color scale.
        """
        print('You can choose the camp in following range: ', plt.colormaps())
        
        # 从传入的data_names中生成dataframes
        dataframes = [
                eval(name) if isinstance(name, str) else name
                for name in data_names
            ]
        dataframe_names = []
        name_count = {}
        for name in data_names:
            if '.' in name:  # 保留原逻辑
                base_name = name.split('.')[-1].split('[')[0]
                if base_name in name_count:
                    name_count[base_name] += 1
                    unique_name = f"{base_name}_{name_count[base_name]}"
                else:
                    name_count[base_name] = 1
                    unique_name = base_name
                dataframe_names.append(unique_name)
            else:  # 修改只在变量处理部分
                # 获取变量的名称
                var_name = [key for key, value in globals().items() if value is name]
                if var_name:
                    base_name = var_name[0]  # 提取变量名
                    if base_name in name_count:
                        name_count[base_name] += 1
                        unique_name = f"{base_name}_{name_count[base_name]}"
                    else:
                        name_count[base_name] = 1
                        unique_name = base_name
                    dataframe_names.append(unique_name)
                else:
                    dataframe_names.append("unknown")  # 如果变量名未知，标记为未知
    
        # 初始化存储统计量的字典，之后转换为DataFrame
        stats_data = {name: [] for name in dataframe_names}
    
        # 遍历每个列名和相应的统计方法
        if not bool(set(statistical_methods) & set(['both', 'unique'])):
            # 生成示例数据
            data1 = pd.DataFrame(dataframes[0][columns].value_counts()).iloc[:10,:]
            data2 = pd.DataFrame(dataframes[1][columns].value_counts()).iloc[:10,:]
            
            # 计算颜色条的统一范围
            vmin = min(data1.min().values[0], data2.min().values[0])
            vmax = max(data1.max().values[0], data2.max().values[0])
            # 创建图形和子图
            fig, axes = plt.subplots(1, 2, figsize=(6, 6), gridspec_kw={'width_ratios': [1, 1], 'wspace': 0})
            
            # 设置两个子图的刻度字体和轴标题
            for ax in axes:
                # x 轴刻度
                ax.tick_params(
                    axis='x',
                    labelrotation=90,
                    labelsize=xaxis_label_font_size,
                    labelcolor=xaxis_label_color,
                    width=splitline_width
                )
                # y 轴刻度
                ax.tick_params(
                    axis='y',
                    labelsize=yaxis_label_font_size,
                    labelcolor=yaxis_label_color,
                    width=splitline_width
                )
                # （可选）如果你想给每个子图单独加 axis title
                if xaxis_title:
                    ax.set_xlabel(xaxis_title, fontsize=xaxis_title_font_size, fontfamily='Arial')
                if yaxis_title:
                    ax.set_ylabel(yaxis_title, fontsize=yaxis_title_font_size, fontfamily='Arial')

            # 绘制第一个热图
            sns.heatmap(data1, ax=axes[0], cmap="coolwarm", cbar=False, vmin=vmin, vmax=vmax,
                        annot=True, fmt=".0f", annot_kws={"size": annotation_font_size, 'color':'black', 'fontfamily':'Arial'})
            axes[0].set_title(dataframe_names[0], fontsize=annotation_font_size)
            axes[0].tick_params(axis='y', labelrotation=0)
            # 绘制第二个热图
            sns.heatmap(data2, ax=axes[1], cmap="coolwarm", cbar=False, vmin=vmin, vmax=vmax,
                        annot=True, fmt=".0f", annot_kws={"size": annotation_font_size, 'color':'black', 'fontfamily':'Arial'})
            axes[1].set_title(dataframe_names[1], fontsize=annotation_font_size)
            axes[1].tick_params(axis='y', labelrotation=0, labelright=True, labelleft=False)
            axes[1].tick_params(left=False)
            axes[1].set_ylabel("")  # 去掉 y 轴标签
            # 在底部添加统一的颜色条
            cbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.02])  # 调整位置 [left, bottom, width, height]
            sns.heatmap(data1, cbar_ax=cbar_ax, cbar=False, cmap="coolwarm", vmin=vmin, vmax=vmax, 
                        cbar_kws={'orientation': 'horizontal', 'label': 'Value'})
            # 调整布局
            plt.tight_layout(rect=[0, 0.1, 1, 1])  # 留出颜色条的空间
            cbar_ax.remove()
            # plt.show()

        if bool(set(statistical_methods) & set(['both', 'unique'])):
            for column, method in zip(columns, statistical_methods):  
                for idx, df in enumerate(dataframes):
                    if '+' in column:
                        # 如果列名包含 '+'，则需要组合多列
                        col_parts = column.split('+')
                        # 确认每个部分在数据集中都存在
                        if all(part in df.columns for part in col_parts):
                            # 创建一个临时 DataFrame
                            temp_data = df[col_parts].astype(str).copy()
                            # 对每一列，若含有 ';'，则分割并扩展
                            explode_cols = []
                            for part in col_parts:
                                if temp_data[part].str.contains(';').any():
                                    temp_data[part] = temp_data[part].str.split(';')
                                    explode_cols.append(part)
                            if explode_cols:
                                temp_data = temp_data.explode(explode_cols)
                            # 合并多个列为一个字符串列
                            combined_values = temp_data[col_parts[0]].astype(str)
                            for part in col_parts[1:]:
                                combined_values += temp_data[part].astype(str)
                            # 根据方法统计
                            if method == 'unique':
                                count = combined_values.dropna().nunique()
                            elif method == 'both':
                                count = combined_values.dropna().shape[0]
                        else:
                            # 如果列不存在，返回 NaN
                            count = float('nan')
                    else:
                        # 处理没有 '+' 的单列
                        if column in df.columns:
                            # 创建一个临时 DataFrame
                            temp_data = df[column].astype(str).copy()
                            # 对每一列，若含有 ';'，则分割并扩展
                            if temp_data.str.contains(';').any():  # 检查是否包含分号
                                temp_data = temp_data.str.split(';')
                                temp_data = temp_data.explode()
                            if method == 'both':
                                # 计算非空值的行数
                                count = temp_data.dropna().shape[0]
                            elif method == 'unique':
                                # 计算唯一值的个数
                                count = temp_data.dropna().nunique()
                            else:
                                count = float('nan')
                        else:
                            count = float('nan')
                    
                    # 将计算结果加入对应的dataframe列下
                    stats_data[dataframe_names[idx]].append(count)
        
            # 将统计结果字典转换为DataFrame
            df = pd.DataFrame(stats_data, index=columns)
        
            minvalue = df.min().min()
            maxvalue = df.max().max()
            
            # 如果颜色未设置，提供默认颜色
            if colors is None:
                colors = 'coolwarm'
    
            # 根据聚类选项进行聚类处理
            if cluster == 'row':
                # 仅对行进行层次聚类
                row_linkage = linkage(df, method=cluster_method)
                g = sns.clustermap(df, row_linkage=row_linkage, cmap=colors, vmin=minvalue, vmax=maxvalue, 
                               xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,
                               linewidths=splitline_width, linecolor=splitline_color)
            elif cluster == 'col':
                # 仅对列进行层次聚类
                col_linkage = linkage(df.T, method=cluster_method)
                g = sns.clustermap(df, col_linkage=col_linkage, cmap=colors, vmin=minvalue, vmax=maxvalue, 
                               xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,
                               linewidths=splitline_width, linecolor=splitline_color)
            elif cluster == 'both':
                # 对行和列都进行层次聚类
                row_linkage = linkage(df, method=cluster_method)
                col_linkage = linkage(df.T, method=cluster_method)
                g = sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage, cmap=colors, 
                               vmin=minvalue, vmax=maxvalue, xticklabels=xaxis_label_show, 
                               yticklabels=yaxis_label_show,
                               linewidths=splitline_width, linecolor=splitline_color)
            else:
                # 不进行聚类
                g = sns.heatmap(df, cmap=colors, vmin=minvalue, vmax=maxvalue, annot=True, fmt='d',
                            xticklabels=xaxis_label_show, yticklabels=yaxis_label_show,
                            linewidths=splitline_width, linecolor=splitline_color,
                            annot_kws={"size": annotation_font_size, 'fontfamily':'Arial'})
            if hasattr(g, 'ax_heatmap'):
                ax_heatmap = g.ax_heatmap
            else:
                ax_heatmap = g    
            # 1) 坐标轴标题
            ax_heatmap.set_xlabel(
                xaxis_title, 
                fontsize=xaxis_title_font_size, 
                fontfamily='Arial'
            )
            ax_heatmap.set_ylabel(
                yaxis_title, 
                fontsize=yaxis_title_font_size, 
                fontfamily='Arial'
            )
        
            # 2) 刻度标签
            ax_heatmap.tick_params(
                axis='x',
                labelrotation=90,
                labelsize=xaxis_label_font_size,
                labelcolor=xaxis_label_color
            )
            ax_heatmap.tick_params(
                axis='y',
                labelsize=yaxis_label_font_size,
                labelcolor=yaxis_label_color
            )

            # # 设置x轴标签样式
            # plt.xticks(rotation=45, fontsize=xaxis_label_font_size, color=xaxis_label_color, 
            #            weight=xaxis_label_font_weight)
    
            # # 设置y轴标签样式
            # plt.yticks(fontsize=yaxis_label_font_size, color=yaxis_label_color, 
            #            weight=yaxis_label_font_weight)
            
            ax = plt.gca()       
            # ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
            # ax.set_ylabel(yaxis_title,fontproperties = 'Arial', size = yaxis_title_font_size) 
            ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
            
            if cluster != 'none':
                for a in g.ax_row_dendrogram.collections:
                    a.set_linewidth(clusterline_width)  
                for a in g.ax_col_dendrogram.collections:
                    a.set_linewidth(clusterline_width) 
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_heatmap_multi_data_{timestamp}.pdf"), dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_heatmap_multi_data_{timestamp}.png")
        plt.savefig(png_file, dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'heatmap_multi_data', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'heatmap_multi_data',f"{columns}_heatmap_multi_data.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_heatmap_multi_data_{timestamp}.png"), 
                'legend': figure_description}
    
    def complexheatmap(self, data, columns, 
                       row_annotation_data = None,
                       row_annotation_data_log2 = None, 
                       row_annotation_plot_type = None,
                       col_annotation_data = None,
                       col_annotation_data_log2 = None, 
                       col_annotation_plot_type = None,
                       z_score = 0,
                       log2 = False,
                       col_cluster = True,
                       row_cluster = True,
                       col_cluster_method = 'ward',
                       row_cluster_method = 'ward',
                       show_colnames = True,
                       show_rownames = True, 
                       col_dendrogram = True,
                       row_dendrogram = True,
                       col_split = 2, 
                       row_split = 2,
                       tree_line_cmap = 'Dark2',
                       tree_line_width = 3,
                       col_split_gap = 2, 
                       row_split_gap = 2,
                       linewidths = 0.1,
                       linecolor = 'white',
                       cmap = 'coolwarm', 
                       subfolder='plot',
                       filename = None,
                       figure_description = 'Complex heatmap',
                       ):
        """
        Create a complex clustered heatmap with optional row/column annotations.

        This function uses **pyComplexHeatmap (pch)** to generate advanced heatmaps 
        supporting row/column dendrograms, clustering, splitting, custom annotations, 
        and flexible styling. It is suitable for high-dimensional omics data visualization.

        Args:
            data (pd.DataFrame): Input dataframe.
            columns (list[str]): Columns to include in the heatmap.
            row_annotation_data (list[pd.Series] or None): Optional row annotation data.
            row_annotation_data_log2 (list[bool] or None): Whether to log2-transform each row annotation.
            row_annotation_plot_type (list[str] or None): Plot type for each row annotation.
                Options: "bar", "scatter".
            col_annotation_data (list[pd.Series] or None): Optional column annotation data.
            col_annotation_data_log2 (list[bool] or None): Whether to log2-transform each column annotation.
            col_annotation_plot_type (list[str] or None): Plot type for each column annotation.
                Options: "bar", "scatter".
            z_score (int): Apply z-score normalization (0=row-wise, 1=col-wise, None=off).
            log2 (bool): Apply log2 transformation to the main data.
            col_cluster (bool): Whether to cluster columns.
            row_cluster (bool): Whether to cluster rows.
            col_cluster_method (str): Linkage method for column clustering. Default 'ward'.
            row_cluster_method (str): Linkage method for row clustering. Default 'ward'.
            show_colnames (bool): Whether to display column names.
            show_rownames (bool): Whether to display row names.
            col_dendrogram (bool): Whether to display column dendrogram.
            row_dendrogram (bool): Whether to display row dendrogram.
            col_split (int): Number of groups to split columns into. Default 2.
            row_split (int): Number of groups to split rows into. Default 2.
            tree_line_cmap (str): Colormap for dendrogram tree lines. Default "Dark2".
            tree_line_width (float): Width of dendrogram tree lines. Default 3.
            col_split_gap (int): Gap between column split groups.
            row_split_gap (int): Gap between row split groups.
            linewidths (float): Width of grid lines in the heatmap. Default 0.1.
            linecolor (str): Color of grid lines in the heatmap. Default "white".
            cmap (str): Colormap for the heatmap. Default "coolwarm".
            subfolder (str): Output subfolder under `./plot/`. Default "plot".
            filename (str): Base name for output files. Timestamp appended automatically.
            figure_description (str): Short description of the figure for logging.

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description.

        Side Effects:
            - Saves the heatmap as PDF and PNG in `./plot/<subfolder>/`.
            - Trims whitespace around the PNG.
            - Logs parameters and outputs via `data_manager`.

        Notes:
            - Row and column annotations can be multiple. 
            Each is dynamically named (Bar1, Bar2, ...).
            - Supports log2-transformation separately for main data and annotation data.
            - Dendrogram styling controlled by `tree_line_cmap` and `tree_line_width`.
            - Splitting (`row_split`, `col_split`) is useful for grouped datasets.
        """
        print('You can choose the camp in following range: ', plt.colormaps())
        
        # 主要数据处理
        data1 = data[columns]
        data1 = data1.apply(pd.to_numeric, errors='coerce')
        if log2:
            data1 = np.log2(data1)

        # 处理注释数据
        if row_annotation_data is None:
            row_annotation_data = []
        if row_annotation_plot_type is None:
            row_annotation_plot_type = ['bar'] * len(row_annotation_data)  
        if row_annotation_data_log2 is None:
            row_annotation_data_log2 = [False] * len(row_annotation_data)  
            
        if col_annotation_data is None:
            col_annotation_data = []
        if col_annotation_plot_type is None:
            col_annotation_plot_type = ['bar'] * len(col_annotation_data)  
        if col_annotation_data_log2 is None:
            col_annotation_data_log2 = [False] * len(col_annotation_data)  

        # 动态生成 annotation
        row_annotation_dict = {}
        for i, (anno_data, plot_type, log2) in enumerate(zip(row_annotation_data, row_annotation_plot_type, row_annotation_data_log2)):
            column_name = f'Bar{i + 1}'  # 动态生成Bar的名称
            if log2:
                anno_data = np.log2(anno_data)  # 进行 log2 转换
            if plot_type == 'bar':
                row_annotation_dict[column_name] = pch.anno_barplot(anno_data, height=18, cmap='winter', legend=True)
            elif plot_type == 'scatter':
                row_annotation_dict[column_name] = pch.anno_scatterplot(anno_data, height=18, cmap='cool', legend=True)
        col_annotation_dict = {}
        for i, (anno_data, plot_type, log2) in enumerate(zip(col_annotation_data, col_annotation_plot_type, col_annotation_data_log2)):
            column_name = f'Bar{i + 1}'  # 动态生成Bar的名称
            if log2:
                anno_data = np.log2(anno_data)  # 进行 log2 转换
            if plot_type == 'bar':
                col_annotation_dict[column_name] = pch.anno_barplot(anno_data, height=18, cmap='winter', legend=True)
            elif plot_type == 'scatter':
                col_annotation_dict[column_name] = pch.anno_scatterplot(anno_data, height=18, cmap='cool', legend=True)

        # 生成 HeatmapAnnotation
        # row_ha = pch.HeatmapAnnotation(
        #     **row_annotation_dict,  # 动态生成的注释数据
        #     label_kws={'rotation': 15, 'horizontalalignment': 'left', 'verticalalignment': 'bottom'},
        #     axis=0,
        #     verbose=0
        # )
        
        if row_annotation_dict:
            row_ha = pch.HeatmapAnnotation(
                **row_annotation_dict,  
                label_kws={'rotation': 15, 'horizontalalignment': 'left', 'verticalalignment': 'bottom'},
                axis=0,wgap=2,
                verbose=0
            )
        else:
            row_ha = None  
            
        if col_annotation_dict:
            col_ha = pch.HeatmapAnnotation(
                **col_annotation_dict,  
                label_kws={'rotation': 15, 'horizontalalignment': 'left', 'verticalalignment': 'bottom'},
                axis=1,
                verbose=0,hgap=2
            )
        else:
            col_ha = None  
        
        # 绘制 ClusterMap
        plt.figure(figsize=(9, 6))
        cm = pch.ClusterMapPlotter(data=data1, z_score = z_score,
                                   right_annotation=row_ha if row_ha else None, 
                                   top_annotation=col_ha if col_ha else None,
                                   col_cluster = col_cluster, row_cluster = row_cluster, 
                                   col_cluster_method = col_cluster_method, row_cluster_method = row_cluster_method,
                                   show_rownames = show_rownames, row_names_side='left',
                                   show_colnames = show_colnames, col_names_side='bottom',
                                   col_dendrogram = col_dendrogram, row_dendrogram = row_dendrogram,
                                   col_split = col_split, row_split = row_split, 
                                   tree_kws = {'row_cmap': tree_line_cmap, 'linewidth': tree_line_width},
                                   col_split_gap = col_split_gap, row_split_gap = row_split_gap, 
                                   cmap = cmap, rasterized=True, legend=True,
                                   xticklabels_kws={'labelrotation': -90, 'labelcolor': 'black'},
                                   yticklabels_kws={'labelcolor': 'black'},
                                   linewidths = linewidths, linecolor = linecolor,
                                   )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_complexheatmap_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_complexheatmap_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'complexheatmap', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'complexheatmap',f"{columns}_complexheatmap.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_complexheatmap_{timestamp}.png"), 
                'legend': figure_description}
    
    def line(self, data, y_column, x_columns, colors = None, 
              if_stack = False,
              is_smooth = False,
              area_opacity = 0.5,
              line_width = 2,
              symbol = 'circle',
              symbol_size = 15,
              is_step = False,
              label_show = False,
              label_color = 'black',
              label_font_weight = 'normal',
              label_font_size = 12,
              bar_width = "10%",
              category_gap = "40%",
              xaxis_label_show = True,
              xaxis_label_color = 'black',
              xaxis_label_font_weight = 'normal',
              xaxis_label_font_size = 9,
              xaxis_label_rotate = -15,
              xaxis_label_margin = 15,
              xaxis_label_text_split = 0,
              xaxis_line_show = True,
              xaxis_line_width = 2,
              xaxis_line_color = 'black',
              xaxis_tick_show = True,
              xaxis_tick_length = 8,
              xaxis_tick_width = 2,
              xaxis_tick_color = 'black',
              xaxis_splitline_show = True,
              xaxis_splitline_width = 0.5,
              xaxis_splitline_color = 'grey',
              top_xaxis_line_show = True,
              yaxis_label_show = True,
              yaxis_label_color = 'black',
              yaxis_label_font_weight = 'normal',
              yaxis_label_font_size = 14,
              yaxis_label_margin = 15,
              yaxis_line_show = True,
              yaxis_line_width = 2,
              yaxis_line_color = 'black',
              yaxis_tick_show = True,
              yaxis_tick_length = 8,
              yaxis_tick_width = 2,
              yaxis_tick_color = 'black',
              yaxis_splitline_show = True,
              yaxis_splitline_width = 0.5,
              yaxis_splitline_color = 'grey',
              right_yaxis_line_show = True,
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              xaxis_title = None,
              xaxis_title_gap = 25,
              yaxis_title = None,
              yaxis_title_gap = 40,
              subfolder='plot',
              filename = None,
              figure_description = 'Line chart',
              ):
        """
        Create a customizable line chart with multiple series.

        This function uses **pyecharts Line** to draw line plots supporting:
        - Multiple series defined by `y_column` groups
        - Optional stacking (area stacked line chart)
        - Smooth curves or step lines
        - Customizable markers, labels, and axis styles
        - Export to HTML, SVG, PDF, PNG with white border trimming

        Args:
            data (pd.DataFrame): Input data.
            y_column (str): Column defining different series (categories).
            x_columns (list[str]): Columns used as x-axis points.
            colors (list[str] or None): Custom colors for series. Default color palette if None.
            if_stack (bool): Whether to stack the series. Default False.
            is_smooth (bool): Whether to smooth the lines. Default False.
            area_opacity (float): Opacity of the filled area under the line. Default 0.5.
            line_width (int): Line thickness. Default 2.
            symbol (str): Marker symbol type. Options: 'circle', 'rect', 'roundRect', 
                        'triangle', 'diamond', 'pin', 'arrow', 'none'.
            symbol_size (int): Marker size. Default 15.
            is_step (bool): Whether to draw step line. Default False.
            label_show (bool): Whether to display labels on points. Default False.
            label_color (str): Label text color.
            label_font_weight (str): Label font weight. Options: 'normal', 'bold', 'bolder', 'lighter'.
            label_font_size (int): Label font size.
            xaxis_* / yaxis_* (various): Control label, line, tick, and splitline styles.
            legend (list[str] or None): Custom legend names. Defaults to series names.
            legend_font_size (int): Legend text size. Default 12.
            plot_title (str): Plot title.
            xaxis_title (str): X-axis title.
            yaxis_title (str): Y-axis title.
            subfolder (str): Output subfolder under `./plot/`. Default "plot".
            filename (str): Base filename for output files.
            figure_description (str): Figure description for logging.

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description.

        Side Effects:
            - Saves interactive chart as HTML, static images as SVG/PDF/PNG.
            - PNG image is post-processed to remove white border.
            - Parameters and outputs are logged via `data_manager`.

        Notes:
            - Supports multi-series plotting by grouping on `y_column`.
            - Use `if_stack=True` to create stacked area line charts.
            - Use `is_smooth=True` for smooth curves, `is_step=True` for step lines.
            - X-axis labels can be wrapped with `xaxis_label_text_split`.
        """
        data = data.dropna(axis=1, how='all')
        data = data.replace(np.nan, 0)
        
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(data[y_column].unique())]
        
        if if_stack in [True, False]:
            if if_stack:
                stack = 'stack1'
                bar_width = None
            else:
                stack = None
                bar_width = bar_width
                
        if is_smooth not in [True, False]:
            is_smooth = False
            
        if symbol not in ['circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none']:
            symbol = 'circle'
                
        if label_show not in [True, False]:
            label_show = False
        
        if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            label_font_weight = 'normal' 
            
        if is_step not in [True, False]:
            is_step = False

        # x_columns = list(set(x_columns) & set(data.columns))

        def generate_all_lists(data, y_column, x_columns):
            # 获取唯一的 Core_structure
            unique_core_structures = data[y_column].tolist()
            # 初始化大列表
            all_lists = []
            for core_structure in unique_core_structures:
                all_lists.append(list(data[data[y_column] == core_structure][x_columns].values[0]))
        
            return all_lists

        all_lists = generate_all_lists(data, y_column, x_columns)
        
        c = Line(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff',
            ))

        if xaxis_label_text_split:
            xaxis_labels = [
                '\n'.join([label[i:i+xaxis_label_text_split] 
                           for i in range(0, len(label), xaxis_label_text_split)])
                for label in list(x_columns)
            ]
        else:
            xaxis_labels = list(x_columns)
        c.add_xaxis(xaxis_labels)
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        if legend is None:
            legend = list(data[y_column])
        
        # 动态添加 y 轴数据
        for idx, glycan_list in enumerate(all_lists):
            c.add_yaxis(legend[idx], 
                        glycan_list, 
                        stack = stack,
                        is_smooth = is_smooth,
                        areastyle_opts = opts.AreaStyleOpts(
                                            opacity = area_opacity,
                                         ),
                        linestyle_opts = opts.LineStyleOpts(
                            width = line_width, 
                            ),
                        symbol = symbol,
                        symbol_size = symbol_size,
                        is_step = is_step,
                        )
        c.set_series_opts(
            label_opts=opts.LabelOpts(
                is_show=label_show,
                color=label_color,
                font_family='Arial',
                font_weight=label_font_weight,
                font_size=label_font_size,
                position="right",
                formatter=JsCode(
                    "function(x){return Number(x.data.percent * 100).toFixed(2) + '%';}"
                ),
            )
        )
        
        c.extend_axis(list(x_columns), 
            xaxis=opts.AxisOpts(
                type_="category",
                position='top',
                axisline_opts=opts.AxisLineOpts(
                    is_show = top_xaxis_line_show,
                    is_on_zero = False, 
                    linestyle_opts = opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                ),
                axislabel_opts=opts.LabelOpts(
                      is_show=False),
                axistick_opts = opts.AxisTickOpts(
                    is_show = False,
                    ),
                splitline_opts = opts.SplitLineOpts(is_show=False),
            )
        )
        
        c.extend_axis(yaxis=opts.AxisOpts(position="right",
                                          axisline_opts=opts.AxisLineOpts(
                                                is_show = right_yaxis_line_show,
                                                linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                                            ),
                                          axislabel_opts = opts.LabelOpts(
                                                is_show=False,
                                                ),
                                          splitline_opts = opts.SplitLineOpts(is_show=False),
                                          )
                      )

        c.set_colors(colors)
        
        c.set_global_opts(
            xaxis_opts=opts.AxisOpts(
                name = xaxis_title,
                name_location='center',
                name_gap = xaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = xaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = xaxis_label_font_weight,
                                                    font_size = xaxis_label_font_size),
                position='bottom',
                axislabel_opts = opts.LabelOpts(
                    is_show = xaxis_label_show,
                    font_size = xaxis_label_font_size,  
                    color = xaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = xaxis_label_font_weight,
                    rotate = xaxis_label_rotate,
                    margin = xaxis_label_margin,
                    interval = 0,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = xaxis_line_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_line_width, 
                                                      color = xaxis_line_color) 
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = xaxis_tick_show,
                    length = xaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_tick_width , 
                                                      color = xaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = xaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_splitline_width, 
                                                      color = xaxis_splitline_color)
                    )
            ),
            yaxis_opts=opts.AxisOpts(
                name = yaxis_title,
                name_location='center',
                name_gap = yaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                position='left',
                axislabel_opts=opts.LabelOpts(
                    is_show = yaxis_label_show,
                    font_size = yaxis_label_font_size,  
                    color = yaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = yaxis_label_font_weight,
                    margin = yaxis_label_margin,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = yaxis_line_show,
                    linestyle_opts = opts.LineStyleOpts(
                                                      width = yaxis_line_width, 
                                                      color = yaxis_line_color)  
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = yaxis_tick_show,
                    length =  yaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_tick_width, 
                                                      color = yaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = yaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_splitline_width, 
                                                      color = yaxis_splitline_color)
                    )
            ),
            title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
            ),
            legend_opts=opts.LegendOpts(
                pos_right="center",  # 图例靠右
                pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
            ),
        )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 文件名保存和渲染
        file_name = os.path.join(output_dir, f"{filename}_line_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_line_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        # 将 SVG 转换为 PDF
        pdf_file = os.path.join(output_dir, f"{filename}_line_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  # 读取 SVG 文件
        renderPDF.drawToFile(drawing, pdf_file)  # 保存为 PDF 文件
        
        png_file = os.path.join(output_dir, f"{filename}_line_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'line', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'line',f"{y_column}_line.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def scatter(self, data, group_column, x_column, y_column, colors = None, 
              symbols = ['circle', 'rect', 'triangle'],
              symbol_size = 15,
              is_step = False,
              label_show = False,
              label_color = 'black',
              label_font_weight = 'normal',
              label_font_size = 12,
              bar_width = "10%",
              category_gap = "40%",
              xaxis_label_show = True,
              xaxis_label_color = 'black',
              xaxis_label_font_weight = 'normal',
              xaxis_label_font_size = 9,
              xaxis_label_rotate = -15,
              xaxis_label_margin = 15,
              xaxis_label_text_split = 0,
              xaxis_line_show = True,
              xaxis_line_width = 2,
              xaxis_line_color = 'black',
              xaxis_tick_show = True,
              xaxis_tick_length = 8,
              xaxis_tick_width = 2,
              xaxis_tick_color = 'black',
              xaxis_splitline_show = True,
              xaxis_splitline_width = 0.5,
              xaxis_splitline_color = 'grey',
              top_xaxis_line_show = True,
              yaxis_label_show = True,
              yaxis_label_color = 'black',
              yaxis_label_font_weight = 'normal',
              yaxis_label_font_size = 14,
              yaxis_label_margin = 15,
              yaxis_line_show = True,
              yaxis_line_width = 2,
              yaxis_line_color = 'black',
              yaxis_tick_show = True,
              yaxis_tick_length = 8,
              yaxis_tick_width = 2,
              yaxis_tick_color = 'black',
              yaxis_splitline_show = True,
              yaxis_splitline_width = 0.5,
              yaxis_splitline_color = 'grey',
              right_yaxis_line_show = True,
              legend = None,
              legend_font_size = 12,
              plot_title = None,
              xaxis_title = None,
              xaxis_title_gap = 25,
              yaxis_title = None,
              yaxis_title_gap = 40,
              subfolder='plot',
              filename = None,
              figure_description = 'Scatter plot',
              ):
        """
        Create a customizable scatter plot grouped by a categorical column.

        This function uses **pyecharts Scatter** to plot scatter points grouped by 
        `group_column`, with support for multiple marker shapes, colors, and 
        comprehensive axis/legend styling. It can export plots to HTML, SVG, PDF, and PNG.

        Args:
            data (pd.DataFrame): Input dataset.
            group_column (str): Column to group data (defines separate scatter series).
            x_column (str): Column for x-axis values.
            y_column (str): Column for y-axis values.
            colors (list[str] or None): Custom colors per group. Defaults to a predefined palette.
            symbols (list[str]): Marker symbols per group. 
                                Options include: 'circle', 'rect', 'roundRect', 'triangle', 
                                'diamond', 'pin', 'arrow', 'none'.
            symbol_size (int): Marker size. Default 15.
            label_show (bool): Whether to show labels on scatter points.
            label_color (str): Label text color.
            label_font_weight (str): Font weight for labels. Options: 'normal', 'bold', 'bolder', 'lighter'.
            label_font_size (int): Font size for labels.
            xaxis_* / yaxis_* (various): Parameters to control axis lines, ticks, gridlines, and labels.
            legend (list[str] or None): Custom legend labels. Defaults to unique values in `group_column`.
            legend_font_size (int): Font size of the legend text.
            plot_title (str): Title of the plot.
            xaxis_title (str): Label for x-axis.
            yaxis_title (str): Label for y-axis.
            subfolder (str): Output subfolder under `./plot/`. Default "plot".
            filename (str): Base filename for saving output files.
            figure_description (str): Textual description of the figure, recorded in logs.

        Returns:
            dict:
                - "file_path": Path to the saved PNG file.
                - "legend": Figure description.

        Side Effects:
            - Exports interactive chart as HTML, and static images as SVG/PDF/PNG.
            - PNG image is post-processed to trim white borders.
            - Logs parameters and outputs using `data_manager`.

        Notes:
            - Each unique value in `group_column` is plotted as a separate scatter series.
            - Use `symbols` and `colors` to customize group-wise appearance.
            - Axis labels can be wrapped using `xaxis_label_text_split`.
        """
        data = data.replace(np.nan, 0)
        
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"][:len(data[y_column].unique())]
        
        print("You can choose symbol in: ['circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none']")
            
        if label_show not in [True, False]:
            label_show = False
        
        if label_font_weight not in ['normal', 'bold', 'bolder', 'lighter']:
            label_font_weight = 'normal' 
        
        c = Scatter(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        if legend is None:
            legend = list(data[group_column].unique())
        
        # 按 group_column（如 sex）进行分组并绘制散点图
        for idx, group_value in enumerate(data[group_column].unique()):
            subset = data[data[group_column] == group_value]
            c.add_xaxis(list(subset[x_column]))  # 只在第一次添加 X 轴数据
            c.add_yaxis(
                str(legend[idx]),
                list(subset[y_column]),
                symbol=symbols[idx],
                symbol_size = symbol_size,
                itemstyle_opts=opts.ItemStyleOpts(color=colors[idx]),
            )
        c.set_series_opts(
            label_opts=opts.LabelOpts(
                is_show=label_show,
                color=label_color,
                font_family='Arial',
                font_weight=label_font_weight,
                font_size=label_font_size,
                ),
        )
        
        c.extend_axis(list(x_column), 
            xaxis=opts.AxisOpts(
                type_="value",
                position='top',
                axisline_opts=opts.AxisLineOpts(
                    is_show = top_xaxis_line_show,
                    is_on_zero = False, 
                    linestyle_opts = opts.LineStyleOpts(width=xaxis_line_width, color=xaxis_line_color)
                ),
                axislabel_opts=opts.LabelOpts(
                      is_show=False),
                axistick_opts = opts.AxisTickOpts(
                    is_show = False,
                    ),
                splitline_opts = opts.SplitLineOpts(is_show=False),
            )
        )
        
        c.extend_axis(yaxis=opts.AxisOpts(position="right",
                                          axisline_opts=opts.AxisLineOpts(
                                                is_show = right_yaxis_line_show,
                                                linestyle_opts=opts.LineStyleOpts(width=yaxis_line_width, color=yaxis_line_color)
                                            ),
                                          axislabel_opts = opts.LabelOpts(
                                                is_show=False,
                                                ),
                                          splitline_opts = opts.SplitLineOpts(is_show=False),
                                          )
                      )

        c.set_colors(colors)
        
        if xaxis_label_text_split != 0:
            wrapped_labels = ['\n'.join([label[i:i+xaxis_label_text_split] for i in range(0, len(label), xaxis_label_text_split)]) 
                      for label in list(x_columns)]
            c.add_xaxis(wrapped_labels)
        
        c.set_global_opts(
            xaxis_opts=opts.AxisOpts(
                name = xaxis_title,
                name_location='center',
                name_gap = xaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = xaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = xaxis_label_font_weight,
                                                    font_size = xaxis_label_font_size),
                type_ = 'value',
                # min_ = 30, 
                # max_ = 90,
                position='bottom',
                axislabel_opts = opts.LabelOpts(
                    is_show = xaxis_label_show,
                    font_size = xaxis_label_font_size,  
                    color = xaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = xaxis_label_font_weight,
                    rotate = xaxis_label_rotate,
                    margin = xaxis_label_margin,
                    interval = 0,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = xaxis_line_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_line_width, 
                                                      color = xaxis_line_color) 
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = xaxis_tick_show,
                    length = xaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_tick_width , 
                                                      color = xaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = xaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = xaxis_splitline_width, 
                                                      color = xaxis_splitline_color)
                    )
            ),
            yaxis_opts=opts.AxisOpts(
                name = yaxis_title,
                name_location='center',
                name_gap = yaxis_title_gap,
                name_textstyle_opts = opts.TextStyleOpts(color = yaxis_label_color, 
                                                    font_family = 'Arial',
                                                    font_weight = yaxis_label_font_weight,
                                                    font_size = yaxis_label_font_size),
                type_ = 'value',
                # min_ = 150, 
                # max_ = 190,
                position='left',
                axislabel_opts=opts.LabelOpts(
                    is_show = yaxis_label_show,
                    font_size = yaxis_label_font_size,  
                    color = yaxis_label_color,  
                    font_family = 'Arial',
                    font_weight = yaxis_label_font_weight,
                    margin = yaxis_label_margin,
                ),
                axisline_opts=opts.AxisLineOpts(
                    is_show = yaxis_line_show,
                    linestyle_opts = opts.LineStyleOpts(
                                                      width = yaxis_line_width, 
                                                      color = yaxis_line_color)  
                ),
                axistick_opts = opts.AxisTickOpts(
                    is_show = yaxis_tick_show,
                    length =  yaxis_tick_length,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_tick_width, 
                                                      color = yaxis_tick_color)
                    ),
                splitline_opts = opts.SplitLineOpts(
                    is_show = yaxis_splitline_show,
                    linestyle_opts=opts.LineStyleOpts(
                                                      width = yaxis_splitline_width, 
                                                      color = yaxis_splitline_color)
                    )
            ),
            title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="center",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
            ),
            legend_opts=opts.LegendOpts(
                pos_right="center",  # 图例靠右
                pos_top="top",
                    textstyle_opts = opts.TextStyleOpts(
                                                    font_family = 'Arial',
                                                    font_weight = 'normal',
                                                    font_size = legend_font_size),
            ),

        )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)   
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 文件名保存和渲染
        file_name = os.path.join(output_dir, f"{filename}_scatter_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_scatter_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        # 将 SVG 转换为 PDF
        pdf_file = os.path.join(output_dir, f"{filename}_scatter_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  # 读取 SVG 文件
        renderPDF.drawToFile(drawing, pdf_file)  # 保存为 PDF 文件
        
        png_file = os.path.join(output_dir, f"{filename}_scatter_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'scatter', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'scatter',f"{y_column}_scatter.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def up_down_scatter(self, datasets, labels, 
                           fc_threshold=1, 
                           rectangle_width=0.55,
                           show_xaxis=False, 
                           spine_width=1, 
                           ytick_labelsize=10, 
                           scatter_size=15, 
                           scatter_edgecolor='k', 
                           up_color='red', 
                           down_color='blue',
                           rectangle_colors=None, 
                           bbox_facecolor='orange', 
                           bbox_textsize=11,
                           xaxis_title = None,
                           xaxis_title_font_size = 20,
                           yaxis_title = None,
                           yaxis_title_font_size = 20,
                           plot_title = None,
                           plot_title_font_size = 20,
                           legend = None,
                           legend_fontsize = 12,
                           subfolder='plot',
                           filename = None,
                           figure_description = 'Up and down scatter plot',
                           ):
        """
        Create an up/down scatter plot to visualize fold-change distributions across datasets.

        Each dataset is plotted along the x-axis with log2 fold-change (FC) values on the y-axis.
        Points above `fc_threshold` are considered "up" (colored `up_color`), while points below 
        are "down" (colored `down_color`). For each dataset, the up/down counts are displayed 
        at the top/bottom, and a label is drawn in the center with a bounding box.

        Args:
            datasets (list[pd.DataFrame]): List of datasets, each containing a 'fc' column.
            labels (list[str]): Labels corresponding to each dataset.
            fc_threshold (float): Fold-change cutoff for up/down classification. Default 1.
            rectangle_width (float): Width of the background rectangle for each dataset.
            show_xaxis (bool): Whether to show the x-axis line and ticks.
            spine_width (int): Width of plot spines (borders).
            ytick_labelsize (int): Font size of y-axis tick labels.
            scatter_size (int): Size of scatter points.
            scatter_edgecolor (str): Edge color of scatter points.
            up_color (str): Color for upregulated points.
            down_color (str): Color for downregulated points.
            rectangle_colors (list[str] or None): Background rectangle colors per dataset.
            bbox_facecolor (str): Background color of dataset labels.
            bbox_textsize (int): Font size of dataset labels inside bounding box.
            xaxis_title (str): Label for the x-axis.
            xaxis_title_font_size (int): Font size for x-axis title.
            yaxis_title (str): Label for the y-axis.
            yaxis_title_font_size (int): Font size for y-axis title.
            plot_title (str): Title of the plot.
            plot_title_font_size (int): Font size of the plot title.
            legend (list[str] or None): Custom legend labels. Defaults to ["Up", "Down"].
            legend_fontsize (int): Font size for legend text.
            subfolder (str): Subfolder under `./plot/` to save output files.
            filename (str): Base name for saved files.
            figure_description (str): Description of the figure for logging.

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Side Effects:
            - Saves the plot as PDF and PNG (trimmed to remove white borders).
            - Logs parameters and output using `data_manager`.

        Notes:
            - Up counts are displayed above each dataset, down counts below.
            - Dataset label is centered inside a rounded box (`bbox_facecolor`).
            - Jittering is applied in the x-direction to avoid overplotting.
        """
        def jitter(center, size, width=0.6):
            return np.random.uniform(center - width / 2.3, center + width / 2.3, size=size)
        
        def generate_centers(num_datasets, start=0.9, interval=0.65):
            return [start + i * interval for i in range(num_datasets)]
    
        if rectangle_colors is None:
            rectangle_colors = ['orange'] * len(datasets)
    
        fig, ax = plt.subplots(figsize=(6, 5))
    
        centers = generate_centers(len(datasets))
        colors = []  
        x_all, y_all = [], []  
        up_counts, down_counts = [], []  
    
        for i, (data, label) in enumerate(zip(datasets, labels)):
            center = centers[i]
            data_up = np.log2(data[data['fc'] > fc_threshold]['fc'])
            data_down = np.log2(data[data['fc'] < fc_threshold]['fc'])
            x_up = jitter(center, len(data_up))
            x_down = jitter(center, len(data_down))
            x_all.extend(x_up)
            y_all.extend(data_up)
            colors.extend([up_color] * len(data_up))
            x_all.extend(x_down)
            y_all.extend(data_down)
            colors.extend([down_color] * len(data_down))
            
            ax.add_patch(plt.Rectangle((center - rectangle_width / 2, -3.5), rectangle_width, 7,
                                       color=rectangle_colors[i], alpha=0.1, zorder=0))
            
            ax.text(center, 3.2, str(len(data_up)), ha='center', fontsize=10, fontweight='bold')
            ax.text(center, -3.2, str(len(data_down)), ha='center', fontsize=10, fontweight='bold')
            
            ax.text(center, 0, label, ha='center', va='center',
                    fontsize=bbox_textsize, fontweight='bold', color='black',
                    bbox=dict(facecolor=bbox_facecolor, edgecolor='black',
                              boxstyle='round,pad=0.3', lw=1.5))
    
        ax.scatter(x_all, y_all, c=colors, alpha=0.7, edgecolor=scatter_edgecolor, s=scatter_size)
    
        legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=up_color, markersize=8, label='Up'),
                          plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=down_color, markersize=8, label='Down')]
        ax.legend(handles=legend_handles, loc='best', frameon=False)
    
        ax.set_xlim(centers[0] - 0.5, centers[-1] + 0.5)
        ax.set_ylim(-3.5, 3.5)
        if not show_xaxis:
            ax.set_xticks([])
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
        ax.tick_params(axis='y', width=spine_width, labelsize=ytick_labelsize)
    
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(spine_width)
        ax.spines['bottom'].set_linewidth(spine_width)
   
        ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
        ax.set_ylabel(yaxis_title,fontproperties = 'Arial', size = yaxis_title_font_size) 
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)

        plt.tight_layout()
        # plt.show()
    
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        plt.savefig(os.path.join(output_dir, f"{filename}_up_down_scatter_plot_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_up_down_scatter_plot_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        params = locals()
        params.pop('self')
        self.data_manager.log_params('StrucGAP_DataVisualization', 'up_down_scatter', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'up_down_scatter',f"{labels}_up_down_scatter_plot.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_up_down_scatter_plot_{timestamp}.png"), 
                'legend': figure_description}
    
    def tree(self, data, node1_column, node2_column, screen_index_value, screen_index = None, colors = None, 
             rootname = 'pathway',
             layout = 'radial',
             symbol = 'emptyCircle',
             symbol_size = 15,
             collapse_interval = 2,
             legend_labels = None,  # 新增参数，用于标注简写的全名
             tree_pos_top = '30%',
             tree_pos_right = '50%',
             label_show = True,
             label_color = 'black',
             label_font_weight = 'normal',
             label_font_size = 10,
             annotation_pos_left = '60%',
             annotation_pos_top = '0%',
             annotation_font_size = 10,
             plot_title = None,
             subfolder='plot',
             filename = None,
             figure_description = 'Tree plot',
             ):
        """
        Create a hierarchical tree plot from pathway or enrichment results.

        This function visualizes terms (e.g., pathways) as nodes and expands them 
        into gene members. Each term is abbreviated with a short index (e.g., 0, 1, 2...) 
        while a legend annotation on the side shows the mapping from abbreviations 
        to full names.

        Args:
            data (pd.DataFrame): Input table containing term and gene information.
                                Must include a `Term` column and a `Genes` column 
                                (semicolon-separated gene list).
            node1_column (str): Placeholder, reserved for hierarchical data inputs (not used directly).
            node2_column (str): Placeholder, reserved for hierarchical data inputs (not used directly).
            screen_index_value (float): Threshold for filtering terms when `screen_index` is provided.
            screen_index (str, optional): Column name used to filter terms (e.g., `p.adjust`).
            colors (list[str], optional): Custom node colors (currently unused).
            rootname (str): Name of the root node of the tree. Default "pathway".
            layout (str): Tree layout, one of {"radial", "orthogonal"}. Default "radial".
            symbol (str): Node shape. Options: {"emptyCircle", "circle", "rect", 
                            "roundRect", "triangle", "diamond", "pin", "arrow", "none"}.
            symbol_size (int): Size of node symbols. Default 15.
            collapse_interval (int): Interval for collapsing tree levels.
            legend_labels (dict, optional): Mapping of short term IDs to full names.
                                            Defaults to auto-generated from `Term`.
            tree_pos_top (str): Vertical position of tree plot (e.g., "30%").
            tree_pos_right (str): Horizontal position of tree plot (e.g., "50%").
            label_show (bool): Whether to display node labels.
            label_color (str): Label font color.
            label_font_weight (str): Label font weight, e.g., {"normal", "bold"}.
            label_font_size (int): Label font size.
            annotation_pos_left (str): Left offset for legend annotation text box.
            annotation_pos_top (str): Top offset for legend annotation text box.
            annotation_font_size (int): Font size for legend annotation.
            plot_title (str, optional): Title of the tree plot.
            subfolder (str): Subfolder under "./plot/" to save results.
            filename (str): Base filename for saved outputs.
            figure_description (str): Description of the figure for logging.

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Side Effects:
            - Saves tree visualization as HTML, SVG, PDF, and PNG.
            - Trims white borders in the final PNG.
            - Logs parameters and outputs via `data_manager`.

        Notes:
            - Terms are automatically abbreviated to numeric IDs for clarity.
            - The right-hand side annotation shows mapping between short IDs and full term names.
            - Genes under each term are displayed as leaf nodes.
        """
        if layout not in ['radial', 'orthogonal']:
            layout = 'radial'
            
        if symbol not in ['emptyCircle', 'circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none']:
            symbol = 'emptyCircle'
        
        if screen_index:
            df = pd.DataFrame(data[data[screen_index] < screen_index_value]).copy()  # 创建副本
        if screen_index is None:
            df = pd.DataFrame(data).copy()  # 创建副本

        # 去除括号内的内容
        df['Term'] = df['Term'].apply(lambda x: re.sub(r'\s*\(.*\)', '', x))
        
        # 用 df.index 代替生成字母进行映射
        term_mapping = {term: str(idx) for idx, term in zip(df.index, df['Term'])}
        # 将 df['Term'] 按照 term_mapping 进行映射，生成简化后的 'Short_Term' 列
        df['Short_Term'] = df['Term'].map(term_mapping)

        # 定义函数生成嵌套字典
        def generate_dict(df, first_name):
            result = {'name': first_name, 'children': []}
            for _, row in df.iterrows():
                term = row['Short_Term']  # 使用简短的名称
                genes = row['Genes'].split(';')
                genes = [gene for gene in genes if gene]  # 去除空值
                gene_children = [{'name': gene, 'value': 1} for gene in genes]
                result['children'].append({'name': term, 'children': gene_children})
            return result
        
        # 生成嵌套字典
        pathway_dict = generate_dict(df, rootname)
        
        # 图例标注
        if legend_labels is None:
            legend_labels = {v: k for k, v in term_mapping.items()}  # 默认映射
        
        # 构造图例标注的文字内容
        legend_text = "\n".join([f"{abbr}: {full}" for abbr, full in legend_labels.items()])
        
        c = (
            Tree(init_opts=opts.InitOpts(
                renderer=RenderType.SVG,bg_color='#fff'
            ))
            .add("", 
                 [pathway_dict], 
                 layout = layout,
                 symbol = symbol,
                 collapse_interval = collapse_interval,
                 pos_top = tree_pos_top,
                 pos_right = tree_pos_right,
                 label_opts = opts.LabelOpts(
                    is_show = label_show,
                    color = label_color,
                    font_family = 'Arial',
                    font_weight = label_font_weight,
                    font_size = label_font_size,
                 ),
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(
                    title = plot_title,
                    pos_left="3%",  # 标题居中
                    pos_top="0%"       # 标题位置靠上
                ),
                graphic_opts=[
                    opts.GraphicGroup(
                        graphic_item=opts.GraphicItem(
                            # 控制整体的位置
                            left = annotation_pos_left,  # 调整文本框位置到右侧空白处
                            top = annotation_pos_top,
                        ),
                        children=[
                            # 添加文字
                            opts.GraphicText(
                                graphic_item=opts.GraphicItem(
                                    left="center",
                                    top="middle",
                                    z=100,
                                ),
                                graphic_textstyle_opts=opts.GraphicTextStyleOpts(
                                    text=legend_text,  # 图例标注内容
                                    font_size = annotation_font_size,
                                    graphic_basicstyle_opts=opts.GraphicBasicStyleOpts(
                                        fill="#333"
                                    )
                                )
                            )
                        ]
                    )
                ],
            )
            # .render("tree_layout.html")
        )
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 文件名保存和渲染
        file_name = os.path.join(output_dir, f"{filename}_tree_{timestamp}.html")
        c.render(file_name)
        svg_file = os.path.join(output_dir, f"{filename}_tree_{timestamp}.svg")
        make_snapshot(snapshot, file_name, svg_file)
        
        # 将 SVG 转换为 PDF
        pdf_file = os.path.join(output_dir, f"{filename}_tree_{timestamp}.pdf")
        drawing = svg2rlg(svg_file)  # 读取 SVG 文件
        renderPDF.drawToFile(drawing, pdf_file)  # 保存为 PDF 文件
        
        png_file = os.path.join(output_dir, f"{filename}_tree_{timestamp}.png")
        pdf_document = fitz.open(pdf_file)
        page = pdf_document.load_page(0)
        pix = page.get_pixmap(matrix=fitz.Matrix(8,8))
        pix.save(png_file)
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'tree', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'tree',f"{node1_column}_tree.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def venn_diagram(self, *data_names, 
                     colors = 'Tropic',
                     legend_loc='best', 
                     font_name='Arial',
                     subfolder='plot',
                     plot_title = None,
                     plot_title_font_size = 20,
                     number_fontsize = 20,
                     legend = None,
                     legend_fontsize = 12,
                     filename = None,
                     figure_description = 'Venn diagram',
                     ):
        """
        Draw a Venn diagram from multiple datasets.

        This function creates a Venn diagram visualization for 2–5 sets, 
        using either colormaps or user-provided color lists. 
        It supports automatic or custom legends and saves the plot in multiple formats.

        Args:
            *data_names (list of str or pd.Series): Datasets or variable names to visualize. 
                Each dataset should be a pandas Series or iterable of elements.
            colors (str or list, optional): 
                - If a colormap name (e.g., 'viridis', 'coolwarm'), colors are sampled automatically.
                - If a list of color strings (e.g., ['red', 'blue']), they are applied directly.
                Default: "Tropic".
            legend_loc (str): Location of the legend. One of:
                {"best","upper right","upper left","lower left","lower right",
                "right","center left","center right","lower center","upper center","center"}.
            font_name (str): Font family used for labels (default: "Arial").
            subfolder (str): Subfolder under `./plot/` to save outputs.
            plot_title (str, optional): Title for the Venn diagram.
            plot_title_font_size (int): Font size of the title.
            number_fontsize (int): Font size for numbers inside diagram regions.
            legend (list[str] or str, optional): Custom legend labels. 
                If None, dataframe variable names are used.
            legend_fontsize (int): Font size of the legend labels.
            filename (str, optional): Base filename for saving files.
            figure_description (str): Description for logging purposes.

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Side Effects:
            - Saves the Venn diagram as PDF and PNG files.
            - Logs parameters and outputs via `data_manager`.

        Notes:
            - Supports 2–5 sets (depending on backend `matplotlib-venn` / `venn` package).
            - If dataset variable names are not accessible, defaults to "unknown".
            - Uses `venn` package for visualization, which provides flexible styling.

        Example:
            >>> set1 = pd.Series(['A', 'B', 'C', 'D'])
            >>> set2 = pd.Series(['C', 'D', 'E', 'F'])
            >>> set3 = pd.Series(['A', 'E', 'G'])
            >>> viz.venn_diagram(set1, set2, set3, colors=['red','blue','green'],
            ...                  plot_title="3-way Venn", filename="example")
        """
        print('You can choose colors in following range: ', plt.colormaps())
        print("\n or enter a color list like ['red', 'blue']")
        
        if legend_loc not in ['best', 'upper right', 'upper left','lower left','lower right',
                              'right','center left','center right','lower center','upper center','center']:
            legend_loc = 'best'
    
        plt.rcParams['font.family'] = font_name
        dataframes = [
            eval(name) if isinstance(name, str) else name
            for name in data_names
        ]
        dataframe_names = []
        name_count = {}
        for name in data_names:
            if '[' in name:  # 保留原逻辑
                base_name = name.split('.')[-1].split('[')[0]
                if base_name in name_count:
                    name_count[base_name] += 1
                    unique_name = f"{base_name}_{name_count[base_name]}"
                else:
                    name_count[base_name] = 1
                    unique_name = base_name
                dataframe_names.append(unique_name)
            else:  # 修改只在变量处理部分
                # 获取变量的名称
                var_name = [key for key, value in globals().items() if value is name]
                if var_name:
                    base_name = var_name[0]  # 提取变量名
                    if base_name in name_count:
                        name_count[base_name] += 1
                        unique_name = f"{base_name}_{name_count[base_name]}"
                    else:
                        name_count[base_name] = 1
                        unique_name = base_name
                    dataframe_names.append(unique_name)
                else:
                    dataframe_names.append("unknown")  # 如果变量名未知，标记为未知

        sets = {name: set(series) for name, series in zip(dataframe_names, dataframes)}
        n_sets = len(sets)
        if not (2 <= n_sets <= 6):
            raise ValueError(f"Venn only supports 2–6 sets, got {n_sets}")
        
        if isinstance(colors, (list, tuple)):
            # 你这个 venn 版本不支持 colors=，那就只能退一步：
            # 把颜色列表做成一个 ListedColormap，再走 cmap=
            from matplotlib.colors import ListedColormap
            cmap_obj = ListedColormap(list(colors))
        else:
            cmap_obj = plt.get_cmap(colors)  # 关键：返回 colormap 对象，不是字符串
        
        venn_diagram = venn(sets, cmap=cmap_obj, fontsize=number_fontsize)
        #
        handles = [plt.Line2D([0], [0], color=venn_diagram.patches[i].get_facecolor(), lw=4) 
                   for i in range(len(dataframe_names))]
        
        if legend is not None:
            if isinstance(legend, str):
                legend = [legend]
        
        if legend is None:
            plt.legend(handles, dataframe_names, loc=legend_loc, fontsize=legend_fontsize)
        else:
            plt.legend(handles, legend, loc=legend_loc, fontsize=legend_fontsize)
        
        ax = plt.gca()       
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_venn_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_venn_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight', pad_inches=0)
        plt.close('all')
        
        # png_file = png_file
        # image = Image.open(png_file)
        # trimmed_image = self.trim_white_border(image)
        # trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'venn', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'venn',f"{len(data_names)}_data_venn.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_venn_{timestamp}.png"), 
                'legend': figure_description}
    
    def upset_plot(self, *data_names, 
                   colors = 'navyblue',
                   font_name='Arial', 
                   legend_loc='best', 
                   legend_fontsize=12,
                   subfolder='plot',
                   filename = None,
                   figure_description = 'Upset plot',
                   ):
        """
        Draw an UpSet plot for multiple datasets.

        This function visualizes set intersections across multiple datasets 
        using an UpSet plot (a scalable alternative to Venn diagrams). 
        It supports 2 or more sets and outputs publication-ready figures.

        Args:
            *data_names (list of str or pd.Series): 
                Datasets (Series, list, or set) or variable names to include in the UpSet plot.
            colors (str, optional): Color for intersection bars. 
                Accepts matplotlib color names or hex codes. Default: 'navyblue'.
            font_name (str, optional): Font family for all text. Default: 'Arial'.
            legend_loc (str, optional): Legend location (not heavily used in UpSet plots). 
                One of {"best","upper right","upper left","lower left","lower right",
                "right","center left","center right","lower center","upper center","center"}.
            legend_fontsize (int, optional): Font size of the legend labels. Default: 12.
            subfolder (str, optional): Subfolder under `./plot/` to save outputs. Default: "plot".
            filename (str, optional): Base filename for saving files. Default: None.
            figure_description (str, optional): Description for logging. Default: "Upset plot".

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Side Effects:
            - Saves the UpSet plot as both PDF and PNG files.
            - Logs parameters and outputs via `data_manager`.

        Notes:
            - Works best when there are more than 3 sets (where Venn diagrams become unreadable).
            - Uses the `UpSetPlot` library for visualization.
            - `show_counts` and `show_percentages` are enabled by default.

        Example:
            >>> set1 = pd.Series(['A','B','C','D'])
            >>> set2 = pd.Series(['C','D','E'])
            >>> set3 = pd.Series(['B','E','F','G'])
            >>> viz.upset_plot(set1, set2, set3, colors='darkred', filename="example_upset")
        """
        plt.rcParams['font.family'] = font_name
        dataframes = [
                eval(name) if isinstance(name, str) else name
                for name in data_names
            ]
        dataframe_names = []
        name_count = {}
        for name in data_names:
            if '[' in name:  # 保留原逻辑
                base_name = name.split('.')[-1].split('[')[0]
                if base_name in name_count:
                    name_count[base_name] += 1
                    unique_name = f"{base_name}_{name_count[base_name]}"
                else:
                    name_count[base_name] = 1
                    unique_name = base_name
                dataframe_names.append(unique_name)
            else:  # 修改只在变量处理部分
                # 获取变量的名称
                var_name = [key for key, value in globals().items() if value is name]
                if var_name:
                    base_name = var_name[0]  # 提取变量名
                    if base_name in name_count:
                        name_count[base_name] += 1
                        unique_name = f"{base_name}_{name_count[base_name]}"
                    else:
                        name_count[base_name] = 1
                        unique_name = base_name
                    dataframe_names.append(unique_name)
                else:
                    dataframe_names.append("unknown")  # 如果变量名未知，标记为未知
        sets = {name: set(series) for name, series in zip(dataframe_names, dataframes)}
        all_elements = set.union(*sets.values())
        intersection_data = {}
        for bool_combination in product([False, True], repeat=len(sets)):
            index = tuple(bool_combination)
            included_sets = [name for name, include in zip(sets.keys(), bool_combination) if include]
            excluded_sets = [name for name, include in zip(sets.keys(), bool_combination) if not include]
            if included_sets:
                current_intersection = set.intersection(*[sets[set_name] for set_name in included_sets])
            else:
                current_intersection = all_elements
            if excluded_sets:
                current_intersection -= set.union(*[sets[set_name] for set_name in excluded_sets])
            intersection_data[index] = len(current_intersection)
        index_names = [name for name in sets.keys()]
        multi_index = pd.MultiIndex.from_tuples(intersection_data.keys(), names=index_names)
        result_series = pd.Series(intersection_data, index=multi_index)
        UpSet(result_series, show_counts="%d", show_percentages="{:.2%}", facecolor=colors).plot()
        # plt.show()
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_upset_plot_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_upset_plot_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'upset_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'upset_plot',f"{len(data_names)}_upset_plot.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_upset_plot_{timestamp}.png"), 
                'legend': figure_description}
    
    def violin_plot(self, data, item_column, item_name, group1_columns, group2_columns, 
                    p_data, p_column, colors=None, 
                    showmeans = False,
                    showmedians = True,
                    violin_width = 0.5,
                    xaxis_label_color='black', 
                    xaxis_label_font_weight='normal', 
                    xaxis_label_font_size=9, 
                    yaxis_label_color='black', 
                    yaxis_label_font_weight='normal', 
                    yaxis_label_font_size=14, 
                    axis_line_width=2, 
                    axis_tick_length=5, 
                    axis_tick_width=2,  
                    add_group_lines=True,  
                    group_line_color='black', 
                    group_line_style='-', 
                    group_line_width=2,  
                    add_x_grid=False,  
                    add_y_grid=True,  
                    x_grid_line_color='grey', 
                    x_grid_line_style='-', 
                    x_grid_line_width=0.5,
                    y_grid_line_color='grey', 
                    y_grid_line_style='-', 
                    y_grid_line_width=0.5,
                    p_text_offset=0.05,  
                    scatter_alpha=0.7,  # 新增参数：散点透明度
                    scatter_color='black',  # 新增参数：散点颜色
                    line_color='black',  # 新增参数：小提琴线和中位数颜色
                    xaxis_title = None,
                    xaxis_title_font_size = 20,
                    yaxis_title = None,
                    yaxis_title_font_size = 20,
                    plot_title = None,
                    plot_title_font_size = 20,
                    legend = None,
                    legend_fontsize = 12,
                    subfolder='plot',
                    filename = None,
                    figure_description = 'Violin plot',
                    ):  
        """
        Draw a violin plot with two groups per item, scatter overlay, and p-value annotations.

        Args:
            data (pd.DataFrame): Input data.
            item_column (str): Column name specifying the items (e.g., features).
            item_name (list of str): List of item names to plot.
            group1_columns (list of str): Column names for group 1 values.
            group2_columns (list of str): Column names for group 2 values.
            p_data (pd.DataFrame): DataFrame containing p-values for each item.
            p_column (str): Column name in `p_data` with p-values.
            colors (list, optional): Two colors for group1 and group2 violins.
            showmeans (bool, optional): Whether to show the mean line. Default False.
            showmedians (bool, optional): Whether to show the median line. Default True.
            violin_width (float, optional): Width of violins. Default 0.5.
            scatter_alpha (float, optional): Transparency of scatter points. Default 0.7.
            scatter_color (str, optional): Color of scatter points. Default 'black'.
            line_color (str, optional): Color of violin borders and median line. Default 'black'.
            add_group_lines (bool, optional): Whether to draw separator lines between groups. Default True.
            p_text_offset (float, optional): Vertical offset for p-value text. Default 0.05.
            xaxis_title, yaxis_title, plot_title (str, optional): Axis and plot titles.
            filename (str, optional): Base name for output files.
            subfolder (str, optional): Folder under ./plot/ to save outputs. Default 'plot'.
            figure_description (str, optional): Figure description for logging. Default 'Violin plot'.

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Notes:
            - Each item has two violins (group1 and group2).
            - Scatter points are overlaid for raw values.
            - p-values from `p_data` are annotated above violins.
            - Outputs are saved as PDF and PNG files, with white border automatically trimmed.
        """
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]

        # 提取group1和group2的数据
        v1 = []
        v2 = []
        for i in item_name:
            v1.append(list(data[data[item_column] == i][group1_columns].values[0]))
            v2.append(list(data[data[item_column] == i][group2_columns].values[0]))

        # 将数据合并
        all_data = []
        for d1, d2 in zip(v1, v2):
            all_data.append(d1)
            all_data.append(d2)
        
        # 创建图形
        plt.figure(figsize=(10, 6))
        parts = plt.violinplot(all_data, 
                               showmeans = showmeans, 
                               showmedians = showmedians,
                               widths = violin_width
                               )

        # 修改小提琴图中的线条颜色为黑色
        for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
            vp = parts[partname]
            vp.set_edgecolor(line_color)
            vp.set_linewidth(1.5)
        
        # 设置小提琴颜色
        for idx, pc in enumerate(parts['bodies']):
            if idx % 2 == 0:
                pc.set_facecolor(colors[0])  # 奇数index为第一组颜色
            else:
                pc.set_facecolor(colors[1])  # 偶数index为第二组颜色
            pc.set_edgecolor('black')
            pc.set_alpha(0.8)
        
        # 添加散点图
        for i in range(len(item_name)):
            # 分别绘制两组的散点
            group1_data = v1[i]
            group2_data = v2[i]
            
            # A组散点
            plt.scatter([i * 2 + 1] * len(group1_data), group1_data, 
                        color=scatter_color, alpha=scatter_alpha)
            
            # B组散点
            plt.scatter([i * 2 + 2] * len(group2_data), group2_data, 
                        color=scatter_color, alpha=scatter_alpha)

        # 设置 x 轴标签，确保每个 item 只在中间显示
        group_labels = [f"{name}" for name in item_name]
        mid_positions = [i * 2 + 1.5 for i in range(len(item_name))]  
        plt.xticks(mid_positions, group_labels, rotation=0, fontsize=xaxis_label_font_size, color=xaxis_label_color)
        plt.xlabel('', fontsize=xaxis_label_font_size, color=xaxis_label_color, weight=xaxis_label_font_weight)
        
        # 设置 y 轴标签
        plt.yticks(fontsize=yaxis_label_font_size, color=yaxis_label_color, weight=yaxis_label_font_weight)

        # 只添加分组之间的垂直分割线
        if add_group_lines:
            for pos in mid_positions[:-1]:  
                plt.axvline(x=pos + 1, color=group_line_color, linestyle=group_line_style, linewidth=group_line_width)
        
        # 显示 p 值
        p_list = []
        for i in item_name:
            p_value = round(p_data[p_data.index == i][p_column].values[0], 5)
            p_list.append(p_value)

        # 将 p 值注释向上偏移
        for idx, p_value in enumerate(p_list):
            plt.text(mid_positions[idx], plt.ylim()[1] + p_text_offset, f"p = {p_value}", ha='center', fontsize=xaxis_label_font_size)

        # 设置坐标轴的线条粗细及标签的长度和粗细
        ax = plt.gca()
        ax.spines['top'].set_linewidth(axis_line_width)
        ax.spines['right'].set_linewidth(axis_line_width)
        ax.spines['left'].set_linewidth(axis_line_width)
        ax.spines['bottom'].set_linewidth(axis_line_width)
        
        ax.set_xlabel(xaxis_title,fontproperties = 'Arial', size = xaxis_title_font_size)
        ax.set_ylabel(yaxis_title,fontproperties = 'Arial', size = yaxis_title_font_size) 
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        if legend:
            plt.legend(legend, loc='best', fontsize=legend_fontsize)
        
        ax.tick_params(axis='both', which='major', length=axis_tick_length, width=axis_tick_width)
        
        # 添加x轴网格线
        if add_x_grid:
            ax.xaxis.grid(True, color=x_grid_line_color, linestyle=x_grid_line_style, linewidth=x_grid_line_width)

        # 添加y轴网格线
        if add_y_grid:
            ax.yaxis.grid(True, color=y_grid_line_color, linestyle=y_grid_line_style, linewidth=y_grid_line_width)

        # 显示图像
        plt.tight_layout()
        # plt.show()
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_violin_plot_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_violin_plot_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'violin_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'violin_plot',f"{item_column}_item_{p_column}_violin_plot.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_violin_plot_{timestamp}.png"), 
                'legend': figure_description}
    
    def volcano_plot(self, data, fc_column, p_column, colors=None, 
                    fc = 2,
                    p_value = 0.05,
                    color_na = 'gray',
                    opacity_na = 0.5,
                    size_na = 20,
                    color_middle = 'orange',
                    opacity_middle = 1,
                    size_middle = 20,
                    color_right = '#870303',
                    opacity_right = 1,
                    size_right = 40,
                    color_left = '#006BAC',
                    opacity_left = 1,
                    size_left = 40,
                    p_line_color = 'black',
                    p_line_linestyle = '--',
                    p_line_linewidth = 2,
                    fc_line_right_color = '#870303',
                    fc_line_right_linestyle = '--',
                    fc_line_right_linewidth = 2,
                    fc_line_left_color = '#006BAC',
                    fc_line_left_linestyle = '--',
                    fc_line_left_linewidth = 2,
                    plot_title = None,
                    plot_title_font_size = 20,
                    subfolder='plot',
                    filename = None,
                    figure_description = 'Volcano plot',
                    ):  
        """
        Generate a volcano plot for fold-change and p-value visualization.

        Args:
            data (pd.DataFrame): Input data.
            fc_column (str): Column name for fold change values.
            p_column (str): Column name for p-values.
            fc (float, optional): Base for log transformation of fold change. Default = 2.
            p_value (float, optional): Threshold for significance. Default = 0.05.
            color_na, opacity_na, size_na: Appearance of non-significant points.
            color_middle, opacity_middle, size_middle: Appearance of significant but moderate FC points.
            color_right, opacity_right, size_right: Appearance of significantly up-regulated points.
            color_left, opacity_left, size_left: Appearance of significantly down-regulated points.
            p_line_color, p_line_linestyle, p_line_linewidth: Style of p-value threshold line.
            fc_line_right_color, fc_line_right_linestyle, fc_line_right_linewidth: Style of right FC threshold line.
            fc_line_left_color, fc_line_left_linestyle, fc_line_left_linewidth: Style of left FC threshold line.
            plot_title (str, optional): Plot title.
            plot_title_font_size (int, optional): Font size of plot title. Default = 20.
            subfolder (str, optional): Output subfolder under ./plot/. Default = 'plot'.
            filename (str, optional): Base filename for output files.
            figure_description (str, optional): Description used for logging. Default = 'Volcano plot'.

        Returns:
            dict:
                - "file_path" (str): Path to saved PNG file.
                - "legend" (str): Figure description.

        Notes:
            - Fold change values are log-transformed to base `fc`.
            - P-values are transformed to -log10(p).
            - Data points are categorized into four groups:
                * Non-significant (p >= threshold)
                * Significant with moderate FC
                * Significantly up-regulated
                * Significantly down-regulated
            - Threshold lines for p-value and FC are drawn automatically.
            - Output includes PDF and PNG with trimmed borders.
        """
        data_g = data[[fc_column,p_column]]
        data_g[fc_column] = np.log(data_g[fc_column])/np.log(fc)
        data_g['MinusLog10PValue'] = -np.log10(data_g[p_column])
        
        # 动态确定横坐标范围
        fc_max = data_g[fc_column].abs().max()
        x_limit = round(fc_max+1)  # 四舍五入到最近的整数
        # x_limit = max(x_limit, 5)  # 确保横坐标范围至少为-5到5
        # 设置Fold Change的阈值和P-value的阈值
        fc_threshold_up = 1  # 上调FC阈值
        fc_threshold_down = -1  # 下调FC阈值
        pvalue_threshold = p_value  # P-value阈值
        
        plt.figure(figsize=(10, 6))
        plt.scatter(data_g[(data_g[p_column] >= pvalue_threshold)][fc_column], 
                    data_g[(data_g[p_column] >= pvalue_threshold)]['MinusLog10PValue'], 
                    color = color_na, alpha = opacity_na, s = size_na)
        plt.scatter(data_g[(data_g[p_column] < pvalue_threshold) & (data_g[fc_column] <= fc_threshold_up) & (data_g[fc_column] >= fc_threshold_down)][fc_column],
                    data_g[(data_g[p_column] < pvalue_threshold) & (data_g[fc_column] <= fc_threshold_up) & (data_g[fc_column] >= fc_threshold_down)]['MinusLog10PValue'],
                    color = color_middle, alpha = opacity_middle, s = size_middle,
                    label = f'Significant, moderate change (n={len(data_g[(data_g[p_column] < pvalue_threshold) & (data_g[fc_column] <= fc_threshold_up) & (data_g[fc_column] >= fc_threshold_down)])})')
        plt.scatter(data_g[(data_g[fc_column] > fc_threshold_up) & (data_g[p_column] < pvalue_threshold)][fc_column],
                    data_g[(data_g[fc_column] > fc_threshold_up) & (data_g[p_column] < pvalue_threshold)]['MinusLog10PValue'],
                    color = color_right, alpha = opacity_right, s = size_right,
                    label = f'Up-regulated (n={len(data_g[(data_g[fc_column] > fc_threshold_up) & (data_g[p_column] < pvalue_threshold)])})')
        plt.scatter(data_g[(data_g[fc_column] < fc_threshold_down) & (data_g[p_column] < pvalue_threshold)][fc_column],
                    data_g[(data_g[fc_column] < fc_threshold_down) & (data_g[p_column] < pvalue_threshold)]['MinusLog10PValue'],
                    color = color_left, alpha = opacity_left, s = size_left,
                    label = f'Down-regulated (n={len(data_g[(data_g[fc_column] < fc_threshold_down) & (data_g[p_column] < pvalue_threshold)])})',)
        plt.xlim([-x_limit, x_limit])
        plt.axhline(-np.log10(pvalue_threshold), color = p_line_color, linestyle = p_line_linestyle, linewidth = p_line_linewidth)
        plt.axvline(fc_threshold_up, color = fc_line_right_color, linestyle = fc_line_right_linestyle, linewidth = fc_line_right_linewidth)
        plt.axvline(fc_threshold_down, color = fc_line_left_color, linestyle = fc_line_left_linestyle, linewidth = fc_line_left_linewidth)
        
        plt.tick_params(width=8,length=5,color='k')
        plt.legend(prop={'family' : 'Arial', 'size'   : 20},loc='best')
        ax = plt.gca()
        ax.spines['bottom'].set_linewidth('3')
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_linewidth('3')
        ax.spines['top'].set_color('k')
        ax.spines['left'].set_linewidth('3')
        ax.spines['left'].set_color('k')
        ax.spines['right'].set_linewidth('3')
        ax.spines['right'].set_color('k')
        plt.yticks(fontproperties = 'Arial', size = 25)
        plt.xticks(fontproperties = 'Arial', size = 25)
        plt.xlabel(f'Log{fc}(Fold Change)',fontproperties = 'Arial', size = 25)
        plt.ylabel('-log10(P-value)',fontproperties = 'Arial', size = 25)
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        # plt.show()
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_volcano_plot_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_volcano_plot_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'volcano_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'volcano_plot',f"{fc_column}_{p_column}_volcano_plot.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_volcano_plot_{timestamp}.png"), 
                'legend': figure_description}
    
    def dotplot_col(self, data, category, p_column, top, term, dot_color_column, dot_size_column,
                annotation_cmap = 'Pastel2', col_split=True, row_cluster=False,col_cluster=True,
                dot_cmap = 'RdYlGn',
                xaxis_font_size = 20,
                yaxis_font_size = 20,
                subfolder='plot',
                filename = None,
                figure_description = 'Dot plot',
                ):
        """
        Generate a clustered dot plot with column annotations.

        Args:
            data (pd.DataFrame): Input data containing enrichment results or gene set statistics.
            category (str): Column name representing dataset categories (e.g., conditions, groups).
            p_column (str): Column name containing p-values for sorting significance.
            top (int): Number of top terms per category to include (lowest p-values).
            term (str): Column name containing term names (e.g., pathways, GO terms).
            dot_color_column (str): Column name to map dot color (e.g., enrichment score).
            dot_size_column (str): Column name to map dot size (e.g., count or hits/total).
            annotation_cmap (str, optional): Colormap for top annotation (categories). Default = "Pastel2".
            col_split (bool, optional): Whether to split columns by dataset category. Default = True.
            row_cluster (bool, optional): Whether to cluster rows. Default = False.
            col_cluster (bool, optional): Whether to cluster columns. Default = True.
            dot_cmap (str, optional): Colormap for dot colors. Default = "RdYlGn".
            xaxis_font_size (int, optional): Font size for x-axis labels. Default = 20.
            yaxis_font_size (int, optional): Font size for y-axis labels. Default = 20.
            subfolder (str, optional): Output subfolder under ./plot/. Default = "plot".
            filename (str, optional): Base filename for saved output.
            figure_description (str, optional): Description string for logging. Default = "Dot plot".

        Returns:
            dict:
                - "file_path" (str): Path to saved PNG file.
                - "legend" (str): Figure description.

        Notes:
            - For each category, the top `top` terms with lowest p-values are selected.
            - `dot_size_column` values are automatically converted from "hits/total" format
            to integer counts (only numerator is used).
            - Column annotations show dataset categories with colors from `annotation_cmap`.
            - Supports optional clustering of rows/columns and splitting by category.
            - Output files include both PDF and PNG versions with trimmed borders.
        """
        gene_set_categories = data[category].unique()

        top10_list = []
        for i in gene_set_categories:
            df_category = data[data[category] == i]
            top10_category = df_category.sort_values(by=p_column, ascending=True).head(top)
            top10_list.append(top10_category)
        top10_combined_df = pd.concat(top10_list)
        data = top10_combined_df
        
        df_col = pd.DataFrame(index = gene_set_categories)
        df_col['Dataset'] = list(gene_set_categories)
        
        col_ha = HeatmapAnnotation(
                            Dataset=anno_simple(df_col.Dataset,cmap=annotation_cmap,
                                                add_text=True,legend=True,text_kws={'fontsize':xaxis_font_size},),
                            verbose=0,label_side='left',label_kws={'horizontalalignment':'right','fontsize':xaxis_font_size})

        matplotlib.use('Qt5Agg')  # 或者 'TkAgg'，取决于你系统的设置
        
        data[dot_size_column] = data[dot_size_column].apply(lambda x: int(x.split('/')[0]))
        
        if col_split:
            col_split=df_col.Dataset
        
        plt.figure(figsize=(3, 4.5))
        cm = DotClustermapPlotter(data=data, x=category, y=term, value=p_column, 
                                  c = dot_color_column,
                                  s = dot_size_column,
                                  cmap = dot_cmap,
                                  # hue='EnrichType', 
                                  row_cluster=row_cluster,col_cluster=col_cluster,
                                  # cmap={'Enrich':'RdYlGn_r','Depletion':'coolwarm_r'},
                                  # colors={'Enrich':'red','Depletion':'blue'},
                                  #marker={'Enrich':'^','Depletion':'v'},
                                  top_annotation=col_ha,
                                  # right_annotation=row_ha,
                                  # col_split=df_col.Dataset,row_split=df_row.Category, 
                                  col_split_gap=0.5,row_split_gap=1,
                                  show_rownames=True,show_colnames=False,row_dendrogram=False,
                                  verbose=1,legend_gap=7,spines=True, col_split=col_split,
                                  xticklabels_kws={'labelsize': xaxis_font_size},
                                  yticklabels_kws={'labelsize': yaxis_font_size}) 
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_dotplot_col_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_dotplot_col_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'dotplot_col', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'dotplot_col',f"{term}_dotplot_col.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_dotplot_col_{timestamp}.png"), 
                'legend': figure_description}
    
    def dotplot_row(self, data, category, p_column, top, term, dot_color_column, dot_size_column,
                annotation_cmap = 'Pastel2',row_split=True,row_cluster=True,col_cluster=False,
                dot_cmap = 'RdYlGn',
                xaxis_font_size = 20,
                yaxis_font_size = 20,
                subfolder='plot',
                filename = None,
                figure_description = 'Dot plot',
                ):
        """
        Generate a clustered dot plot with row annotations.

        Args:
            data (pd.DataFrame): Input data containing enrichment results or gene set statistics.
            category (str): Column name representing dataset categories (e.g., experimental groups).
            p_column (str): Column name containing p-values for sorting significance.
            top (int): Number of top terms per category to include (lowest p-values).
            term (str): Column name containing term names (e.g., pathways, GO terms).
            dot_color_column (str): Column name mapped to dot colors (e.g., enrichment score).
            dot_size_column (str): Column name mapped to dot sizes (e.g., "hits/total" format, converted to ratio).
            annotation_cmap (str, optional): Colormap for row annotations. Default = "Pastel2".
            row_split (bool, optional): Whether to split rows by category. Default = True.
            row_cluster (bool, optional): Whether to cluster rows. Default = True.
            col_cluster (bool, optional): Whether to cluster columns. Default = False.
            dot_cmap (str, optional): Colormap for dot colors. Default = "RdYlGn".
            xaxis_font_size (int, optional): Font size for x-axis labels. Default = 20.
            yaxis_font_size (int, optional): Font size for y-axis labels. Default = 20.
            subfolder (str, optional): Output subfolder under ./plot/. Default = "plot".
            filename (str, optional): Base filename for saved output.
            figure_description (str, optional): Description string for logging. Default = "Dot plot".

        Returns:
            dict:
                - "file_path" (str): Path to saved PNG file.
                - "legend" (str): Figure description.

        Notes:
            - For each category, the top `top` terms with lowest p-values are selected.
            - `dot_size_column` values are expected in "hits/total" format and are converted
            to the ratio (hits ÷ total).
            - Row annotations (`annotation_cmap`) indicate the category of each term.
            - Supports optional clustering of rows/columns and row splitting.
            - Output files include both PDF and PNG versions with trimmed borders.
        """
        print('You can choose the camp in following range: ', plt.colormaps())
        
        gene_set_categories = data[category].unique()

        top10_list = []
        for i in gene_set_categories:
            df_category = data[data[category] == i]
            top10_category = df_category.sort_values(by=p_column, ascending=True).head(top)
            top10_list.append(top10_category)
        top10_combined_df = pd.concat(top10_list)
        data = top10_combined_df
        
        df_col = data[[term, category]]  
        df_col.set_index(term, inplace=True, drop=True)
                
        row_ha = HeatmapAnnotation(
                                    Dataset=anno_simple(df_col[category],cmap=annotation_cmap,
                                                        add_text=False,legend=True,text_kws={'fontsize':yaxis_font_size}),
                                    verbose=0,axis=0,label_kws={'horizontalalignment':'right','fontsize':yaxis_font_size})
        
        matplotlib.use('Qt5Agg')  # 或者 'TkAgg'，取决于你系统的设置
        
        data[dot_size_column] = data[dot_size_column].apply(lambda x: int(x.split('/')[0]) / int(x.split('/')[1]))
        
        if row_split:
            row_split=df_row[row_category]
        
        plt.figure(figsize=(3, 4.5))
        cm = DotClustermapPlotter(data=data, x=category, y=term, value=p_column, 
                                  c = dot_color_column,
                                  s = dot_size_column,
                                  cmap = dot_cmap,
                                  # hue='EnrichType', 
                                  row_cluster=row_cluster,col_cluster=col_cluster,
                                  # cmap={'Enrich':'RdYlGn_r','Depletion':'coolwarm_r'},
                                  # colors={'Enrich':'red','Depletion':'blue'},
                                  #marker={'Enrich':'^','Depletion':'v'},
                                  # top_annotation=col_ha,
                                  right_annotation=row_ha,
                                  # col_split=df_col.Dataset,row_split=df_row.Category, 
                                  col_split_gap=0.5,row_split_gap=1,
                                  show_rownames=True,show_colnames=False,row_dendrogram=False,
                                  verbose=1,legend_gap=7,spines=True, row_split=row_split, 
                                  xticklabels_kws={'labelsize': xaxis_font_size},
                                  yticklabels_kws={'labelsize': yaxis_font_size},) 
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_dotplot_row_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_dotplot_row_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'dotplot_row', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'dotplot_row',f"{term}_dotplot_row.pdf")

        return {'file_path': os.path.join(output_dir, f"{filename}_dotplot_row_{timestamp}.png"), 
                'legend': figure_description}
    
    def dotplot_row_col(self, data, row_category, col_category, p_column, top, term, 
                        dot_color_column, 
                        dot_size_column,
                        annotate_rows=False, 
                        annotate_cols=False,
                        col_split=True, 
                        row_split=True,
                        annotation_cmap='Pastel2', 
                        dot_cmap='RdYlGn', 
                        xaxis_font_size = 20,
                        yaxis_font_size = 20,
                        subfolder='plot',
                        filename = None,
                        figure_description = 'Dot plot',
                        ):
        """
        Generate a two-dimensional dot plot (dotplot) with both row and column categories.

        This function creates a clustered dot plot where:
            - The x-axis represents datasets (from `col_category`).
            - The y-axis represents terms or gene sets (from `row_category`).
            - Dot size encodes values from `dot_size_column` (e.g., hits/total ratio).
            - Dot color encodes values from `dot_color_column` (e.g., enrichment score or fold change).
            - Rows/columns can optionally include category annotations and be split by categories.

        Args:
            data (pd.DataFrame): Input dataframe containing enrichment results or gene set statistics.
            row_category (str): Column name defining row categories (e.g., gene set group).
            col_category (str): Column name defining column categories (e.g., datasets).
            p_column (str): Column name containing p-values used for ranking.
            top (int): Number of top terms (lowest p-values) per category to include.
            term (str): Column name containing the term or pathway names.
            dot_color_column (str): Column name mapped to dot color scale.
            dot_size_column (str): Column name mapped to dot sizes. If formatted as "hits/total",
                                it will be converted into a ratio.
            annotate_rows (bool, optional): Whether to show row annotation for row categories. Default = False.
            annotate_cols (bool, optional): Whether to show column annotation for column categories. Default = False.
            col_split (bool, optional): Whether to split columns by category. Default = True.
            row_split (bool, optional): Whether to split rows by category. Default = True.
            annotation_cmap (str, optional): Colormap for annotations. Default = "Pastel2".
            dot_cmap (str, optional): Colormap for dot colors. Default = "RdYlGn".
            xaxis_font_size (int, optional): Font size of x-axis labels. Default = 20.
            yaxis_font_size (int, optional): Font size of y-axis labels. Default = 20.
            subfolder (str, optional): Output subfolder under `./plot/`. Default = "plot".
            filename (str, optional): Base filename for saved outputs. Default = None.
            figure_description (str, optional): Description string for logging. Default = "Dot plot".

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Notes:
            - The function selects the top `top` terms with the lowest p-values for each combination
            of `row_category` × `col_category`.
            - If `dot_size_column` contains values like "5/100", they are converted to ratios (e.g., 0.05).
            - Annotations (`annotate_rows` / `annotate_cols`) add extra categorical information
            alongside the heatmap axes.
            - Supports row/column splitting for improved visualization of grouped categories.
            - The output includes both PDF and PNG formats, with white borders trimmed.
        """
        data_set_categories = data[col_category].unique()
        gene_set_categories = data[row_category].unique()
    
        # 筛选每类中的前 top 项
        result_list = []
        for i in data_set_categories:  
            df_category = data[data[col_category] == i]
            top_list = []
            for j in gene_set_categories:  
                gf_category = df_category[df_category[row_category] == j]
                top_category = gf_category.sort_values(by=p_column, ascending=True).head(top)
                top_list.append(top_category)
            if top_list:  
                top_combined_df = pd.concat(top_list)
                result_list.append(top_combined_df)
        if result_list:  
            data = pd.concat(result_list)
        else:
            data = pd.DataFrame()  
    
        # 处理点大小
        if '/' in data[dot_size_column][0]:
            data[dot_size_column] = data[dot_size_column].apply(
                lambda x: int(x.split('/')[0]) / int(x.split('/')[1])
            )
        
        data['id'] = data[col_category]
        
        data = data.sort_values(by=row_category, ascending=False)
    
        # 创建行和列的标注
        col_ha = None
        if annotate_cols:
            df_col = data[['id', col_category]]
            df_col.set_index('id', inplace=True, drop=True)
            df_col = df_col[~df_col.index.duplicated()]
            col_ha = HeatmapAnnotation(
                Dataset=anno_simple(
                    df_col[col_category],
                    cmap=annotation_cmap,
                    add_text=True,
                    legend=True,
                    text_kws={'fontsize': xaxis_font_size}
                ),
                verbose=0,
                label_side='left',
                label_kws={'horizontalalignment': 'right','fontsize':xaxis_font_size}
            )
    
        row_ha = None
        if annotate_rows:
            df_row = data[[term, row_category]].copy()
            df_row.set_index(term, inplace=True, drop=True)
            df_row = df_row[~df_row.index.duplicated()]
            row_ha = HeatmapAnnotation(
                Dataset=anno_simple(
                    df_row[row_category],
                    cmap=annotation_cmap,
                    add_text=False,
                    legend=True,
                    text_kws={'fontsize': yaxis_font_size}
                ),
                verbose=0,
                axis=0,
                label_kws={'horizontalalignment': 'right','fontsize':yaxis_font_size}
            )
            
        if col_split:
            col_split=df_col[col_category]
        if row_split:
            row_split=df_row[row_category]
    
        # 绘制 Dotplot  data1[term]
        plt.figure(figsize=(3, 4.5))
        cm = DotClustermapPlotter(
            data=data,
            x='id',
            y=term,
            value=p_column,
            c=dot_color_column,
            s=dot_size_column,
            cmap=dot_cmap,
            row_cluster=row_cluster,
            col_cluster=col_cluster,
            top_annotation=col_ha,
            right_annotation=row_ha,
            col_split_gap=0.5,
            row_split_gap=1,
            show_rownames=True,
            show_colnames=False,
            row_dendrogram=False,
            verbose=1,
            legend_gap=7,
            spines=True,
            xticklabels_kws={'labelsize': xaxis_font_size},
            yticklabels_kws={'labelsize': yaxis_font_size},
            col_split=col_split,
            row_split=row_split, 
        )
    
        # 输出目录和文件
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
    
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        plt.savefig(
            os.path.join(output_dir, f"{filename}_dotplot_{timestamp}.pdf"),
            dpi=900,
            bbox_inches='tight'
        )
        png_file = os.path.join(output_dir, f"{filename}_dotplot_{timestamp}.png")
        plt.savefig(
            png_file,
            dpi=900,
            bbox_inches='tight'
        )
        plt.close('all')
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
    
        # 自动记录参数
        params = locals()
        params.pop('self')
        self.data_manager.log_params('StrucGAP_DataVisualization', 'dotplot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'dotplot', f"{term}_dotplot.pdf")
    
        return {'file_path': os.path.join(output_dir, f"{filename}_dotplot_{timestamp}.png"), 
                'legend': figure_description}
    
    def sankey_dotplot(self, data, category, p_column, top, term, dot_color_column, dot_size_column, genes,
                colors_term = None,
                color_node = None,
                annotation_cmap = 'mrybm',
                dot_cmap = 'rainbow',
                sankey_node_fontsize = 6,
                sankey_node_padding = 1.5,
                sankey_node_width = 50,
                dot_xaxis_fontsize = 20,
                dot_size_min = 10,
                dot_opacity = 1,
                axis_line_color = 'black',
                axis_line_width = 3,
                grid_line_color = 'grey',
                grid_line_width = 2,
                subfolder='plot',
                filename = None,
                figure_description = 'Sankey & dot plot',
                ):
        """
        Generate a combined Sankey diagram and dot plot for enriched terms and their associated genes.

        This function integrates two visualization styles:
            1. Sankey diagram: Shows relationships between genes (left nodes) and enriched terms (right nodes),
            where links represent gene-term associations.
            2. Dot plot: Highlights term-level statistics (e.g., significance, overlap ratio), where
            dot color encodes values from a specified column, and dot size encodes gene overlap counts.

        Args:
            data (pd.DataFrame): Input dataframe containing enrichment results.
            category (str): Column name specifying the dataset or category grouping terms.
            p_column (str): Column name containing p-values, used to rank top terms.
            top (int): Number of top terms (lowest p-values) to select per category.
            term (str): Column name containing term names (e.g., pathways).
            dot_color_column (str): Column mapped to dot colors (e.g., -log10(p-value)).
            dot_size_column (str): Column mapped to dot size. Values like "5/100" will be converted to ratios.
            genes (str): Column containing gene lists (semicolon-separated) for each term.

            colors_term (list, optional): List of colors to assign categories in the Sankey diagram.
            color_node (list, optional): Additional node color settings. Default = None.
            annotation_cmap (str, optional): Colormap for annotations. Default = "mrybm".
            dot_cmap (str, optional): Colormap for dot plot coloring. Default = "rainbow".

            sankey_node_fontsize (int, optional): Font size for Sankey node labels. Default = 6.
            sankey_node_padding (float, optional): Padding between Sankey nodes. Default = 1.5.
            sankey_node_width (int, optional): Node width in Sankey diagram. Default = 50.

            dot_xaxis_fontsize (int, optional): Font size of x-axis labels in dot plot. Default = 20.
            dot_size_min (int, optional): Minimum dot size in dot plot. Default = 10.
            dot_opacity (float, optional): Transparency of dots (0–1). Default = 1.

            axis_line_color (str, optional): Color of axis lines in dot plot. Default = "black".
            axis_line_width (float, optional): Line width for axes. Default = 3.
            grid_line_color (str, optional): Gridline color. Default = "grey".
            grid_line_width (float, optional): Gridline width. Default = 2.

            subfolder (str, optional): Output subfolder under `./plot/`. Default = "plot".
            filename (str, optional): Base filename for saved outputs. Default = None.
            figure_description (str, optional): Description string for logging. Default = "Sankey & dot plot".

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Notes:
            - The Sankey diagram links individual genes (sources) to enriched terms (targets).
            - Dot plot encodes overlap ratio on the x-axis, with dot size proportional to hit counts.
            - Only the top `top` terms (lowest p-values) per category are included.
            - Colors for terms are derived from `colors_term`, while gene nodes are assigned colors automatically.
            - Outputs are saved as both `.pdf` and `.png` with trimmed white borders.
        """
        # print('You can choose the cmap in following range: ', plt.colormaps())
        
        if colors_term is None:
            colors_term = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]

        gene_set_categories = data[category].unique()
        top10_list = []
        for i in gene_set_categories:
            df_category = data[data[category] == i]
            top10_category = df_category.sort_values(by=p_column, ascending=True).head(top)
            top10_list.append(top10_category)
        top10_combined_df = pd.concat(top10_list)
        data = top10_combined_df
        
        def calculate_ratio(overlap_str):
            num, denom = overlap_str.split('/')  
            return int(num) / int(denom)  
        
        # 准备数据
        df = data[[category, dot_size_column, term, p_column, genes]].copy()  # 使用 .copy() 明确表示你在创建一个副本
        df.loc[:, 'Overlap_Ratio'] = df[dot_size_column].apply(calculate_ratio)  # 使用 .loc[] 确保赋值行为明确
        df.loc[:, genes] = df[genes].str.split(';')  # 使用 .loc[] 来确保不会产生副本赋值的警告
        df[term] = [x[0] for x in df[term].str.split(' \(')]
        df[term] = pd.Categorical(df[term], list(df[term].values))
        df.sort_values([term], inplace = True)
        df = df.reset_index(drop=True)
        
        # 创建桑基图数据   type(data[dot_size_column][0])
        sources = []
        targets = []
        values = []
        
        for idx, row in df.iterrows():
            for gene in row[genes]:
                sources.append(gene)
                targets.append(row[term])
                values.append(1)
        
        # 根据Gene_set为右侧Term节点着色
        unique_categories = df[category].unique()
        # 手动输入颜色
        input_colors = colors_term
        # 生成 color_mapping 字典
        color_mapping = dict(zip(unique_categories, input_colors))

        term_colors = [color_mapping[row[category]] for _, row in df.iterrows()]
        
        # 创建桑基图节点标签（基因 + Term）
        node_labels = []
        for target in list(df[term]):
            node_labels.append(target)
            for source in list(df[df[term] == target][genes])[0]:
                if source not in node_labels:
                        node_labels.append(source)
        
        # 这样node_labels中将按照sources和targets的原始顺序排列
        
        # 获取基因的唯一列表
        unique_genes = list(set(sources))
        # 使用 cmap 自动为基因节点生成颜色，并转换为 hex 格式
        cmap = plt.get_cmap('random100')  # 可以选择不同的cmap，例如 'Set3', 'tab10', 'Paired'
        gene_colors = [mcolors.rgb2hex(cmap(i / len(unique_genes))) for i in range(len(unique_genes))]
        # 创建基因颜色映射
        gene_color_mapping = dict(zip(unique_genes, gene_colors))
        # 为节点分配颜色，基因节点使用 gene_color_mapping 中的颜色
        node_colors = [
            gene_color_mapping[lbl] if lbl in sources else color_mapping[df[df[term] == lbl][category].values[0]] 
            for lbl in node_labels
        ]

        # 将源和目标映射到节点索引
        source_indices = [node_labels.index(gene) for gene in sources]
        target_indices = [node_labels.index(term) for term in targets]
        
        # 创建桑基图
        sankey_trace = go.Sankey(
            arrangement = 'perpendicular',
            textfont=dict(family='Arial', size=sankey_node_fontsize),
            node=dict(
                pad = sankey_node_padding,
                thickness = sankey_node_width,
                label=node_labels,
                color=node_colors,
                y=[0.1, 0.2, 0.5, 0.7, 0.9],
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values
            )
        )
        
        # 创建气泡图
        df = df.sort_index(ascending=False)
        # 计算每个term对应的节点高度，根据基因数量
        df['Gene_count'] = df[genes].apply(lambda x: len(x))  # 计算每个term的基因数量
        # 基于基因数量的比例来调整气泡图中的y坐标
        first_half = df['Gene_count'].iloc[0] / 2
        last_half = df['Gene_count'].iloc[-1] / 2
        middle_sum = df['Gene_count'].iloc[1:-1].sum()
        max_height = first_half + middle_sum + last_half
        df['y_scaled'] = df['Gene_count'] / max_height
        df['y_processed'] = 0.0
        for i in range(len(df)-1, -1, -1):  # 从 len(df)-1 递减到 0
            if i == len(df)-1:
                df.at[i, 'y_processed'] = 0  
            elif i == 0:
                df.at[i, 'y_processed'] = 1 
            else:
                df.at[i, 'y_processed'] = (df.at[len(df)-1, 'y_scaled'] / 2) + df['y_scaled'][1:len(df)-i-1].sum() + (df.at[i, 'y_scaled'] / 2)
        
        dot_color = ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
                     'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
                     'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
                     'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
                     'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
                     'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
                     'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
                     'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
                     'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
                     'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
                     'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
                     'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
                     'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
                     'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr',
                     'ylorrd']
        
        # print('You can choose dot color in the following list: ', dot_color)
        
        # 定义最大、中间、最小气泡的值
        max_size_value = df[dot_size_column].apply(lambda x: int(x.split('/')[0])).max()
        min_size_value = df[dot_size_column].apply(lambda x: int(x.split('/')[0])).min()
        median_size_value = df[dot_size_column].apply(lambda x: int(x.split('/')[0])).median()

        bubble_trace = go.Scatter(
            x=df['Overlap_Ratio'],
            y=df['y_processed'],
            mode='markers',
            marker=dict(
                size=df[dot_size_column].apply(lambda x: int(x.split('/')[0])),
                sizemin = dot_size_min,
                sizemode='diameter',
                color=df[dot_color_column],
                opacity = dot_opacity,
                showscale=True,
                colorscale=dot_cmap,
                colorbar=dict(
                    title='P-value',
                    len=0.5,  # 减小颜色条长度
                    x=1.05
                )
            ),
            text=df[term]
        )
                
        # 创建子图布局
        fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4], specs=[[{"type": "sankey"}, {"type": "scatter"}]])
        
        # 添加桑基图和气泡图
        fig.add_trace(sankey_trace, row=1, col=1)
        fig.add_trace(bubble_trace, row=1, col=2)
        
        # 自定义气泡大小图例
        fig.add_trace(go.Scatter(
            x=[df['Overlap_Ratio'].max()*2], y=[0.95],
            mode="markers+text",
            marker=dict(size=[max_size_value], color="black"),
            text=[f'{max_size_value}'],
            textposition="middle right",
            showlegend=False
        ))
        
        fig.add_trace(go.Scatter(
            x=[df['Overlap_Ratio'].max()*2], y=[0.85],
            mode="markers+text",
            marker=dict(size=[median_size_value], color="black"),
            text=[f'{median_size_value}'],
            textposition="middle right",
            showlegend=False
        ))
        
        fig.add_trace(go.Scatter(
            x=[df['Overlap_Ratio'].max()*2], y=[0.75],
            mode="markers+text",
            marker=dict(size=[min_size_value], color="black"),
            text=[f'{min_size_value}'],
            textposition="middle right",
            showlegend=False
        ))

        # 隐藏气泡图的y轴标签
        fig.update_yaxes(showticklabels=False, row=1, col=2)
        fig.update_xaxes(tickfont=dict(size=dot_xaxis_fontsize))
        
        # 设置布局，应用背景颜色
        fig.update_layout(
            # 只修改气泡图的背景（假设气泡图在第2列）
            xaxis=dict(
                showgrid = True,
                gridcolor = grid_line_color,
                gridwidth = grid_line_width,
                zeroline = False,
                zerolinecolor = 'blue',
                zerolinewidth = 0,
                linecolor = axis_line_color,
                linewidth = axis_line_width,
                mirror = True,
            ),
            yaxis=dict(
                showgrid = True,
                gridcolor = grid_line_color,
                gridwidth = grid_line_width,
                zeroline = False,
                zerolinecolor = 'blue',
                zerolinewidth = 0,
                linecolor = axis_line_color,
                linewidth = axis_line_width,
                mirror = True,
            ),
            template='none',
            # plot_bgcolor='rgba(0,0,0,0)'  
        )
        fig.update_layout(width=1280, height=720)
        # 显示图表
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # 保存到 plot 文件夹中的路径
        output_path = os.path.join(output_dir, "{filename}_sankey_dotplot_{timestamp}.pdf")
        pio.write_image(fig, output_path, engine="orca")
        png_file = os.path.join(output_dir, "{filename}_sankey_dotplot_{timestamp}.png")
        pio.write_image(fig, png_file, engine="orca")
        # fig.show()
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # # 显示图像
        # plt.savefig(f"{term}_dotplot_row.pdf",dpi=900, bbox_inches='tight')
        # plt.close('all')

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'sankey_dotplot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'sankey_dotplot',"sankey_dotplot.pdf")

        return {'file_path': png_file, 'legend': figure_description}
    
    def dimension_reduction(self, data, data_columns, sample_group, filter_data, p_column, p_value, fc, 
                            method='tsne',
                            dimension_number = 3,
                            colors = None,
                            marker_shape = ['o', 'v'],
                            marker_size = [100, 100],
                            marker_alpha = [1, 1],
                            ellipse_alpha = 0.3,
                            ellipse_linestyle = '-',
                            ellipse_linewidth = 4,
                            axis_line_width = 3,
                            axis_line_color = 'black',
                            subfolder = 'plot',
                            random_state = None,
                            show_labels = False,
                            plot_title = None,
                            plot_title_font_size = 20,
                            filename = None,
                            figure_description = 'Dimension reduction plot',
                            ):
        """
        Perform dimensionality reduction (PCA, t-SNE, or UMAP) and visualize samples in 2D with group annotations.

        This function reduces high-dimensional data (e.g., expression or intensity values) into
        two dimensions using PCA, t-SNE, or UMAP. Samples are plotted as scatter points, optionally
        grouped by experimental conditions with distinct markers, colors, and confidence ellipses.

        Args:
            data (pd.DataFrame): Input data matrix with features as rows and samples as columns.
            data_columns (list[str]): List of columns in `data` to include in dimensionality reduction.
            sample_group (pd.DataFrame): DataFrame mapping sample IDs (index) to group labels 
                (must contain a column named 'group').
            filter_data (pd.DataFrame or None): Optional filtering DataFrame (e.g., differential expression results). 
                If provided, only samples/genes passing thresholds are kept.
            p_column (str): Column name in `filter_data` containing p-values for filtering.
            p_value (float): P-value threshold; only entries below this threshold are retained.
            fc (float): Fold-change threshold; entries outside [1/fc, fc] are retained.
            method (str, optional): Dimensionality reduction method: 'pca', 'tsne', or 'umap'. Default = 'tsne'.
            dimension_number (int, optional): Number of dimensions to reduce to (usually 2 or 3). Default = 3.
            colors (list, optional): List of colors for groups. Default uses preset color palette.
            marker_shape (list, optional): List of matplotlib marker styles for groups. Default = ['o', 'v'].
            marker_size (list, optional): Marker sizes for groups. Default = [100, 100].
            marker_alpha (list, optional): Transparency of markers (0–1). Default = [1, 1].
            ellipse_alpha (float, optional): Transparency of group confidence ellipses. Default = 0.3.
            ellipse_linestyle (str, optional): Line style for ellipses ('-', '--', '-.', ':'). Default = '-'.
            ellipse_linewidth (int, optional): Line width of ellipses. Default = 4.
            axis_line_width (int, optional): Width of axis lines. Default = 3.
            axis_line_color (str, optional): Color of axis lines. Default = "black".
            subfolder (str, optional): Subdirectory for saving results. Default = 'plot'.
            random_state (int or None, optional): Random seed for reproducibility. Default = None.
            show_labels (bool, optional): Whether to annotate sample points with labels. Default = False.
            plot_title (str, optional): Title of the plot. Default = None.
            plot_title_font_size (int, optional): Font size of plot title. Default = 20.
            filename (str, optional): Base filename for saving results. Default = None.
            figure_description (str, optional): Description string for logging. Default = "Dimension reduction plot".

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG plot file.
                - "legend" (str): Figure description for downstream reporting.

        Notes:
            - PCA is deterministic; t-SNE and UMAP are stochastic, so results vary unless `random_state` is fixed.
            - If `filter_data` is provided, only significant and fold-change-filtered features are used.
            - Groups in `sample_group` are assigned colors and markers in plotting.
            - Confidence ellipses are drawn using group covariance to approximate data spread.
            - Outputs include `.pdf` and `.png` files with trimmed borders.
        """
        if colors is None:
            colors = ["#bd221f", "#099eda", "#fee301", "#abb7bd", "#A07EBA", "#293a6e", 
                      "#d6c223", "#6ebb53", "#d75d73", "#e63b29", "#e0592b", "#58b7b3"]
            
        if ellipse_linestyle not in ['-', '--', '-.', ':', '']:
            ellipse_linestyle = '-'
        
        for i in marker_shape:
            if i not in ['.', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']:
                marker_shape = ['o']*len(sample_group['group'].unique())

        sample_group.index = sample_group.index.astype(str)
        
        # 检查降维方法
        if method not in ['pca', 'tsne', 'umap']:
            method = 'tsne'
            print('Your input method is not available, the method was set to tsne.')
        
        # 筛选数据
        data = data[data_columns]
        if filter_data is not None:
            data = data.loc[data.index.isin(filter_data[filter_data[p_column] < p_value].index)]
            if fc is not None:
                data = data.loc[data.index.isin(filter_data[(filter_data['fc'] < 1/fc)|(filter_data['fc'] > fc)].index)]
        
        # 选择降维方法
        if method == 'pca':
            reducer = PCA(n_components=dimension_number, random_state = random_state)
            reduced_data = reducer.fit_transform(data.T)
            x_label, y_label = "PCA1", "PCA2"
        
        elif method == 'tsne':
            reducer = TSNE(n_components=dimension_number, random_state = random_state)
            reduced_data = reducer.fit_transform(data.T)
            x_label, y_label = "tSNE1", "tSNE2"
        
        elif method == 'umap':
            reducer = umap.UMAP(n_components=dimension_number, random_state = random_state)
            reduced_data = reducer.fit_transform(data.T)
            x_label, y_label = "UMAP1", "UMAP2"
        
        # 创建降维后的 DataFrame
        reduced_df = pd.DataFrame(reduced_data[:, :2], columns=[x_label, y_label], index=data_columns)
        reduced_df = pd.concat([reduced_df, sample_group], axis=1)
        
        # 可视化
        plt.figure(figsize=(8, 6))
        groups = reduced_df['group'].unique()
        
        for i, group in enumerate(groups):
            subset = reduced_df[reduced_df['group'] == group]
            plt.scatter(subset[x_label], subset[y_label], label=group,
                        marker=marker_shape[i],
                        color=colors[i % len(colors)],
                        s=marker_size[i],
                        alpha=marker_alpha[i])
            
            # 在绘制散点后，添加标签（根据 show_labels 参数）
            if show_labels:
                for _, row in subset.iterrows():
                    plt.text(row[x_label], row[y_label], row.name,  # 使用样本索引作为标签
                             fontsize=10, ha='center', va='center', color='black')

            # 绘制椭圆包围区域
            x_mean, y_mean = subset[x_label].mean(), subset[y_label].mean()
            cov = np.cov(subset[[x_label, y_label]].T)
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
            width, height = 4 * np.sqrt(eigenvalues)
            ellipse = Ellipse((x_mean, y_mean), width, height, angle, 
                              edgecolor=colors[i % len(colors)], 
                              facecolor=colors[i % len(colors)], 
                              alpha = ellipse_alpha,
                              linestyle = ellipse_linestyle,
                              linewidth = ellipse_linewidth,
                              )
            plt.gca().add_patch(ellipse)
            
        plt.tick_params(width=axis_line_width,length=5,color=axis_line_color)  
        ax = plt.gca()
        ax.spines['bottom'].set_linewidth(axis_line_width)
        ax.spines['bottom'].set_color(axis_line_color)
        ax.spines['top'].set_linewidth(axis_line_width)
        ax.spines['top'].set_color(axis_line_color)
        ax.spines['left'].set_linewidth(axis_line_width)
        ax.spines['left'].set_color(axis_line_color)
        ax.spines['right'].set_linewidth(axis_line_width)
        ax.spines['right'].set_color(axis_line_color)

        plt.xlabel(x_label,fontproperties = 'Arial', size = 25)
        plt.ylabel(y_label,fontproperties = 'Arial', size = 25)
        plt.yticks(fontproperties = 'Arial', size = 25)
        plt.xticks(fontproperties = 'Arial', size = 25)
        plt.legend(prop={'family' : 'Arial', 'size'   : 20},loc='best')
        ax.set_title(plot_title,fontproperties = 'Arial', size = plot_title_font_size)
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        plt.savefig(os.path.join(output_dir, f'{filename}_{method}_{timestamp}.pdf'),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f'{filename}_{method}_{timestamp}.png')
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        # plt.show()
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'dimension_reduction', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'dimension_reduction',f'{data_columns}_{method}.pdf')

        return {'file_path': os.path.join(output_dir, f'{filename}_{method}_{timestamp}.png'), 
                'legend': figure_description}
    
    def network_plot(self, data, term, gene, p_column, p,
                     species=10090,
                     pathway_size=1500,
                     pathway_font_size=0,
                     protein_size=800,
                     protein_font_size=0,
                     pathway_protein_weight=2,
                     subfolder='plot',
                     filename = None,
                     figure_description = 'Network plot',
                     ):
        """
        Construct and visualize a pathway–protein interaction network using STRING database.

        This function builds a bipartite network of pathways and proteins, integrating STRING 
        protein–protein interaction (PPI) data to visualize both pathway–protein associations 
        and intra-protein interactions. Nodes representing pathways and proteins are 
        drawn with customizable sizes, and multi-pathway proteins are shown as pie-chart 
        nodes with multiple colors.

        Args:
            data (pd.DataFrame): Input table containing enrichment results or associations, 
                with at least columns for terms, genes, and p-values.
            term (str): Column name in `data` containing pathway/term identifiers.
            gene (str): Column name in `data` containing associated genes (semicolon-separated).
            p_column (str): Column name containing p-values.
            p (float): P-value threshold; only terms below this cutoff are used.
            species (int, optional): STRING taxonomy ID for the species (default = 10090 for mouse).
            pathway_size (int, optional): Node size for pathway nodes. Default = 1500.
            pathway_font_size (int, optional): Font size for pathway labels. Default = 0 (hidden).
            protein_size (int, optional): Node size for protein nodes. Default = 800.
            protein_font_size (int, optional): Font size for protein labels. Default = 0 (hidden).
            pathway_protein_weight (float, optional): Edge width for pathway–protein connections. Default = 2.
            subfolder (str, optional): Subdirectory under `./plot/` to save results. Default = 'plot'.
            filename (str, optional): Base filename for output files. Default = None.
            figure_description (str, optional): Text description of the figure. Default = "Network plot".

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file of the network plot.
                - "legend" (str): Figure description (for downstream reporting).

        Notes:
            - This function uses the STRING database API (version 12.0) to retrieve protein–protein 
            interactions (PPIs). An internet connection is required.
            - Pathways are represented as colored nodes; proteins connected to multiple pathways 
            are drawn as pie-chart nodes with multiple colors.
            - Both PDF and PNG versions of the plot are saved in the specified subfolder.
            - Edges:
                * Pathway–protein edges: fixed width (`pathway_protein_weight`).
                * Protein–protein edges: weighted by STRING combined score.
            - Outputs are trimmed to remove white margins.

        Example:
            >>> vis.network_plot(data=df, 
                                term='Pathway', 
                                gene='Genes', 
                                p_column='pvalue', 
                                p=0.05, 
                                species=9606,
                                filename='KEGG_network')
            # Returns a dict with file path and description, and saves PDF/PNG plots.
        """
        data = data[data[p_column] < p]
        data = data[[term, gene]]
        genes_split = data[gene].str.split(';', expand=True)
        data = pd.concat([data[[term]], genes_split], axis=1)
        data.set_index(term, inplace=True, drop=True)
        
        def get_string_interactions(proteins, species=species, required_score=400, caller_identity="my_tool"):
            identifiers = "%0D".join(proteins)
            # 使用特定版本的STRING URL
            versioned_url = "https://version-12-0.string-db.org/api/json/network"
            url = f"{versioned_url}?identifiers={identifiers}&species={species}&required_score={required_score}&caller_identity={caller_identity}"
            
            try:
                response = requests.get(url)
                response.raise_for_status()  # 检查HTTP请求是否成功
                data = response.json()
                
                # 检查返回的数据是否为空或格式异常
                if not isinstance(data, list):
                    print(f"Warning: Unexpected data format received for proteins {proteins}")
                    return []  # 返回空列表，表示没有有效的交互信息
                
                return data
        
            except requests.exceptions.RequestException as e:
                print(f"Network error when fetching interactions for {proteins}: {e}")
                return []  # 返回空列表，表示请求失败
        
            except ValueError as e:
                print(f"JSON decoding error when fetching interactions for {proteins}: {e}")
                return []  # 返回空列表，表示解析失败
            
        def map_to_string_ids(proteins, species=species, caller_identity="my_tool"):
            # 构建STRING的映射API URL
            url = "https://version-12-0.string-db.org/api/json/get_string_ids"
            
            # 请求的参数，包括蛋白质列表、物种和调用者身份
            params = {
                "identifiers": "%0D".join(proteins),  # 用 "%0D" 连接所有标识符，符合API要求
                "species": species,                   # 物种编号，如9606代表人类
                "caller_identity": caller_identity    # 标识调用者身份
            }
            
            try:
                # 发送请求
                response = requests.get(url, params=params)
                response.raise_for_status()
                
                # 将返回的结果解析为JSON格式
                data = response.json()
                
                # 检查数据结构并提取STRING ID
                mapped_ids = {}
                for entry in data:
                    original_id = entry["queryItem"]       # 原始输入的蛋白质名称
                    string_id = entry["stringId"]          # 对应的STRING ID
                    mapped_ids[original_id] = string_id    # 将STRING ID保存到字典
                
                return mapped_ids
        
            except requests.exceptions.RequestException as e:
                print(f"Network error during mapping: {e}")
                return {}
            except ValueError as e:
                print(f"JSON decoding error during mapping: {e}")
                return {}

        string_dict = {}
        for i in data.index:  # i='Cell-Matrix Adhesion (GO:0007160)'
            proteins = data.loc[i]
            proteins = list(proteins.dropna())
            proteins = map_to_string_ids(proteins)
            interactions = get_string_interactions(list(proteins.values()))
            # 创建反向字典，将STRING ID作为键，原始蛋白质名称作为值
            reverse_proteins = {v: k for k, v in proteins.items()}
            # 构建DataFrame并还原原始蛋白质名称
            # 构建DataFrame并还原原始蛋白质名称
            string = pd.DataFrame([
                {
                    'node1': reverse_proteins[interaction["stringId_A"]],
                    'node2': reverse_proteins[interaction["stringId_B"]],
                    'combined_score': interaction["score"]
                }
                for interaction in interactions
                if interaction["stringId_A"] in reverse_proteins and interaction["stringId_B"] in reverse_proteins
            ])

            string_dict[i] = string

        G = nx.Graph()

        pathway_list = list(string_dict.keys())
        cmap = plt.get_cmap('random100')
        colors = [mcolors.rgb2hex(cmap(i / len(pathway_list))) for i in range(len(pathway_list))]
        pathway_colors = dict(zip(pathway_list, colors))
        node_colors = {}
        node_sizes = {}
        node_labels = {}
        edge_colors = []
        edge_widths = []

        gene_pathways = {}

        for idx, pname in enumerate(pathway_list):
            G.add_node(pname, type='pathway')
            node_colors[pname] = [colors[idx]]
            node_sizes[pname] = pathway_size
            node_labels[pname] = pname

            unique_proteins = [protein for protein in list(data.loc[pname]) if protein is not None]

            for gene in unique_proteins:
                G.add_node(gene, type='protein')
                if gene not in gene_pathways:
                    gene_pathways[gene] = []
                gene_pathways[gene].append(pname)
                node_sizes[gene] = protein_size
                node_labels[gene] = gene

                G.add_edge(pname, gene, color=colors[idx], weight=pathway_protein_weight)
                edge_colors.append(colors[idx])
                edge_widths.append(1)

            if not string_dict[pname].empty:
                for i, j, k in zip(string_dict[pname]['node1'], string_dict[pname]['node2'], string_dict[pname]['combined_score']):
                    G.add_edge(i, j, color=colors[idx], weight=k * 5)
                    edge_colors.append(colors[idx])
                    edge_widths.append(k * 5)
                    
        # 为每个基因节点设置颜色
        for gene in gene_pathways:
            pathways = gene_pathways[gene]
            if len(pathways) == 1:
                node_colors[gene] = [pathway_colors[pathways[0]]]
            else:
                node_colors[gene] = [pathway_colors[p] for p in pathways]

        # 计算多色节点的半径，使其面积与单色节点相等
        def calculate_total_size(protein_size, scaling_factor=0.0015):
            return np.sqrt(protein_size / np.pi) * scaling_factor
        
        total_radius = calculate_total_size(protein_size)

        def draw_piechart_node(ax, x, y, sizes, colors, total_radius):
            start = 0.
            for frac, color in zip(sizes, colors):
                end = start + frac
                wedges = Wedge(center=(x, y), r=total_radius, theta1=start * 360, theta2=end * 360,
                               facecolor=color, lw=0.5)
                ax.add_patch(wedges)
                start = end

        pos = nx.spring_layout(G, k=0.3, iterations=50, seed=42)

        plt.figure(figsize=(20, 20))
        ax = plt.gca()

        # 绘制边
        for u, v in G.edges():
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)],
                                   width=G[u][v]['weight'],
                                   edge_color=G[u][v]['color'],
                                   ax=ax)

        # 绘制节点
        for node in G.nodes():
            x, y = pos[node]
            sizes = node_sizes[node]
            colors = node_colors[node]
            if len(colors) == 1:
                # 单一颜色节点
                ax.scatter(x, y, s=sizes, c=colors[0], zorder=3)
            else:
                # 多颜色节点，绘制饼图
                fractions = [1 / len(colors)] * len(colors)
                draw_piechart_node(ax, x, y, fractions, colors, total_radius)  # 使用计算出的半径

        # 绘制节点标签
        for node in G.nodes():
            x, y = pos[node]
            ax.text(x, y, node_labels[node], fontsize=protein_font_size, ha='center', va='center', zorder=4)

        # 创建图例句柄
        legend_handles = []
        for pname, color in pathway_colors.items():
            # 为图例创建标记
            legend_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=pname,
                                             markersize=10, markerfacecolor=color))
        
        # 添加图例
        plt.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.05, 0.5), title="Pathways", fontsize=10, title_fontsize=12)

        ax.set_aspect('equal')
        plt.axis('off')
        plt.tight_layout()
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        plt.savefig(os.path.join(output_dir, '{filename}_network_plot_{timestamp}.pdf'), format='pdf')
        png_file = os.path.join(output_dir, '{filename}_network_plot_{timestamp}.png')
        plt.savefig(png_file, format='png')
        plt.close()
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)
        
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'network_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'network_plot','network_plot.pdf')

        return {'file_path': os.path.join(output_dir, '{filename}_network_plot_{timestamp}.png'), 
                'legend': figure_description}
    
    def glyconetwork(self, network_dict, regulation_type, 
                                  feature, 
                                  fc_scale=300, 
                                  overlap_scale=15, 
                                  center_scale=3000,
                                  center_node_fontsize = 20,
                                  other_node_fontsize = 18,
                                  subfolder='plot',
                                  filename = None,
                                  figure_description = 'Network',):
        """
        Visualize a glycan-centered regulatory network.

        This function constructs and plots a directed network for a selected glycan feature, 
        integrating enzymes, glycan-binding proteins (GBPs), and downstream pathways. 
        The glycan center serves as the hub, with upstream regulators (enzymes/GBPs) 
        connected by correlation-weighted edges, and downstream pathways connected 
        according to enrichment overlaps. Node sizes are scaled by fold change or overlap, 
        and edges are colored according to regulatory effect.

        Args:
            network_dict (dict): Nested dictionary containing network data in the format:
                network_dict[regulation_type][feature] = {
                    'glycan_center': {'id': ..., 'fc': ...},
                    'enzymes': [{'id': ..., 'type': 'enzyme', 'fc': ..., 'edge_weight': ...}, ...],
                    'gbps': [{'id': ..., 'type': 'GBP', 'fc': ..., 'edge_weight': ...}, ...],
                    'pathway': [{'term': ..., 'overlap': ...}, ...]
                }
            regulation_type (str): The regulation category (e.g., "up", "down").
            feature (str): The glycan feature to visualize as the network center.
            fc_scale (float, optional): Scaling factor for enzyme/GBP node sizes by fold change. Default = 300.
            overlap_scale (float, optional): Scaling factor for downstream pathway node sizes by overlap. Default = 15.
            center_scale (float, optional): Fixed size for the glycan center node. Default = 3000.
            center_node_fontsize (int, optional): Font size for the center node label. Default = 20.
            other_node_fontsize (int, optional): Font size for other node labels. Default = 18.
            subfolder (str, optional): Subdirectory under `./plot/` for saving outputs. Default = 'plot'.
            filename (str, optional): Base filename for outputs. Default = None.
            figure_description (str, optional): Description of the figure for logging. Default = "Network".

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Visualization details:
            - **Nodes**:
                * Center (glycan feature): Purple (#5C0A98), fixed size.
                * Enzymes: Orange (#E99E75), size proportional to |fold change| × `fc_scale`.
                * GBPs: Brown (#776483), size proportional to |fold change| × `fc_scale`.
                * Pathways: Gray, size proportional to overlap × `overlap_scale`.
            - **Edges**:
                * Enzyme/GBP → glycan: Blue (#3558AE) for positive correlation, 
                Pink (#B64074) for negative correlation.
                * Glycan → pathway: Green (#1F5F5B), weighted by overlap.
            - Legends for node types and edge types are automatically added.

        Output:
            - PDF and PNG plots are saved in `./plot/{subfolder}/` with time-stamped filenames.
            - White margins are trimmed from the PNG.
            - Parameters and outputs are logged via `self.data_manager`.

        Example:
            >>> vis.glyconetwork(network_dict=my_network,
                                regulation_type='up',
                                feature='HexNAc',
                                filename='glyco_up_HexNAc')
            # Saves PDF/PNG network plot and returns file path + description.
        """
        net = network_dict[regulation_type][feature]
        G = nx.DiGraph()

        # 1. 添加中心节点
        glycan_center_id = net['glycan_center']['id']
        glycan_center_fc = net['glycan_center']['fc']
        G.add_node(glycan_center_id, label=feature, type='center', size=center_scale)

        # 2. 添加上游节点及边
        for node in net['enzymes'] + net['gbps']:
            G.add_node(node['id'], label=node['id'], type=node['type'], size=abs(node['fc']*fc_scale))
            G.add_edge(node['id'], glycan_center_id, color='#3558AE' if node['edge_weight'] > 0 else '#B64074',
                       weight=max(0.5, abs(node['edge_weight'])*8), label=f"{node['edge_weight']:.2f}")

        # 3. 添加下游节点及边
        for node in net['pathway']:
            term = node['term']
            G.add_node(term, label=term, type='pathway', size=(node['overlap'] or 1)*overlap_scale)
            G.add_edge(glycan_center_id, term, color='#1F5F5B', weight=max(0.5, (node['overlap'] or 1)/2))

        # ==== 定义手工布局 ====
        pos = {}
        # 假设中心节点为(0, 0)，上游在左，下游在右
        center_pos = np.array([0, 0])
        r1 = 5
        
        # 上游节点在左（π到2π）
        upstream = [n for n in G.nodes if G.nodes[n]['type'] in ('enzyme', 'GBP')]
        # 上游节点均匀分布在左侧
        n_up = len(upstream)
        for i, n in enumerate(upstream):
            # θ: 135°~225°（π*3/4~π*5/4），均匀分布，保证全部在左侧
            theta = (3*np.pi/4) + (np.pi/2) * (i+1)/(n_up+1)
            pos[n] = np.array([
                center_pos[0] + r1 * np.cos(theta),
                center_pos[1] + r1 * np.sin(theta)
            ])        
        # 下游节点在右（0到π）
        downstream = [n for n in G.nodes if G.nodes[n]['type']=='pathway']
        n_down = len(downstream)
        for i, n in enumerate(downstream):
            # θ: -45°~+45°（-π/4~+π/4），全部在右侧
            theta = (-np.pi/4) + (np.pi/2) * (i+1)/(n_down+1)
            pos[n] = np.array([
                center_pos[0] + r1 * np.cos(theta),
                center_pos[1] + r1 * np.sin(theta)
            ])
            
        # 中心节点
        pos[glycan_center_id] = center_pos

        # ==== 绘图 ====
        node_sizes = [G.nodes[n]['size'] for n in G.nodes]
        node_colors = [
            '#5C0A98' if G.nodes[n]['type']=='center' else
            '#E99E75' if G.nodes[n]['type']=='enzyme' else
            '#776483' if G.nodes[n]['type']=='GBP' else
            'gray'
            for n in G.nodes
        ]
        edge_colors = [G[u][v]['color'] for u,v in G.edges]
        edge_widths = [G[u][v]['weight'] for u,v in G.edges]

        plt.figure(figsize=(13,10))
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, alpha=0.95)
        
        def get_radius(node_size, k=1600):
            return np.sqrt(node_size/np.pi)/k

        ax = plt.gca()
        r_dict = {n: get_radius(G.nodes[n]['size']) for n in G.nodes}
        
        for u, v in G.edges:
            x1, y1 = pos[u]
            x2, y2 = pos[v]
            r_src = r_dict[u]
            r_dst = r_dict[v]
            dx, dy = x2-x1, y2-y1
            dist = np.hypot(dx, dy)
            dx, dy = dx/dist, dy/dist
        
            # 判断是上游->中心还是中心->下游
            if v == glycan_center_id:
                # 上游→中心，缩短终点（到中心节点外圈）
                x1_new = x1 + dx * r_src * 40
                y1_new = y1 + dy * r_src * 40
                x2_new = x2 - dx * r_dst * 30  # 稍微大于1，防止被盖住
                y2_new = y2 - dy * r_dst * 30
            elif u == glycan_center_id:
                # 中心→下游，缩短起点（从中心节点外圈出发）
                x1_new = x1 + dx * r_src * 35
                y1_new = y1 + dy * r_src * 35
                x2_new = x2 - dx * r_dst * 35
                y2_new = y2 - dy * r_dst * 35
            else:
                # 其它边（默认）
                x1_new = x1 + dx * r_src
                y1_new = y1 + dy * r_src
                x2_new = x2 - dx * r_dst
                y2_new = y2 - dy * r_dst
        
            color = G[u][v]['color']
            width = G[u][v]['weight']
            arrow = FancyArrowPatch(
                (x1_new, y1_new), (x2_new, y2_new),
                arrowstyle='-|>',
                color=color,
                linewidth=width,
                mutation_scale=28,
                alpha=0.9,
                zorder=1
            )
            ax.add_patch(arrow)
        
        for n in G.nodes:
            x, y = pos[n]
            label = G.nodes[n]['label']
            if n == glycan_center_id:  # 中心节点
                plt.text(x, y + 0.15*r1, label, ha='center', va='bottom', fontsize=14, fontweight='bold')
                if len(G.nodes) == 1:
                    margin = 2
                    plt.xlim(x - margin, x + margin)
                    plt.ylim(y - margin, y + margin)
            elif n in upstream:  # 上游
                plt.text(x - 0.05*r1, y, label, ha='right', va='center', fontsize=12)
            elif n in downstream:  # 下游
                plt.text(x + 0.05*r1, y, label, ha='left', va='center', fontsize=12)
            else:
                plt.text(x, y, label, ha='center', va='center', fontsize=12)
        # nx.draw_networkx_labels(G, pos, labels={n:G.nodes[n]['label'] for n in G.nodes}, font_size=12)
        # plt.title(f"Glyco Network for {feature}", fontsize=20)
        plt.axis('off')
        plt.tight_layout()
        
        if len(G.nodes) != 1:
            # 获取中心节点坐标
            x_center, y_center = pos[glycan_center_id]
            r_center = get_radius(G.nodes[glycan_center_id]['size'])
            # ---- 1. 节点图例（上方） ----
            node_legend_handles = [
                mpatches.Circle((0, 0), radius=0.1, color='#E99E75', label='Enzyme'),
                mpatches.Circle((0, 0), radius=0.1, color='#776483', label='Glycan-binding protein'),
                mpatches.Circle((0, 0), radius=0.1, color='gray', label='Downstream pathway')
            ]
            legend_y_gap = 2.9  # 与中心节点的距离（可以微调）
            legend_y = y_center + legend_y_gap
            legend_x = x_center
            circle_colors = ['#E99E75', '#776483', 'gray']
            circle_labels = ['Enzyme', 'Glycan-binding protein', 'Downstream pathway']
            for i, (c, lab) in enumerate(zip(circle_colors, circle_labels)):
                plt.scatter([legend_x - 1.3], [legend_y - i*0.3], s=450, color=c, zorder=10)
                plt.text(legend_x - 0.9, legend_y - i*0.3, lab, va='center', ha='left', fontsize=14)
            # ---- 2. 连线图例（下方） ----
            edge_legend_handles = [
                mlines.Line2D([], [], color='#3558AE', linewidth=5, label='Positive correlation'),
                mlines.Line2D([], [], color='#B64074', linewidth=5, label='Negative correlation'),
                mlines.Line2D([], [], color='#1F5F5B', linewidth=5, label='Pathway overlap')
            ]
            edge_legend_y_gap = 2.3  # 下方距离
            edge_legend_y = y_center - edge_legend_y_gap
            edge_legend_x = x_center
            for i, handle in enumerate(edge_legend_handles):
                plt.plot([edge_legend_x - 1.5, edge_legend_x -1],
                         [edge_legend_y - i*0.3, edge_legend_y - i*0.3],
                         color=handle.get_color(), linewidth=handle.get_linewidth())
                plt.text(edge_legend_x - 0.9, edge_legend_y - i*0.3, handle.get_label(), va='center', ha='left', fontsize=14)
        
        # plt.show()
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        plt.savefig(os.path.join(output_dir, f'{regulation_type}_{feature}_network_{timestamp}.pdf'), format='pdf')
        png_file = os.path.join(output_dir, f'{regulation_type}_{feature}_network_{timestamp}.png')
        plt.savefig(png_file, format='png')
        plt.close()
        
        png_file = png_file
        image = Image.open(png_file)
        trimmed_image = self.trim_white_border(image)
        trimmed_image.save(png_file)

        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'glyconetwork', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'glyconetwork',f'{regulation_type}_{feature}_network.pdf')

        return {'file_path': os.path.join(output_dir, f'{regulation_type}_{feature}_network_{timestamp}.png'), 
                'legend': figure_description}
    
    def draw_glycans(self, chain_list, linewidth=0.2,
                     node_size = 12,
                     subfolder='plot',
                     figure_description = 'Glycans plot',
                     filename = None):
        """
        Visualize glycan chain structures as symbolic tree-like diagrams.

        This function parses simplified string representations of glycan chains
        into hierarchical tree structures, applies color and layout rules,
        and renders each glycan in a separate subplot. Nodes are drawn with
        specific shapes and colors depending on monosaccharide type, and 
        branches are connected with lines. Optional fucose (F5) residues are 
        drawn as special side branches.

        Args:
            chain_list (list of str): 
                List of glycan chain strings to visualize. 
                Each chain string encodes branching with uppercase letters for opening 
                residues (with numeric types, e.g., 'A1', 'B2', 'F5') and lowercase letters 
                for branch closures.
                Example: "A1B2b" (root A1 with child B2, then branch closed).
            linewidth (float, optional): 
                Line width for edges and node borders. Default = 0.2.
            node_size (int, optional):
                Monosaccharide node size (if your matplotlib version > 3.5.3, please select 12 as input, otherwise select 400 for lower version).
            subfolder (str, optional): 
                Subdirectory under `./plot/` for saving outputs. Default = 'plot'.
            figure_description (str, optional): 
                Text description of the figure, logged in metadata. Default = 'Glycans plot'.
            filename (str, optional): 
                Base filename for saving figures (without extension). Default = None.

        Returns:
            dict:
                - "file_path" (str): Path to the saved PNG file.
                - "legend" (str): Figure description.

        Visualization details:
            - **Node types (num_type → color & marker):**
                * 1 → Green circle (`#00C832`)
                * 2 → Blue square (`#0000FA`)
                * 3 → Purple diamond (`#C800C8`)
                * 4 → White diamond
                * 5 → Red triangle (`#FA0000`, typically fucose)
            - **Branching logic:**
                * Uppercase + digit starts a branch (e.g., 'A1').
                * Lowercase closes the matching branch (e.g., 'a' closes 'A').
                * F5 nodes (fucose) attach laterally as `fucose_child`.
            - **Layout:**
                * X-coordinates assigned to center branches evenly.
                * Y-coordinates increase by `level_gap` for each level.
                * Fucose (F5) residues are placed at a side offset.
            - **Labels:**
                * Each chain string is annotated below its corresponding diagram.

        Output:
            - One subplot per glycan chain.
            - PDF and PNG files saved in `./plot/{subfolder}/` with a timestamp.
            - Parameters and outputs logged using `self.data_manager`.

        Example:
            >>> vis.draw_glycans(["A1B2b", "C1D1dE2e"])
            # Generates and saves glycan diagrams for both chains.
        """
        # 内部 Node 类（每个节点记录字母、类型、颜色、坐标等信息）
        class Node:
            def __init__(self, letter, num_type):
                self.letter = letter
                self.num_type = num_type
                self.children = []
                self.fucose_child = None  # F5 特殊分支（红色三角形）
                self.x = self.y = 0
                self.width = 1  # 默认宽度
                if num_type == 1:
                    self.color = "#00C832"
                elif num_type == 2:
                    self.color = "#0000FA"
                elif num_type == 3:
                    self.color = "#C800C8"
                elif num_type == 4:
                    self.color = "white"
                elif num_type == 5:
                    self.color = "#FA0000"
        # 解析单个糖链字符串，构建糖链树
        def parse_chain(chain_str):
            stack = []
            root = None
            i = 0
            while i < len(chain_str):
                ch = chain_str[i]
                if ch.isupper():  # 新分支开始
                    letter = ch
                    if i+1 >= len(chain_str) or not chain_str[i+1].isdigit():
                        raise ValueError(f"Invalid format at position {i}")
                    num = int(chain_str[i+1])
                    node = Node(letter, num)
                    if stack:
                        parent = stack[-1]
                        if num == 5:
                            parent.fucose_child = node
                        else:
                            parent.children.append(node)
                    else:
                        root = node
                    stack.append(node)
                    i += 1  
                elif ch.islower():  # 结束分支
                    if not stack:
                        raise ValueError(f"Unmatched closing at position {i}")
                    node = stack.pop()
                    if node.letter.lower() != ch:
                        raise ValueError(f"Mismatched branch closing '{ch}' for '{node.letter}'")
                else:
                    raise ValueError(f"Invalid character '{ch}' in chain string")
                i += 1
            if stack:
                raise ValueError("Unclosed branches in chain string")
            return root
        # 调整颜色规则（可根据需要修改）
        def apply_color_rules(node):
            if node.num_type == 1:
                parent = getattr(node, 'parent', None)
                # 这里可根据需求添加或修改规则
                if (parent and parent.letter == 'E' and parent.num_type == 2 and node.letter == 'F') or (node.letter >= 'G'):
                    node.color = "#FFFF00"
            if node.num_type == 2:
                parent = getattr(node, 'parent', None)
                if parent and parent.letter == 'E' and parent.num_type == 2 and node.letter == 'F':
                    node.color = "#FFFF00"
            if node.fucose_child:
                node.fucose_child.parent = node
                apply_color_rules(node.fucose_child)
            for child in node.children:
                child.parent = node
                apply_color_rules(child)
        # 计算节点的子树宽度，并考虑 F5 额外间距以及特殊调整
        def compute_width(node):
            if not node.children:
                node.width = 1.3   # 末端分支间隔
            else:
                total_width = sum(compute_width(c) for c in node.children)
                node.width = total_width
            # 调整：让 C1 下面的 D1 们对称
            if node.letter == 'C' and node.children and all(c.letter == 'D' for c in node.children):
                max_width = max(c.width for c in node.children)
                for c in node.children:
                    c.width = max_width * 0.9  # D层分支间隔
                node.width = len(node.children) * max_width
            if node.fucose_child:
                node.width += 2.5  # F5 分支占用额外宽度
            return node.width
        # 分配坐标，x_center 为当前糖链的中心位置，y_level 为起始纵坐标
        def assign_coords(node, x_center, y_level=0, level_gap=1.7, side_offset=1.5):
            node.x, node.y = x_center, y_level
            if node.children:
                total_width = sum(c.width for c in node.children)
                start_x = x_center - (total_width - 1) / 2.0
                current_x = start_x
                for c in node.children:
                    center_x = current_x + (c.width - 1) / 2.0
                    assign_coords(c, center_x, y_level + level_gap, level_gap, side_offset)
                    current_x += c.width
            if node.fucose_child:
                fuc = node.fucose_child
                fuc.x = x_center + side_offset * 1.7  # 调整 F5 分支的横向位置
                fuc.y = y_level
                if hasattr(fuc, "children") and fuc.children:
                    assign_coords(fuc, fuc.x, y_level + level_gap, level_gap, side_offset)
        # 获取糖链中最大/最小 y 值
        def get_max_y(node):
            max_y = node.y
            if node.fucose_child:
                max_y = max(max_y, get_max_y(node.fucose_child))
            for child in node.children:
                max_y = max(max_y, get_max_y(child))
            return max_y
        def get_min_y(node):
            min_y = node.y
            if node.fucose_child:
                min_y = min(min_y, get_min_y(node.fucose_child))
            for child in node.children:
                min_y = min(min_y, get_min_y(child))
            return min_y
        # 递归绘制节点和连线
        def draw_node(node, ax):
            if hasattr(node, "parent") and node.parent:
                px, py = node.parent.x, node.parent.y
                ax.plot([px, node.x], [py, node.y], color='black', linewidth=linewidth, zorder=1)
            # 根据 num_type 设置不同形状
            if version.parse(matplotlib.__version__) > version.parse("3.5.3"):
                if node.num_type == 1:
                    marker = 'o'
                    size = int(node_size)
                elif node.num_type == 2:
                    marker = 's'
                    size = int(node_size)
                elif node.num_type in (3, 4):
                    marker = 'D'
                    size = int(node_size) / 2
                elif node.num_type == 5:
                    marker = '<'
                    size = int(node_size)
                else:
                    marker = 'o'
                    size = int(node_size)
            else:
                if node.num_type == 1:
                    marker = 'o'
                    size = int(node_size)
                elif node.num_type == 2:
                    marker = 's'
                    size = int(node_size)
                elif node.num_type in (3, 4):
                    marker = 'D'
                    size = int(node_size) / 2
                elif node.num_type == 5:
                    marker = '<'
                    size = int(node_size)
                else:
                    marker = 'o'
                    size = int(node_size)
            ax.scatter(node.x, node.y, s=size, marker=marker, 
                       facecolor=node.color, edgecolor='black', linewidths=linewidth, zorder=2)
            if node.fucose_child:
                draw_node(node.fucose_child, ax)
            for c in node.children:
                draw_node(c, ax)
        # 为了让每条糖链单独绘制在一个 subplot 上，这里创建多个子图
        fig, axs = plt.subplots(
            nrows=len(chain_list),
            ncols=1,
            figsize=(10, 4*len(chain_list))
        )
        # 如果只有一条糖链，axs 不是列表，需要特殊处理
        if len(chain_list) == 1:
            axs = [axs]
        # 遍历糖链列表，逐个绘制到各自的子图
        for ax, chain_str in zip(axs, chain_list):
            ax.set_aspect('equal')
            ax.axis('off')
            # 构建糖链树
            root = parse_chain(chain_str)
            apply_color_rules(root)
            compute_width(root)
            # 在当前子图内，将糖链居中绘制
            glycan_center = (root.width - 1) / 2.0
            assign_coords(root, glycan_center, y_level=0, level_gap=1.5, side_offset=0.8)
            # 绘制节点
            draw_node(root, ax)
            # 获取最大最小 y 值，用于设定显示范围
            max_y = get_max_y(root)
            min_y = get_min_y(root)
            # 在糖链下方添加标注（y 值比最小值再低一点）
            ax.text(glycan_center, min_y - 0.8, chain_str,
                    ha='center', va='top', fontsize=10)
            # 根据糖链大小动态设置坐标范围
            ax.set_xlim(-1, root.width + 1)
            ax.set_ylim(min_y - 2, max_y + 2)
        plt.tight_layout(pad=2.0)
        
        output_dir = os.path.join('./plot', subfolder)
        os.makedirs(output_dir, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        # 显示图像
        plt.savefig(os.path.join(output_dir, f"{filename}_glycans_plot_{timestamp}.pdf"),dpi=900, bbox_inches='tight')
        png_file = os.path.join(output_dir, f"{filename}_glycans_plot_{timestamp}.png")
        plt.savefig(png_file,dpi=900, bbox_inches='tight')
        plt.close('all')
        # png_file = png_file
        # image = Image.open(png_file)
        # trimmed_image = self.trim_white_border(image)
        # trimmed_image.save(png_file)
        # 自动获取所有参数
        params = locals()
        # 去掉不需要记录的局部变量 'self'
        params.pop('self')
        # 使用data_manager记录这些参数
        self.data_manager.log_params('StrucGAP_DataVisualization', 'glycans_plot', params)
        self.data_manager.log_output('StrucGAP_DataVisualization', 'glycans_plot',f"glycans_plot.pdf")
        return {'file_path': os.path.join(output_dir, f"{filename}_glycans_plot_{timestamp}.png"), 
                'legend': figure_description}
        fig.savefig(output_pdf)
    
    def add_figure(self, figure_meta: dict, figure_name: str = "default"):
        """
        Adds a figure to a specified figure collection.
    
        Parameters:
            figure_meta (dict): The metadata dictionary returned by the plotting function.
            figure_name (str): The name of the figure collection (optional, default is "default").
    
        Returns:
            None
            
        """
        if figure_name not in self.figure_collections:
            self.figure_collections[figure_name] = []
        self.figure_collections[figure_name].append(figure_meta)
        self.current_figure = figure_name
        
    def _calculate_layout(self, num_plots: int) -> dict:
        """
        Calculates the layout for arranging figures on an A4 page.
        
        Parameters:
            num_plots (int): The number of plots to arrange.
    
        Returns:
            dict: The layout dictionary with page dimensions, row and column configurations, and image regions.
        
        """
        page_width, page_height = A4
        margins = 10 * mm  # 可以适当加大页边距以免过于贴边
    
        # 固定 5 行，每行高度相同
        total_rows = 5
        row_height = (page_height - 2 * margins) / total_rows
    
        # 每行固定 3 列
        cols = 3
        col_width = (page_width - 2 * margins) / cols
    
        # 实际需要的行数(仅用于放图片的) = min(4, ceil(num_plots / 3))
        # 最多支持 12 张图(4 行×3 列)
        image_rows = min(4, ceil(num_plots / 3))
    
        # 如果超过 12 张图，这里可酌情处理(报错或者另开新页等)
        if num_plots > 12:
            raise ValueError("当前布局最多只支持 12 张图，请自行扩展逻辑。")
    
        layout = {
            "page_width": page_width,
            "page_height": page_height,
            "margins": margins,
            "row_height": row_height,
            "cols": cols,
            "col_width": col_width,
            "image_rows": image_rows
        }
        return layout
    
    def _number_to_letter(self, num: int) -> str:
        """
        Converts a number to a corresponding lowercase letter.
    
        Parameters:
            num (int): The number to convert.
    
        Returns:
            str: The corresponding lowercase letter.
            
        """
        return chr(ord('a') + num - 1)
    
    def compose_figures(self, output_path: str, figure_name: str = None, custom_sizes: list = None,
                       pdf_description: str = None, pdf_description_bg_color='white'):
        """
        Composes a combined PDF of figures, allowing for custom image area assignments.
    
        Parameters:
            output_path (str): The file path for the output PDF.
            figure_name (str): The name of the target figure collection (optional).
            custom_sizes (list): A list specifying custom areas for each image (optional).
    
        Returns:
            None
            
        """
        # 确定目标 figure
        target_figure = figure_name or self.current_figure
        if not target_figure or target_figure not in self.figure_collections:
            raise ValueError(f"未找到指定 figure 集合: {target_figure}")
        figures = self.figure_collections[target_figure]
        if not figures:
            print("警告：当前 figure 集合为空")
            return

        # 检查 custom_sizes
        if custom_sizes and len(custom_sizes) != len(figures):
            raise ValueError("custom_sizes 长度必须与图片数量一致")
        if custom_sizes:
            used_positions = set()
            for positions in custom_sizes:
                for pos in positions:
                    if pos < 1 or pos > 12:
                        raise ValueError(f"区域编号 {pos} 超出范围（1-12）")
                    if pos in used_positions:
                        raise ValueError(f"区域编号 {pos} 被重复使用")
                    used_positions.add(pos)

        # 计算布局
        layout = self._calculate_layout(len(figures) if not custom_sizes else 12)  # 默认占满12个区域
        page_width = layout["page_width"]
        page_height = layout["page_height"]
        margins = layout["margins"]
        row_height = layout["row_height"]
        cols = layout["cols"]
        col_width = layout["col_width"]
        image_rows = layout["image_rows"]

        # 创建画布
        c = canvas.Canvas(output_path, pagesize=A4)
        # -------- 添加标题区域（不影响图像布局） --------
        if pdf_description:
            title_font_size = 10
            title_padding = 4  # 上下 padding，适配10号字体
            title_height = 21
            c.setFillColor(HexColor(pdf_description_bg_color))  
            c.rect(
                0, 
                page_height - title_height, 
                page_width, 
                title_height, 
                stroke=0, 
                fill=1
            )
            c.setFillColor(HexColor("#000000"))  
            c.setFont("Helvetica-Bold", title_font_size)
            c.drawString(
                2 * mm, 
                page_height - title_height + 6,  
                pdf_description
            )

        # -------- 绘制图片部分 --------
        max_row = 0  # 记录图片占用的最大行数（从0开始）
        for idx, fig_meta in enumerate(figures):
            if custom_sizes:
                # 使用 custom_sizes 指定位置
                positions = custom_sizes[idx]  # 该图片占据的区域编号列表
            else:
                # 默认布局：按顺序填充
                positions = [idx + 1]

            # 计算占用区域的边界
            min_pos = min(positions)
            max_pos = max(positions)
            row_start = (min_pos - 1) // cols  # 起始行（从0开始）
            row_end = (max_pos - 1) // cols    # 结束行
            col_start = (min_pos - 1) % cols   # 起始列
            col_end = (max_pos - 1) % cols     # 结束列

            # 更新最大行数
            max_row = max(max_row, row_end)

            # 计算图片的坐标和大小
            x = margins + col_start * col_width
            row_top_y = page_height - margins - row_start * row_height
            width = (col_end - col_start + 1) * col_width
            height = (row_end - row_start + 1) * row_height

            # 留出顶部编号空间和底部边距
            label_space = 8
            img_y = row_top_y - height + 5
            img_h = height - label_space - 10

            # 绘制图片
            img = ImageReader(fig_meta["file_path"])
            c.drawImage(
                img,
                x,
                img_y,
                width=width - 5,
                height=img_h,
                preserveAspectRatio=True
            )

            # 绘制“Fig.x”编号
            c.setFont("Helvetica", 10)
            label_x = x + 2
            label_y = row_top_y - 2
            c.drawString(label_x, label_y, self._number_to_letter(idx + 1))

        # -------- 绘制文字说明 --------
        # 文字区域位于图片区域的下一行（max_row + 1），最大为第4行（索引3）
        text_row_index = min(max_row + 1, 4)  # 确保不超过第4行（索引3）
        text_top_y = page_height - margins - text_row_index * row_height
        current_y = text_top_y - 15
        styles = getSampleStyleSheet()
        normal_style = styles["Normal"]
        normal_style.fontName = "Helvetica"
        normal_style.fontSize = 10
        normal_style.leading = 12
        max_text_width = page_width - 2 * margins
        for i, fig_meta in enumerate(figures):
            legend_str = fig_meta.get('legend', '')
            full_text = f"{self._number_to_letter(i + 1)}: {legend_str}"
            lines = simpleSplit(full_text, normal_style.fontName, normal_style.fontSize, max_text_width)
            for line in lines:
                c.drawString(margins, current_y, line)
                current_y -= normal_style.leading

        # 保存并清理
        c.save()
        del self.figure_collections[target_figure]
        


