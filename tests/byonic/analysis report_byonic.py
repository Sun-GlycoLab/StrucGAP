# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 11:38:48 2025

@author: 28051
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 12 17:08:44 2025

@author: 28051
"""
from PyPDF2 import PdfMerger
import random
def generate_random_colors(n):
    colors = []
    for _ in range(n):
        # 随机生成一个 RGB 颜色并转为十六进制字符串
        color = "#{:02X}{:02X}{:02X}".format(
            random.randint(0, 255),
            random.randint(0, 255),
            random.randint(0, 255)
        )
        colors.append(color)
    return colors
# 示例：生成10个颜色
random_colors = generate_random_colors(10)
print(random_colors)

import matplotlib.spines as msp
_original_set_linewidth = msp.Spine.set_linewidth
def _patched_set_linewidth(self, lw):
    if isinstance(lw, str):
        lw = float(lw)
    return _original_set_linewidth(self, lw)
msp.Spine.set_linewidth = _patched_set_linewidth

# StrucGAP_GlycanStructure_1
fig1 = module7.bar(data = module2.GlycanComposition_rank.copy(),
            top = 10,
            y_column='GlycanComposition',
            y_column_value='GlycanComposition_count',
            xaxis_label_rotate = 45,
            xaxis_label_margin = 20,
            transform_ratio = False,
            subfolder='StrucGAP_GlycanStructure_1',
            colors = '#0078FF',
            y_max = None,
            yaxis_splitline_show = False,
            legend = 'Glycan composition number',
            filename='Isomer counts',
            xaxis_title = 'Top 10 glycan compositions',
            xaxis_title_gap = 110,
            yaxis_title = 'Number of distinct glycan compositions',
            yaxis_title_gap = 60,
            xaxis_label_font_size = 30,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            figure_description = 'Number of top 10 glycan compositions.',
            )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_1")
module7.compose_figures("StrucGAP_GlycanStructure_1.pdf", figure_name="StrucGAP_GlycanStructure_1",
                        pdf_description = "StrucGAP_GlycanStructure_1: Overview of dominant glycan structures and substructure features.",
                        pdf_description_bg_color = "#89C7CB",
                        custom_sizes=[[1]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoSite_1
plot_data = pd.DataFrame(index=['Preprocessed data'],
                         columns=['PSMs','Unique IGPs','Glycoproteins','Glycosites',
                                 'Peptides','Glycan compositions'])
plot_data['PSMs'] = module1.data_peptide_fdr_data.shape[0]
plot_data['Unique IGPs'] = module1.data_psm_filtered.shape[0]
temp_data = module1.data_psm_filtered.copy()
temp_data = pd.DataFrame(temp_data[['ProteinID', 'Glycosite_Position']])
temp_data = temp_data.dropna(subset=['ProteinID', 'Glycosite_Position'])
temp_data = temp_data[(temp_data['ProteinID'].str.len() > 0) & (temp_data['Glycosite_Position'].str.len() > 0)]
temp_data['ProteinID'] = temp_data['ProteinID'].str.split(';')
temp_data['Glycosite_Position'] = temp_data['Glycosite_Position'].str.split(';')
temp_data = temp_data.explode(['ProteinID', 'Glycosite_Position'])
temp_data['glycosite'] = temp_data['ProteinID'] + temp_data['Glycosite_Position']
plot_data['Glycoproteins'] = temp_data['ProteinID'].nunique()
plot_data['Glycosites'] = temp_data['glycosite'].nunique()
plot_data['Peptides'] = module1.data_psm_filtered['PeptideSequence'].nunique()
plot_data['Glycan compositions'] = module1.data_psm_filtered['GlycanComposition'].nunique()
# plot_data['Glycan structures'] = module1.data_psm_filtered['structure_coding'].nunique()
fig1 = module7.heatmap2(data = plot_data.T,
                 colors = 'tab10',
                 columns=['Preprocessed data'],
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None, 
                 z_score = None,
                 minvalue = plot_data.min().min(), 
                 centervalue = 3000,
                 maxvalue = plot_data.max().max(),
                 splitline_width = 1,
                 yaxis_label_show=True, 
                 text_annotation = True,
                 text_size = 30,
                 xaxis_label_font_size=30, 
                 yaxis_label_font_size=30, 
                 xaxis_label_rotate = 0,
                 text_color = 'white',
                 subfolder='StrucGAP_GlycoSite_1',
                 figure_description = 'Overall profiling of preprocessed glycoproteomics data.',
                 filename='overview heatmap',
                 )

fig2 = module7.nested_pie(data = module3.glycoprotein_glycosite_count, 
                   item_column = 'glycoprotein', 
                   number_column = 'glycosite_count',
                   value_counts_column = 'glycosite_count',
                   subfolder='StrucGAP_GlycoSite_1',
                   inner_data_first=False,
                   label_font_size = 15,
                   legend_font_size = 15,
                   split=9,
                   filename = 'glycosite_count',
                   figure_description = 'Number of glycoproteins grouped by the count of glycosylation sites.',
                   )

fig3 = module7.bar(data = module3.glycoprotein_glycan_count,
            top = 10,
            y_column='gene_name',
            y_column_value='glycan_count',
            subfolder='StrucGAP_GlycoSite_1',
            colors = '#159947',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Count of glycans',
            xaxis_title = 'Glycoprotein',
            xaxis_title_gap = 140,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Number of glycans',
            figure_description = 'Top 10 glycoproteins ranked by the number of glycan across all glycosylation sites.',
            yaxis_title_gap = 70,
            xaxis_label_font_size = 30,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='glycoprotein_glycan_type'
            )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoSite_1")
# module7.add_figure(fig4, figure_name="StrucGAP_GlycoSite_1")
# module7.add_figure(fig5, figure_name="StrucGAP_GlycoSite_1")
# module7.add_figure(fig6, figure_name="StrucGAP_GlycoSite_1")
module7.compose_figures("StrucGAP_GlycoSite_1.pdf", figure_name="StrucGAP_GlycoSite_1",
                        pdf_description = "StrucGAP_GlycoSite_1: Comprehensive overview of the dataset and glycoprotein glycosite diversity at the site level.",
                        pdf_description_bg_color = "#89C7CB",
                        custom_sizes=[[1], [2,3,5,6], [4]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_both_glycoproteins_1
module5 = StrucGAP_FunctionAnnotation(module1, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  

plot_data = module5.ora_no_background_both_proteins_result.copy()
def clean_term(term):
    term = re.sub(r'\s*\(.*?\)', '', term)  
    term = term.title()  
    return term
plot_data['Term'] = plot_data['Term'].apply(clean_term)
fig1 = module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'Gene_set',
    p_column = 'P-value',
    top = 10,
    xaxis_font_size = 10,
    yaxis_font_size = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    filename="both proteins dotplot",
    figure_description = 'Enrichment results of identified glycoproteins based on GO enrichment.',
    col_cluster=False,
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.compose_figures("StrucGAP_FunctionAnnotation_both_glycoproteins_1.pdf", figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1",
                        pdf_description = "StrucGAP_FunctionAnnotation_1: GO enrichment of glycoproteins and association analysis with substructure features.",
                        pdf_description_bg_color = "#FDAE86",
                        custom_sizes=[[1,2,3]])  # 生成后自动清理figure1数据


## key insights report

# if you want to get the combined pdf categoried by module, you can ...
from PyPDF2 import PdfMerger

merger = PdfMerger()
merger.append("StrucGAP_GlycoSite_1.pdf")
merger.append("StrucGAP_GlycanStructure_1.pdf")
merger.append("StrucGAP_FunctionAnnotation_both_glycoproteins_1.pdf")
merger.append("Glycan structures coding rule.pdf")
merger.write("glycodecipher_analysis_report.pdf")
merger.close()


