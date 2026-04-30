# -*- coding: utf-8 -*-
# 初始化数据管理系统
from strucgap.preprocess import StrucGAP_Preprocess
from strucgap.glycanstructure import StrucGAP_GlycanStructure
from strucgap.glycosite import StrucGAP_GlycoSite
from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
from strucgap.functionannotation import StrucGAP_FunctionAnnotation
from strucgap.glyconetwork import StrucGAP_GlycoNetwork
from strucgap.datavisualization import StrucGAP_DataVisualization
from strucgap.insighttracker import StrucGAP_InsightTracker
import os
data_manager = StrucGAP_InsightTracker()
# data_manager.module_records.keys()
# data_manager.analysis_params
# data_manager.outputs
# data_manager.key_information_extraction()

# data_manager.read_pickle()

os.chdir('D:\\doctor\\analysisys\\StrucGAP')

# 
module1 = StrucGAP_Preprocess(data_dir="D:\\doctor\\StrucGAP\\tests\\mouse uterus.xlsx",
                      data_sheet_name = '1 PSM',
                      sample_group_data_dir = 'D:\\doctor\\analysisys\\data\\sample_group.xlsx',
                      branch_list_dir = "D:\\doctor\\wyq\\branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager)
module1.data_cleaning(data_type='tmt')
module1.cv_raw(threshold='no')
# 1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464
module1.fdr(feature_type='no')
module1.outliers(abundance_ratio=[1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464],
                 samplewise_normalization = False,
                 total_intensity_normalization=True,
                 total_intensity_method='mean')
module1.cv(threshold = 'no')
module1.psm(psm_number = 'no')
# Using glytoucan = True and biosynthetic_pathways = True is a very time-consuming task, due to the limitations of the GlyTouCan and KEGG APIs. Please be patient when enabling these two annotations. If you prefer faster execution, set both options to False.
module1.annotation(glytoucan = True, glytoucan_structure = True, glytoucan_wurcs_file = "D:\\doctor\\StrucGAP\\tests\\glycosmos_glycans_wurcs.csv", biosynthetic_pathways = True, glycobiology_filter = True)
module1.output() 

# 
module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = False) 
module2.structure_statistics()
module2.lacdinac()
module2.cor()
module2.isoforms()
module2.output()

#
module3 = StrucGAP_GlycoSite(module1, data_manager=data_manager)
module3.glycoprotein_site()
module3.glycopeptide_site()
module3.specific_site()
module3.output()

#
module4 = StrucGAP_GlycoPeptideQuant(module1, data_type = 'psm_filtered', data_manager=data_manager)
module4.statistics()
module4.statistics_index()
module4.differential_analysis(pvalue_type='pvalue_ttest', fc = 1.7)
module4.threshold_variation_analysis(pvalue_type='pvalue_ttest',statistic_index='fc')
module4.glycopeptide_glycosite_glycan_variation(fc = 1.7)
module4.glycoprotein_glycosite_glycan_variation(fc = 1.7)
module4.output()

data_manager.key_information_extraction(module='StrucGAP_GlycoPeptideQuant')

#  
module5 = StrucGAP_FunctionAnnotation(module1, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
# module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
#                                  data_manager=data_manager)  
# module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
# module5.kegg_function_structure(function_data = 'ora_no_background_both_proteins_result')  
# module5.output()

module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', enrich_feature='glycopeptide',
            background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest') # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_up_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

module5.go_function_structure(function_data = 'ora_no_background_down_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

# module5.go_function_structure(function_data = 'gsea_result')  
# module5.output()

# module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 119
# module5.gsea()
# module5.kegg_function_structure(function_data = 'ora_no_background_up_result')
# module5.output()
# module5.kegg_function_structure(function_data = 'ora_no_background_down_result')
# module5.output()
# module5.kegg_function_structure(function_data = 'gsea_result')
# module5.output()

module5 = StrucGAP_FunctionAnnotation(module6, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
# module5.gsea()
module5.go_function_structure(function_data = 'ora_no_background_up_result')  
module5.output()
module5.go_function_structure(function_data = 'ora_no_background_down_result')  
module5.output()
# module5.go_function_structure(function_data = 'gsea_result')  
# module5.output()

#
module6 = StrucGAP_GlycoNetwork(module4, data_manager=data_manager)
module6.proteomic(protein_data_dir="D:\\doctor\\wyq\\WYQ_Mus_uterus_global.xlsx",
                  data_sheet_name = '1 Proteins',
                  cv = 'no', psm = 1, 
                  normalization_samplewise_method = 'no',
                  normalization_featurewise_method = 'no',
                  )
module6.phosphorylation(phospho_data_dir="D:\doctor\zzd\_Rat_Phospho_mixThymus_TMT6c_.xlsx",
                        data_sheet_name='PeptideGroups')
module6.glycosyltransferases(glycosyltransferases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", 
                             data_sheet_name="glycosyltransferases")
module6.glycosidases(glycosidases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", 
                     data_sheet_name='glycosidases')
module6.sialyltransferases()
module6.fucosyltransferase()
module6.glycan_binding_protein()
module6.output()

data_manager.key_information_extraction(module='StrucGAP_GlycoNetwork')
#
data_manager.output_pickle()

#
module7 = StrucGAP_DataVisualization(data_manager=data_manager)


# overview
module7.venn_diagram(
    'module1.data["ProteinID"]',
    'module6.protein_raw_data["Accession"]',
    colors = 'Tropic',
    subfolder='overview',
    legend = ['Glycoproteins', 'Proteins'],
    filename='glyco and proteomic protein'
)

# module7.upset_plot(
#     'module1.data["ProteinID"]',
#     'module6.protein_raw_data["Accession"]',
#     colors = 'blue',
#     subfolder='overview',
#     filename='glyco and proteomic protein'
# )

# bothfdr = module1.data
# module1 = StrucGAP_Preprocess(data_dir="D:\\doctor\\wyq\\uterus\\RatOvary_IGP_onlypeptides.xlsx",
#                       data_sheet_name = 'PSM',
#                       sample_group_data_dir = 'D:\\doctor\\analysisys\\data\\sample_group.xlsx',
#                       branch_list_dir = "D:\\doctor\\wyq\\branch_structures_18_mice uterus.0240401.xlsx",
#                       data_manager=data_manager)
# module1.data_cleaning()
# peptidefdr = module1.data

# # 从保存的文件中加载变量
# with open('all_variables.pkl', 'rb') as f:
#     loaded_variables = pickle.load(f)
# # 将变量重新加载到当前全局命名空间中
# globals().update(loaded_variables)

module7.heatmap_multi_data('module1.data',  
              'module1.data_psm_filtered',   
              columns=['MS2Scan', 'ProteinID', 'PeptideSequence', 'structure_coding', 
                       'GlycanComposition', 'ProteinID+Glycosite_Position'],
              statistical_methods = ['both', 'unique', 'unique', 'unique', 
                                     'unique', 'unique'],
              subfolder='overview',
              colors = 'coolwarm',
              filename='cleaned and psm filtered data'
              )

# glycosite
plot_data = pd.DataFrame(module3.glycoprotein_glycosite_count['glycosite_count'].value_counts())
plot_data['glycosite_type'] = plot_data.index
plot_data = plot_data.sort_values('glycosite_type', ascending=True)
plot_data['ratio'] = plot_data['glycosite_count']/plot_data['glycosite_count'].sum()

module7.pie(data = plot_data, 
            top = 9,
              item_column = 'glycosite_type', 
              number_column = 'glycosite_count',
              radius = ['50%', '70%'],
              rosetype = None,
              subfolder='glycosite',
              filename='top 9 glycosite_count',
              )

module7.pie(data = plot_data, 
            top = None,
            end = 5,
              item_column = 'glycosite_type', 
              number_column = 'glycosite_count',
              radius = ['10%', '70%'],
              rosetype = 'radius',
              subfolder='glycosite',
              filename='end 5 glycosite_count',
              colors = ['#F8ED70','#EB8B5B','#B53656','#632A69','#48A6B5']
              )

plot_data = pd.DataFrame(module3.protein_glycosite_glycan_count)
# plot_data['glycosite_type'] = plot_data.index
# plot_data = plot_data.sort_values('glycosite_type', ascending=True)
# plot_data['ratio'] = plot_data['glycosite_count']/plot_data['glycosite_count'].sum()

module7.pie(data = plot_data, 
            top = 3050,
              item_column = 'gene_name', 
              number_column = 'glycan_type_count',
              radius = ['0%', '100%'],
              rosetype = None,
              subfolder='glycosite',
              filename='glycan_type_count',
              )

# glycanstructure
module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = False) 
data1 = module2.structure_coding_rank
module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = True)
data2 = module2.structure_coding_rank 
module7.butterfly_plot(data1,data2,
                       item_column='Structure_coding',
                       count_column='Structure_coding_count',
                       colors = ['lightblue', 'darkblue'],
                       xaxis_title = 'Number of Glycan',
                       plot_title = 'Top 10 glycans',
                       legend = ['With oligo mannose', 'Without oligo mannose'],
                       subfolder='glycanstructure',
                       filename='top 10 structure coding'
                       )

module7.draw_glycans(data1['Structure_coding'][:10], linewidth=0.2, filename="withom")
module7.draw_glycans(data2['Structure_coding'][:10], linewidth=0.2, filename="withoutom")

module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = False) 
module2.structure_statistics()
module2.lacdinac()
module2.cor()
module2.isoforms()

plot_data = pd.DataFrame(module2.glycan_composition_isoforms['GlycanComposition'].value_counts())
plot_data = pd.DataFrame(plot_data['GlycanComposition'].value_counts())
plot_data['isomer_count'] = plot_data.index
plot_data = plot_data.sort_values('isomer_count', ascending=False)
module7.polar1('plot_data', 
              columns=['isomer_count'],
              number_column = 'GlycanComposition',
              subfolder='glycanstructure',
              if_unique = False,
              radiusaxis_label_show = True,
              colors = ['#81B29A'],
              plot_title = None,
              filename='isomer',
              )

plot_data = pd.DataFrame(module2.glycan_composition_isoforms)
plot_data = plot_data[plot_data['GlycanComposition']=='N4H5F1G1']
module7.draw_glycans(plot_data['structure_coding'], linewidth=0.2, filename="isomers")

module7.pie(data = module2.core_structure, 
              item_column = 'Core_structure', 
              number_column = 'Core_structure_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='glycanstructure',
              filename='corestructure',
              colors = ['#F8ED70','#EB8B5B','#B53656','#632A69']
              )

module7.pie(data = module2.glycan_type, 
              item_column = 'Glycan_type', 
              number_column = 'Glycan_type_count',
              radius = ['0%', '70%'],
              rosetype = 'area',
              subfolder='glycanstructure',
              filename='glycantype',
              colors = ['#F8ED70','#EB8B5B','#B53656','#632A69']
              )

module7.bar(data = module2.branches_structure,
            top = 20,
            y_column='Branches',
            y_column_value='Branches_count',
            xaxis_label_font_size = 6,
            subfolder='glycanstructure',
            colors = '#8FBFB8',
            transform_ratio = True,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Branches structure',
            xaxis_title = 'Branches structure',
            xaxis_title_gap = 35,
            yaxis_title = 'Percentage (%)',
            plot_title = 'Branches structure',
            filename='branches structure'
            )

module7.bar(data = module2.branches_structure,
            top = None,
            end = 8,
            y_column='Branches',
            y_column_value='Branches_count',
            xaxis_label_font_size = 6,
            subfolder='glycanstructure',
            colors = '#8FBFB8',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Branches structure',
            xaxis_title = 'Branches structure',
            xaxis_title_gap = 35,
            yaxis_title = 'Count',
            plot_title = 'Branches structure',
            filename='end8 branches structure'
            )

module7.bar(data = module2.branches_count,
            top = 10,
            y_column='BranchNumber',
            y_column_value='BranchNumber_count',
            xaxis_label_rotate = 0,
            xaxis_label_margin = 20,
            transform_ratio = True,
            subfolder='glycanstructure',
            colors = '#FFD804',
            y_max = None,
            yaxis_splitline_show = False,
            legend = 'Branch number',
            filename='branch number',
            plot_title = 'Branch number',
            xaxis_title = 'Branch number of glycans',
            xaxis_title_gap = 35,
            yaxis_title = 'Percentage of IGPs',
            )

data1 = module1.data_psm_filtered
data1 = data1[(data1['fucosylated type']=='core fucosylated')|(data1['fucosylated type']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data1 = data1.iloc[:,0]
data2 = module1.data_psm_filtered
data2 = data2[(data2['fucosylated type']=='antenna fucosylated')|(data2['fucosylated type']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data2 = data2.iloc[:,0]
module7.venn_diagram(
    data1,
    data2,
    colors = 'Tropic',
    subfolder='glycanstructure',
    legend = ['core fucosylated', 'antenna fucosylated'],
    filename='fucosylated type'
)

data1 = module1.data_psm_filtered
data1 = data1[(data1['Ac/Gc']=='Ac')|(data1['Ac/Gc']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data1 = data1.iloc[:,0]
data2 = module1.data_psm_filtered
data2 = data2[(data2['Ac/Gc']=='Gc')|(data2['Ac/Gc']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data2 = data2.iloc[:,0]
module7.venn_diagram(
    data1,
    data2,
    colors = 'Tropic',
    subfolder='glycanstructure',
    legend = ['Ac', 'Gc'],
    filename='acgc'
)

plot_data = module2.branches_structure_core_structure
plot_data = plot_data[(plot_data['Branches']=='E2F1G4gfF4fe')|
                      (plot_data['Branches']=='E2F1G3gfe')|
                      (plot_data['Branches']=='E2F1G3gfF5fe')]
module7.sunburst(data = plot_data,
                 root_column = 'Branches',
                 child_column = 'Core_structure',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'braches_corestructure',
                 subfolder='glycanstructure',
                 )

module7.pie(data = module2.lacdinac_count, 
              item_column = 'lacdinac', 
              number_column = 'lacdinac_count',
              radius = ['10%', '70%'],
              rosetype = None,
              subfolder='glycanstructure',
              filename='lacdinac',
              colors = ['#4E2A5C','#C54E83','#4EBABC']
              )

module7.pie(data = module2.fsg, 
              item_column = 'FSG', 
              number_column = 'FSG_count',
              radius = ['10%', '70%'],
              rosetype = 'area',
              subfolder='glycanstructure',
              filename='FSG',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9']
              )

# quantification
module7.volcano_plot(data = module4.fc_result,
                    fc_column = 'fc',
                    p_column = 'pvalue_ttest',
                    fc = 1.7,
                    p_value = 0.05,
                    subfolder='quantification',
                    )

module7.dimension_reduction(
    data=module4.data_quant,
    data_columns=['126.1277', '127.1248', '127.1311', '128.1281', '128.1344',
                  '129.1315', '129.1378', '130.1348', '130.1411', '131.1382'],
    sample_group=module4.sample_group,
    filter_data=module4.fc_result, # filter_data=module4.fc_result,
    p_column='pvalue_ttest',
    p_value=0.05,
    fc = 1.7,
    method='pca',
    dimension_number = 2,
    random_state = 0,
    colors = ['#3558AE', '#B64074'],
    subfolder='quantification',
)

module7.heatmap2(data = module4.data_quant,
                 columns=['126.1277','127.1248','127.1311','128.1281','128.1344', 
                          '129.1315','129.1378','130.1348','130.1411','131.1382'],
                 filter_data = module4.fc_result,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.7, 0.05],
                 log = False,
                 z_score = 0,
                 splitline_width = 0.0000001,
                 subfolder='quantification',
                 cluster_method='complete', 
                 )

up_data = module4.up_data[module4.up_data['Glycan_type']!='Oligo mannose']
down_data = module4.down_data[module4.down_data['Glycan_type']!='Oligo mannose']
module7.heatmap_multi_data(up_data,  
              down_data,   
              columns=['structure_coding'],
              statistical_methods = ['count'],
              subfolder='quantification',
              colors = 'coolwarm',
              filename = 'top 10 differential glycan',
              )

module7.draw_glycans(up_data['structure_coding'].value_counts()[:10].index, 
                     linewidth=0.2, 
                     subfolder='quantification',
                     filename="top 10 up")
module7.draw_glycans(down_data['structure_coding'].value_counts()[:10].index, 
                     linewidth=0.2, 
                     subfolder='quantification',
                     filename="top 10 down")

module7.bar_up_down(data = module4.differential_analysis_branches_count,
                 x_column = 'BranchNumber',
                 up_column = 'Up_data_BranchNumber_count',
                 down_column = 'Down_data_BranchNumber_count',
                 subfolder='quantification',
                 colors = ['#FF4359', '#0078FF'],
                 filename="da branches count",
                 xaxis_label_text_split = 0,
                 if_stack = False,
                 )

module7.bar_up_down(data = module4.differential_analysis_branches_count,
                 x_column = 'BranchNumber',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='quantification',
                 colors = ['#FF4359', '#0078FF'],
                 filename="da branches count",
                 xaxis_label_text_split = 0,
                 )

module7.bar_up_down_ratio(feature='glycan_type', 
                          colors=['blue', 'green', 'red'],
                          subfolder='quantification',
                          filename="glycan type",
                          )

plot_data = module4.differential_analysis_data
plot_data = plot_data[(plot_data['ProteinID']=='A2AX52')&(plot_data['Glycosite_Position']=='986')]
module7.heatmap2(data = plot_data,
                 columns=['126.1277','127.1248','127.1311','128.1281','128.1344', 
                          '129.1315','129.1378','130.1348','130.1411','131.1382'],
                 filter_data = module4.fc_result,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.7, 0.05],
                 log = False,
                 z_score = 0,
                 splitline_width = 0.0000001,
                 filename="Site-specific analyses",
                 cluster = 'row',
                 subfolder='quantification',
                 )


# overview function annotation
module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  

plot_data = module5.ora_no_background_both_proteins_result
# plot_data['Numerator'] = plot_data['Overlap'].str.split('/').str[0].astype(int)
module7.dotplot_col(
    data = plot_data,  
    category = 'Gene_set',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    subfolder='functionannotation',
)

plot_data = pd.read_excel("D:\\doctor\\wyq\\uterus_20241230\\analysis_result\\StrucGAP_FunctionAnnotation_GO__ora_no_background_both_proteins_result.xlsx",sheet_name='cc_core_structure_ratio')
plot_data = plot_data.iloc[:,1:]
module7.line(data = plot_data,
             colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='functionannotation',
             symbol_size = 5,
             plot_title = None,
             xaxis_title_gap = 35,
             yaxis_title = 'Percentage (%)',
             filename='gocc core structure'
             )

plot_data = module5.cc_glycan_type
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='functionannotation',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='gocc glycan type'
)

plot_data = module5.cc_branches_structure
module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='functionannotation',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    filename='gocc branch structure'
)

plot_data = module1.data_psm_filtered.copy()
plot_data1 = plot_data[plot_data['fucosylated type']=='core fucosylated']
module5 = StrucGAP_FunctionAnnotation(plot_data1, 
                                  data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
# module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
plot_data1 = module5.ora_no_background_both_proteins_result
plot_data1 = plot_data1[plot_data1['Gene_set']=='GO:BP']
plot_data1 = plot_data1[plot_data1['P-value']<0.05]
plot_data1['type'] = 'core fucosylated'
plot_data1 = plot_data1.iloc[:10,:]
plot_data2 = plot_data[plot_data['fucosylated type']=='antenna fucosylated']
module5 = StrucGAP_FunctionAnnotation(plot_data2,
                                  data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
# module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
plot_data2 = module5.ora_no_background_both_proteins_result
plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:BP']
plot_data2 = plot_data2[plot_data2['P-value']<0.05]
plot_data2['type'] = 'antenna fucosylated'
plot_data2 = plot_data2.iloc[:10,:]
plot_data = pd.concat([plot_data1, plot_data2], axis=0)
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'type',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='fucosylation gobp',
    subfolder='functionannotation',
)

plot_data = module1.data_psm_filtered.copy()
plot_data1 = plot_data[plot_data['Ac/Gc']=='Ac']
module5 = StrucGAP_FunctionAnnotation(plot_data1, 
                                  data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
plot_data1 = module5.ora_no_background_both_proteins_result
plot_data1 = plot_data1[plot_data1['Gene_set']=='GO:BP']
plot_data1 = plot_data1[plot_data1['P-value']<0.05]
plot_data1['type'] = 'Ac'
plot_data1 = plot_data1.iloc[:10,:]
plot_data2 = plot_data[plot_data['Ac/Gc']=='Gc']
module5 = StrucGAP_FunctionAnnotation(plot_data2, 
                                  data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
plot_data2 = module5.ora_no_background_both_proteins_result
plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:BP']
plot_data2 = plot_data2[plot_data2['P-value']<0.05]
plot_data2['type'] = 'Gc'
plot_data2 = plot_data2.iloc[:10,:]
plot_data = pd.concat([plot_data1, plot_data2], axis=0)
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'type',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='AcGc gobp',
    subfolder='functionannotation',
)

plot_data = module5.bp_glycan_type
plot_data = plot_data[['Glycan_type','MACROMOLECULE GLYCOSYLATION','PROTEIN GLYCOSYLATION', 
                         'GLYCOSYLATION']]
plot_data.set_index('Glycan_type',inplace=True)
plot_data = plot_data.T
plot_data.reset_index(inplace=True)
module7.radar('plot_data',  
              columns = list(plot_data.columns[1:]),
              text_split = 10,
              subfolder='functionannotation',
              screen_column = 'Glycan_type',
              screen_values = ['Oligo mannose', 'Hybrid', 'Complex'],
              text_font_size = 6,
              filename='gobp glycan type'
              )

# differential function annotation
module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest',
            selected_terms=['GO:MF'], enrich_feature = 'glycopeptide') # 65,72,79
plot_data1 = module5.ora_no_background_up_result
# plot_data1 = pd.read_excel("D:\\doctor\\analysisys\\StrucGAP\\analysis_result\\StrucGAP_FunctionAnnotation_GO__ora_no_background_up_result.xlsx",sheet_name='GO_ora_no_bg_up_result')
plot_data1 = plot_data1[plot_data1['Gene_set']=='GO:MF']
plot_data1 = plot_data1[plot_data1['P-value']<0.05]
plot_data1['type'] = 'up'
plot_data2 = module5.ora_no_background_down_result
# plot_data2 = pd.read_excel("D:\\doctor\\analysisys\\StrucGAP\\analysis_result\\StrucGAP_FunctionAnnotation_GO__ora_no_background_up_result.xlsx",sheet_name='GO_ora_no_bg_down_result')
plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:MF']
plot_data2 = plot_data2[plot_data2['P-value']<0.05]
plot_data2['type'] = 'down'
plot_data = pd.concat([plot_data1, plot_data2], axis=0)

# plot_data_down_copy = plot_data.copy()
# plot_data_down_copy['type'] = 'down'
# plot_data = pd.concat([plot_data, plot_data_down_copy], axis=0).reset_index(drop=True)

module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'type',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='differential peptide',
    subfolder='functionannotation',
    col_cluster=False,
)

module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest',
            selected_terms=['GO:BP'], enrich_feature = 'glycopeptide') # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_up_result')
plot_data = module5.bp_core_structure
module7.line(data = plot_data,
             colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='functionannotation',
             symbol_size = 5,
             plot_title = None,
             xaxis_title_gap = 35,
             yaxis_title = 'Percentage (%)',
             filename='up gobp core structure'
             )

plot_data = module5.bp_glycan_type
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='functionannotation',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='up gobp glycan type'
)

plot_data = module5.bp_branches_structure
module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='functionannotation',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    filename='up gobp branch structure'
)

module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest',
            selected_terms=['GO:BP'], enrich_feature = 'glycopeptide') # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_down_result')
plot_data = module5.bp_core_structure
module7.line(data = plot_data,
             colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='functionannotation',
             symbol_size = 5,
             plot_title = None,
             xaxis_title_gap = 35,
             yaxis_title = 'Percentage (%)',
             filename='down gobp core structure'
             )

plot_data = module5.bp_glycan_type
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='functionannotation',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='down gobp glycan type'
)

plot_data = module5.bp_branches_structure
module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='functionannotation',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    filename='down gobp branch structure'
)

# StrucGAP_GlycoNetwork  plot_data1=plot_data1[plot_data1['']]
plot_data1 = module6.protein_no_glyco_up[['fc_g','fc_p','normalized_fc_g']]
plot_data1[plot_data1['normalized_fc_g']>1.7].shape[0]
plot_data1['group'] = 'up'
plot_data2 = module6.protein_no_glyco_down[['fc_g','fc_p','normalized_fc_g']]
plot_data2[plot_data2['normalized_fc_g']<0.59].shape[0]
plot_data2['group'] = 'down'
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data.reset_index(inplace=True,drop=True)
module7.scatter(data = plot_data,
             group_column = 'group',
             x_column = 'fc_p',
             y_column = 'fc_g',
             subfolder='StrucGAP_GlycoNetwork',
             top_xaxis_line_show = False,
             right_yaxis_line_show = False,
             xaxis_splitline_show = False,
             yaxis_splitline_show = False,
             )

plot_data = module6.protein_up_glyco_down[['fc_g','fc_p','normalized_fc_g']]
plot_data.reset_index(inplace=True)
module7.scatter(data = plot_data,
             group_column = 'index',
             x_column = 'fc_p',
             y_column = 'fc_g',
             subfolder='StrucGAP_GlycoNetwork',
             top_xaxis_line_show = False,
             right_yaxis_line_show = False,
             xaxis_splitline_show = False,
             yaxis_splitline_show = False,
             filename = 'protein_up_glyco_down'
             )
module7.draw_glycans(module6.protein_up_glyco_down['structure_coding'], 
                     linewidth=0.2, 
                     subfolder='StrucGAP_GlycoNetwork',
                     filename="protein_up_glyco_down")

plot_data = module6.pg_fc[(module6.pg_fc['ProteinID']=='P09055')&(module6.pg_fc['Glycosite_Position']=='406')].copy()
plot_data = plot_data[plot_data['p_g']<0.05]
plot_data['fc'] = plot_data['fc_g']
module7.up_down_scatter(
    [plot_data],
    ['P09055+406'],
    fc_threshold=1,
    show_xaxis=False,
    spine_width=1.2,
    ytick_labelsize=10,
    scatter_size=200,
    scatter_edgecolor='k',
    up_color='red',
    down_color='blue',
    bbox_facecolor='yellow',
    bbox_textsize=11,
    subfolder='StrucGAP_GlycoNetwork',
)

plot_data = module6.pg_fc[(module6.pg_fc['ProteinID']=='P09055')&(module6.pg_fc['Glycosite_Position']=='406')].copy()
plot_data = plot_data[plot_data['p_g']<0.05]
plot_data['fc'] = plot_data['fc_g'] / plot_data['fc_p']
plot_data = plot_data[['structure_coding','fc_g', 'fc_p', 'fc']]
plot_data = plot_data.set_index('structure_coding',drop=True)
plot_data = plot_data.sort_values('fc_g', ascending=False)
module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 centervalue = 1,
                 minvalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 filename = 'P09055_406',
                 subfolder='StrucGAP_GlycoNetwork',
                 )
module7.draw_glycans(plot_data.index[:10], linewidth=1, filename="0_10",subfolder='StrucGAP_GlycoNetwork')
module7.draw_glycans(plot_data.index[10:20], linewidth=1, filename="10_20",subfolder='StrucGAP_GlycoNetwork')
module7.draw_glycans(plot_data.index[20:30], linewidth=1, filename="20_30",subfolder='StrucGAP_GlycoNetwork')
module7.draw_glycans(plot_data.index[30:40], linewidth=1, filename="30_40",subfolder='StrucGAP_GlycoNetwork')


plot_data = module6.glycosyltransferases
plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
plot_data.reset_index(inplace=True)
plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
plot_data.set_index('gene_id',inplace=True)
module7.heatmap2(data = plot_data,
                 colors='Spectral',
                 columns=plot_data.columns[1:11],
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 yaxis_label_show = True,
                 z_score = None,
                 splitline_width = 0.5,
                 filename = 'glycosyltransferases',
                 subfolder='StrucGAP_GlycoNetwork',
                 )

plot_data = module6.cv_filter_data.reset_index()
module7.boxplot(data = plot_data,
                  item_column = 'Accession',
                  item_name = ['Q8BJT9'],
                  group1_columns = plot_data.columns[1:6],
                  group2_columns = plot_data.columns[6:11],
                  p_data = module6.proteomic_fc,
                  p_column = 'pvalue_ttest',
                  filename = 'Q8BJT9',
                  subfolder='StrucGAP_GlycoNetwork',
                  )
module7.boxplot(data = plot_data,
                  item_column = 'Accession',
                  item_name = ['Q6YGZ1'],
                  group1_columns = plot_data.columns[1:6],
                  group2_columns = plot_data.columns[6:11],
                  p_data = module6.proteomic_fc,
                  p_column = 'pvalue_ttest',
                  filename = 'Q6YGZ1',
                  subfolder='StrucGAP_GlycoNetwork',
                  )

plot_data1 = module6.pg_fc
plot_data1 = pd.DataFrame(plot_data1[plot_data1['structure_coding'].str.contains('3')])
plot_data1 = plot_data1[['fc_g', 'fc_p']]
plot_data1['glycosylation'] = plot_data1['fc_g'] / plot_data1['fc_p']
# plot_data1 = plot_data1.T
greater_than_1_5_percent = (plot_data1 > 1.5).sum(axis=0) / plot_data1.shape[0] * 100
plot_data2 = pd.DataFrame({
    "fc>1.5": greater_than_1_5_percent
})

plot_data1 = module6.cv_filter_data[module6.cv_filter_data.index.isin(module6.sialyltransferases.index)]
plot_data2 = module6.sialyltransferases[['fc']]
plot_data3 = module6.sialyltransferases[['pvalue_ttest']]
module7.complexheatmap(data = plot_data1, 
                       columns = plot_data1.columns[:10], 
                       row_annotation_data = [plot_data2, plot_data3], 
                       row_annotation_data_log2 = [False, False],
                       row_annotation_plot_type = ['bar', 'scatter'],
                       col_annotation_data = None, 
                       col_annotation_data_log2 = None,
                       col_annotation_plot_type = None,
                       log2 = True,
                       col_cluster = False,
                       row_cluster = False,
                       subfolder='StrucGAP_GlycoNetwork',
                       col_split=None,
                       cmap = 'Blues', 
                       z_score = 0,
                       show_rownames = True,
                       row_split = None,
                       filename = 'sialyltransferases',
                       linewidths = 3,
                       )

plot_data1 = module6.cv_filter_data[module6.cv_filter_data.index.isin(module6.fucosyltransferase.index)]
plot_data2 = module6.fucosyltransferase[['fc']]
plot_data3 = module6.fucosyltransferase[['pvalue_ttest']]
module7.complexheatmap(data = plot_data1, 
                       columns = plot_data1.columns[:10], 
                       row_annotation_data = [plot_data2, plot_data3], 
                       row_annotation_data_log2 = [False, False],
                       row_annotation_plot_type = ['bar', 'scatter'],
                       col_annotation_data = None, 
                       col_annotation_data_log2 = None,
                       col_annotation_plot_type = None,
                       log2 = True,
                       col_cluster = False,
                       row_cluster = False,
                       subfolder='StrucGAP_GlycoNetwork',
                       col_split=None,
                       cmap = 'Blues', 
                       z_score = 0,
                       show_rownames = True,
                       row_split = None,
                       filename = 'fucosyltransferase',
                       linewidths = 3,
                       )

plot_data = module6.glycan_binding_protein
plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
plot_data.reset_index(inplace=True)
plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
plot_data.set_index('gene_id',inplace=True,drop=False)
module7.violin_plot(data=plot_data,
                    item_column='gene_id',
                    item_name=['Emcn',
                               'Fuom',
                               'Klra4',
                               'Klra7'],
                    group1_columns=['Abundances (Normalized): F1: 126, Control',
                                    'Abundances (Normalized): F1: 127N, Control',
                                    'Abundances (Normalized): F1: 127C, Control',
                                    'Abundances (Normalized): F1: 128N, Control',
                                    'Abundances (Normalized): F1: 128C, Control',],
                    group2_columns=['Abundances (Normalized): F1: 129N, Sample',
                                    'Abundances (Normalized): F1: 129C, Sample',
                                    'Abundances (Normalized): F1: 130N, Sample',
                                    'Abundances (Normalized): F1: 130C, Sample',
                                    'Abundances (Normalized): F1: 131, Sample',],
                    p_data=plot_data,
                    p_column='pvalue_ttest',
                    subfolder='StrucGAP_GlycoNetwork',
                    filename = 'glycan binding protins',
                    ) 

plot_data1 = module6.glycan_binding_protein
plot_data1 = plot_data1[(plot_data1.index=='Q9R0H2')]
plot_data1 = plot_data1.iloc[:,:10]
plot_data2 = module6.sialyltransferases
plot_data2 = plot_data2[(plot_data2.index=='O88829')|(plot_data2.index=='Q11204')]
plot_data2 = plot_data2.iloc[:,:10]
plot_data3 = module4.up_data
plot_data3 = plot_data3[~plot_data3['Ac/Gc'].isnull()]
plot_data3 = plot_data3[['126.1277','127.1248','127.1311','128.1281','128.1344','129.1315','129.1378','130.1348','130.1411','131.1382']]
plot_data3 = pd.DataFrame([plot_data3.mean(axis=0)])
plot_circle_data = [plot_data1, plot_data2, plot_data3]
module7.correlation_heatmap(plot_circle_data,
                            colors = 'Tropic',
                            minvalue = 0.6, 
                            centervalue = 0.8,
                            maxvalue = 1,
                            filename = 'sia correlation',
                            subfolder='StrucGAP_GlycoNetwork',
                            )

plot_data1 = module6.glycan_binding_protein
plot_data1 = plot_data1[(plot_data1.index=='Q8R2K1')]
plot_data1 = plot_data1.iloc[:,:10]
plot_data2 = module6.fucosyltransferase
plot_data2 = plot_data2[(plot_data2.index=='Q9WTS2')]
plot_data2 = plot_data2.iloc[:,:10]
plot_data3 = module4.up_data
plot_data3 = plot_data3[~plot_data3['fucosylated type'].isnull()]
plot_data3 = plot_data3[['126.1277','127.1248','127.1311','128.1281','128.1344','129.1315','129.1378','130.1348','130.1411','131.1382']]
plot_data3 = pd.DataFrame([plot_data3.mean(axis=0)])
plot_circle_data = [plot_data1, plot_data2, plot_data3]
module7.correlation_heatmap(plot_circle_data,
                            colors = 'Tropic',
                            minvalue = 0.7, 
                            centervalue = 0.85,
                            maxvalue = 1,
                            filename = 'fucos correlation',
                            subfolder='StrucGAP_GlycoNetwork',
                            )

glyco = module4.fc_result[module4.fc_result['pvalue_ttest'] < 0.05]
globaldata = module6.proteomic_fc[module6.proteomic_fc['pvalue_ttest'] < 0.05]
module7.up_down_scatter(
    [glyco, globaldata],
    ['Glyco proteome', 'Global proteome'],
    fc_threshold=1,
    show_xaxis=False,
    spine_width=1.2,
    ytick_labelsize=10,
    # ylabel_fontsize=12,
    scatter_size=20,
    scatter_edgecolor='k',
    up_color='red',
    down_color='blue',
    bbox_facecolor='yellow',
    bbox_textsize=11,
    subfolder='StrucGAP_GlycoNetwork',
)

plot_data1 = module4.differential_analysis_data.copy()
plot_data1 = plot_data1[(plot_data1['fc']>1.7)|(plot_data1['fc']<1/1.7)]
plot_data2 = module6.proteomic_fc[module6.proteomic_fc['pvalue_ttest'] < 0.05].copy()
plot_data2 = plot_data2[(plot_data2['fc']>1.5)|(plot_data2['fc']<1/1.5)]
module7.venn_diagram(
    plot_data1["ProteinID"].tolist(),
    plot_data2.index.tolist(),
    colors = 'Tropic',
    subfolder='StrucGAP_GlycoNetwork',
    legend = ['Glycoproteins', 'Proteins'],
    filename='glyco and proteomic protein'
)

plot_data1 = module4.differential_analysis_data.copy()
plot_data1 = plot_data1[(plot_data1['fc']>1.7)|(plot_data1['fc']<1/1.7)]
plot_data1 = plot_data1.iloc[:,:-6]
module5 = StrucGAP_FunctionAnnotation(plot_data1, 
                                  data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
plot_data1 = module5.ora_no_background_both_proteins_result
plot_data1 = plot_data1[plot_data1['Gene_set']=='GO:CC']
plot_data1 = plot_data1[plot_data1['P-value']<0.05]
plot_data1['type'] = 'igp'
plot_data1 = plot_data1.iloc[:10,:]
plot_data2 = module6.proteomic_fc[module6.proteomic_fc['pvalue_ttest'] < 0.05].copy()
plot_data2 = plot_data2[(plot_data2['fc']>1.5)|(plot_data2['fc']<1/1.5)]
plot_data2['Accession'] = plot_data2.index
plot_data2 = module6.convert_accession_to_gene(plot_data2, 'Accession')
plot_data2 = plot_data2.rename(columns={'gene_id':'GeneName'})
plot_data2 = plot_data2.drop(columns=['pvalue_ttest_mannwhitneyu'])
module5 = StrucGAP_FunctionAnnotation(plot_data2, 
                                  data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
plot_data2 = module5.ora_no_background_both_proteins_result
plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:CC']
plot_data2 = plot_data2[plot_data2['P-value']<0.05]
plot_data2['type'] = 'protein'
plot_data2 = plot_data2.iloc[:10,:]
plot_data = pd.concat([plot_data1, plot_data2], axis=0)
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'type',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='igp protein gocc',
    subfolder='StrucGAP_GlycoNetwork',
    col_cluster=False,
)

igps = module4.differential_analysis_data.copy()
igps = igps[igps['pvalue_ttest']<0.05]
igps = igps[(igps['fc']>1.7)|(igps['fc']<1/1.7)]
igps = pd.DataFrame(igps['GeneName'])
proteins = module6.proteomic_fc.copy()
proteins = proteins[proteins['pvalue_ttest']<0.05]
proteins = proteins[(proteins['fc']>1.5)|(proteins['fc']<1/1.5)]
proteins['GeneName'] = proteins.index
proteins = pd.DataFrame(proteins['GeneName'])
plot_data = pd.concat([igps, proteins],axis=0)
module5 = StrucGAP_FunctionAnnotation(plot_data, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', selected_terms=['KEGG'],
            background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest')
plot_data = module5.ora_no_background_both_proteins_result.copy()
plot_data = plot_data[plot_data['P-value']<0.05]
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'Gene_set',
    p_column = 'P-value',
    top = 20,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='igp protein gocc',
    subfolder='StrucGAP_GlycoNetwork',
    col_cluster=False,
    col_split=False
    )



# StrucGAP_InsightTracker
module7.bar_up_down_ratio(feature='core_structure', 
                          colors=['#F3A400', '#D11547', '#274192', '#302A40'],
                          screen_feature = ['A2B2C1D1dD1dcbB5'],
                          subfolder='StrucGAP_InsightTracker',
                          filename="core structure",
                          )

module4.data = module4.data[module4.data['core_structure']=='A2B2C1D1dD1dcbB5']
module5 = StrucGAP_FunctionAnnotation(module4,
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest',
            enrich_feature='glycopeptide') # 65,72,79
plot_data1 = module5.ora_no_background_up_result
plot_data1 = plot_data1[plot_data1['Gene_set']=='GO:MF']
plot_data1['type'] = 'up'
plot_data2 = module5.ora_no_background_down_result
plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:MF']
plot_data2['type'] = 'down'
plot_data = pd.concat([plot_data1, plot_data2], axis=0)
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'type',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='core2',
    subfolder='StrucGAP_InsightTracker',
    col_cluster=False,
)

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

import random
def generate_random_colors(n=20, hex_format=True):
    colors = []
    for _ in range(n):
        r = random.randint(0, 255)
        g = random.randint(0, 255)
        b = random.randint(0, 255)
        if hex_format:
            colors.append('#{:02X}{:02X}{:02X}'.format(r, g, b))
        else:
            colors.append((r, g, b))
    return colors
# 示例：生成20个十六进制颜色
color_list = generate_random_colors(20)

module7.bar_up_down_ratio(feature='branches_structure', 
                          colors=color_list,
                          screen_feature = ['E2F1G5gfF5fe','E2F1G3gfF5fe'],
                          subfolder='StrucGAP_InsightTracker',
                          filename="branches structure",
                          )

module7.bar_up_down_ratio(feature='glycan_type', 
                          colors=generate_random_colors(1),
                          screen_feature = ['Hybrid'],
                          subfolder='StrucGAP_InsightTracker',
                          filename="glycan type",
                          )

module7.bar_up_down_ratio(feature='fucosylated_type', 
                          colors=generate_random_colors(3),
                          subfolder='StrucGAP_InsightTracker',
                          filename="fucosylated type",
                          )

module7.bar_up_down_ratio(feature='acgc', 
                          colors=generate_random_colors(1),
                          screen_feature = ['Ac'],
                          subfolder='StrucGAP_InsightTracker',
                          filename="acgc",
                          )

list1 = []
list2 = []
for i in [1.2, 1.5, 2, 2.5, 3]:
    module2 = StrucGAP_GlycanStructure(gs_data=module4, data_manager=data_manager, 
                              fc = i, pvalue_type='pvalue_ttest', differential_data_type='up')
    module2.statistics(remove_oligo_mannose = False) 
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    export_data = module2.glycan_type_branches_structure.copy()
    export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['Branches_structure']=='E2F1G5gfF5fe')&
                      (export_data['Glycan_type']=='Hybrid')]['Ratio'])
    list1.append(export_data[(export_data['Branches_structure']=='E2F1G5gfF5fe')&
                      (export_data['Glycan_type']=='Hybrid')]['Ratio'].values[0])
    export_data = module2.core_structure_branches_structure.copy()
    export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['Branches_structure']=='E2F1G5gfF5fe')&
                      (export_data['Core_structure']=='A2B2C1D1dD1dcbB5')]['Ratio'])
    list2.append(export_data[(export_data['Branches_structure']=='E2F1G5gfF5fe')&
                      (export_data['Core_structure']=='A2B2C1D1dD1dcbB5')]['Ratio'].values[0])
plot_data = pd.DataFrame()
plot_data['Hybrid_E2F1G5gfF5fe'] = list1
plot_data['A2B2C1D1dD1dcbB5_E2F1G5gfF5fe'] = list2
plot_data = plot_data.T
plot_data.reset_index(inplace=True)
module7.line(data = plot_data,
             colors = ['#E07A5F','#3D405B','#81B29A','#F2CC8F'],
             y_column = 'index',
             x_columns = plot_data.columns[1:],
             subfolder='StrucGAP_InsightTracker',
             symbol_size = 5,
             plot_title = None,
             xaxis_title_gap = 35,
             yaxis_title = 'Percentage (%)',
             filename='modest'
             )

list1 = []
list2 = []
for i in [1.2, 1.5, 2, 2.5, 3]:
    module2 = StrucGAP_GlycanStructure(gs_data=module4, data_manager=data_manager, 
                              fc = i, pvalue_type='pvalue_ttest', differential_data_type='up')
    module2.statistics(remove_oligo_mannose = False) 
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    export_data = module2.fucosylated_type_branches_structure.copy()
    export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['Branches']=='E2F1G5gfF5fe')&
                      (export_data['FucosylatedType']=='dual')]['Ratio'])
    list1.append(export_data[(export_data['Branches']=='E2F1G5gfF5fe')&
                      (export_data['FucosylatedType']=='dual')]['Ratio'].values[0])
    export_data = module2.fucosylated_type_branches_structure.copy()
    export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['FucosylatedType']=='dual')&
                      (export_data['Branches']=='E2F1G3gfF5fe')]['Ratio'])
    list2.append(export_data[(export_data['FucosylatedType']=='dual')&
                      (export_data['Branches']=='E2F1G3gfF5fe')]['Ratio'].values[0])
plot_data = pd.DataFrame()
plot_data['dual fucosylation_E2F1G5gfF5fe'] = list1
plot_data['dual fucosylation_E2F1G3gfF5fe'] = list2
plot_data = plot_data.T
plot_data.reset_index(inplace=True)
module7.line(data = plot_data,
             colors = ['#E07A5F','#3D405B','#81B29A','#F2CC8F'],
             y_column = 'index',
             x_columns = plot_data.columns[1:],
             subfolder='StrucGAP_InsightTracker',
             symbol_size = 5,
             plot_title = None,
             xaxis_title_gap = 35,
             yaxis_title = 'Percentage (%)',
             filename='expected'
             )

list1 = []
list2 = []
list3 = []
list4 = []
for i in [1.2, 1.5, 2, 2.5, 3]:
    module2 = StrucGAP_GlycanStructure(gs_data=module4, data_manager=data_manager, 
                              fc = i, pvalue_type='pvalue_ttest', differential_data_type='up')
    module2.statistics(remove_oligo_mannose = False) 
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    # export_data = module2.acgc_fucosylated_type.copy() 
    export_data = module2.fucosylated_type_acgc.copy()
    export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['FucosylatedType']=='dual')&
                      (export_data['Ac/Gc']=='Ac')]['Ratio'])
    list1.append(export_data[(export_data['FucosylatedType']=='dual')&
                      (export_data['Ac/Gc']=='Ac')]['Ratio'].values[0])
    export_data = module2.core_structure_glycan_type.copy()
    export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['Core_structure']=='A2B2C1D1dD1dcbB5')&
                      (export_data['Glycan_type']=='Hybrid')]['Ratio'])
    list2.append(export_data[(export_data['Core_structure']=='A2B2C1D1dD1dcbB5')&
                      (export_data['Glycan_type']=='Hybrid')]['Ratio'].values[0])
    export_data = module2.glycan_type_core_structure.copy()
    export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['Core_structure']=='A2B2C1D1dD1dcbB5')&
                      (export_data['Glycan_type']=='Hybrid')]['Ratio'])
    list3.append(export_data[(export_data['Core_structure']=='A2B2C1D1dD1dcbB5')&
                      (export_data['Glycan_type']=='Hybrid')]['Ratio'].values[0])
    export_data = module2.core_structure_fucosylated_type.copy()
    export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
    print(export_data[(export_data['Core_structure']=='A2B2C1D1dD1dcbB5')&
                      (export_data['fucosylated type']=='dual')]['Ratio'])
    list4.append(export_data[(export_data['Core_structure']=='A2B2C1D1dD1dcbB5')&
                      (export_data['fucosylated type']=='dual')]['Ratio'].values[0])
    
plot_data = pd.DataFrame()
plot_data['ac_fucosylation'] = list1
plot_data['core2_hybrid'] = list2
plot_data['hybrid_core2'] = list3
plot_data['core2_fucosylation'] = list4
plot_data = plot_data.T
plot_data.reset_index(inplace=True)
module7.line(data = plot_data,
             colors = ['#E07A5F','#3D405B','#81B29A','#F2CC8F'],
             y_column = 'index',
             x_columns = plot_data.columns[1:],
             subfolder='StrucGAP_InsightTracker',
             symbol_size = 5,
             plot_title = None,
             xaxis_title_gap = 35,
             yaxis_title = 'Percentage (%)',
             filename='stronger'
             )


plot_data = pd.read_excel("D:\\doctor\\analysisys\\GAP_manuscript\\Supplementary Tables\\analysis_result\\Supplementary Table 10_StrucGAP_FunctionAnnotation_GO_ora_no_background_up_result_key_information.xlsx",
                          sheet_name = 'bp_core_structure')
plot_data = plot_data.iloc[6:,1:]
plot_data.columns = plot_data.iloc[0]
plot_data = plot_data[1:].reset_index(drop=True)
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
    y_column = "Core_structure",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_InsightTracker',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='bp_core_structure'
)

plot_data = pd.read_excel("D:\\doctor\\analysisys\\GAP_manuscript\\Supplementary Tables\\analysis_result\\Supplementary Table 10_StrucGAP_FunctionAnnotation_GO_ora_no_background_up_result_key_information.xlsx",
                          sheet_name = 'bp_glycan_type')
plot_data = plot_data.iloc[5:,1:]
plot_data.columns = plot_data.iloc[0]
plot_data = plot_data[1:].reset_index(drop=True)
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#F8ED70','#EB8B5B','#B53656','#632A69'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_InsightTracker',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='bp_glycan_type'
)

plot_data = pd.read_excel("D:\\doctor\\analysisys\\GAP_manuscript\\Supplementary Tables\\analysis_result\\Supplementary Table 10_StrucGAP_FunctionAnnotation_GO_ora_no_background_up_result_key_information.xlsx",
                          sheet_name = 'bp_fucosylated_type')
plot_data = plot_data.iloc[5:,1:]
plot_data.columns = plot_data.iloc[0]
plot_data = plot_data[1:].reset_index(drop=True)
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#F8ED70','#EB8B5B','#B53656','#632A69'],
    y_column = "Fucosylated_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_InsightTracker',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='bp_fucosylated_type'
)

plot_data = pd.read_excel("D:\\doctor\\analysisys\\GAP_manuscript\\Supplementary Tables\\analysis_result\\Supplementary Table 10_StrucGAP_FunctionAnnotation_GO_ora_no_background_up_result_key_information.xlsx",
                          sheet_name = 'bp_acgc')
plot_data = plot_data.iloc[5:,1:]
plot_data.columns = plot_data.iloc[0]
plot_data = plot_data[1:].reset_index(drop=True)
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#F8ED70','#EB8B5B','#B53656','#632A69'],
    y_column = "Ac/Gc",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_InsightTracker',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='bp_acgc'
)

plot_data = pd.read_excel("D:\\doctor\\analysisys\\GAP_manuscript\\Supplementary Tables\\analysis_result\\Supplementary Table 11_StrucGAP_FunctionAnnotation_GO_ora_no_background_down_result_key_information.xlsx",
                          sheet_name = 'bp_core_structure')
plot_data = plot_data.iloc[6:,1:]
plot_data.columns = plot_data.iloc[0]
plot_data = plot_data[1:].reset_index(drop=True)
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
    y_column = "Core_structure",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_InsightTracker',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='down bp_core_structure'
)

plot_data = pd.read_excel("D:\\doctor\\analysisys\\GAP_manuscript\\Supplementary Tables\\analysis_result\\Supplementary Table 11_StrucGAP_FunctionAnnotation_GO_ora_no_background_down_result_key_information.xlsx",
                          sheet_name = 'bp_acgc')
plot_data = plot_data.iloc[5:,1:]
plot_data.columns = plot_data.iloc[0]
plot_data = plot_data[1:].reset_index(drop=True)
module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#F8ED70','#EB8B5B','#B53656','#632A69'],
    y_column = "Ac/Gc",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_InsightTracker',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 4,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage (%)',
    filename='down bp_fucosylated_type'
)

plot_data = module4.up_data.copy()
plot_data = plot_data[(plot_data['Ac/Gc']=='Ac')&
                                  (plot_data['fucosylated type']=='dual')]
module5 = StrucGAP_FunctionAnnotation(plot_data, data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.2,pvalue_type='pvalue_ttest') # 65,72,79
plot_data = module5.ora_no_background_up_result
plot_data = plot_data[plot_data['Gene_set']=='GO:BP']
plot_data = plot_data.iloc[:15,:]
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'Gene_set',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='Ac dual fucosylated',
    subfolder='StrucGAP_InsightTracker',
    col_cluster=False,
    col_split = False
)

plot_data = module4.up_data.copy()
plot_data = plot_data[(plot_data['core_structure']=='A2B2C1D1dD1dcbB5')&
                                  (plot_data['Glycan_type']=='Hybrid')]
module5 = StrucGAP_FunctionAnnotation(plot_data, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.2,pvalue_type='pvalue_ttest') # 65,72,79
plot_data = module5.ora_no_background_up_result
plot_data = plot_data[plot_data['Gene_set']=='GO:BP']
plot_data = plot_data.iloc[:15,:]
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'Gene_set',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='core 2 hybrid',
    subfolder='StrucGAP_InsightTracker',
    col_cluster=False,
    col_split = False
)

plot_data = module4.up_data.copy()
plot_data = plot_data[(plot_data['core_structure']=='A2B2C1D1dD1dcbB5')&
                                  (plot_data['fucosylated type']=='dual')]
module5 = StrucGAP_FunctionAnnotation(plot_data, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.2,pvalue_type='pvalue_ttest') # 65,72,79
plot_data = module5.ora_no_background_up_result
plot_data = plot_data[plot_data['Gene_set']=='GO:BP']
plot_data = plot_data.iloc[:15,:]
module7.dotplot_col(
    data = plot_data,  
    dot_cmap = 'Hawaii',
    category = 'Gene_set',
    p_column = 'P-value',
    top = 10,
    term = 'Term',
    dot_color_column = 'P-value',
    dot_size_column = 'Overlap',
    filename='core 2 dual fucosylated',
    subfolder='StrucGAP_InsightTracker',
    col_cluster=False,
    col_split = False
)

network = data_manager.StrucGAP_GlycoNetwork_key_information.copy()
module7.glyconetwork(network, 
                     'Up', 
                     'Sialylation',
                     figure_description = 'Sialylation network',
                     )




