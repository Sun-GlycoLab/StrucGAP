# -*- coding: utf-8 -*-

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

# StrucGAP_GlycanStructure_1
plot_data = pd.DataFrame(module2.glycan_composition_isoforms['GlycanComposition'].value_counts())
plot_data = pd.DataFrame(plot_data['GlycanComposition'].value_counts())
plot_data['isomer_count'] = plot_data.index
plot_data = plot_data.sort_values('isomer_count', ascending=False)
fig1 = module7.polar1('plot_data', 
              columns=['isomer_count'],
              number_column = 'GlycanComposition',
              subfolder='StrucGAP_GlycanStructure_1',
              if_unique = False,
              radiusaxis_label_show = True,
              colors = ['#81B29A'],
              plot_title = None,
              filename='isomer',
              radiusaxis_label_font_size = 20,
              angleaxis_label_font_size = 20,
              legend_font_size = 20,
              legend = 'Isomer counts',
              figure_description = 'Overveiw of glycan structure isomers identified.',
              )

module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = False) 
data1 = module2.structure_coding_rank
module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = True)
data2 = module2.structure_coding_rank 
fig2 = module7.butterfly_plot(data1,data2,
                       item_column='Structure_coding',
                       count_column='Structure_coding_count',
                       colors = ['lightblue', 'darkblue'],
                       xaxis_title = 'Number of Glycan',
                       plot_title = 'Top 10 glycans',
                       legend = ['With oligo mannose', 'Without oligo mannose'],
                       subfolder='StrucGAP_GlycanStructure_1',
                       filename='top 10 structure coding',
                       label_font_size = 30,
                       xaxis_title_font_size = 30,
                       plot_title_font_size = 0,
                       legend_fontsize = 30,
                       legend_loc=[0,-0.4], 
                       figure_description = 'Top 10 glycan structures identified based on the number of their modified N-glycosites.',
                       )

module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
module2.statistics(remove_oligo_mannose = False) 
module2.structure_statistics()
module2.lacdinac()
module2.cor()
module2.isoforms()

plot_data = module2.core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig3 = module7.pie(data = plot_data, 
              item_column = 'Core_structure', 
              number_column = 'Core_structure_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_1',
              filename='corestructure',
              label_font_size = 20,
              legend_font_size = 20,
              figure_description = 'Distribution of core structures among all unique IGPs, with Core-I (N2H3), Core-II (core-fucosylated), Core-III (bisected), and Core-IV (dual features) quantified.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

fig4 = module7.pie(data = module2.glycan_type, 
              item_column = 'Glycan_type', 
              number_column = 'Glycan_type_count',
              radius = ['0%', '70%'],
              rosetype = 'area',
              subfolder='StrucGAP_GlycanStructure_1',
              filename='glycantype',
              label_font_size = 20,
              legend_font_size = 20,
              figure_description = 'Classification of glycan types (complex, hybrid, and oligo-mannose) across all unique IGPs.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

fig5 = module7.bar(data = module2.branches_structure,
            top = 10,
            y_column='Branches',
            y_column_value='Branches_count',
            subfolder='StrucGAP_GlycanStructure_1',
            colors = '#8FBFB8',
            transform_ratio = True,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Branches structure',
            xaxis_title = 'Branches structure',
            xaxis_title_gap = 150,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Percentage of total branch structures',
            figure_description = 'Relative abundance of top 10 detected glycan branch structures.',
            yaxis_title_gap = 70,
            xaxis_label_font_size = 30,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='branches structure'
            )

fig6 = module7.bar(data = module2.branches_structure,
            top = None,
            end = 8,
            y_column='Branches',
            y_column_value='Branches_count',
            subfolder='StrucGAP_GlycanStructure_1',
            colors = '#8FBFB8',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Branches structure',
            xaxis_title = 'Branches structure',
            xaxis_title_gap = 180,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Count of IGPs',
            figure_description = 'Absolute counts for the eight least abundant structures.',
            yaxis_title_gap = 50,
            xaxis_label_font_size = 30,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='end8 branches structure'
            )

fig7 = module7.bar(data = module2.branches_count,
            top = 10,
            y_column='BranchNumber',
            y_column_value='BranchNumber_count',
            xaxis_label_rotate = 0,
            xaxis_label_margin = 20,
            transform_ratio = True,
            subfolder='StrucGAP_GlycanStructure_1',
            colors = '#FFD804',
            y_max = None,
            yaxis_splitline_show = False,
            legend = 'Branch number',
            filename='branch number',
            xaxis_title = 'Branch number per glycan',
            xaxis_title_gap = 50,
            yaxis_title = 'Ratio of total IGPs by glycan branch number',
            yaxis_title_gap = 60,
            xaxis_label_font_size = 30,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            figure_description = 'Branch count distribution among identified glycans, with frequencies shown per branch number category.',
            )

fig8 = module7.pie(data = module2.sialicacid_count, 
              item_column = 'Sialicacid_count', 
              number_column = 'Number',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_1',
              filename='sialicacid_count',
              label_font_size = 20,
              legend_font_size = 20,
              figure_description = 'The distribution and proportion of IGPs containing different numbers of sialic acid.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

fig9 = module7.pie(data = module2.fucose_count, 
              item_column = 'Fucose_count', 
              number_column = 'Number',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_1',
              filename='fucose_count',
              label_font_size = 20,
              legend_font_size = 20,
              figure_description = 'The distribution and proportion of IGPs containing different numbers of fucose.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

fig10 = module7.pie(data = module2.acgc, 
              item_column = 'Ac/Gc', 
              number_column = 'Ac/Gc_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_1',
              filename='acgc',
              label_font_size = 20,
              legend_font_size = 20,
              figure_description = 'The distribution and proportion of IGPs containing different sialylated types.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

fig11 = module7.pie(data = module2.lacdinac_count, 
              item_column = 'lacdinac', 
              number_column = 'lacdinac_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_1',
              filename='lacdinac',
              label_font_size = 20,
              legend_font_size = 20,
              figure_description = 'The distribution and proportion of IGPs containing different LacdiNAc.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig2, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig3, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig4, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig5, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig6, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_1")
module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_1")
module7.compose_figures("StrucGAP_GlycanStructure_1.pdf", figure_name="StrucGAP_GlycanStructure_1",
                        custom_sizes=[[1], [2,3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycanStructure_2
fig1 = module7.pie(data = module2.fucosylated_type, 
              item_column = 'FucosylatedType', 
              number_column = 'FucosylatedType_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_2',
              filename='fucosylated_type',
              label_font_size = 20,
              legend_font_size = 12,
              figure_description = 'The distribution and proportion of IGPs containing different fucosylated types.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

fig2 = module7.pie(data = module2.fsg, 
              item_column = 'FSG', 
              number_column = 'FSG_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='StrucGAP_GlycanStructure_2',
              filename='fsg',
              label_font_size = 20,
              legend_font_size = 12,
              figure_description = 'The distribution and proportion of IGPs containing different fucosylation/sialylation patterns.',
              colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF']
              )

data1 = module1.data_psm_filtered
data1 = data1[(data1['fucosylated type']=='core fucosylated')|(data1['fucosylated type']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data1 = data1.iloc[:,0]
data2 = module1.data_psm_filtered
data2 = data2[(data2['fucosylated type']=='antenna fucosylated')|(data2['fucosylated type']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data2 = data2.iloc[:,0]
fig3 = module7.venn_diagram(
    data1,
    data2,
    colors = 'Tropic',
    subfolder='StrucGAP_GlycanStructure_2',
    legend = ['Core fucosylated', 'Antenna fucosylated'],
    filename='fucosylated type',
    figure_description = 'Composition of fucosylated type.',
    plot_title_font_size = 25,
    legend_fontsize = 20,
    number_fontsize = 25,
    legend_loc='lower center', 
)

data1 = module1.data_psm_filtered
data1 = data1[(data1['Ac/Gc']=='Ac')|(data1['Ac/Gc']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data1 = data1.iloc[:,0]
data2 = module1.data_psm_filtered
data2 = data2[(data2['Ac/Gc']=='Gc')|(data2['Ac/Gc']=='dual')]['PeptideSequence+structure_coding+ProteinID']
data2 = data2.iloc[:,0]
fig4 = module7.venn_diagram(
    data1,
    data2,
    colors = 'Tropic',
    subfolder='StrucGAP_GlycanStructure_2',
    legend = ['Neu5Ac', 'Neu5Gc'],
    filename='acgc',
    figure_description = 'Composition of sialylated type.',
    plot_title_font_size = 25,
    legend_fontsize = 20,
    number_fontsize = 25,
    legend_loc='lower center', 
)

plot_data = module2.lacdinac_core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig5 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'Core_structure',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_core_structure',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and core structures.',
                 )

plot_data = module2.lacdinac_glycan_type
fig6 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'Glycan_type',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_glycan_type',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and glycan types.',
                 )

plot_data = module2.lacdinac_branches_structure
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig7 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_2',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and branch structures.',
                 )
# fig7 = module7.sunburst(data = plot_data,
#                  root_column = 'Lacdinac',
#                  child_column = 'Branches',
#                  child_column_value = 'Count',
#                  hide_root_label=False,
#                  filename = 'lacdinac_branches_structure',
#                  subfolder='StrucGAP_GlycanStructure_2',
#                  colors = generate_random_colors(13),
#                  figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and branch structures.',
#                  )

plot_data = module2.lacdinac_branches_count
fig8 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'Branch_Number',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_branches_count',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and branch counts.',
                 )

plot_data = module2.lacdinac_sialicacid_count
fig9 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'Sialicacid_count',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_sialicacid_count',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and sialic acid counts.',
                 )

plot_data = module2.lacdinac_fucose_count
fig10 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'Fucose_count',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_fucose_count',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and fucose counts.',
                 )

plot_data = module2.lacdinac_fucosylated_type
fig11 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'fucosylated type',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_fucosylated_type',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and fucosylated types.',
                 )

plot_data = module2.lacdinac_acgc
fig12 = module7.sunburst(data = plot_data,
                 root_column = 'Lacdinac',
                 child_column = 'Ac/Gc',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'lacdinac_acgc',
                 subfolder='StrucGAP_GlycanStructure_2',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between LacdiNAc and sialylated types.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig2, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig3, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig4, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig5, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig6, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_2")
module7.add_figure(fig12, figure_name="StrucGAP_GlycanStructure_2")
module7.compose_figures("StrucGAP_GlycanStructure_2.pdf", figure_name="StrucGAP_GlycanStructure_2",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycanStructure_3
plot_data = module2.core_structure_lacdinac
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig1 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and LacdiNAc.',
                 )

plot_data = module2.core_structure_glycan_type
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig2 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and glycan types.',
                 )

plot_data = module2.core_structure_branches_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig3 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and branch structures.',
                 )

plot_data = module2.core_structure_branches_count
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig4 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and branch counts.',
                 )

plot_data = module2.core_structure_sialicacid_count
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig5 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and sialic acid counts.',
                 )

plot_data = module2.core_structure_fucose_count
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig6 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and fucose counts.',
                 )

plot_data = module2.core_structure_fucosylated_type
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig7 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and fucosylated types.',
                 )

plot_data = module2.core_structure_acgc
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig8 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between core structures and sialylated types.',
                 )

plot_data = module2.glycan_type_lacdinac
fig9 = module7.sunburst(data = plot_data,
                 root_column = 'Glycan_type',
                 child_column = 'Lacdinac',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'glycan_type_lacdinac',
                 subfolder='StrucGAP_GlycanStructure_3',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and LacdiNAc.',
                 )

plot_data = module2.glycan_type_core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig10 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and core structures.',
                 )

plot_data = module2.glycan_type_branches_structure
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig11 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and branch structures.',
                 )

plot_data = module2.glycan_type_branches_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig12 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_3',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and branch counts.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig2, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig3, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig4, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig5, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig6, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_3")
module7.add_figure(fig12, figure_name="StrucGAP_GlycanStructure_3")
module7.compose_figures("StrucGAP_GlycanStructure_3.pdf", figure_name="StrucGAP_GlycanStructure_3",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycanStructure_4
plot_data = module2.glycan_type_sialicacid_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig1 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and sialic acid counts.',
                 )

plot_data = module2.glycan_type_fucose_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig2 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and fucose counts.',
                 )

plot_data = module2.glycan_type_fucosylated_type
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig3 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and fucosylated types.',
                 )

plot_data = module2.glycan_type_acgc
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig4 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between glycan types and sialylated types.',
                 )

random_colors = generate_random_colors(30)

plot_data = module2.branches_structure_core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig5 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch structures and core structures.',
                 )

plot_data = module2.branches_structure_glycan_type
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig6 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch structures and glycan types.',
                 )

plot_data = module2.branches_structure_fucosylated_type
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig7 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch structures and fucosylated types.',
                 )

plot_data = module2.branches_structure_acgc
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig8 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch structures and sialylated types.',
                 )

plot_data = module2.branches_structure_sialicacid_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig9 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch structures and sialic acid counts.',
                 )

plot_data = module2.branches_structure_fucose_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig10 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch structures and fucose counts.',
                 )

plot_data = module2.branches_count_sialicacid_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig11 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch counts and sialic acid counts.',
                 )

plot_data = module2.branches_count_fucose_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig12 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_4',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between branch counts and fucose counts.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig2, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig3, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig4, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig5, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig6, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_4")
module7.add_figure(fig12, figure_name="StrucGAP_GlycanStructure_4")
module7.compose_figures("StrucGAP_GlycanStructure_4.pdf", figure_name="StrucGAP_GlycanStructure_4",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycanStructure_5
plot_data = module2.fucosylated_type_lacdinac
fig1 = module7.sunburst(data = plot_data,
                 root_column = 'FucosylatedType',
                 child_column = 'Lacdinac',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'fucosylated_type_lacdinac',
                 subfolder='StrucGAP_GlycanStructure_5',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and LacdiNAc.',
                 )

plot_data = module2.fucosylated_type_core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig2 = module7.sunburst(data = plot_data,
                 root_column = 'FucosylatedType',
                 child_column = 'Core_structure',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'fucosylated_type_core_structure',
                 subfolder='StrucGAP_GlycanStructure_5',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and core structures.',
                 )

plot_data = module2.fucosylated_type_glycan_type
fig3 = module7.sunburst(data = plot_data,
                 root_column = 'FucosylatedType',
                 child_column = 'Glycan_type',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'fucosylated_type_glycan_type',
                 subfolder='StrucGAP_GlycanStructure_5',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and glycan types.',
                 )

plot_data = module2.fucosylated_type_branches_structure
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig4 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_5',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and branch structures.',
                 )

plot_data = module2.fucosylated_type_branches_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig5 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_5',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and branch counts.',
                 )

plot_data = module2.fucosylated_type_sialicacid_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig6 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_5',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and sialic acid counts.',
                 )

plot_data = module2.fucosylated_type_fucose_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig7 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_5',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and fucose counts.',
                 )

plot_data = module2.fucosylated_type_acgc
fig8 = module7.sunburst(data = plot_data,
                 root_column = 'FucosylatedType',
                 child_column = 'Ac/Gc',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'fucosylated_type_acgc',
                 subfolder='StrucGAP_GlycanStructure_5',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between fucosylated types and sialylated types.',
                 )

plot_data = module2.acgc_lacdinac
fig9 = module7.sunburst(data = plot_data,
                 root_column = 'Ac/Gc',
                 child_column = 'Lacdinac',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'acgc_lacdinac',
                 subfolder='StrucGAP_GlycanStructure_5',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and LacdiNAc.',
                 )

plot_data = module2.acgc_core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig10 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_5',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and core structures.',
                 )

plot_data = module2.acgc_glycan_type
fig11 = module7.sunburst(data = plot_data,
                 root_column = 'Ac/Gc',
                 child_column = 'Glycan_type',
                 child_column_value = 'Count',
                 hide_root_label=False,
                 filename = 'acgc_glycan_type',
                 subfolder='StrucGAP_GlycanStructure_5',
                 colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and glycan types.',
                 )

plot_data = module2.acgc_branches_structure
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig12 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_5',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and branch structures.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig2, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig3, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig4, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig5, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig6, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_5")
module7.add_figure(fig12, figure_name="StrucGAP_GlycanStructure_5")
module7.compose_figures("StrucGAP_GlycanStructure_5.pdf", figure_name="StrucGAP_GlycanStructure_5",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycanStructure_6
plot_data = module2.acgc_branches_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig1 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_6',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and branch structures.',
                 )

plot_data = module2.acgc_sialicacid_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig2 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_6',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and sialic acid counts.',
                 )

plot_data = module2.acgc_fucose_count
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig3 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_6',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and fucose counts.',
                 )

plot_data = module2.acgc_fucosylated_type
plot_data = plot_data.pivot_table(index=plot_data.columns[0], columns=plot_data.columns[1], values=plot_data.columns[2], fill_value=0).T
fig4 = module7.heatmap2(data = plot_data,
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 yaxis_label_show=True, 
                 splitline_width = 1,
                 subfolder='StrucGAP_GlycanStructure_6',
                 figure_description = 'Co-occurrence analysis revealing preferential pairing between sialylated types and fucosylated types.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycanStructure_6")
module7.add_figure(fig2, figure_name="StrucGAP_GlycanStructure_6")
module7.add_figure(fig3, figure_name="StrucGAP_GlycanStructure_6")
module7.add_figure(fig4, figure_name="StrucGAP_GlycanStructure_6")
# module7.add_figure(fig5, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig6, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig12, figure_name="StrucGAP_GlycanStructure_5")
module7.compose_figures("StrucGAP_GlycanStructure_6.pdf", figure_name="StrucGAP_GlycanStructure_6",
                        custom_sizes=[[1], [2], [3], [4]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoSite_1
fig1 = module7.nested_pie(data = module3.glycoprotein_glycosite_count, 
                   item_column = 'glycoprotein', 
                   number_column = 'glycosite_count',
                   value_counts_column = 'glycosite_count',
                   subfolder='StrucGAP_GlycoSite_1',
                   inner_data_first=False,
                   label_font_size = 15,
                   legend_font_size = 15,
                   split=9,
                   filename = 'glycosite_count',
                   figure_description = 'Quantification of glycosylation site distributions across glycoproteins, illustrating overall site occupancy.',
                   )

fig2 = module7.funnel(data = module3.glycoprotein_glycan_count, 
              item_column = 'gene_name', 
              number_column = 'glycan_count',
              top = 10,
              label_font_size = 15,
              legend_font_size = 15,
              subfolder='StrucGAP_GlycoSite_1',
              filename = 'glycan type',
              figure_description = 'Glycan count per glycoprotein based on glycosite.',
              )

fig3 = module7.bar(data = module3.glycoprotein_glycan_type,
            top = 10,
            y_column='gene_name',
            y_column_value='glycan_type',
            subfolder='StrucGAP_GlycoSite_1',
            colors = '#159947',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Count of different glycan types',
            xaxis_title = 'Glycoprotein',
            xaxis_title_gap = 140,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Glycan type number',
            figure_description = 'Glycan type per glycoprotein based on glycosite.',
            yaxis_title_gap = 70,
            xaxis_label_font_size = 12,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='glycoprotein_glycan_type'
            )

plot_data = module3.protein_glycosite_glycan_count.copy()
plot_data['Protein+glycosite'] = plot_data['gene_name'] + "+" + plot_data['Glycosite_Position']
fig4 = module7.bar(data = plot_data,
            top = 10,
            y_column='Protein+glycosite',
            y_column_value='glycan_type_count',
            subfolder='StrucGAP_GlycoSite_1',
            colors = '#159947',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Count of different glycan types',
            xaxis_title = 'Glycoprotein+glycosite',
            xaxis_title_gap = 140,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Glycan type number',
            figure_description = 'Glycan type count per glycosite based on glycoprotein.',
            yaxis_title_gap = 70,
            xaxis_label_font_size = 12,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='protein_glycosite_glycan_count'
            )

plot_data = module3.protein_glycosite_glycan_composition_count.copy()
plot_data['Protein+glycosite+composition'] = plot_data['gene_name'] + "+" + plot_data['Glycosite_Position'] + "+" + plot_data['GlycanComposition']
fig5 = module7.bar(data = plot_data,
            top = 10,
            y_column='Protein+glycosite+composition',
            y_column_value='count',
            subfolder='StrucGAP_GlycoSite_1',
            colors = '#159947',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Count of different glycan compositions',
            xaxis_title = 'Glycoprotein+glycosite+glycan composition',
            xaxis_title_gap = 140,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Glycan composition number',
            figure_description = 'Glycan composition count per glycosite based on glycoprotein.',
            yaxis_title_gap = 70,
            xaxis_label_font_size = 12,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='protein_glycosite_glycan_composition_count'
            )

plot_data = module3.protein_glycosite_isoforms_count.copy()
plot_data['Protein+glycosite+composition+isomers'] = plot_data['gene_name'] + "+" + plot_data['Glycosite_Position'] + "+" + plot_data['GlycanComposition'] + "+" + plot_data['isoforms']
fig6 = module7.bar(data = plot_data,
            top = 10,
            y_column='Protein+glycosite+composition+isomers',
            y_column_value='count',
            subfolder='StrucGAP_GlycoSite_1',
            colors = '#159947',
            transform_ratio = False,
            y_max = None,
            yaxis_splitline_show = False,
            xaxis_splitline_show = False,
            legend = 'Count of isomers',
            xaxis_title = 'Glycoprotein+glycosite+glycan composition+isomers',
            xaxis_title_gap = 140,
            xaxis_label_rotate = -45,
            xaxis_label_margin = 25,
            xaxis_label_text_split = 20,
            yaxis_title = 'Isomer number',
            figure_description = 'Glycan isomer count per glycosite based on glycoprotein.',
            yaxis_title_gap = 70,
            xaxis_label_font_size = 12,
            yaxis_label_font_size = 30,
            legend_font_size = 30,
            filename='protein_glycosite_isoforms_count'
            )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoSite_1")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoSite_1")
# module7.add_figure(fig7, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig8, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig9, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig10, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig11, figure_name="StrucGAP_GlycanStructure_5")
# module7.add_figure(fig12, figure_name="StrucGAP_GlycanStructure_5")
module7.compose_figures("StrucGAP_GlycoSite_1.pdf", figure_name="StrucGAP_GlycoSite_1",
                        custom_sizes=[[1], [2], [3], [4], [5], [6]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoPeptideQuant_1
fig1 = module7.volcano_plot(data = module4.fc_result,
                    fc_column = 'fc',
                    p_column = 'pvalue_ttest',
                    fc = 1.5,
                    p_value = 0.05,
                    subfolder='StrucGAP_GlycoPeptideQuant_1',
                    figure_description = 'Volcano plot of differentially expressed glycopeptides in the dataset, showing significant up- and downregulated IGPs (FC > 1.5 or < 0.67, P value < 0.05).',
                    )

fig2 = module7.dimension_reduction(
    data=module4.data_quant,
    data_columns=['126.1277', '127.1248', '127.1311', '128.1281', '128.1344',
                  '129.1315', '129.1378', '130.1348', '130.1411', '131.1382'],
    sample_group=module4.sample_group,
    filter_data=module4.fc_result, # filter_data=module4.fc_result,
    p_column='pvalue_ttest',
    p_value=0.05,
    fc = 1.5,
    method='pca',
    dimension_number = 2,
    random_state = 0,
    colors = ['#3558AE', '#B64074'],
    subfolder='StrucGAP_GlycoPeptideQuant_1',
    figure_description = 'PCA plot based on glycopeptide expression profiles.',
)

fig3 = module7.heatmap2(data = module4.data_quant,
                 columns=['126.1277','127.1248','127.1311','128.1281','128.1344', 
                          '129.1315','129.1378','130.1348','130.1411','131.1382'],
                 filter_data = module4.fc_result,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = 0,
                 splitline_width = 0.0000001,
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Heatmap of differentially expressed glycopeptides, revealing global expression trends and group-wise clustering.',
                 )

up_data = module4.up_data[module4.up_data['Glycan_type']!='Oligo mannose']
down_data = module4.down_data[module4.down_data['Glycan_type']!='Oligo mannose']
fig4 = module7.heatmap_multi_data(up_data,  
              down_data,   
              columns=['structure_coding'],
              statistical_methods = ['count'],
              subfolder='StrucGAP_GlycoPeptideQuant_1',
              colors = 'coolwarm',
              filename = 'top 10 differential glycan',
              figure_description = 'Top 10 most frequently observed upregulated and downregulated glycan structures, excluding oligo-mannose types, highlighting distinct compositional shifts.',
              )

plot_data = module4.roc_result
plot_data = plot_data[plot_data['auc_pvalue']<0.05]
plot_data = plot_data.sort_values('auc_pvalue', ascending=True)
plot_data = plot_data.sort_values('auc', ascending=False)
plot_data = plot_data.iloc[:10,:]
fig5 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 minvalue=0, 
                 centervalue = 0.5,
                 maxvalue=1,
                 splitline_width = 0.5,
                 filename = 'top 10 differential igps auc',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Top 10 differential IGPs based on auc.',
                 )

plot_data = module4.ml_result
plot_data = plot_data.sort_values('randomforest_features_importance_means', ascending=False)
plot_data = pd.DataFrame(plot_data.iloc[:10,0])
fig6 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 filename = 'top 10 differential igps randomforest',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Top 10 differential IGPs based on random forest-based features importance score.',
                 )

plot_data = module4.ml_result
plot_data = plot_data.sort_values('xgbclassifier_features_importance_means', ascending=False)
plot_data = pd.DataFrame(plot_data.iloc[:10,0])
fig7 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 filename = 'top 10 differential igps xgbclassifier',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Top 10 differential IGPs based on XGBClassifier-based features importance score.',
                 )

plot_data = module4.pca_result
plot_data = plot_data.sort_values('pca_features_importance', ascending=False)
plot_data = pd.DataFrame(plot_data.iloc[:10,0])
fig8 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 filename = 'top 10 differential igps pca',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Top 10 differential IGPs based on pca-based features importance score.',
                 )

plot_data = module4.anova_result
plot_data = plot_data[plot_data['anova_pvalue']<0.05]
plot_data = plot_data.sort_values('anova_pvalue', ascending=True)
plot_data = plot_data.sort_values('f score', ascending=False)
plot_data = plot_data.iloc[:10,:]
fig9 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 filename = 'top 10 differential igps anova',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Top 10 differential IGPs based on anova-based f score.',
                 )

plot_data = module4.chi2_result
plot_data = plot_data[plot_data['chi2_pvalue']<0.05]
plot_data = plot_data.sort_values('chi2_pvalue', ascending=True)
plot_data = plot_data.sort_values('chi2 score', ascending=False)
plot_data = plot_data.iloc[:10,:]
fig10 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 cluster = None,
                 yaxis_label_show = True,
                 z_score = None,
                 minvalue = None, 
                 centervalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 filename = 'top 10 differential igps chi2',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Top 10 differential IGPs based on chi2-based chi2 score.',
                 )

fig11 = module7.bar_up_down(data = module4.differential_analysis_overview,
                 x_column = 'item',
                 up_column = 'Up_data_item_count',
                 down_column = 'Down_data_item_count',
                 xaxis_label_margin = 25,
                 colors = ['#FF4359', '#0078FF'],
                 filename = 'differential_analysis_overview',
                 subfolder='StrucGAP_GlycoPeptideQuant_1',
                 figure_description = 'Overview of identified features from differential analysis.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig10, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.add_figure(fig11, figure_name="StrucGAP_GlycoPeptideQuant_1")
# module7.add_figure(fig12, figure_name="StrucGAP_GlycoPeptideQuant_1")
module7.compose_figures("StrucGAP_GlycoPeptideQuant_1.pdf", figure_name="StrucGAP_GlycoPeptideQuant_1",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoPeptideQuant_2
plot_data = module4.differential_analysis_core_structure
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig1 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'Core_structure',
                 up_column = 'Up_data_Core_structure_count',
                 down_column = 'Down_data_Core_structure_count',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_core_structure",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of core structures',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs counts',
                 figure_description = 'Comparison of core structures in up- versus downregulated IGPs.',
                 )

fig2 = module7.bar_up_down(data = plot_data,
                 x_column = 'Core_structure',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_core_structure ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of core structures',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of core structures in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_glycan_type
fig3 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'Glycan_type',
                 up_column = 'Up_data_Glycan_type_count',
                 down_column = 'Down_data_Glycan_type_count',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_glycan_type",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of glycan types',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs counts',
                 figure_description = 'Comparison of glycan types in up- versus downregulated IGPs.',
                 )

fig4 = module7.bar_up_down(data = plot_data,
                 x_column = 'Glycan_type',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_glycan_type ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of glycan types',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of glycan types in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_branches_structure
fig5 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'Branches',
                 up_column = 'Up_data_Branches_count',
                 down_column = 'Down_data_Branches_count',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_branches_structure",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of branch structures',
                 xaxis_label_margin = 4,
                 xaxis_title_gap = 45,
                 xaxis_label_rotate = -45,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 8,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs counts',
                 figure_description = 'Comparison of branch structures in up- versus downregulated IGPs.',
                 )

fig6 = module7.bar_up_down(data = plot_data,
                 x_column = 'Branches',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_branches_structure ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of branch structures',
                 xaxis_label_margin = 4,
                 xaxis_title_gap = 45,
                 xaxis_label_rotate = -45,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 8,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of branch structures in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_branches_count
fig7 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'BranchNumber',
                 up_column = 'Up_data_BranchNumber_count',
                 down_column = 'Down_data_BranchNumber_count',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_branches_count",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'Number of branches per glycan',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs count',
                 figure_description = 'Comparison of branch number in up- versus downregulated IGPs.',
                 )

fig8 = module7.bar_up_down(data = plot_data,
                 x_column = 'BranchNumber',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_branches_count ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'Number of branches per glycan',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of branch number in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_sialicacid_count
fig9 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'Sialicacid_count',
                 up_column = 'Up_data_Number',
                 down_column = 'Down_data_Number',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_sialicacid_count",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'Number of sialic acids per glycan',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs count',
                 figure_description = 'Comparison of sialic acid count in up- versus downregulated IGPs.',
                 )

fig10 = module7.bar_up_down(data = plot_data,
                 x_column = 'Sialicacid_count',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_sialicacid_count ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'Number of sialic acids per glycan',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of sialic acid count in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_fucose_count
fig11 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'Fucose_count',
                 up_column = 'Up_data_Number',
                 down_column = 'Down_data_Number',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_fucose_count",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'Number of fucoses per glycan',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs count',
                 figure_description = 'Comparison of fucose count in up- versus downregulated IGPs.',
                 )

fig12 = module7.bar_up_down(data = plot_data,
                 x_column = 'Fucose_count',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_2',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_fucose_count ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'Number of fucoses per glycan',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of fucose count in up- versus downregulated IGPs based on total number of related glycan.',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig10, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig11, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.add_figure(fig12, figure_name="StrucGAP_GlycoPeptideQuant_2")
module7.compose_figures("StrucGAP_GlycoPeptideQuant_2.pdf", figure_name="StrucGAP_GlycoPeptideQuant_2",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoPeptideQuant_3
plot_data = module4.differential_analysis_lacdinac
fig1 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'lacdinac',
                 up_column = 'Up_data_lacdinac',
                 down_column = 'Down_data_lacdinac',
                 subfolder='StrucGAP_GlycoPeptideQuant_3',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_lacdinac",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of LacdiNAc',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs count',
                 figure_description = 'Comparison of LacdiNAc in up- versus downregulated IGPs.',
                 )

fig2 = module7.bar_up_down(data = plot_data,
                 x_column = 'lacdinac',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_3',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_lacdinac ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of LacdiNAc',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of LacdiNAc in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_fucosylated_type
fig3 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'fucosylated type',
                 up_column = 'Up_data_fucosylated_type',
                 down_column = 'Down_data_fucosylated_type',
                 subfolder='StrucGAP_GlycoPeptideQuant_3',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_fucosylated_type",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of fucosylation',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs count',
                 figure_description = 'Comparison of fucosylated type in up- versus downregulated IGPs.',
                 )

fig4 = module7.bar_up_down(data = plot_data,
                 x_column = 'fucosylated type',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_3',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_fucosylated_type ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of fucosylation',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of fucosylated type in up- versus downregulated IGPs based on total number of related glycan.',
                 )

plot_data = module4.differential_analysis_acgc
fig5 = module7.bar_up_down(data = plot_data,
                           if_stack = False,
                 x_column = 'Ac/Gc',
                 up_column = 'Up_data_Ac/Gc',
                 down_column = 'Down_data_Ac/Gc',
                 subfolder='StrucGAP_GlycoPeptideQuant_3',
                 colors = ['#F9C3D7', '#3558AE'],
                 filename="differential_analysis_acgc",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of sialylation',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs count',
                 figure_description = 'Comparison of sialylated type in up- versus downregulated IGPs.',
                 )

fig6 = module7.bar_up_down(data = plot_data,
                 x_column = 'Ac/Gc',
                 up_column = 'Up_ratio',
                 down_column = 'Down_ratio',
                 subfolder='StrucGAP_GlycoPeptideQuant_3',
                 colors = ['#B64074', '#2A255C'],
                 filename="differential_analysis_acgc ratio",
                 xaxis_label_text_split = 0,
                 xaxis_title = 'IGPs containing different types of sialylatione',
                 xaxis_title_gap = 40,
                 xaxis_label_rotate = 0,
                 yaxis_title_gap = 50,
                 xaxis_label_font_size = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title = 'Up and downregulated(-) IGPs ratio',
                 figure_description = 'Ratio of sialylated type in up- versus downregulated IGPs based on total number of related glycan.',
                 )

fig7 = module7.bar_up_down_ratio(feature='core_structure', 
                          colors=['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                          subfolder='StrucGAP_GlycoPeptideQuant_3',
                          filename="core_structure",
                          xaxis_font_size = 25,
                          yaxis_font_size = 25,
                          legend_fontsize = 20,
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for core structures.',
                          )

fig8 = module7.bar_up_down_ratio(feature='glycan_type', 
                          colors=['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
                          subfolder='StrucGAP_GlycoPeptideQuant_3',
                          filename="glycan_type",
                          xaxis_font_size = 25,
                          yaxis_font_size = 25,
                          legend_fontsize = 20,
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for glycan types.',
                          )

fig9 = module7.bar_up_down_ratio(feature='branches_count', 
                          colors=random_colors,
                          subfolder='StrucGAP_GlycoPeptideQuant_3',
                          filename="branches_count",
                          xaxis_font_size = 25,
                          yaxis_font_size = 25,
                          legend_fontsize = 20,
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for branch counts.',
                          )

fig10 = module7.bar_up_down_ratio(feature='branches_structure', 
                          colors=random_colors,
                          subfolder='StrucGAP_GlycoPeptideQuant_3',
                          filename="branches_structure",
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for branch structures.',
                          )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.add_figure(fig10, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig11, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig12, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.compose_figures("StrucGAP_GlycoPeptideQuant_3.pdf", figure_name="StrucGAP_GlycoPeptideQuant_3",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10,11,12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoPeptideQuant_4
fig1 = module7.bar_up_down_ratio(feature='lacdinac', 
                          colors=random_colors,
                          subfolder='StrucGAP_GlycoPeptideQuant_4',
                          filename="lacdinac",
                          xaxis_font_size = 25,
                          yaxis_font_size = 25,
                          legend_fontsize = 20,
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for LacdiNAc.',
                          )

fig2 = module7.bar_up_down_ratio(feature='fucosylated_type', 
                          colors=random_colors,
                          subfolder='StrucGAP_GlycoPeptideQuant_4',
                          filename="fucosylated_type",
                          xaxis_font_size = 25,
                          yaxis_font_size = 25,
                          legend_fontsize = 20,
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for fucosylated type.',
                          )

fig3 = module7.bar_up_down_ratio(feature='acgc', 
                          colors=random_colors,
                          subfolder='StrucGAP_GlycoPeptideQuant_4',
                          filename="acgc",
                          xaxis_font_size = 25,
                          yaxis_font_size = 25,
                          legend_fontsize = 20,
                          figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends for sialylated type.',
                          )

plot_data = module4.result_glycoprotein_glycosite_glycan_variation.copy()
plot_data['sum'] = plot_data['Up count'] + plot_data['Down count']
plot_data = plot_data.sort_values('sum', ascending=False)
plot_data = plot_data[plot_data['Up count']!=0]
plot_data = plot_data[plot_data['Down count']!=0]
protein = list(plot_data['Glycoprotein'])[0]
site = list(plot_data['Glycosite_Position'])[0]
plot_data = module6.pg_fc[(module6.pg_fc['ProteinID']==protein)&(module6.pg_fc['Glycosite_Position']==site)].copy()
plot_data = plot_data[plot_data['p_g']<0.05]
plot_data['fc'] = plot_data['fc_g']
fig4 = module7.up_down_scatter(
    [plot_data],
    [protein + "+" + site],
    fc_threshold=1,
    show_xaxis=False,
    spine_width=1.2,
    ytick_labelsize=30,
    yaxis_title = 'Log2(FC) of IGPs',
    yaxis_title_font_size = 30,
    scatter_size=200,
    scatter_edgecolor='k',
    up_color='red',
    down_color='blue',
    bbox_facecolor='yellow',
    bbox_textsize=20,
    legend_fontsize = 20,
    subfolder='StrucGAP_GlycoPeptideQuant_4',
    filename=protein + "+" + site,
    figure_description = f'Highly variation of glycosite Asn-{site} on glycoprotein {protein} showed in glycosite-level analysis.',
)

module7.add_figure(fig1, figure_name="StrucGAP_GlycoPeptideQuant_4")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoPeptideQuant_4")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoPeptideQuant_4")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoPeptideQuant_4")
# module7.add_figure(fig5, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig6, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig7, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig8, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig9, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig10, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig11, figure_name="StrucGAP_GlycoPeptideQuant_3")
# module7.add_figure(fig12, figure_name="StrucGAP_GlycoPeptideQuant_3")
module7.compose_figures("StrucGAP_GlycoPeptideQuant_4.pdf", figure_name="StrucGAP_GlycoPeptideQuant_4",
                        custom_sizes=[[1], [2], [3], [4]])  # 生成后自动清理figure1数据

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
)

plot_data = module5.bp_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig2 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='bp_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:BP-enriched terms.',
             )

plot_data = module5.mf_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig3 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='mf_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:MF-enriched terms.',
             )

plot_data = module5.cc_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig4 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='cc_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:CC-enriched terms.',
             )

plot_data = module5.bp_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig10 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.compose_figures("StrucGAP_FunctionAnnotation_both_glycoproteins_1.pdf", figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1",
                        custom_sizes=[[1,2,3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_both_glycoproteins_2
plot_data = module5.bp_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='bp_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='mf_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='cc_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig4 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.compose_figures("StrucGAP_FunctionAnnotation_both_glycoproteins_2.pdf", figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2",
                        custom_sizes=[[1,4], [2,5], [3,6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_both_glycoproteins_3
plot_data = module5.bp_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig4 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig8 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig9 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_both_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:CC-enriched terms.',
              )

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.compose_figures("StrucGAP_FunctionAnnotation_both_glycoproteins_3.pdf", figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_up_glycoproteins_1
module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_up_result')  

plot_data = module5.ora_no_background_up_result.copy()
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
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    filename="both proteins dotplot",
    figure_description = 'Enrichment results of upregulated glycoproteins based on GO enrichment.',
)

plot_data = module5.bp_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig2 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='bp_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:BP-enriched terms.',
             )

plot_data = module5.mf_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig3 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='mf_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:MF-enriched terms.',
             )

plot_data = module5.cc_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig4 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='cc_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:CC-enriched terms.',
             )

plot_data = module5.bp_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig10 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.compose_figures("StrucGAP_FunctionAnnotation_up_glycoproteins_1.pdf", figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_1",
                        custom_sizes=[[1,2,3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_up_glycoproteins_2
plot_data = module5.bp_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='bp_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='mf_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='cc_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig4 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.compose_figures("StrucGAP_FunctionAnnotation_up_glycoproteins_2.pdf", figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_2",
                        custom_sizes=[[1,4], [2,5], [3,6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_up_glycoproteins_3
plot_data = module5.bp_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig4 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig8 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig9 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_up_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:CC-enriched terms.',
              )

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.compose_figures("StrucGAP_FunctionAnnotation_up_glycoproteins_3.pdf", figure_name="StrucGAP_FunctionAnnotation_up_glycoproteins_3",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_down_glycoproteins_1
module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_down_result')  

plot_data = module5.ora_no_background_down_result.copy()
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
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    filename="both proteins dotplot",
    figure_description = 'Enrichment results of downregulated glycoproteins based on GO enrichment',
)

plot_data = module5.bp_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig2 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='bp_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:BP-enriched terms.',
             )

plot_data = module5.mf_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig3 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='mf_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:MF-enriched terms.',
             )

plot_data = module5.cc_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig4 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='cc_core_structure',
             figure_description = 'Distribution of different core structures across the top 10 GO:CC-enriched terms.',
             )

plot_data = module5.bp_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_glycan_type',
    figure_description = 'Distribution of different glycan types across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig10 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_1',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_branches_count',
    figure_description = 'Distribution of different branch counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.compose_figures("StrucGAP_FunctionAnnotation_down_glycoproteins_1.pdf", figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_1",
                        custom_sizes=[[1,2,3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_down_glycoproteins_2
plot_data = module5.bp_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='bp_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='mf_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='cc_branches_structure',
    figure_description = 'Distribution of different branch structures across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig4 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_sialicacid_count',
    figure_description = 'Distribution of different sialic acid counts across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_2',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_fucose_count',
    figure_description = 'Distribution of different fucose counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.compose_figures("StrucGAP_FunctionAnnotation_down_glycoproteins_2.pdf", figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_2",
                        custom_sizes=[[1,4], [2,5], [3,6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_FunctionAnnotation_down_glycoproteins_3
plot_data = module5.bp_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_lacdinac',
              figure_description = 'Distribution of different LacdiNAc across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig4 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_fucosylated_type',
              figure_description = 'Distribution of different fucosylated types across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig8 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig9 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_FunctionAnnotation_down_glycoproteins_3',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_acgc',
              figure_description = 'Distribution of different sialylated types across the top 10 GO:CC-enriched terms.',
              )

module7.add_figure(fig1, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig2, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig3, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig4, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig5, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig6, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig7, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig8, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
module7.add_figure(fig9, figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.compose_figures("StrucGAP_FunctionAnnotation_down_glycoproteins_3.pdf", figure_name="StrucGAP_FunctionAnnotation_down_glycoproteins_3",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_1
fig1 = module7.volcano_plot(data = module6.proteomic_fc,
                    fc_column = 'fc',
                    p_column = 'pvalue_ttest',
                    fc = 1.5,
                    p_value = 0.05,
                    subfolder='StrucGAP_GlycoNetwork_1',
                    filename='proteomic volcano',
                    figure_description = 'Volcano plot of differentially expressed proteins.',
                    )

plot_data = module6.proteomic_fc.copy()
cols = list(plot_data.columns)
cols[:10] = ['126.1277', '127.1248', '127.1311', '128.1281', '128.1344',
                  '129.1315', '129.1378', '130.1348', '130.1411', '131.1382']
plot_data.columns = cols
fig2 = module7.dimension_reduction(
    data = plot_data,
    data_columns = ['126.1277', '127.1248', '127.1311', '128.1281', '128.1344',
                  '129.1315', '129.1378', '130.1348', '130.1411', '131.1382'],
    sample_group = module4.sample_group,
    filter_data = plot_data, # filter_data=module4.fc_result,
    p_column = 'pvalue_ttest',
    p_value = 0.05,
    fc = 1.5,
    method='pca',
    dimension_number = 2,
    random_state = 0,
    colors = ['#3558AE', '#B64074'],
    subfolder='StrucGAP_GlycoNetwork_1',
    filename='proteomic pca',
    figure_description = 'PCA plot based on differential proteins expression profiles.',
)

fig3 = module7.heatmap2(data = plot_data,
                 columns=['126.1277','127.1248','127.1311','128.1281','128.1344', 
                          '129.1315','129.1378','130.1348','130.1411','131.1382'],
                 filter_data = plot_data,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 z_score = 0,
                 splitline_width = 0.0000001,
                 subfolder='StrucGAP_GlycoNetwork_1',
                 filename='proteomic heatmap',
                 figure_description = 'Heatmap of differentially expressed proteins.',
                 )

plot_data = module6.proteomic_protein_glycosite_count.copy()
plot_data['sum'] = plot_data['Up_glycopeptide_count'] + plot_data['Down_glycopeptide_count']
plot_data = plot_data[plot_data['Protein_type']=='None significance']
plot_data = plot_data.sort_values('sum', ascending=False)
plot_data = plot_data[plot_data['Up_glycopeptide_count']!=0]
plot_data = plot_data[plot_data['Down_glycopeptide_count']!=0]
fc_data = module6.proteomic_protein_glycosite_value.copy()
fc_data = fc_data[
    (fc_data['Up_glycopeptide_fc'] > 1.5) &
    (fc_data['Down_glycopeptide_fc'] < 0.667)
]
merged = pd.merge(
    plot_data,
    fc_data,
    on=['ProteinID', 'Glycosite_Position'],
    suffixes=('_count', '_fc')
)
plot_data = merged.sort_values('sum', ascending=False)
protein = list(plot_data['ProteinID'])[0]
site = list(plot_data['Glycosite_Position'])[0]
plot_data = module6.pg_fc[(module6.pg_fc['ProteinID']==protein)&(module6.pg_fc['Glycosite_Position']==site)].copy()
plot_data = plot_data[plot_data['p_g']<0.05]
plot_data['fc'] = plot_data['fc_g']
fig4 = module7.up_down_scatter(
    [plot_data],
    [protein + "+" + site],
    fc_threshold=1,
    show_xaxis=False,
    spine_width=1.2,
    ytick_labelsize=30,
    yaxis_title = 'Log2(FC) of IGPs',
    yaxis_title_font_size = 30,
    scatter_size=200,
    scatter_edgecolor='k',
    up_color='red',
    down_color='blue',
    bbox_facecolor='yellow',
    bbox_textsize=11,
    legend_fontsize = 20,
    subfolder='StrucGAP_GlycoNetwork_1',
    filename=protein + "+" + site,
    figure_description = f'Glycosites Asn-{site} on protein {protein} exhibited mixed glycan regulation despite stable protein levels.',
)

plot_data = module6.proteomic_protein_glycosite_same_direction.copy()
plot_data['Protein+glycosite'] = plot_data['GeneName'] + '+' + plot_data['Glycosite_Position']
plot_data = plot_data[['Protein+glycosite', 'fc_p', 'fc_g']]
plot_data = plot_data.sort_values('Protein+glycosite', ascending=True)
plot_data = plot_data.set_index('Protein+glycosite')
fig5 = module7.heatmap2(data = plot_data,
                 columns=['fc_p', 'fc_g'],
                 filter_data = None,
                 cluster = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = True,
                 z_score = None,
                 minvalue=None, 
                 centervalue = 0,
                 maxvalue=None,
                 splitline_width = 0.1,
                 yaxis_label_show=True,
                 yaxis_label_font_size=5,
                 subfolder='StrucGAP_GlycoNetwork_1',
                 filename='proteomic_protein_glycosite_same_direction',
                 figure_description = 'Glycopeptide showed same regulation direction based on the same glycosite.',
                 )

plot_data = module6.proteomic_protein_glycosite_different_direction.copy()
plot_data['Protein+glycosite'] = plot_data['GeneName'] + '+' + plot_data['Glycosite_Position']
plot_data = plot_data[['Protein+glycosite', 'fc_p', 'fc_g']]
plot_data = plot_data.sort_values('Protein+glycosite', ascending=True)
plot_data = plot_data.set_index('Protein+glycosite')
fig6 = module7.heatmap2(data = plot_data,
                 columns=['fc_p', 'fc_g'],
                 filter_data = None,
                 cluster = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = True,
                 z_score = None,
                 minvalue=None, 
                 centervalue = 0,
                 maxvalue=None,
                 splitline_width = 0.00001,
                 yaxis_label_show=False,
                 yaxis_label_font_size=5,
                 subfolder='StrucGAP_GlycoNetwork_1',
                 filename='proteomic_protein_glycosite_different_direction',
                 figure_description = 'Glycopeptide showed different regulation direction based on the same glycosite.',
                 )

plot_data1 = module6.protein_no_glyco_up[['fc_g','fc_p','normalized_fc_g']]
plot_data1['group'] = 'up'
plot_data2 = module6.protein_no_glyco_down[['fc_g','fc_p','normalized_fc_g']]
plot_data2['group'] = 'down'
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data.reset_index(inplace=True,drop=True)
fig7 = module7.scatter(data = plot_data,
             group_column = 'group',
             x_column = 'fc_p',
             y_column = 'fc_g',
             subfolder='StrucGAP_GlycoNetwork_1',
             top_xaxis_line_show = False,
             right_yaxis_line_show = False,
             xaxis_splitline_show = False,
             yaxis_splitline_show = False,
             xaxis_label_rotate = 0,
             xaxis_label_font_size = 30,
             yaxis_label_font_size = 30,
             legend_font_size = 30,
             xaxis_title = 'FC of proteins',
             yaxis_title = 'FC of IGPs',
             xaxis_title_gap = 40,
             yaxis_title_gap = 30,
             filename='protein_no_glyco_up_down',
             figure_description = 'Identification of N-glycan features that were altered solely at the glycopeptide level, independent of protein-level changes.',
             )

plot_data1 = module6.protein_up_glyco_up[['fc_g','fc_p','normalized_fc_g']]
plot_data1['group'] = 'up'
plot_data2 = module6.protein_up_glyco_down[['fc_g','fc_p','normalized_fc_g']]
plot_data2['group'] = 'down'
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data.reset_index(inplace=True,drop=True)
fig8 = module7.scatter(data = plot_data,
             group_column = 'group',
             x_column = 'fc_p',
             y_column = 'fc_g',
             subfolder='StrucGAP_GlycoNetwork_1',
             top_xaxis_line_show = False,
             right_yaxis_line_show = False,
             xaxis_splitline_show = False,
             yaxis_splitline_show = False,
             xaxis_label_rotate = 0,
             xaxis_label_font_size = 30,
             yaxis_label_font_size = 30,
             legend_font_size = 30,
             xaxis_title = 'FC of proteins',
             yaxis_title = 'FC of IGPs',
             xaxis_title_gap = 40,
             yaxis_title_gap = 30,
             filename='protein_up_glyco_up_down',
             figure_description = 'Identification of N-glycan features that were altered both at the glycopeptide level and protein-level (upregulated).',
             )

plot_data1 = module6.protein_down_glyco_up[['fc_g','fc_p','normalized_fc_g']]
plot_data1['group'] = 'up'
plot_data2 = module6.protein_down_glyco_down[['fc_g','fc_p','normalized_fc_g']]
plot_data2['group'] = 'down'
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data.reset_index(inplace=True,drop=True)
fig9 = module7.scatter(data = plot_data,
             group_column = 'group',
             x_column = 'fc_p',
             y_column = 'fc_g',
             subfolder='StrucGAP_GlycoNetwork_1',
             top_xaxis_line_show = False,
             right_yaxis_line_show = False,
             xaxis_splitline_show = False,
             yaxis_splitline_show = False,
             xaxis_label_rotate = 0,
             xaxis_label_font_size = 30,
             yaxis_label_font_size = 30,
             legend_font_size = 30,
             xaxis_title = 'FC of proteins',
             yaxis_title = 'FC of IGPs',
             xaxis_title_gap = 40,
             yaxis_title_gap = 50,
             filename='protein_no_glyco_up_down',
             figure_description = 'Identification of N-glycan features that were altered both at the glycopeptide level and protein-level (downregulated).',
             )

plot_data1 = module6.protein_no_glyco_up.copy()
plot_data2 = module6.protein_no_glyco_down.copy()
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data['fc'] = plot_data['fc_g'] / plot_data['fc_p']
plot_data = plot_data[['structure_coding','fc_g', 'fc_p', 'fc']]
plot_data = plot_data.set_index('structure_coding',drop=True)
plot_data = plot_data.sort_values('fc_g', ascending=False)
fig10 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = True,
                 cluster = None,
                 yaxis_label_show = False,
                 z_score = None,
                 centervalue = 0,
                 minvalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 subfolder='StrucGAP_GlycoNetwork_1',
                 filename='protein_no_glyco_up_down',
                 figure_description = 'Identification of N-glycan features that were altered solely at the glycopeptide level, independent of protein-level changes.',
                 )

plot_data1 = module6.protein_up_glyco_up.copy()
plot_data2 = module6.protein_up_glyco_down.copy()
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data['fc'] = plot_data['fc_g'] / plot_data['fc_p']
plot_data = plot_data[['structure_coding','fc_g', 'fc_p', 'fc']]
plot_data = plot_data.set_index('structure_coding',drop=True)
plot_data = plot_data.sort_values('fc_g', ascending=False)
fig11 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = True,
                 cluster = None,
                 yaxis_label_show = False,
                 z_score = None,
                 centervalue = 0,
                 minvalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 subfolder='StrucGAP_GlycoNetwork_1',
                 filename='protein_up_glyco_up_down',
                 figure_description = 'Identification of N-glycan features that were altered both at the glycopeptide level and protein-level (upregulated).',
                 )

plot_data1 = module6.protein_down_glyco_up.copy()
plot_data2 = module6.protein_down_glyco_down.copy()
plot_data = pd.concat([plot_data1,plot_data2],axis=0)
plot_data['fc'] = plot_data['fc_g'] / plot_data['fc_p']
plot_data = plot_data[['structure_coding','fc_g', 'fc_p', 'fc']]
plot_data = plot_data.set_index('structure_coding',drop=True)
plot_data = plot_data.sort_values('fc_g', ascending=False)
fig12 = module7.heatmap2(data = plot_data,
                 colors='seismic',
                 columns=plot_data.columns,
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = True,
                 cluster = None,
                 yaxis_label_show = False,
                 z_score = None,
                 centervalue = 0,
                 minvalue = None,
                 maxvalue = None,
                 splitline_width = 0.5,
                 subfolder='StrucGAP_GlycoNetwork_1',
                 filename='protein_down_glyco_up_down',
                 figure_description = 'Identification of N-glycan features that were altered both at the glycopeptide level and protein-level (downregulated).',
                 )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig10, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig11, figure_name="StrucGAP_GlycoNetwork_1")
module7.add_figure(fig12, figure_name="StrucGAP_GlycoNetwork_1")
module7.compose_figures("StrucGAP_GlycoNetwork_1.pdf", figure_name="StrucGAP_GlycoNetwork_1",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_2
plot_data = module6.glycosyltransferases.copy()
plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
plot_data.reset_index(inplace=True)
plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
plot_data.set_index('gene_id',inplace=True)
fig1 = module7.heatmap2(data = plot_data,
                 colors='Spectral',
                 columns=plot_data.columns[1:11],
                 filter_data = None,
                 filter_columns = ['fc', 'pvalue_ttest'],
                 filter_values = [1.5, 0.05],
                 log = False,
                 yaxis_label_show = True,
                 xaxis_label_show = False,
                 z_score = None,
                 splitline_width = 0.5,
                 xaxis_title_font_size=25, 
                 yaxis_title_font_size=25,
                 xaxis_label_font_size=15, 
                 yaxis_label_font_size=10, 
                 filename = 'glycosyltransferases',
                 subfolder='StrucGAP_GlycoNetwork_2',
                 figure_description = 'Quantitative profiling of altered glycosyltransferases (P value < 0.05).',
                 )

plot_data = module6.glycosidases.copy()
plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
plot_data = plot_data.sort_values(by='fc', ascending=False)
top2 = list(plot_data.index)[:2]
plot_data = module6.cv_filter_data.reset_index().copy()
plot_data = plot_data[plot_data['Accession'].isin(top2)]
plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
plot_data.set_index('gene_id',inplace=True)
fig2 = module7.violin_plot(data=plot_data,
                    item_column = 'Accession',
                    item_name = [top2[0]],
                    group1_columns = plot_data.columns[1:6],
                    group2_columns = plot_data.columns[6:11],
                    p_data = module6.proteomic_fc,
                    p_column='pvalue_ttest',
                    subfolder='StrucGAP_GlycoNetwork_2',
                    filename = plot_data.index[0],
                    xaxis_label_font_size = 30,
                    yaxis_label_font_size = 30,
                    legend_fontsize = 30,
                    p_text_offset=30,  
                    yaxis_title = 'Quantification value',
                    yaxis_title_font_size = 30,
                    figure_description = f'Significantly upregulated glycosidases {plot_data.index[0]}.',
                  ) 

fig3 = module7.violin_plot(data=plot_data,
                    item_column = 'Accession',
                    item_name = [top2[1]],
                    group1_columns = plot_data.columns[1:6],
                    group2_columns = plot_data.columns[6:11],
                    p_data = module6.proteomic_fc,
                    p_column='pvalue_ttest',
                    subfolder='StrucGAP_GlycoNetwork_2',
                    filename = plot_data.index[1],
                    xaxis_label_font_size = 30,
                    yaxis_label_font_size = 30,
                    legend_fontsize = 30,
                    p_text_offset=2, 
                    yaxis_title = 'Quantification value',
                    yaxis_title_font_size = 30,
                    figure_description = f'Significantly upregulated glycosidases {plot_data.index[1]}.',
                  )

plot_data1 = module6.cv_filter_data[module6.cv_filter_data.index.isin(module6.sialyltransferases.index)]
plot_data2 = module6.sialyltransferases[['fc']]
plot_data3 = module6.sialyltransferases[['pvalue_ttest']]
plot_data1.columns = ['Control 1','Control 2','Control 3','Control 4','Control 5',
                      'Sample 1','Sample 2','Sample 3','Sample 4','Sample 5']
fig4 = module7.complexheatmap(data = plot_data1, 
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
                       subfolder='StrucGAP_GlycoNetwork_2',
                       col_split=None,
                       cmap = 'Blues', 
                       z_score = 0,
                       show_rownames = True,
                       row_split = None,
                       filename = 'sialyltransferases',
                       linewidths = 3,
                       figure_description = 'Heatmap of identified sialyltransferases.',
                       )

plot_data1 = module6.cv_filter_data[module6.cv_filter_data.index.isin(module6.fucosyltransferase.index)]
plot_data2 = module6.fucosyltransferase[['fc']]
plot_data3 = module6.fucosyltransferase[['pvalue_ttest']]
plot_data1.columns = ['Control 1','Control 2','Control 3','Control 4','Control 5',
                      'Sample 1','Sample 2','Sample 3','Sample 4','Sample 5']
fig5 = module7.complexheatmap(data = plot_data1, 
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
                       subfolder='StrucGAP_GlycoNetwork_2',
                       col_split=None,
                       cmap = 'Blues', 
                       z_score = 0,
                       show_rownames = True,
                       row_split = None,
                       filename = 'fucosyltransferase',
                       linewidths = 3,
                       figure_description = 'Heatmap of identified fucosyltransferases.',
                       )

plot_data = module6.glycan_binding_protein.copy()
plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
plot_data = plot_data[(plot_data['fc']>1.5)|(plot_data['fc']<1/1.5)]
plot_data.reset_index(inplace=True)
plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
plot_data.set_index('gene_id',inplace=True,drop=False)
fig6 = module7.violin_plot(data=plot_data,
                    item_column='gene_id',
                    item_name=list(plot_data['gene_id']),
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
                    subfolder='StrucGAP_GlycoNetwork_2',
                    filename = 'glycan binding protins',
                    xaxis_label_font_size = 30,
                    yaxis_label_font_size = 30,
                    legend_fontsize = 30,
                    yaxis_title = 'Quantification value',
                    yaxis_title_font_size = 30,
                    figure_description = 'Expression patterns of significantly altered glycan-binding proteins (P value < 0.05, FC > 1.5 or < 0.67).',
                    ) 

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_2")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_2")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_2")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_2")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_2")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_2")
# module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_2")
# module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_2")
# module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_2")
# module7.add_figure(fig10, figure_name="StrucGAP_GlycoNetwork_2")
# module7.add_figure(fig11, figure_name="StrucGAP_GlycoNetwork_2")
# module7.add_figure(fig12, figure_name="StrucGAP_GlycoNetwork_2")
module7.compose_figures("StrucGAP_GlycoNetwork_2.pdf", figure_name="StrucGAP_GlycoNetwork_2",
                        custom_sizes=[[1], [2], [3], [4], [5], [6]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_3
module5 = StrucGAP_FunctionAnnotation(module6, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
module5.go_function_structure(function_data = 'ora_no_background_up_result') 
 
plot_data = module5.ora_no_background_up_result.copy()
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
    subfolder='StrucGAP_GlycoNetwork_3',
    filename="both proteins dotplot (protein_no_glyco_up)",
    figure_description = 'Enrichment results of upregulated glycoproteins based on GO enrichment (upregulated glycopeptides with stable protein levels; this applies to all panels below).',
)

plot_data = module5.bp_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig2 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_GlycoNetwork_3',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='bp_core_structure',
             figure_description = 'Distribution of core structures across the top 10 GO:BP-enriched terms.',
             )

plot_data = module5.mf_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig3 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_GlycoNetwork_3',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='mf_core_structure',
             figure_description = 'Distribution of core structures across the top 10 GO:MF-enriched terms.',
             )

plot_data = module5.cc_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig4 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_GlycoNetwork_3',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='cc_core_structure',
             figure_description = 'Distribution of core structures across the top 10 GO:CC-enriched terms.',
             )

plot_data = module5.bp_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_3',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_glycan_type',
    figure_description = 'Distribution of glycan types across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_3',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_glycan_type',
    figure_description = 'Distribution of glycan types across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_3',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_glycan_type',
    figure_description = 'Distribution of glycan types across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_3',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_branches_count',
    figure_description = 'Distribution of branch counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_3',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_branches_count',
    figure_description = 'Distribution of branch counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig10 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_3',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_branches_count',
    figure_description = 'Distribution of branch counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_3")
module7.add_figure(fig10, figure_name="StrucGAP_GlycoNetwork_3")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.compose_figures("StrucGAP_GlycoNetwork_3.pdf", figure_name="StrucGAP_GlycoNetwork_3",
                        custom_sizes=[[1,2,3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_4
plot_data = module5.bp_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='bp_branches_structure',
    figure_description = 'Distribution of branch structures across the top 10 GO:BP-enriched terms (upregulated glycopeptides with stable protein levels; this applies to all panels below).',
)

plot_data = module5.mf_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='mf_branches_structure',
    figure_description = 'Distribution of branch structures across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='cc_branches_structure',
    figure_description = 'Distribution of branch structures across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig4 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_sialicacid_count',
    figure_description = 'Distribution of sialic acid counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_sialicacid_count',
    figure_description = 'Distribution of sialic acid counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_sialicacid_count',
    figure_description = 'Distribution of sialic acid counts across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_fucose_count',
    figure_description = 'Distribution of fucose counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_fucose_count',
    figure_description = 'Distribution of fucose counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_4',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_fucose_count',
    figure_description = 'Distribution of fucose counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_4")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_4")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.compose_figures("StrucGAP_GlycoNetwork_4.pdf", figure_name="StrucGAP_GlycoNetwork_4",
                        custom_sizes=[[1,4], [2,5], [3,6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_5
plot_data = module5.bp_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_lacdinac',
              figure_description = 'Distribution of LacdiNAc across the top 10 GO:BP-enriched terms (upregulated glycopeptides with stable protein levels; this applies to all panels below).',
              )

plot_data = module5.mf_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_lacdinac',
              figure_description = 'Distribution of LacdiNAc across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_lacdinac.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_lacdinac',
              figure_description = 'Distribution of LacdiNAc across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig4 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_fucosylated_type',
              figure_description = 'Distribution of fucosylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_fucosylated_type',
              figure_description = 'Distribution of fucosylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_fucosylated_type',
              figure_description = 'Distribution of fucosylated types across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_acgc',
              figure_description = 'Distribution of sialylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig8 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_acgc',
              figure_description = 'Distribution of sialylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig9 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_5',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_acgc',
              figure_description = 'Distribution of sialylated types across the top 10 GO:CC-enriched terms.',
              )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_5")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_5")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.compose_figures("StrucGAP_GlycoNetwork_5.pdf", figure_name="StrucGAP_GlycoNetwork_5",
                        custom_sizes=[[1], [2], [3], [4], [5], [6], [7], [8], [9]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_6
module5 = StrucGAP_FunctionAnnotation(module6, 
                                 data_manager=data_manager, data_type = 'protein_no_glyco_down')  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
module5.go_function_structure(function_data = 'ora_no_background_down_result') 
 
plot_data = module5.ora_no_background_down_result.copy()
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
    subfolder='StrucGAP_GlycoNetwork_6',
    filename="both proteins dotplot (protein_no_glyco_down)",
    figure_description = 'Enrichment results of downregulated glycoproteins based on GO enrichment (downregulated glycopeptides with stable protein levels; this applies to all panels below).',
)

plot_data = module5.bp_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig2 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_GlycoNetwork_6',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='bp_core_structure',
             figure_description = 'Distribution of core structures across the top 10 GO:BP-enriched terms.',
             )

plot_data = module5.mf_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig3 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_GlycoNetwork_6',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='mf_core_structure',
             figure_description = 'Distribution of core structures across the top 10 GO:MF-enriched terms.',
             )

plot_data = module5.cc_core_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
plot_data = plot_data.replace('A2B2C1D1dD2dD1','Core-III')
plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
fig4 = module7.line(data = plot_data,
             colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
             y_column = 'Core_structure',
             x_columns = plot_data.columns[1:11],
             subfolder='StrucGAP_GlycoNetwork_6',
             symbol_size = 5,
             plot_title = None,
             xaxis_label_rotate = -15,
             xaxis_title_gap = 35,
             xaxis_label_font_size = 10,
             xaxis_label_text_split = 20,
             yaxis_label_font_size = 20,
             legend_font_size = 20,
             yaxis_title_gap = 40,
             yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
             filename='cc_core_structure',
             figure_description = 'Distribution of core structures across the top 10 GO:CC-enriched terms.',
             )

plot_data = module5.bp_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_6',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_glycan_type',
    figure_description = 'Distribution of glycan types across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_6',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_glycan_type',
    figure_description = 'Distribution of glycan types across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_glycan_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Glycan_type",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_6',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_glycan_type',
    figure_description = 'Distribution of glycan types across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_6',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_branches_count',
    figure_description = 'Distribution of branch counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_6',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_branches_count',
    figure_description = 'Distribution of branch counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['BranchNumber'] = plot_data['BranchNumber'].astype(str)
fig10 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "BranchNumber",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_6',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_branches_count',
    figure_description = 'Distribution of branch counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_6")
module7.add_figure(fig10, figure_name="StrucGAP_GlycoNetwork_6")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_1")
module7.compose_figures("StrucGAP_GlycoNetwork_6.pdf", figure_name="StrucGAP_GlycoNetwork_6",
                        custom_sizes=[[1,2,3], [4], [5], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_7
plot_data = module5.bp_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='bp_branches_structure',
    figure_description = 'Distribution of branch structures across the top 10 GO:BP-enriched terms (downregulated glycopeptides with stable protein levels; this applies to all panels below).',
)

plot_data = module5.mf_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='mf_branches_structure',
    figure_description = 'Distribution of branch structures across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_branches_structure.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.multi_bar(
    data = plot_data,  
    x_column = "Branches",
    y_column = plot_data.columns[:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    bar_width="50%",
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 10,
    yaxis_label_font_size = 10,
    yaxis_label_margin=35,
    filename='cc_branches_structure',
    figure_description = 'Distribution of branch structures across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig4 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_sialicacid_count',
    figure_description = 'Distribution of sialic acid counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig5 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_sialicacid_count',
    figure_description = 'Distribution of sialic acid counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_sialicacid_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Sialicacid_count'] = plot_data['Sialicacid_count'].astype(str)
fig6 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Sialicacid_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_sialicacid_count',
    figure_description = 'Distribution of sialic acid counts across the top 10 GO:CC-enriched terms.',
)

plot_data = module5.bp_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig7 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='bp_fucose_count',
    figure_description = 'Distribution of fucose counts across the top 10 GO:BP-enriched terms.',
)

plot_data = module5.mf_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig8 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='mf_fucose_count',
    figure_description = 'Distribution of fucose counts across the top 10 GO:MF-enriched terms.',
)

plot_data = module5.cc_fucose_count.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
plot_data['Fucose_count'] = plot_data['Fucose_count'].astype(str)
fig9 = module7.bar_multi_columns(
    data = plot_data,  
    colors = ['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                        '#003C71','#8FC1D9','#521887','#F4A2B9',
                        '#FF7D5B','#FBFTEF'],
    y_column = "Fucose_count",
    x_columns = plot_data.columns[1:11], 
    subfolder='StrucGAP_GlycoNetwork_7',
    xaxis_splitline_show = False,
    yaxis_splitline_show = False,
    xaxis_label_font_size = 10,
    xaxis_label_text_split = 20,
    yaxis_label_font_size = 20,
    legend_font_size = 20,
    y_max=1,
    xaxis_title_gap = 35,
    yaxis_title = 'Percentage of IGPs carrying\neach core structures (per term)',
    filename='cc_fucose_count',
    figure_description = 'Distribution of fucose counts across the top 10 GO:CC-enriched terms.',
)

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_7")
module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_7")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_2")
module7.compose_figures("StrucGAP_GlycoNetwork_7.pdf", figure_name="StrucGAP_GlycoNetwork_7",
                        custom_sizes=[[1,4], [2,5], [3,6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据

# StrucGAP_GlycoNetwork_8
# plot_data = module5.bp_lacdinac.copy()
# plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
# if not plot_data.empty:
#     fig1 = module7.radar('plot_data',  
#                   columns = list(plot_data.columns[1:11]),
#                   text_font_size = 15,
#                   legend_font_size = 15,
#                   text_split = 15,
#                   subfolder='StrucGAP_GlycoNetwork_8',
#                   screen_column = plot_data.columns[0],
#                   screen_values = list(plot_data.iloc[:,0]),
#                   filename='bp_lacdinac',
#                   figure_description = 'Distribution of lacdinac across the top 10 GO:BP-enriched terms (protein_no_glyco_down)',
#                   )

# plot_data = module5.mf_lacdinac.copy()
# plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
# if not plot_data.empty:
#     fig2 = module7.radar('plot_data',  
#                   columns = list(plot_data.columns[1:11]),
#                   text_font_size = 15,
#                   legend_font_size = 15,
#                   text_split = 15,
#                   subfolder='StrucGAP_GlycoNetwork_8',
#                   screen_column = plot_data.columns[0],
#                   screen_values = list(plot_data.iloc[:,0]),
#                   filename='mf_lacdinac',
#                   figure_description = 'Distribution of lacdinac across the top 10 GO:MF-enriched terms (protein_no_glyco_down)',
#                   )

# plot_data = module5.cc_lacdinac.copy()
# plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
# if not plot_data.empty:
#     fig3 = module7.radar('plot_data',  
#                   columns = list(plot_data.columns[1:11]),
#                   text_font_size = 15,
#                   legend_font_size = 15,
#                   text_split = 15,
#                   subfolder='StrucGAP_GlycoNetwork_8',
#                   screen_column = plot_data.columns[0],
#                   screen_values = list(plot_data.iloc[:,0]),
#                   filename='cc_lacdinac',
#                   figure_description = 'Distribution of lacdinac across the top 10 GO:CC-enriched terms (protein_no_glyco_down)',
#                   )

plot_data = module5.bp_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig1 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_8',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_fucosylated_type',
              figure_description = 'Distribution of fucosylated types across the top 10 GO:BP-enriched terms (downregulated glycopeptides with stable protein levels; this applies to all panels below).',
              )

plot_data = module5.mf_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig2 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_8',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_fucosylated_type',
              figure_description = 'Distribution of fucosylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_fucosylated_type.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig3 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_8',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_fucosylated_type',
              figure_description = 'Distribution of fucosylated types across the top 10 GO:CC-enriched terms.',
              )

plot_data = module5.bp_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig4 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_8',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='bp_acgc',
              figure_description = 'Distribution of sialylated types across the top 10 GO:BP-enriched terms.',
              )

plot_data = module5.mf_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig5 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_8',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='mf_acgc',
              figure_description = 'Distribution of sialylated types across the top 10 GO:MF-enriched terms.',
              )

plot_data = module5.cc_acgc.copy()
plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
fig6 = module7.radar('plot_data',  
              columns = list(plot_data.columns[1:11]),
              text_font_size = 15,
              legend_font_size = 15,
              text_split = 15,
              subfolder='StrucGAP_GlycoNetwork_8',
              screen_column = plot_data.columns[0],
              screen_values = list(plot_data.iloc[:,0]),
              filename='cc_acgc',
              figure_description = 'Distribution of sialylated types across the top 10 GO:CC-enriched terms.',
              )

module7.add_figure(fig1, figure_name="StrucGAP_GlycoNetwork_8")
module7.add_figure(fig2, figure_name="StrucGAP_GlycoNetwork_8")
module7.add_figure(fig3, figure_name="StrucGAP_GlycoNetwork_8")
module7.add_figure(fig4, figure_name="StrucGAP_GlycoNetwork_8")
module7.add_figure(fig5, figure_name="StrucGAP_GlycoNetwork_8")
module7.add_figure(fig6, figure_name="StrucGAP_GlycoNetwork_8")
# module7.add_figure(fig7, figure_name="StrucGAP_GlycoNetwork_8")
# module7.add_figure(fig8, figure_name="StrucGAP_GlycoNetwork_8")
# module7.add_figure(fig9, figure_name="StrucGAP_GlycoNetwork_8")
# module7.add_figure(fig10, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig11, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
# module7.add_figure(fig12, figure_name="StrucGAP_FunctionAnnotation_both_glycoproteins_3")
module7.compose_figures("StrucGAP_GlycoNetwork_8.pdf", figure_name="StrucGAP_GlycoNetwork_8",
                        custom_sizes=[[1], [2], [3], [4], [5], [6]])  # 生成后自动清理figure1数据


# key insights StrucGAP_GlycoPeptideQuant
pdf_index = 1
figs = []
current_figs = []
for i in ['core_structure','glycan_type','branches_structure','branches_count',
          'lacdinac','fucosylated_type','acgc']:
    generated_figs = []  # 本轮新生成的fig，存储顺序：fig1, fig2, fig3
    if i == 'core_structure':
        j = 'core structures'
    elif i == 'glycan_type':
        j = 'glycan types'
    elif i == 'branches_structure':
        j = 'branch structures'    
    elif i == 'branches_count':
        j = 'branch counts'
    elif i == 'lacdinac':
        j = 'LacdiNAc'
    elif i == 'fucosylated type':
        j = 'fucosylated types'
    elif i == 'acgc':
        j = 'sialylated types'
    # ...生成fig1/fig2/fig3的逻辑...
    # 下面是假定逻辑，根据你的实际生成fig的代码复制粘贴即可
    data = pd.read_excel("test/analysis_result/StrucGAP_GlycoPeptideQuant_key_information.xlsx", sheet_name=i)
    if not data.empty:
        data = data.replace('A2B2C1D1dD1','Core-I') \
                          .replace('A2B2C1D1dD1dcbB5','Core-II') \
                          .replace('A2B2C1D1dD2dD1','Core-III') \
                          .replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
        fig1 = module7.bar_up_down_ratio(feature = i, 
                              colors=['#4E2A5C','#C54E83','#4EBABC','#FFC599',
                                    '#003C71','#8FC1D9','#521887','#F4A2B9',
                                    '#FF7D5B','#FBFTEF'],
                              screen_feature = [x for x in list(data[0].unique()) if not (isinstance(x, float) and np.isnan(x))][1:],
                              subfolder = f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}",
                              filename = i,
                              figure_description = f'Key information from substructure: {j.replace("_", " ")}.',
                  )
        generated_figs.append(fig1)

    da_data = pd.read_excel("test/analysis_result/StrucGAP_GlycoPeptideQuant_key_information.xlsx", sheet_name=f'da_{i}')
    if not da_data.empty:
        attr_name = f'differential_analysis_{i}'
        plot_data = getattr(module4, attr_name)
        plot_data = plot_data.replace('A2B2C1D1dD1','Core-I') \
                          .replace('A2B2C1D1dD1dcbB5','Core-II') \
                          .replace('A2B2C1D1dD2dD1','Core-III') \
                          .replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
        da_data[da_data.columns[1]]
        
        fig2 = module7.bar_up_down(data = plot_data,
                                   if_stack = False,
                         x_column = plot_data.columns[0],
                         up_column = plot_data.columns[1],
                         down_column = plot_data.columns[2],
                         subfolder=f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}",
                         colors = ['#F9C3D7', '#3558AE'],
                         filename = f"differential_analysis_{i}",
                         xaxis_label_text_split = 0,
                         xaxis_title = f'IGPs containing different types of {i.replace("_", " ")}s',
                         xaxis_title_gap = 40,
                         xaxis_label_rotate = 0,
                         yaxis_title_gap = 50,
                         xaxis_label_font_size = 20,
                         yaxis_label_font_size = 20,
                         legend_font_size = 20,
                         yaxis_title = 'Up and downregulated(-) IGPs counts',
                         figure_description = f'Comparison of {j.replace("_", " ")}s in up- versus downregulated IGPs.',
                         )
        generated_figs.append(fig2)
        fig3 = module7.bar_up_down(data = plot_data,
                         x_column = plot_data.columns[0],
                         up_column = plot_data.columns[3],
                         down_column = plot_data.columns[4],
                         subfolder=f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}",
                         colors = ['#B64074', '#2A255C'],
                         filename = f"differential_analysis_{i}_ratio",
                         xaxis_label_text_split = 0,
                         xaxis_title = f'IGPs containing different types of {i.replace("_", " ")}s',
                         xaxis_title_gap = 40,
                         xaxis_label_rotate = 0,
                         yaxis_title_gap = 50,
                         xaxis_label_font_size = 20,
                         yaxis_label_font_size = 20,
                         legend_font_size = 20,
                         yaxis_title = 'Up and downregulated(-) IGPs ratio',
                         figure_description = f'Ratio of {j.replace("_", " ")}s in up- versus downregulated IGPs based on total number of related glycan.',
                         )
        generated_figs.append(fig3)

    # 检查加上本轮新fig后是否超出12张，如果是，先拼接之前的
    if len(figs) + len(generated_figs) > 12:
        # 拼接当前figs为一个pdf
        figure_name = f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}"
        pdf_file = f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}.pdf"
        for idx, f in enumerate(figs):
            module7.add_figure(f, figure_name=figure_name)
        custom_sizes = [[j+1] for j in range(len(figs))]
        module7.compose_figures(pdf_file, figure_name=figure_name, custom_sizes=custom_sizes)
        # 清空figs，pdf_index+1
        figs = []
        pdf_index += 1

    # 本轮生成的fig加入figs
    figs.extend(generated_figs)

# 循环结束后，把最后不足12张的也拼接保存
if figs:
    figure_name = f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}"
    pdf_file = f"StrucGAP_GlycoPeptideQuant_key_information_{pdf_index}.pdf"
    for idx, f in enumerate(figs):
        module7.add_figure(f, figure_name=figure_name)
    custom_sizes = [[j+1] for j in range(len(figs))]
    module7.compose_figures(pdf_file, figure_name=figure_name, custom_sizes=custom_sizes)
        
# key insights StrucGAP_FunctionAnnotation
sheet_names = pd.ExcelFile("test/analysis_result/StrucGAP_FunctionAnnotation_GO_ora_no_background_up_result_key_information.xlsx").sheet_names
figs = []
pdf_index = 1  
for sheet in sheet_names:
    df = pd.read_excel("test/analysis_result/StrucGAP_FunctionAnnotation_GO_ora_no_background_up_result_key_information.xlsx", sheet_name=sheet, header=None)
    df = df.iloc[:,1:]
    df = df.replace('A2B2C1D1dD1','Core-I') \
                          .replace('A2B2C1D1dD1dcbB5','Core-II') \
                          .replace('A2B2C1D1dD2dD1','Core-III') \
                          .replace('A2B2C1D1dD2dD1dcbB5','Core-IV')
    if df.shape[1] > 3:
        empty_row_idx = df[df.isnull().all(axis=1)].index[0]
        df = df.iloc[empty_row_idx+1:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df[1:].reset_index(drop=True)
        df.columns = [df.columns[0]] + [col.capitalize() for col in df.columns[1:]]
        parts = sheet.split('_', 1)
        prefix = parts[0]
        main_part = parts[1]
        if prefix == 'bp':
            j = 'GOBP'
        elif prefix == 'cc':
            j = 'GOCC'
        elif prefix == 'mf':
            j = 'GOMF'
        else:
            j = prefix.upper()
        i = main_part.replace('_', ' ')
        if i == 'core structure':
            i = 'core structures'
        elif i == 'glycan type':
            i = 'glycan types'
        elif i == 'branches structure':
            i = 'branch structures' 
        elif i == 'branches count':
            i = 'branch counts'
        elif i == 'lacdinac':
            i = 'LacdiNAc'
        elif i == 'fucosylated type':
            i = 'fucosylated types'
        elif i == 'acgc':
            i = 'sialylated types'
        s = f'Key information on {j} enrichment of the {i} substructure in upregulated IGPs.'
        fig = module7.radar('df',  
                      columns = list(df.columns[1:]),
                      text_font_size = 15,
                      legend_font_size = 12,
                      text_split = 15,
                      subfolder=f'StrucGAP_FunctionAnnotation_key_information_{pdf_index}',
                      screen_column = df.columns[0],
                      screen_values = list(df.iloc[:,0]),
                      filename = sheet,
                      figure_description = s,
                      )
        figs.append(fig)
    
        if len(figs) == 12:
            figure_name = f"StrucGAP_FunctionAnnotation_key_information_{pdf_index}"
            pdf_file = f"StrucGAP_FunctionAnnotation_key_information_{pdf_index}.pdf"
            for f in figs:
                module7.add_figure(f, figure_name=figure_name)
            custom_sizes = [[k+1] for k in range(len(figs))]
            module7.compose_figures(pdf_file, figure_name=figure_name, custom_sizes=custom_sizes)
            figs = []
            pdf_index += 1

if figs:
    figure_name = f"StrucGAP_FunctionAnnotation_key_information_{pdf_index}"
    pdf_file = f"StrucGAP_FunctionAnnotation_key_information_{pdf_index}.pdf"
    for f in figs:
        module7.add_figure(f, figure_name=figure_name)
    custom_sizes = [[k+1] for k in range(len(figs))]
    module7.compose_figures(pdf_file, figure_name=figure_name, custom_sizes=custom_sizes)




