StrucGAP_DataVisualization Module
================================================

Overview
-------------

This module offer 38 customizable plotting options. It also features draw_glycan(), which visually renders glycan structures from glycan coding, aiding structural interpretation. The module supports automated generation of PDF reports summarizing results and visualizations.

How to Use
---------------

Instantiate and use the module as follows:

.. code-block:: python

    from strucgap.datavisualization import StrucGAP_DataVisualization

    module7 = StrucGAP_DataVisualization(data_manager=data_manager)
    # pie chart
    module7.pie(data = module2.core_structure, 
              item_column = 'Core_structure', 
              number_column = 'Core_structure_count',
              radius = ['0%', '70%'],
              rosetype = None,
              subfolder='glycanstructure',
              filename='corestructure',
              colors = ['#F8ED70','#EB8B5B','#B53656','#632A69']
              )
    # draw glycans
    module7.draw_glycans(module6.protein_up_glyco_down['structure_coding'], 
                     linewidth=0.2, 
                     subfolder='StrucGAP_GlycoNetwork',
                     filename="protein_up_glyco_down")
    # analysis report generation

Analytical report generation
--------------------------------------------------------------

We generate two types of analysis reports by default based on the research content. The first type includes results from all analysis modules, while the second type provides a quick overview of key insights. Below, we provide an example of generating such a report. The complete script for report generation can be found in our GitHub repository(https://github.com/Sun-GlycoLab/StrucGAP) under tests/analysis_report.py:

.. code-block:: python

    # Figure 1
    fig1 = module7.heatmap_multi_data('module1.data',  
                  'module1.data_psm_filtered',   
                  columns=['MS2Scan', 'ProteinID', 'PeptideSequence', 
                           'ProteinID+Glycosite_Position','structure_coding', 'GlycanComposition'],
                  statistical_methods = ['both', 'unique', 'unique', 'unique', 
                                         'unique', 'unique'],
                  subfolder='figure1',
                  colors = 'Spectral',
                  figure_description = 'Overall profiling of unique IGPs, glycoproteins, peptides, glycosites, glycan structures and glycan composition identified in cleaned data and fdr filtered data',
                  filename='overview heatmap',
                  annotation_font_size = 25,
                  xaxis_label_font_size = 25,
                  )
    
    fig2 = module7.nested_pie(data = module3.glycoprotein_glycosite_count, 
                       item_column = 'glycoprotein', 
                       number_column = 'glycosite_count',
                       value_counts_column = 'glycosite_count',
                       subfolder='figure1',
                       label_font_size = 20,
                       legend_font_size = 20,
                       figure_description = 'Number of glycosites identified on each glycoprotein',
                       filename='glycosite count',
                       split=10)
    
    module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
    module2.statistics(remove_oligo_mannose = False) 
    data1 = module2.structure_coding_rank
    module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
    module2.statistics(remove_oligo_mannose = True)
    data2 = module2.structure_coding_rank 
    fig3 = module7.butterfly_plot(data1,data2,
                           item_column='Structure_coding',
                           count_column='Structure_coding_count',
                           colors = ['lightblue', 'darkblue'],
                           xaxis_title = 'Number of Glycan',
                           plot_title = 'Top 10 glycans',
                           legend = ['With oligo mannose', 'Without oligo mannose'],
                           subfolder='figure1',
                           filename='top 10 structure coding',
                           label_font_size = 30,
                           xaxis_title_font_size = 30,
                           plot_title_font_size = 0,
                           legend_fontsize = 30,
                           legend_loc=[0,-0.4], 
                           figure_description = 'Top 10 glycan structures identified based on the number of their modified N-glycosites',
                           )
    
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
    fig4 = module7.polar1('plot_data', 
                  columns=['isomer_count'],
                  number_column = 'GlycanComposition',
                  subfolder='figure1',
                  if_unique = False,
                  radiusaxis_label_show = True,
                  colors = ['#81B29A'],
                  plot_title = None,
                  filename='isomer',
                  radiusaxis_label_font_size = 20,
                  angleaxis_label_font_size = 20,
                  legend_font_size = 20,
                  legend = 'Isomer counts',
                  figure_description = 'Overveiw of glycan structure isomers identified',
                  )
    
    plot_data = module2.core_structure
    plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
    plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
    plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5', 'Core-III')
    plot_data = plot_data.replace('A2B2C1D1dD2dD1', 'Core-IV')
    fig5 = module7.pie(data = plot_data, 
                  item_column = 'Core_structure', 
                  number_column = 'Core_structure_count',
                  radius = ['0%', '70%'],
                  rosetype = None,
                  subfolder='figure1',
                  filename='corestructure',
                  label_font_size = 20,
                  legend_font_size = 20,
                  figure_description = 'Porprotion of core structures',
                  colors = ['#F8ED70','#EB8B5B','#B53656','#632A69']
                  )
    
    fig6 = module7.pie(data = module2.glycan_type, 
                  item_column = 'Glycan_type', 
                  number_column = 'Glycan_type_count',
                  radius = ['0%', '70%'],
                  rosetype = 'area',
                  subfolder='figure1',
                  filename='glycantype',
                  label_font_size = 20,
                  legend_font_size = 20,
                  figure_description = 'Porprotion of glycan types',
                  colors = ['#F8ED70','#EB8B5B','#B53656','#632A69']
                  )
    
    fig7 = module7.bar(data = module2.branches_structure,
                top = 10,
                y_column='Branches',
                y_column_value='Branches_count',
                subfolder='figure1',
                colors = '#8FBFB8',
                transform_ratio = True,
                y_max = None,
                yaxis_splitline_show = False,
                xaxis_splitline_show = False,
                legend = 'Branches structure',
                xaxis_title = 'Branches structure',
                xaxis_title_gap = 140,
                xaxis_label_rotate = -45,
                xaxis_label_margin = 25,
                xaxis_label_text_split = 20,
                yaxis_title = 'Percentage (%)',
                figure_description = 'Porprotion of both branch structures',
                yaxis_title_gap = 70,
                xaxis_label_font_size = 30,
                yaxis_label_font_size = 30,
                legend_font_size = 30,
                filename='branches structure'
                )
    
    fig8 = module7.bar(data = module2.branches_structure,
                top = None,
                end = 8,
                y_column='Branches',
                y_column_value='Branches_count',
                subfolder='figure1',
                colors = '#8FBFB8',
                transform_ratio = False,
                y_max = None,
                yaxis_splitline_show = False,
                xaxis_splitline_show = False,
                legend = 'Branches structure',
                xaxis_title = 'Branches structure',
                xaxis_title_gap = 170,
                xaxis_label_rotate = -45,
                xaxis_label_margin = 25,
                xaxis_label_text_split = 20,
                yaxis_title = 'Count',
                figure_description = 'The least count of eight branch structures',
                yaxis_title_gap = 50,
                xaxis_label_font_size = 30,
                yaxis_label_font_size = 30,
                legend_font_size = 30,
                filename='end8 branches structure'
                )
    
    fig9 = module7.bar(data = module2.branches_count,
                top = 10,
                y_column='BranchNumber',
                y_column_value='BranchNumber_count',
                xaxis_label_rotate = 0,
                xaxis_label_margin = 20,
                transform_ratio = True,
                subfolder='figure1',
                colors = '#FFD804',
                y_max = None,
                yaxis_splitline_show = False,
                legend = 'Branch number',
                filename='branch number',
                xaxis_title = 'Branch number of glycans',
                xaxis_title_gap = 50,
                yaxis_title = 'Percentage of IGPs',
                yaxis_title_gap = 60,
                xaxis_label_font_size = 30,
                yaxis_label_font_size = 30,
                legend_font_size = 30,
                figure_description = 'Number of branches per glycan',
                )
    
    data1 = module1.data_psm_filtered
    data1 = data1[(data1['fucosylated type']=='core fucosylated')|(data1['fucosylated type']=='dual')]['PeptideSequence+structure_coding+ProteinID']
    data1 = data1.iloc[:,0]
    data2 = module1.data_psm_filtered
    data2 = data2[(data2['fucosylated type']=='antenna fucosylated')|(data2['fucosylated type']=='dual')]['PeptideSequence+structure_coding+ProteinID']
    data2 = data2.iloc[:,0]
    fig10 = module7.venn_diagram(
        data1,
        data2,
        colors = 'Tropic',
        subfolder='figure1',
        legend = ['Core fucosylated', 'Antenna fucosylated'],
        filename='fucosylated type',
        figure_description = 'Composition of fucosylated type',
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
    fig11 = module7.venn_diagram(
        data1,
        data2,
        colors = 'Tropic',
        subfolder='figure1',
        legend = ['Neu5Ac', 'Neu5Gc'],
        filename='acgc',
        figure_description = 'Composition of sialylated type',
        plot_title_font_size = 25,
        legend_fontsize = 20,
        number_fontsize = 25,
        legend_loc='lower center', 
    )
    
    module7.add_figure(fig1, figure_name="figure1")
    module7.add_figure(fig2, figure_name="figure1")
    module7.add_figure(fig3, figure_name="figure1")
    module7.add_figure(fig4, figure_name="figure1")
    module7.add_figure(fig5, figure_name="figure1")
    module7.add_figure(fig6, figure_name="figure1")
    module7.add_figure(fig7, figure_name="figure1")
    module7.add_figure(fig8, figure_name="figure1")
    module7.add_figure(fig9, figure_name="figure1")
    module7.add_figure(fig10, figure_name="figure1")
    module7.add_figure(fig11, figure_name="figure1")
    module7.compose_figures("figure1.pdf", figure_name="figure1",
                            custom_sizes=[[1], [2], [4,5], [3], [6], [7], [8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据
    
  

API Reference
------------------

.. automodule:: strucgap.datavisualization
    :members:
    :undoc-members:
    :show-inheritance:
