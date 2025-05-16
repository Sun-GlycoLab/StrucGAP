StrucGAP_DataVisualization Module
================================================

Overview
-------------

This module offer 38 customizable plotting options. It also features draw_glycan(), which visually renders glycan structures from glycan coding, aiding structural interpretation. The module supports automated generation of PDF reports summarizing results and visualizations.

How to Use
---------------

Instantiate and use the module as follows:

.. code-block:: python

    from .datavisualization import StrucGAP_DataVisualization

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

By default, we generate a four‑page analysis report here, containing basic glycan substructure feature analysis, functional analysis, quantitative analysis, and upstream/downstream analysis. Moreover, For an analysis report containing all the analysis results, see the analysis report.py file on GitHub. Essentially, this invokes the plotting functions and analysis results of the StrucGAP_DataVisualization module, so users can further adjust it to meet their advanced requirements:

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
    
    
    # Figure 2
    module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,selected_terms=['GO:CC']) # 69,76,83
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    plot_data = module5.ora_no_background_both_proteins_result
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
        subfolder='figure2',
        filename='both glycoproteins gocc',
        figure_description = 'Subcellular localization of identified glycoproteins based on GO:CC enrichment',
    )
    
    module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,selected_terms=['GO:BP']) # 69,76,83
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    plot_data = module5.ora_no_background_both_proteins_result
    def clean_term(term):
        term = re.sub(r'\s*\(.*?\)', '', term)  
        term = term.title()  
        return term
    plot_data['Term'] = plot_data['Term'].apply(clean_term)
    fig2 = module7.dotplot_col(
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
        subfolder='figure2',
        filename='both glycoproteins gobp',
        figure_description = 'Subcellular localization of identified glycoproteins based on GO:BP enrichment',
    )
    
    module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,selected_terms=['GO:MF']) # 69,76,83
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    plot_data = module5.ora_no_background_both_proteins_result
    def clean_term(term):
        term = re.sub(r'\s*\(.*?\)', '', term)  
        term = term.title()  
        return term
    plot_data['Term'] = plot_data['Term'].apply(clean_term)
    fig3 = module7.dotplot_col(
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
        subfolder='figure2',
        filename='both glycoproteins gomf',
        figure_description = 'Subcellular localization of identified glycoproteins based on GO:MF enrichment',
    )
    
    module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,selected_terms=['GO:CC']) # 69,76,83
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    plot_data = module5.cc_core_structure
    plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
    plot_data = plot_data.replace('A2B2C1D1dD1', 'Core-I')
    plot_data = plot_data.replace('A2B2C1D1dD1dcbB5', 'Core-II')
    plot_data = plot_data.replace('A2B2C1D1dD2dD1dcbB5', 'Core-III')
    plot_data = plot_data.replace('A2B2C1D1dD2dD1', 'Core-IV')
    fig4 = module7.line(data = plot_data,
                 colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
                 y_column = 'Core_structure',
                 x_columns = plot_data.columns[1:11],
                 subfolder='figure2',
                 symbol_size = 5,
                 plot_title = None,
                 xaxis_label_rotate = -15,
                 xaxis_title_gap = 35,
                 xaxis_label_font_size = 10,
                 xaxis_label_text_split = 20,
                 yaxis_label_font_size = 20,
                 legend_font_size = 20,
                 yaxis_title_gap = 50,
                 yaxis_title = 'Percentage (%)',
                 filename='gocc core structure',
                 figure_description = 'Distribution of core structures across the top 10 GO:CC-enriched terms',
                 )
    
    plot_data = module5.cc_glycan_type
    plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
    fig5 = module7.bar_multi_columns(
        data = plot_data,  
        colors = ['#4285F4','#F4B400','#DB4437','#OF9D58'],
        y_column = "Glycan_type",
        x_columns = plot_data.columns[1:11], 
        subfolder='figure2',
        xaxis_splitline_show = False,
        yaxis_splitline_show = False,
        xaxis_label_font_size = 10,
        xaxis_label_text_split = 20,
        yaxis_label_font_size = 20,
        legend_font_size = 20,
        y_max=1,
        xaxis_title_gap = 35,
        yaxis_title = 'Percentage (%)',
        filename='gocc glycan type',
        figure_description = 'Distribution of glycan types across the top 10 GO:CC-enriched terms',
    )
    
    plot_data = module5.cc_branches_structure
    plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
    fig6 = module7.multi_bar(
        data = plot_data,  
        x_column = "Branches",
        y_column = plot_data.columns[:11], 
        subfolder='figure2',
        xaxis_splitline_show = False,
        yaxis_splitline_show = False,
        bar_width="50%",
        xaxis_label_font_size = 10,
        xaxis_label_text_split = 10,
        yaxis_label_font_size = 10,
        yaxis_label_margin=35,
        filename='gocc branch structure',
        figure_description = 'Distribution of branch structures across the top 10 GO:CC-enriched terms',
    )
    
    module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,selected_terms=['GO:BP']) # 69,76,83
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    plot_data = module5.bp_fucosylated_type.copy()
    plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
    fig7 = module7.radar('plot_data',  
                  columns = list(plot_data.columns[1:11]),
                  text_font_size = 15,
                  legend_font_size = 15,
                  text_split = 15,
                  subfolder='figure2',
                  screen_column = 'Fucosylated_type',
                  screen_values = ['core fucosylated', 'antenna fucosylated', 'dual'],
                  filename='gobp fucosylated type',
                  figure_description = 'Distribution of fucosylated types across the top 10 GO:BP-enriched terms',
                  )
    
    plot_data = module5.bp_acgc.copy()
    plot_data.columns = [plot_data.columns[0]] + [col.capitalize() for col in plot_data.columns[1:]]
    fig8 = module7.radar('plot_data',  
                  columns = list(plot_data.columns[1:11]),
                  text_font_size = 15,
                  legend_font_size = 15,
                  text_split = 15,
                  subfolder='figure2',
                  screen_column = 'Ac/Gc',
                  screen_values = ['Ac', 'Gc', 'dual'],
                  filename='gobp sialylated type',
                  figure_description = 'Distribution of sialylated types across the top 10 GO:BP-enriched terms',
                  )
    
    module7.add_figure(fig1, figure_name="figure2")
    module7.add_figure(fig2, figure_name="figure2")
    module7.add_figure(fig3, figure_name="figure2")
    module7.add_figure(fig4, figure_name="figure2")
    module7.add_figure(fig5, figure_name="figure2")
    module7.add_figure(fig6, figure_name="figure2")
    module7.add_figure(fig7, figure_name="figure2")
    module7.add_figure(fig8, figure_name="figure2")
    # module7.add_figure(fig9, figure_name="figure2")
    # module7.add_figure(fig10, figure_name="figure2")
    # module7.add_figure(fig11, figure_name="figure2")
    module7.compose_figures("figure2.pdf", figure_name="figure2",
                            custom_sizes=[[1], [2], [3], [4], [5], [6,9,12], [7], [8]])  # 生成后自动清理figure1数据
    
    # Figure 3
    fig1 = module7.volcano_plot(data = module4.fc_result,
                        fc_column = 'fc',
                        p_column = 'pvalue_ttest',
                        fc = 1.5,
                        p_value = 0.05,
                        subfolder='figure3',
                        figure_description = 'Volcano plot of differentially expressed glycopeptides in the aging uterus dataset, showing significant up- and downregulated IGPs (FC > 1.5 or < 0.67, P value < 0.05)',
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
        subfolder='figure3',
        figure_description = 'PCA plot based on glycopeptide expression profiles, demonstrating sample separation between age groups',
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
                     xaxis_title_font_size=25, 
                     yaxis_title_font_size=25,
                     xaxis_label_font_size=25, 
                     yaxis_label_font_size=25, 
                     subfolder='figure3',
                     figure_description = 'Heatmap of differentially expressed glycopeptides, revealing global expression trends and group-wise clustering',
                     )
    
    fig4 = module7.frequency_bar(data = module6.pg_fc, 
                          columns = ['fc_g', 'fc_p'], 
                          log2_transformation = True,
                          ref_lines_value = 1.5, 
                          xaxis_label_font_size=20, 
                          yaxis_label_font_size=20, 
                          legend_fontsize = 20,
                          colors=['red', 'blue'],  
                          xaxis_title = 'Fc',
                          yaxis_title = 'Frequency',
                          plot_title = 'Protein and glycopeptide fc',
                          legend = ['Glycopeptide', 'Protein'],
                          subfolder='figure3',
                          figure_description = 'Frequency distribution of FC for glycopeptides and global proteins',
                          )
    
    up_data = module4.up_data[module4.up_data['Glycan_type']!='Oligo mannose']
    down_data = module4.down_data[module4.down_data['Glycan_type']!='Oligo mannose']
    fig5 = module7.heatmap_multi_data(up_data,  
                  down_data,   
                  columns=['structure_coding'],
                  statistical_methods = ['count'],
                  subfolder='figure3',
                  colors = 'coolwarm',
                  filename = 'top 10 differential glycan',
                  figure_description = 'Top 10 most frequently observed upregulated and downregulated glycan structures, excluding oligo-mannose types, highlighting distinct compositional shifts',
                  )
    
    module5 = StrucGAP_FunctionAnnotation(module4, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,pvalue_type='pvalue_ttest',
                selected_terms=['GO:MF']) # 65,72,79
    plot_data1 = module5.ora_no_background_up_result.copy()
    plot_data1 = plot_data1[plot_data1['Gene_set']=='GO:MF']
    plot_data1 = plot_data1[plot_data1['P-value']<0.05]
    plot_data1['type'] = 'up'
    plot_data2 = module5.ora_no_background_down_result.copy()
    plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:MF']
    plot_data2 = plot_data2[plot_data2['P-value']<0.05]
    plot_data2['type'] = 'down'
    plot_data = pd.concat([plot_data1, plot_data2], axis=0)
    def clean_term(term):
        term = re.sub(r'\s*\(.*?\)', '', term)  
        term = term.title()  
        return term
    plot_data['Term'] = plot_data['Term'].apply(clean_term)
    fig6 = module7.dotplot_col(
        data = plot_data,  
        dot_cmap = 'Hawaii',
        category = 'type',
        p_column = 'P-value',
        top = 10,
        xaxis_font_size = 10,
        yaxis_font_size = 10,
        term = 'Term',
        dot_color_column = 'P-value',
        dot_size_column = 'Overlap',
        filename='differential peptide',
        subfolder='figure3',
        figure_description = 'Differential GO:MF enrichment of upregulated and downregulated glycopeptides (P value < 0.05, FC > 1.5 or < 0.67)',
    )
    
    fig7 = module7.bar_up_down_ratio(feature='core_structure', 
                              colors=['blue', 'green', 'red','yellow'],
                              subfolder='figure3',
                              filename="core structure",
                              figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends of core structures',
                              )
    
    fig8 = module7.bar_up_down_ratio(feature='glycan_type', 
                              colors=['blue', 'green', 'red'],
                              subfolder='figure3',
                              filename="glycan type",
                              figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends of glycan types',
                              )
    
    fig9 = module7.bar_up_down_ratio(feature='fucosylated_type', 
                              colors=['blue', 'green', 'red'],
                              subfolder='figure3',
                              filename="fucosylated type",
                              figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends of fucosylation',
                              )
    
    fig10 = module7.bar_up_down_ratio(feature='acgc', 
                              colors=['blue', 'green', 'red'],
                              subfolder='figure3',
                              filename="acgc",
                              figure_description = 'Threshold variation analysis demonstrating dynamic enrichment trends of core sialylation',
                              )
    
    module7.add_figure(fig1, figure_name="figure3")
    module7.add_figure(fig2, figure_name="figure3")
    module7.add_figure(fig3, figure_name="figure3")
    module7.add_figure(fig4, figure_name="figure3")
    module7.add_figure(fig5, figure_name="figure3")
    module7.add_figure(fig6, figure_name="figure3")
    module7.add_figure(fig7, figure_name="figure3")
    module7.add_figure(fig8, figure_name="figure3")
    module7.add_figure(fig9, figure_name="figure3")
    module7.add_figure(fig10, figure_name="figure3")
    # module7.add_figure(fig11, figure_name="figure2")
    module7.compose_figures("figure3.pdf", figure_name="figure3",
                            custom_sizes=[[1], [2], [3], [4], [5,6], [7,8], [9], [10], [11], [12]])  # 生成后自动清理figure1数据
    
    # Figure 4
    glyco = module4.fc_result[module4.fc_result['pvalue_ttest'] < 0.05]
    globaldata = module6.proteomic_fc[module6.proteomic_fc['pvalue_ttest'] < 0.05]
    fig1 = module7.up_down_scatter(
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
        subfolder='figure4',
        figure_description = 'Scatter plot showing the number of differentially expressed glycopeptides and proteins identified from glycoproteomic and global proteomic datasets, respectively (P value < 0.05)',
    )
    
    plot_data1 = module4.differential_analysis_data.copy()
    plot_data1 = plot_data1[(plot_data1['fc']>1.5)|(plot_data1['fc']<1/1.5)]
    plot_data2 = module6.proteomic_fc[module6.proteomic_fc['pvalue_ttest'] < 0.05].copy()
    plot_data2 = plot_data2[(plot_data2['fc']>1.5)|(plot_data2['fc']<1/1.5)]
    fig2 = module7.venn_diagram(
        list(plot_data1["ProteinID"]),
        plot_data2.index,
        subfolder='figure4',
        legend = ['Glycoproteins', 'Proteins'],
        filename='glyco and proteomic protein',
        figure_description = 'Composition of sialylated type',
        plot_title_font_size = 25,
        legend_fontsize = 20,
        number_fontsize = 25,
        legend_loc='lower center', 
    )
    
    plot_data1 = module4.differential_analysis_data.copy()
    plot_data1 = plot_data1[(plot_data1['fc']>1.5)|(plot_data1['fc']<1/1.5)]
    plot_data1 = plot_data1.iloc[:,:-6]
    module5 = StrucGAP_FunctionAnnotation(plot_data1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                      data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,
                selected_terms=['GO:CC']) # 69,76,83
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
    module5 = StrucGAP_FunctionAnnotation(plot_data2, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
                                      data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5,
                selected_terms=['GO:CC']) # 69,76,83
    plot_data2 = module5.ora_no_background_both_proteins_result
    plot_data2 = plot_data2[plot_data2['Gene_set']=='GO:CC']
    plot_data2 = plot_data2[plot_data2['P-value']<0.05]
    plot_data2['type'] = 'protein'
    plot_data2 = plot_data2.iloc[:10,:]
    plot_data = pd.concat([plot_data1, plot_data2], axis=0)
    def clean_term(term):
        term = re.sub(r'\s*\(.*?\)', '', term)  
        term = term.title()  
        return term
    plot_data['Term'] = plot_data['Term'].apply(clean_term)
    fig3 = module7.dotplot_col(
        data = plot_data,  
        dot_cmap = 'Hawaii',
        category = 'type',
        p_column = 'P-value',
        top = 10,
        xaxis_font_size = 10,
        yaxis_font_size = 10,
        term = 'Term',
        dot_color_column = 'P-value',
        dot_size_column = 'Overlap',
        filename='igp protein gocc',
        subfolder='figure4',
        figure_description = 'GO:CC enrichment analysis comparing subcellular localization of differentially expressed glycoproteins and global proteins',
    )
    
    plot_data1 = module6.protein_no_glyco_up[['fc_g','fc_p','normalized_fc_g']]
    plot_data1['group'] = 'up'
    plot_data2 = module6.protein_no_glyco_down[['fc_g','fc_p','normalized_fc_g']]
    plot_data2['group'] = 'down'
    plot_data = pd.concat([plot_data1,plot_data2],axis=0)
    plot_data.reset_index(inplace=True,drop=True)
    fig4 = module7.scatter(data = plot_data,
                 group_column = 'group',
                 x_column = 'fc_p',
                 y_column = 'fc_g',
                 subfolder='figure4',
                 top_xaxis_line_show = False,
                 right_yaxis_line_show = False,
                 xaxis_splitline_show = False,
                 yaxis_splitline_show = False,
                 xaxis_label_rotate = 0,
                 xaxis_title = 'FC of proteins',
                 xaxis_title_gap = 40,
                 yaxis_title = 'FC of IGPs',
                 yaxis_title_gap = 30,
                 xaxis_label_font_size = 30,
                 yaxis_label_font_size = 30,
                 figure_description = 'Identification of N-glycan features that were altered solely at the glycopeptide level, independent of protein-level changes, suggesting glycosylation-specific regulation',
                 )
    
    plot_data = module6.glycosyltransferases.copy()
    plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
    plot_data.reset_index(inplace=True)
    plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
    plot_data.set_index('gene_id',inplace=True)
    fig5 = module7.heatmap2(data = plot_data,
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
                     subfolder='figure4',
                     figure_description = 'Quantitative profiling of 40 altered glycosyltransferases (P value < 0.05)',
                     )
    
    plot_data = module6.glycosidases.copy()
    plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
    plot_data = plot_data.sort_values(by='fc', ascending=False)
    top2 = list(plot_data.index)[:2]
    plot_data = module6.cv_filter_data.reset_index().copy()
    plot_data = plot_data[plot_data['Accession'].isin(top2)]
    plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
    plot_data.set_index('gene_id',inplace=True)
    fig6 = module7.boxplot(data = plot_data,
                      item_column = 'Accession',
                      item_name = [top2[0]],
                      group1_columns = plot_data.columns[1:6],
                      group2_columns = plot_data.columns[6:11],
                      p_data = module6.proteomic_fc,
                      p_column = 'pvalue_ttest',
                      filename = plot_data.index[0],
                      subfolder='figure4',
                      yaxis_title = 'Quantification value',
                      figure_description = f'Significantly upregulated glycosidases {plot_data.index[0]}',
                      )
    
    fig7 = module7.boxplot(data = plot_data,
                      item_column = 'Accession',
                      item_name = [top2[1]],
                      group1_columns = plot_data.columns[1:6],
                      group2_columns = plot_data.columns[6:11],
                      p_data = module6.proteomic_fc,
                      p_column = 'pvalue_ttest',
                      filename = plot_data.index[1],
                      subfolder='figure4',
                      yaxis_title = 'Quantification value',
                      figure_description = f'Significantly upregulated glycosidases {plot_data.index[1]}',
                      )
    
    plot_data1 = module6.cv_filter_data[module6.cv_filter_data.index.isin(module6.sialyltransferases.index)]
    plot_data2 = module6.sialyltransferases[['fc']]
    plot_data3 = module6.sialyltransferases[['pvalue_ttest']]
    plot_data1.columns = ['Control 1','Control 2','Control 3','Control 4','Control 5',
                          'Sample 1','Sample 2','Sample 3','Sample 4','Sample 5']
    fig8 = module7.complexheatmap(data = plot_data1, 
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
                           subfolder='figure4',
                           col_split=None,
                           cmap = 'Blues', 
                           z_score = 0,
                           show_rownames = True,
                           row_split = None,
                           filename = 'sialyltransferases',
                           linewidths = 3,
                           figure_description = 'Heatmap of identified sialyltransferases',
                           )
    
    plot_data1 = module6.cv_filter_data[module6.cv_filter_data.index.isin(module6.fucosyltransferase.index)]
    plot_data2 = module6.fucosyltransferase[['fc']]
    plot_data3 = module6.fucosyltransferase[['pvalue_ttest']]
    plot_data1.columns = ['Control 1','Control 2','Control 3','Control 4','Control 5',
                          'Sample 1','Sample 2','Sample 3','Sample 4','Sample 5']
    fig9 = module7.complexheatmap(data = plot_data1, 
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
                           subfolder='figure4',
                           col_split=None,
                           cmap = 'Blues', 
                           z_score = 0,
                           show_rownames = True,
                           row_split = None,
                           filename = 'fucosyltransferase',
                           linewidths = 3,
                           figure_description = 'Heatmap of identified fucosyltransferases',
                           )
    
    plot_data = module6.glycan_binding_protein.copy()
    plot_data = plot_data[plot_data['pvalue_ttest']<0.05]
    plot_data = plot_data[(plot_data['fc']>1.5)|(plot_data['fc']<1/1.5)]
    plot_data.reset_index(inplace=True)
    plot_data = module6.convert_accession_to_gene(plot_data, "Accession", species=10090)
    plot_data.set_index('gene_id',inplace=True,drop=False)
    fig10 = module7.violin_plot(data=plot_data,
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
                        subfolder='figure4',
                        filename = 'glycan binding protins',
                        xaxis_label_font_size=20, 
                        yaxis_label_font_size=20, 
                        figure_description = 'Expression patterns of significantly altered glycan-binding proteins (P value < 0.05, FC > 1.5 or < 0.67)',
                        ) 
    
    module7.add_figure(fig1, figure_name="figure4")
    module7.add_figure(fig2, figure_name="figure4")
    module7.add_figure(fig3, figure_name="figure4")
    module7.add_figure(fig4, figure_name="figure4")
    module7.add_figure(fig5, figure_name="figure4")
    module7.add_figure(fig6, figure_name="figure4")
    module7.add_figure(fig7, figure_name="figure4")
    module7.add_figure(fig8, figure_name="figure4")
    module7.add_figure(fig9, figure_name="figure4")
    module7.add_figure(fig10, figure_name="figure4")
    # module7.add_figure(fig11, figure_name="figure2")
    module7.compose_figures("figure4.pdf", figure_name="figure4",
                            custom_sizes=[[1], [2], [4], [5], [3,6], [7], [8], [9], [12], [10,11]])  # 生成后自动清理figure1数据


API Reference
------------------

.. automodule:: strucgap.datavisualization
    :members:
    :undoc-members:
    :show-inheritance:
