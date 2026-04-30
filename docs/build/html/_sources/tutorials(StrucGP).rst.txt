Tutorials(StrucGP)
==============================

A comprehensive analytical pipeline based on glycoproteomic data from our mouse uterus aging (the mass spectrometry data have been deposited to the ProteomeXchange Consortium via the iProX partner repository with the dataset identifier IPX0011747000):

.. code-block:: python

    from strucgap.preprocess import StrucGAP_Preprocess
    from strucgap.glycanstructure import StrucGAP_GlycanStructure
    from strucgap.glycosite import StrucGAP_GlycoSite
    from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
    from strucgap.functionannotation import StrucGAP_FunctionAnnotation
    from strucgap.glyconetwork import StrucGAP_GlycoNetwork
    from strucgap.datavisualization import StrucGAP_DataVisualization
    from strucgap.insighttracker import StrucGAP_InsightTracker

    # Initialization
    data_manager = StrucGAP_InsightTracker()
    # Setting the result storage path (folder)
    os.chdir('tests/')
    # Read it if you've already done the analysis
    data_manager.read_pickle()
    
    # Robust data quality control
    module1 = StrucGAP_Preprocess(data_dir="tests/mouse uterus.xlsx",
                      data_sheet_name = '1 PSM',
                      sample_group_data_dir = 'tests/sample_group.xlsx',
                      branch_list_dir = "tests/branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager)
    module1.data_cleaning(data_type='tmt')
    module1.cv_raw(threshold='no')
    module1.fdr(feature_type='no')
    module1.outliers(abundance_ratio=[1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464],
                 samplewise_normalization = False, total_intensity_normalization=False, total_intensity_method='mean')
    module1.cv(threshold = 'no')
    module1.psm(psm_number = 'no')
    module1.annotation(glytoucan = True, glytoucan_structure = True, glytoucan_wurcs_file = "tests/glycosmos_glycans_wurcs.csv", biosynthetic_pathways = True, glycobiology_filter = True)
    module1.output() 
    
    # Substructural features of site-specific N-glycans 
    module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager, data_type='psm_filtered')
    module2.statistics(remove_oligo_mannose = False) 
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    module2.isoforms()
    module2.output()
    
    # Glycosite information
    module3 = StrucGAP_GlycoSite(module1, data_manager=data_manager)
    module3.glycoprotein_site()
    module3.glycopeptide_site()
    module3.specific_site()
    module3.output()
    
    # Quantification analysis
    module4 = StrucGAP_GlycoPeptideQuant(module1, data_type = 'psm_filtered', data_manager=data_manager)
    module4.statistics()
    module4.statistics_index()
    module4.differential_analysis(pvalue_type='pvalue_ttest', fc = 1.7)
    module4.threshold_variation_analysis(pvalue_type='pvalue_ttest',statistic_index='fc')
    module4.glycopeptide_glycosite_glycan_variation(fc = 1.7)
    module4.glycoprotein_glycosite_glycan_variation(fc = 1.7)
    module4.output()
    # Key informaiton extraction
    data_manager.key_information_extraction(module='StrucGAP_GlycoPeptideQuant')
    
    # Functional enrichment of both glycoproteins
    module5 = StrucGAP_FunctionAnnotation(module1, 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) 
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    module5.output()
    # Key informaiton extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

    # Functional enrichment of differential glycopeptides
    module5 = StrucGAP_FunctionAnnotation(module4, 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', enrich_feature='glycopeptide', background_input=False, up_down_fc_threshold=1.7,pvalue_type='pvalue_ttest') 
    # Upregulated glycopeptides
    module5.go_function_structure(function_data = 'ora_no_background_up_result')  
    module5.output()
    # Key informaiton extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
    # Downregulated glycopeptides
    module5.go_function_structure(function_data = 'ora_no_background_down_result')  
    module5.output()
    # Key informaiton extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
    
    # Upstream and downstream interactions in glycosylation networks
    module6 = StrucGAP_GlycoNetwork(module4, data_manager=data_manager)
    module6.proteomic(protein_data_dir="tests/WYQ_Mus_uterus_global.xlsx",
                      data_sheet_name = '1 Proteins')
    module6.glycosyltransferases(glycosyltransferases_data_dir="tests/enzyme.xlsx", 
                                 data_sheet_name="glycosyltransferases")
    module6.glycosidases(glycosidases_data_dir="tests/enzyme.xlsx", 
                         data_sheet_name='glycosidases')
    module6.sialyltransferases()
    module6.fucosyltransferase()
    module6.glycan_binding_protein()
    module6.output()
    
    # Functional enrichment of global proteins
    module5 = StrucGAP_FunctionAnnotation(module6, 
                                     data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) 
    module5.go_function_structure(function_data = 'ora_no_background_up_result')  
    module5.output()
    module5.go_function_structure(function_data = 'ora_no_background_down_result')  
    module5.output()

    # Output both results as pkl file
    data_manager.output_pickle()
    
    # Plotting
    module7 = StrucGAP_DataVisualization(data_manager=data_manager)
    # Example
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
    # Example
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
    # Example
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
    # ... Other analysis

Please check the 'MODULES' for more details.



