Tutorials(pGlyco3)
==============================

A comprehensive analytical pipeline based on glycoproteomic data from pGlyco3 (data source: DOI: 10.1038/s41467-025-60437-6):

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
    
    # Only for structure data
    # Robust data quality control
    module1 = StrucGAP_Preprocess(data_dir="tests/pd structure.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'tests/sample_group.xlsx',
                      branch_list_dir = "tests/branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager,
                      search_engine = 'pGlyco3')
    module1.data_cleaning(data_type='tmt', fc_recommendation = False)
    module1.cv_raw(threshold='no')
    module1.fdr(feature_type='no')
    module1.outliers(abundance_ratio=[1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464],
                 samplewise_normalization = False)
    module1.cv(threshold = 'no')
    module1.psm(psm_number = 'no', fc_recommendation = False)
    module1.annotation(glytoucan = True, biosynthetic_pathways = True, glycobiology_filter = True)
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


    # Link quantification data with structure data
    # Robust data quality control
    module1 = StrucGAP_Preprocess(data_dir="tests/pd structure.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'tests/sample_group.xlsx',
                      branch_list_dir = "tests/branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager,
                      search_engine = 'pGlyco3')
    module1.data_cleaning(data_type='tmt',quantification_from_no_strucgp = True,
                      quantification_data_dir = "tests/s6.xlsx" , 
                      sheet_name = 'Sheet1', quant_cols = ["Young-SN-1","Young-SN-2","Young-SN-3","PD-SN-1","PD-SN-2","PD-SN-3"])
    module1.cv_raw(threshold='no', fc_recommendation = True)
    module1.fdr(feature_type='no')
    module1.outliers(abundance_ratio=[1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464],
                 samplewise_normalization = False)
    module1.cv(threshold = 'no')
    module1.psm(psm_number = 'no', fc_recommendation = True)
    module1.annotation(glytoucan = True, biosynthetic_pathways = True, glycobiology_filter = True)
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
    module4.differential_analysis(pvalue_type='pvalue_ttest', fc = 4.2)
    module4.threshold_variation_analysis(pvalue_type='pvalue_ttest',statistic_index='fc', fc_range = [3, 6, 10, 50, 100])
    module4.glycopeptide_glycosite_glycan_variation(fc = 4.2)
    module4.glycoprotein_glycosite_glycan_variation(fc = 4.2)
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
    module5.ora(organism='mmusculus', enrich_feature='glycopeptide', enrich_feature='glycopeptide', background_input=False, up_down_fc_threshold=4.2,pvalue_type='pvalue_ttest') 
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
    
Please check the 'MODULES' for more details.



