Tutorials(Glyco-Decipher)
==============================

A comprehensive analytical pipeline based on glycoproteomic data from Glyco-Decipher (data source: PXD031032):

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
    
    # only structure information
    module1 = StrucGAP_Preprocess(data_dir="tests/glycodecipher.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'tests/sample_group.xlsx',
                      branch_list_dir = "tests/branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager,
                      search_engine = 'MSFragger-Glyco')
    module1.data_cleaning(data_type='tmt')
    module1.cv_raw(threshold='no', fc_recommendation = False)
    module1.fdr(feature_type='no')
    module1.outliers(abundance_ratio=[1.240003449, 0, 1.344387558, 0, 1.576533442, 0, 1, 0, 1.956346409, 1.517000766])
    module1.cv(threshold = 'no')
    module1.psm(psm_number = 'no', fc_recommendation = False)
    module1.annotation(glytoucan = True, biosynthetic_pathways = True, glycobiology_filter = False)
    module1.output() 
    
    # Substructural features of site-specific N-glycans 
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
    module5 = StrucGAP_FunctionAnnotation(module1, 
                                 data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    module5.output()
    # Key informaiton extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

Please check the 'MODULES' for more details.



