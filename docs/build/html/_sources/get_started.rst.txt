Get Started
===========

Here is a basic example of how to use StrucGAP. We strongly recommend using the output results from StrucGP as the input for StrucGAP. Although we also provide data processing pipelines for search results from MSFragger-Glyco, pGlyco3, and Glyco-Decipher, the information contained in the outputs of these three search engines is very limited. As a result, only the StrucGAP_GlycanStructure or StrucGAP_GlycoSite modules can be used for subsequent analyses.

.. code-block:: python

    from strucgap.insighttracker import StrucGAP_InsightTracker
    from strucgap.preprocess import StrucGAP_Preprocess

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
    module1.fdr(feature_type='no')
    module1.outliers(abundance_ratio=[1.172277596,1.142983373,1,1.46390136,1.466662624,
    1.449428354,1.109519196,1.387464059,1.291746761,1.487440464])
    module1.cv(threshold = 'no')
    module1.psm()
    module1.output()
    # ... Other analysis

Each module can be instantiated and run independently depending on your workflow.
