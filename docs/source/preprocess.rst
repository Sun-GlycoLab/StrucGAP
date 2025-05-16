StrucGAP_Preprocess Module
==========================

Overview
-----------

This module processes outputs from StrucGP alongside sample metadata and a glycan branch library through sequential correction and filtering steps, harmonizing heterogeneous data into a unified format.

How to Use
--------------

Instantiate and use the module as follows:

.. code-block:: python

    from strucgap.preprocess import StrucGAP_Preprocess

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

API Reference
------------------

.. automodule:: strucgap.preprocess
    :members:
    :undoc-members:
    :show-inheritance:
