StrucGAP_GlycoPeptideQuant Module
==============================================================

Overview
------------

This module is a comprehensive module for differential analysis in glycoproteomics, utilizing various statistical methods and machine learning techniques to identify and quantify dynamic changes in N-glycan substructures and glycosite-specific alterations.

How to Use
--------------

Instantiate and use the module as follows:

.. code-block:: python

    from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant

    module4 = StrucGAP_GlycoPeptideQuant(module1, data_type = 'psm_filtered',   data_manager=data_manager)
    module4.statistics()
    module4.statistics_index()
    module4.differential_analysis(pvalue_type='pvalue_ttest')
    module4.threshold_variation_analysis(pvalue_type='pvalue_ttest', statistic_index='fc')
    module4.glycopeptide_glycosite_glycan_variation()
    module4.glycoprotein_glycosite_glycan_variation()
    module4.output()
    # key information extraction
    data_manager.key_information_extraction(module='StrucGAP_GlycoPeptideQuant')

API Reference
-----------------

.. automodule:: strucgap.glycopeptidequant
    :members:
    :undoc-members:
    :show-inheritance:
