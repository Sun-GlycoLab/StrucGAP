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


API Reference
------------------

.. automodule:: strucgap.datavisualization
    :members:
    :undoc-members:
    :show-inheritance:
