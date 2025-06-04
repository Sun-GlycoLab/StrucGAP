StrucGAP_GlycanStructure Module
=============================================

Overview
------------

This module delineates N-glycan structural features, profiling their distributions and compositions.

How to Use
--------------

Instantiate and use the module as follows:

.. code-block:: python

    from strucgap.glycanstructure import StrucGAP_GlycanStructure

    module2 = StrucGAP_GlycanStructure(gs_data=module1, data_manager=data_manager,  data_type='psm_filtered')
    module2.statistics(remove_oligo_mannose = False) 
    module2.structure_statistics()
    module2.lacdinac()
    module2.cor()
    module2.isoforms()
    module2.output()

API Reference
-----------------

.. automodule:: strucgap.glycanstructure
    :members:
    :undoc-members:
    :show-inheritance:
