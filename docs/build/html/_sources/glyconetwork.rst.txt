StrucGAP_GlycoNetwork Module
==================================================

Overview
------------

This module analyzes upstream and downstream interactions within the glycosylation process, placing site-specific N-glycan alterations in the context of broader biological networks.

How to Use
--------------

Instantiate and use the module as follows:

.. code-block:: python

    from .glyconetwork import StrucGAP_GlycoNetwork

    module6 = StrucGAP_GlycoNetwork(module4, data_manager=data_manager)
    module6.proteomic(protein_data_dir="D:\\doctor\\wyq\\WYQ_Mus_uterus_global.xlsx",
                  data_sheet_name = '1 Proteins')
    module6.glycosyltransferases(glycosyltransferases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", data_sheet_name="glycosyltransferases")
    module6.glycosidases(glycosidases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", data_sheet_name='glycosidases')
    module6.sialyltransferases()
    module6.fucosyltransferase()
    module6.glycan_binding_protein()
    module6.output()

API Reference
-----------------

.. automodule:: strucgap.glyconetwork
    :members:
    :undoc-members:
    :show-inheritance:
