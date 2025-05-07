StrucGAP_GlycoSite Module
=====================================================

Overview
------------

This module mapped site-level heterogeneities.

How to Use
--------------

Instantiate and use the module as follows:

.. code-block:: python

    from .glycosite import StrucGAP_GlycoSite

    module3 = StrucGAP_GlycoSite(module1, data_manager=data_manager)
    module3.glycoprotein_site()
    module3.glycopeptide_site()
    module3.specific_site()
    module3.output()

API Reference
-----------------

.. automodule:: strucgap.glycosite
    :members:
    :undoc-members:
    :show-inheritance:
