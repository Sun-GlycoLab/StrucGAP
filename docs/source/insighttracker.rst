StrucGAP_InsightTracker Module
=================================================

Overview
-----------

This module consolidates outputs from two prior modules and distills them into biologically meaningful patterns.

How to Use
--------------

Instantiate and use the module as follows:

.. code-block:: python

    from .insighttracker import StrucGAP_InsightTracker
    # before the both code
    data_manager = StrucGAP_InsightTracker()
    # output both results in pkl file
    data_manager.output_pickle()
    # input both results
    data_manager.read_pickle()
    # key information extraction
    data_manager.key_information_extraction(module='StrucGAP_GlycoPeptideQuant')
    # key information extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

API Reference
------------------

.. automodule:: strucgap.insighttracker
    :members:
    :undoc-members:
    :show-inheritance:
