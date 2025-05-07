StrucGAP_FunctionAnnotation Module
==================================================

Overview
------------

This module provides functional interpretations of altered glycan traits.

How to Use
----------------

Instantiate and use the module as follows:

.. code-block:: python

    from .functionannotation import StrucGAP_FunctionAnnotation

    # both glycoproteins
    module5 = StrucGAP_FunctionAnnotation(module1, data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5)
    module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
    module5.output()
    # key information extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
    # upregulated glycopeptides
    module5 = StrucGAP_FunctionAnnotation(module4, data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=  1.5,pvalue_type='pvalue_ttest') 
    module5.go_function_structure(function_data = 'ora_no_background_up_result')  
    module5.output()
    # key information extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
    # downregulated glycopeptides
    module5.go_function_structure(function_data = 'ora_no_background_down_result')  
    module5.output()
    # key information extraction
    data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
    # both proteins
    module5 = StrucGAP_FunctionAnnotation(module6,  data_manager=data_manager)  
    module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5)
    module5.go_function_structure(function_data = 'ora_no_background_up_result')  
    module5.output()
    module5.go_function_structure(function_data = 'ora_no_background_down_result')  
    module5.output()

API Reference
------------------

.. automodule:: strucgap.functionannotation
    :members:
    :undoc-members:
    :show-inheritance:
