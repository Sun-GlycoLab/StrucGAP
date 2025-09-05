from strucgap.preprocess import StrucGAP_Preprocess
from strucgap.glycanstructure import StrucGAP_GlycanStructure
from strucgap.glycosite import StrucGAP_GlycoSite
from strucgap.glycopeptidequant import StrucGAP_GlycoPeptideQuant
from strucgap.functionannotation import StrucGAP_FunctionAnnotation
from strucgap.glyconetwork import StrucGAP_GlycoNetwork
from strucgap.datavisualization import StrucGAP_DataVisualization
from strucgap.insighttracker import StrucGAP_InsightTracker
import os
data_manager = StrucGAP_InsightTracker()

os.chdir('D:\\doctor\\analysisys\\StrucGAP\\pglyco3')

# only structure information
module1 = StrucGAP_Preprocess(data_dir="D:\\doctor\\analysisys\\data\\fangpan_nc\\pd structure.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'D:\\doctor\\analysisys\\data\\sample_group.xlsx',
                      branch_list_dir = "D:\\doctor\\wyq\\branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager,
                      search_engine = 'pGlyco3')
module1.data_cleaning(data_type='tmt')
module1.cv_raw(threshold='no', fc_recommendation = False)
# 1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464
module1.fdr(feature_type='no')
module1.outliers(abundance_ratio=[1.240003449, 0, 1.344387558, 0, 1.576533442, 0, 1, 0, 1.956346409, 1.517000766])
module1.cv(threshold = 'no')
module1.psm(psm_number = 'no', fc_recommendation = False)
module1.annotation(glytoucan = True, biosynthetic_pathways = True, glycobiology_filter = True)
module1.output() 

# 
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


# structure information with quantification
module1 = StrucGAP_Preprocess(data_dir="D:\\doctor\\analysisys\\data\\fangpan_nc\\pd structure.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'D:\\doctor\\analysisys\\data\\sample_group.xlsx',
                      branch_list_dir = "D:\\doctor\\wyq\\branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager,
                      search_engine = 'pGlyco3')
module1.data_cleaning(data_type='tmt',quantification_from_no_strucgp = True,
                      quantification_data_dir = "D:\\doctor\\analysisys\\data\\fangpan_nc\\s6.xlsx" , 
                      sheet_name = 'Sheet1', quant_cols = ["Young-SN-1","Young-SN-2","Young-SN-3","PD-SN-1","PD-SN-2","PD-SN-3"])
module1.cv_raw(threshold='no', fc_recommendation = True)
# 1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464
module1.fdr(feature_type='no')
module1.outliers(abundance_ratio=[1.240003449, 0, 1.344387558, 0, 1.576533442, 0, 1, 0, 1.956346409, 1.517000766])
module1.cv(threshold = 'no')
module1.psm(psm_number = 'no', fc_recommendation = True)
module1.annotation(glytoucan = True, biosynthetic_pathways = True, glycobiology_filter = True)
module1.output() 

# 
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
module4 = StrucGAP_GlycoPeptideQuant(module1, data_type = 'psm_filtered', data_manager=data_manager)
module4.statistics()
module4.statistics_index()
module4.differential_analysis(pvalue_type='pvalue_ttest',
                                     fc = 4.2)
# module4.threshold_variation_analysis(pvalue_type='pvalue_ttest',statistic_index='fc',
#                                      fc_range = [3, 6, 9, 12, 15])
module4.threshold_variation_analysis(pvalue_type='pvalue_ttest',statistic_index='fc',
                                     fc_range = [3, 6, 10, 50, 100])
module4.glycopeptide_glycosite_glycan_variation(fc = 4.2)
module4.glycoprotein_glycosite_glycan_variation(fc = 4.2)
module4.output()

data_manager.key_information_extraction(module='StrucGAP_GlycoPeptideQuant')

#  
module5 = StrucGAP_FunctionAnnotation(module1, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')
# module5 = StrucGAP_FunctionAnnotation(module1, GO_data_dir = "D:\\doctor\\analysisys\\GAP\\go-basic.obo", 
#                                  data_manager=data_manager)  
# module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
# module5.kegg_function_structure(function_data = 'ora_no_background_both_proteins_result')  
# module5.output()

module5 = StrucGAP_FunctionAnnotation(module4, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', selected_terms=['GO:MF', 'GO:CC', 'GO:BP'], enrich_feature='glycopeptide',
            background_input=False, up_down_fc_threshold=4.2,pvalue_type='pvalue_ttest') # 65,72,79
# module5.gsea()   test = module5.ora_no_background_up_result
module5.go_function_structure(function_data = 'ora_no_background_up_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

module5.go_function_structure(function_data = 'ora_no_background_down_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

# module5.go_function_structure(function_data = 'gsea_result')  
# module5.output()

# module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 119
# module5.gsea()
# module5.kegg_function_structure(function_data = 'ora_no_background_up_result')
# module5.output()
# module5.kegg_function_structure(function_data = 'ora_no_background_down_result')
# module5.output()
# module5.kegg_function_structure(function_data = 'gsea_result')
# module5.output()

# module5 = StrucGAP_FunctionAnnotation(module6, 
#                                  data_manager=data_manager)  
# module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 69,76,83
# # module5.gsea()
# module5.go_function_structure(function_data = 'ora_no_background_up_result')  
# module5.output()
# module5.go_function_structure(function_data = 'ora_no_background_down_result')  
# module5.output()
# # module5.go_function_structure(function_data = 'gsea_result')  
# # module5.output()

# #
# module6 = StrucGAP_GlycoNetwork(module4, data_manager=data_manager)
# module6.proteomic(protein_data_dir="D:\\doctor\\wyq\\WYQ_Mus_uterus_global.xlsx",
#                   data_sheet_name = '1 Proteins',
#                   cv = 'no', psm = 1, 
#                   # normalization_samplewise_method = 'no',
#                   # normalization_featurewise_method = 'no',
#                   )
# module6.phosphorylation(phospho_data_dir="D:\doctor\zzd\_Rat_Phospho_mixThymus_TMT6c_.xlsx",
#                         data_sheet_name='PeptideGroups')
# module6.glycosyltransferases(glycosyltransferases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", 
#                              data_sheet_name="glycosyltransferases")
# module6.glycosidases(glycosidases_data_dir="D:\\doctor\\analysisys\\GAP\\enzyme.xlsx", 
#                      data_sheet_name='glycosidases')
# module6.sialyltransferases()
# module6.fucosyltransferase()
# module6.glycan_binding_protein()
# module6.output()

# data_manager.key_information_extraction(module='StrucGAP_GlycoNetwork')
#
data_manager.output_pickle()

module7 = StrucGAP_DataVisualization(data_manager=data_manager)




