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

os.chdir('D:\\doctor\\analysisys\\StrucGAP\\msfragger')

# only structure information
module1 = StrucGAP_Preprocess(data_dir="D:\\doctor\\analysisys\\data\\others output\\msfragger.xlsx",
                      data_sheet_name = 'Sheet1',
                      sample_group_data_dir = 'D:\\doctor\\analysisys\\data\\sample_group.xlsx',
                      branch_list_dir = "D:\\doctor\\wyq\\branch_structures_18_mice uterus.0240401.xlsx",
                      data_manager=data_manager,
                      search_engine = 'MSFragger-Glyco')
module1.data_cleaning(data_type='tmt')
module1.cv_raw(threshold='no', fc_recommendation = False)
# 1.172277596,1.142983373,1,1.46390136,1.466662624,1.449428354,1.109519196,1.387464059,1.291746761,1.487440464
module1.fdr(feature_type='no')
module1.outliers(abundance_ratio=[1.240003449, 0, 1.344387558, 0, 1.576533442, 0, 1, 0, 1.956346409, 1.517000766])
module1.cv(threshold = 'no')
module1.psm(psm_number = 'no', fc_recommendation = False)
module1.annotation(glytoucan = True, biosynthetic_pathways = True, glycobiology_filter = False)
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
module5 = StrucGAP_FunctionAnnotation(module1, 
                                 data_manager=data_manager)  
module5.ora(organism='mmusculus', background_input=False, up_down_fc_threshold=1.5) # 65,72,79
module5.go_function_structure(function_data = 'ora_no_background_both_proteins_result')  
module5.output()

data_manager.key_information_extraction(module='StrucGAP_FunctionAnnotation')

#
data_manager.output_pickle()

module7 = StrucGAP_DataVisualization(data_manager=data_manager)




