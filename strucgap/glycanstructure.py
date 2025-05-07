# -*- coding: utf-8 -*-
## 糖链结构分析模块--68
class StrucGAP_GlycanStructure:
    def __init__(self, gs_data, data_manager, data_type = 'psm_filtered',
                 pvalue = 0.05, fc = 1.5, pvalue_type='pvalue_ttest', differential_data_type='both'):
        self.gs_data = gs_data
        self.sample_group = self.gs_data.sample_group
        #
        if gs_data == module1:
            if data_type == 'psm_filtered':
                self.data = self.gs_data.data_psm_filtered
            elif data_type == 'cv_filtered':
                self.data = self.gs_data.data_cv_filtered
            elif data_type == 'outliers_filtered':
                self.data = self.gs_data.data_outliers_filtered
            elif data_type == 'data':
                self.data = self.gs_data.data
            else:
                raise ValueError("Invalid data_type provided.")  
        else:
            if gs_data == module4:
                self.data = pd.concat([self.gs_data.data, self.gs_data.fc_result],axis=1)
                self.data = self.data[self.data[pvalue_type]<pvalue]
                if differential_data_type=='both':
                    self.data = self.data[(self.data['fc']>fc)|(self.data['fc']<(1/fc))]
                elif differential_data_type=='up':
                    self.data = self.data[(self.data['fc']>fc)]
                elif differential_data_type=='down':
                    self.data = self.data[(self.data['fc']<(1/fc))]
        
        self.data_manager = data_manager    
        self.data_manager.register_module('StrucGAP_GlycanStructure', self, {})
        self.data_manager.log_params('StrucGAP_GlycanStructure', 'input_data', {'data_type': data_type})
    
    def statistics(self, remove_oligo_mannose = True): # , rank = None
        """
        Provides a global overview of glycan structures in the dataset.
        
        Parameters:
            remove_oligo_mannose: if exclude high-abundance oligo-mannose glycans to focus on structurally diverse and functionally relevant isomers (True of False).
        
        Returns:
            self.GlycanComposition_rank (The ranking of glycan compositions based on their quantities). 
            self.structure_coding_rank (The ranking of glycan structures based on their quantities).
            
        Return type:
            dataframe
        
        """
        if remove_oligo_mannose == True:
            self.data = self.data[self.data['Glycan_type'] != 'Oligo mannose']
        elif remove_oligo_mannose == False:
            self.data = self.data
        elif remove_oligo_mannose not in [True, False]:
            self.data = self.data[self.data['Glycan_type'] != 'Oligo mannose']
            print("Your input is invalid, the parameter 'remove_oligo_mannose' is set to True! ")
   
        self.GlycanComposition_rank = self.identify_glycan_composition() # .iloc[:self.rank]  
        #
        self.structure_coding_rank = pd.DataFrame(self.data['structure_coding'].value_counts().reset_index()) # .iloc[:rank]  
        self.structure_coding_rank.columns = ['Structure_coding', 'Structure_coding_count']
        
        self.data_manager.log_params('StrucGAP_GlycanStructure', 'statistics', {'remove_oligo_mannose': remove_oligo_mannose})
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'GlycanComposition_rank', self.GlycanComposition_rank)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'structure_coding_rank', self.structure_coding_rank)
        
        return self
    
    def identify_glycan_composition(self):  
        """An auxiliary function called by other functions to identify glycan composition."""
        temp_data = self.data[['GlycanComposition', 'structure_coding']] 
        temp_data = temp_data.replace(np.nan, '')
        glycancomposition = {}
        for i in temp_data['GlycanComposition'].unique(): 
            glycancomposition[i] = temp_data[temp_data['GlycanComposition'] == i]['structure_coding'].value_counts().max()
        GlycanComposition_rank = pd.DataFrame(list(glycancomposition.items()), columns=['GlycanComposition', 'GlycanComposition_count'])
        GlycanComposition_rank = GlycanComposition_rank.sort_values(by='GlycanComposition_count', ascending=False) 
        GlycanComposition_rank = GlycanComposition_rank[GlycanComposition_rank['GlycanComposition']!='']
        return GlycanComposition_rank
    
    def identify_core_structure(self):  
        """An auxiliary function called by other functions to identify core structure."""
        temp_data = self.data[(self.data['GlycanComposition']!='N2H2')&(self.data['GlycanComposition']!='N2H2F1')&(~self.data['structure_coding'].str.contains('A2B2C1D2d'))]
        temp_data = temp_data[['structure_coding','Corefucose', 'Bisection']] 
        core_structure_list = []
        for i, row in temp_data.iterrows():
            if row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD1'
            elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' not in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD1dcbB5'
            elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' not in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD2dD1'
            elif row['structure_coding'][0:8] == 'A2B2C1D1' and 'D1' in row['structure_coding'][9:] and 'D2' in row['structure_coding'][9:] and 'B5' in row['structure_coding'][9:]:
                structure = 'A2B2C1D1dD2dD1dcbB5'
            else:
                continue
            core_structure_list.append([i, structure])
        core_structure = pd.DataFrame(core_structure_list, columns=['PeptideSequence+structure_coding+ProteinID', 'Core_structure']).set_index('PeptideSequence+structure_coding+ProteinID', drop=True)
        return core_structure
    
    def structure_statistics(self):
        """
        Performs independent quantification of all annotated glycan substructural attributes in the dataset.
        
        Parameters:
            None.
        
        Returns:
            self.core_structure (The number of IGPs containing different core structures).
            self.glycan_type (The number of IGPs containing different glycan types).
            self.branches_structure (The number of IGPs containing different branch structures).
            self.branches_count (The number of IGPs containing different branch counts).
            self.sialicacid_count (The number of IGPs containing different sialic acid counts).
            self.fucose_count (The number of IGPs containing different fucose counts).
            self.sialicacid_structure (The number of IGPs containing different sialic acid structures).
            self.fucose_structure (The number of IGPs containing different fucose structures).
            self.fucosylated_type (The number of IGPs containing different fucosylated types).
            self.acgc (The number of IGPs containing different sialylated types).
            self.fsg (The number of IGPs containing different fucosylated/sialylated patterns).
            
        Return type:
            dataframe
        
        """
        self.core_structure = self.identify_core_structure()['Core_structure'].value_counts().reset_index()
        self.core_structure.columns = ['Core_structure', 'Core_structure_count']
        #
        self.glycan_type = pd.DataFrame(self.data['Glycan_type'].value_counts().reset_index())
        self.glycan_type.columns = ['Glycan_type', 'Glycan_type_count']
        #
        self.branches_structure = pd.DataFrame(self.data['Branches'].value_counts().reset_index())
        self.branches_structure.columns = ['Branches', 'Branches_count']
        count_dict = {}
        for index, row in self.branches_structure.iterrows():
            branch_list = literal_eval(row['Branches'])
            count = row['Branches_count']
            for branch in branch_list:
                if branch:
                    if branch in count_dict:
                        count_dict[branch] += count
                    else:
                        count_dict[branch] = count
        self.branches_structure = pd.DataFrame(list(count_dict.items()), columns=['Branches', 'Branches_count'])
        #
        self.branches_count = pd.DataFrame(self.data['BranchNumber'].value_counts().reset_index())
        self.branches_count.columns = ['BranchNumber', 'BranchNumber_count']
        #
        self.sialicacid_count = pd.DataFrame(self.data[~self.data['structure_coding'].isnull()]['structure_coding'].str.count('3').value_counts().reset_index())
        self.sialicacid_count.columns = ['Sialicacid_count', 'Number']
        #
        self.fucose_count = pd.DataFrame(self.data[~self.data['structure_coding'].isnull()]['structure_coding'].str.count('5').value_counts().reset_index())
        self.fucose_count.columns = ['Fucose_count', 'Number']
        #
        self.sialicacid_structure = pd.DataFrame(self.data[~self.data['structure_coding'].isnull() & self.data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().reset_index())
        self.sialicacid_structure.columns = ['sialicacid_structure', 'Number']
        #
        self.fucose_structure = pd.DataFrame(self.data[~self.data['structure_coding'].isnull() & self.data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().reset_index())
        self.fucose_structure.columns = ['fucose_structure', 'Number']
        #
        self.fucosylated_type = pd.DataFrame(self.data['fucosylated type'].value_counts().reset_index())
        self.fucosylated_type = self.fucosylated_type[self.fucosylated_type['index']!=' ']
        self.fucosylated_type.columns = ['FucosylatedType', 'FucosylatedType_count']
        #
        self.acgc = pd.DataFrame(self.data['Ac/Gc'].value_counts().reset_index())
        self.acgc = self.acgc[self.acgc['index']!=' ']
        self.acgc.columns = ['Ac/Gc', 'Ac/Gc_count']
        #
        self.fsg = pd.DataFrame(self.data['FSG'].value_counts().reset_index())
        self.fsg.columns = ['FSG', 'FSG_count']
        #
        self.data_manager.log_params('StrucGAP_GlycanStructure', 'structure_statistics', {})
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure', self.core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type', self.glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure', self.branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_count', self.branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'sialicacid_count', self.sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucose_count', self.fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'sialicacid_structure', self.sialicacid_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucose_structure', self.fucose_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type', self.fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc', self.acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fsg', self.fsg)
        
        return self
    
    def lacdinac(self):
        """An auxiliary function called by other functions to identify LacdiNAc."""
        # def extract_key_strings(s):
        #     matches = re.findall(r'(E2F2.*?fe)', s)
        #     return ', '.join(matches)
        
        # # 应用提取函数到新列
        # self.data['lacdinac'] = self.data['structure_coding'].apply(extract_key_strings)
        self.lacdinac_count = self.data[['PeptideSequence+structure_coding+ProteinID', 'lacdinac']]
        self.lacdinac_count.set_index('PeptideSequence+structure_coding+ProteinID', inplace=True, drop=True)
        self.lacdinac_count = self.lacdinac_count['lacdinac'].str.split(',\s*', expand=True).stack().reset_index()
        self.lacdinac_count = self.lacdinac_count.drop_duplicates(subset=['PeptideSequence+structure_coding+ProteinID', 0])
        self.lacdinac_count = self.lacdinac_count[0]
        self.lacdinac_count = self.lacdinac_count.value_counts().reset_index()
        self.lacdinac_count = self.lacdinac_count[self.lacdinac_count['index']!=' ']
        self.lacdinac_count = self.lacdinac_count.rename(columns={'index': 'lacdinac', 0:'lacdinac_count'})
        
        self.data_manager.log_params('StrucGAP_GlycanStructure', 'lacdinac', {})
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_count', self.lacdinac_count)

        return self
    
    def cor(self):
        """
        Computes pairwise co-occurrence frequencies of glycan substructural features across all IGPs.
        
        Parameters:
            None.
        
        Returns:
            The co-occurrence rate between pairwise glycan substructure features.
            self.struc_cod_lacdinac
            self.struc_cod_core_structure
            self.struc_cod_glycan_type
            self.struc_cod_branches_structure
            self.struc_cod_branches_count
            self.struc_cod_sialicacid_count
            self.struc_cod_fucose_count
            self.struc_cod_fucosylated_type
            self.struc_cod_acgc
            self.lacdinac_core_structure
            self.lacdinac_glycan_type
            self.lacdinac_branches_structure
            self.lacdinac_branches_count
            self.lacdinac_sialicacid_count
            self.lacdinac_fucose_count
            self.lacdinac_fucosylated_type
            self.lacdinac_acgc
            self.core_structure_lacdinac
            self.core_structure_glycan_type
            self.core_structure_branches_structure
            self.core_structure_branches_count
            self.core_structure_sialicacid_count
            self.core_structure_fucose_count
            self.core_structure_fucosylated_type
            self.core_structure_acgc
            self.glycan_type_lacdinac
            self.glycan_type_core_structure
            self.glycan_type_branches_structure
            self.glycan_type_branches_count
            self.glycan_type_sialicacid_count
            self.glycan_type_fucose_count
            self.glycan_type_fucosylated_type
            self.glycan_type_acgc
            self.branches_structure_core_structure
            self.branches_structure_glycan_type
            self.branches_structure_fucosylated_type
            self.branches_structure_acgc
            self.branches_structure_sialicacid_count 
            self.branches_structure_fucose_count
            self.branches_count_sialicacid_count
            self.branches_count_fucose_count
            self.fucosylated_type_lacdinac
            self.fucosylated_type_core_structure
            self.fucosylated_type_glycan_type
            self.fucosylated_type_branches_structure
            self.fucosylated_type_branches_count
            self.fucosylated_type_sialicacid_count
            self.fucosylated_type_fucose_count
            self.fucosylated_type_acgc
            self.acgc_lacdinac
            self.acgc_core_structure
            self.acgc_glycan_type
            self.acgc_branches_structure
            self.acgc_branches_count
            self.acgc_sialicacid_count
            self.acgc_fucose_count
            self.acgc_fucosylated_type
            
        Return type:
            dataframe
        
        """
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['structure_coding'])[~pd.DataFrame(self.data['structure_coding']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['lacdinac']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Lacdinac', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Lacdinac', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.struc_cod_lacdinac = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['structure_coding']),self.identify_core_structure()], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['Core_structure']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Core_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Core_structure', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.struc_cod_core_structure = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['structure_coding']), 
                             pd.DataFrame(self.data['Glycan_type'])[~pd.DataFrame(self.data['Glycan_type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['Glycan_type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Glycan_type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Glycan_type', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['structure_coding']!=' ']
        self.struc_cod_glycan_type = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),pd.DataFrame(self.data['structure_coding'])], axis=1)
        for branch in type_counts.index:
            count_dict = {}
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['Branches']
            item = filtered_data.value_counts().reset_index()
            for index, row in item.iterrows():
                branch_list = literal_eval(row['index'])
                count = row['Branches']
                for i in branch_list:
                    if i:
                        if i in count_dict:
                            count_dict[i] += count
                        else:
                            count_dict[i] = count
            item = pd.Series(count_dict)
            result[branch] = item   
        df = pd.DataFrame(columns=['structure_coding', 'Branches_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branches_structure', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.struc_cod_branches_structure = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['BranchNumber']),pd.DataFrame(self.data['structure_coding'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['BranchNumber']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Branch_number', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branch_number', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.struc_cod_branches_count = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}),pd.DataFrame(self.data['structure_coding'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.struc_cod_sialicacid_count = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')}),pd.DataFrame(self.data['structure_coding'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.struc_cod_fucose_count = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['fucosylated type']),pd.DataFrame(self.data['structure_coding'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['fucosylated type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'fucosylated type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['fucosylated type', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['fucosylated type']!=' ']
        self.struc_cod_fucosylated_type = df
        #
        type_counts = pd.DataFrame(self.data['structure_coding'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Ac/Gc']),pd.DataFrame(self.data['structure_coding'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['structure_coding'] == branch]['Ac/Gc']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['structure_coding', 'Ac/Gc', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Ac/Gc', 'Count']
            temp_df['structure_coding'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.struc_cod_acgc = df
        
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.identify_core_structure())[~pd.DataFrame(self.identify_core_structure()).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['Core_structure']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Core_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Core_structure', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_core_structure = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['Glycan_type'])[~pd.DataFrame(self.data['Glycan_type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['Glycan_type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Glycan_type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Glycan_type', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_glycan_type = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['Branches'])[~pd.DataFrame(self.data['Branches']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            count_dict = {}
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['Branches']
            item = filtered_data.value_counts().reset_index()
            for index, row in item.iterrows():
                branch_list = literal_eval(row['index'])
                count = row['Branches']
                for i in branch_list:
                    if i:
                        if i in count_dict:
                            count_dict[i] += count
                        else:
                            count_dict[i] = count
            item = pd.Series(count_dict)
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Branches', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branches', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_branches_structure = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['BranchNumber'])[~pd.DataFrame(self.data['BranchNumber']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['BranchNumber']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Branch_Number', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branch_Number', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_branches_count = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')})[~pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_sialicacid_count = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')})[~pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('5')}).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_fucose_count = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['fucosylated type'])[~pd.DataFrame(self.data['fucosylated type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['fucosylated type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'fucosylated type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['fucosylated type', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['fucosylated type']!=' ']
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_fucosylated_type = df
        #
        type_counts = self.lacdinac_count['lacdinac']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['Ac/Gc'])[~pd.DataFrame(self.data['Ac/Gc']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['lacdinac'] == branch]['Ac/Gc']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Lacdinac', 'Ac/Gc', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Ac/Gc', 'Count']
            temp_df['Lacdinac'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        df = df[df['Lacdinac']!=' ']
        self.lacdinac_acgc = df
        
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                     pd.DataFrame(self.identify_core_structure())[~pd.DataFrame(self.identify_core_structure()).index.duplicated()], 
                     left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['lacdinac']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Lacdinac', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Lacdinac', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.core_structure_lacdinac = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Glycan_type']),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['Glycan_type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Glycan_type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Glycan_type', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.core_structure_glycan_type = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            count_dict = {}
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['Branches']
            item = filtered_data.value_counts().reset_index()
            for index, row in item.iterrows():
                branch_list = literal_eval(row['index'])
                count = row['Branches']
                for i in branch_list:
                    if i:
                        if i in count_dict:
                            count_dict[i] += count
                        else:
                            count_dict[i] = count
            item = pd.Series(count_dict)
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Branches_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branches_structure', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.core_structure_branches_structure = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['BranchNumber']),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['BranchNumber']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Branch_number', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branch_number', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.core_structure_branches_count = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.core_structure_sialicacid_count = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')}),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.core_structure_fucose_count = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['fucosylated type']),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['fucosylated type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'fucosylated type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['fucosylated type', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['fucosylated type']!=' ']
        self.core_structure_fucosylated_type = df
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Ac/Gc']),self.identify_core_structure()], axis=1)
        for branch in type_counts:
            filtered_data = temp_data[temp_data['Core_structure'] == branch]['Ac/Gc']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Core_structure', 'Ac/Gc', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Ac/Gc', 'Count']
            temp_df['Core_structure'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.core_structure_acgc = df
        
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.merge(self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             pd.DataFrame(self.data['Glycan_type'])[~pd.DataFrame(self.data['Glycan_type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['lacdinac']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'Lacdinac', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Lacdinac', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.glycan_type_lacdinac = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Glycan_type']),self.identify_core_structure()], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['Core_structure']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'Core_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Core_structure', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.glycan_type_core_structure = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        for branch in type_counts.index:
            count_dict = {}
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['Branches']
            item = filtered_data.value_counts().reset_index()
            for index, row in item.iterrows():
                branch_list = literal_eval(row['index'])
                count = row['Branches']
                for i in branch_list:
                    if i:
                        if i in count_dict:
                            count_dict[i] += count
                        else:
                            count_dict[i] = count
            item = pd.Series(count_dict)
            result[branch] = item   
        df = pd.DataFrame(columns=['Glycan_type', 'Branches_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branches_structure', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.glycan_type_branches_structure = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['BranchNumber']),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['BranchNumber']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'Branch_number', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branch_number', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.glycan_type_branches_count = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.glycan_type_sialicacid_count = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')}),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.glycan_type_fucose_count = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['fucosylated type']),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['fucosylated type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'fucosylated type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['fucosylated type', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['fucosylated type']!=' ']
        self.glycan_type_fucosylated_type = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Ac/Gc']),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Glycan_type'] == branch]['Ac/Gc']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Glycan_type', 'Ac/Gc', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Ac/Gc', 'Count']
            temp_df['Glycan_type'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.glycan_type_acgc = df
        
        #
        type_counts = self.core_structure['Core_structure']
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),self.identify_core_structure()], axis=1)
        co_occurrence_count = {}
        unique_branches = set()
        unique_core_structures = set()
        data = temp_data[['Branches', 'Core_structure']]          
        unique_branches.update(data['Branches'].apply(eval).explode().unique())
        unique_core_structures.update(data['Core_structure'].unique())
        for _, row in data.iterrows():
            core = row['Core_structure']
            branches = eval(row['Branches'])  
            for branch in branches:
                pair = (core, branch)
                # print(pair)
                if pair not in co_occurrence_count:
                    co_occurrence_count[pair] = 0
                co_occurrence_count[pair] += 1
        df = pd.DataFrame([
            {
                'Branches': key[1],  
                'Core_structure': key[0],  
                'Count': value   
            }
            for key, value in co_occurrence_count.items()
        ])
        df = df.sort_values(by=['Branches', 'Count'], ascending=[False, False])
        self.branches_structure_core_structure = df
        #
        type_counts = pd.DataFrame(self.data['Glycan_type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),pd.DataFrame(self.data['Glycan_type'])], axis=1)
        co_occurrence_count = {}
        unique_branches = set()
        unique_core_structures = set()
        data = temp_data[['Branches', 'Glycan_type']]          
        unique_branches.update(data['Branches'].apply(eval).explode().unique())
        unique_core_structures.update(data['Glycan_type'].unique())
        for _, row in data.iterrows():
            core = row['Glycan_type']
            branches = eval(row['Branches'])  
            for branch in branches:
                pair = (core, branch)
                # print(pair)
                if pair not in co_occurrence_count:
                    co_occurrence_count[pair] = 0
                co_occurrence_count[pair] += 1
        df = pd.DataFrame([
            {
                'Branches': key[1],  
                'Glycan_type': key[0],  
                'Count': value   
            }
            for key, value in co_occurrence_count.items()
        ])
        df = df.sort_values(by=['Branches', 'Count'], ascending=[False, False])
        self.branches_structure_glycan_type = df
        #
        type_counts = pd.DataFrame(self.data['fucosylated type'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),pd.DataFrame(self.data['fucosylated type'])], axis=1)
        co_occurrence_count = {}
        unique_branches = set()
        unique_core_structures = set()
        data = temp_data[['Branches', 'fucosylated type']]          
        unique_branches.update(data['Branches'].apply(eval).explode().unique())
        unique_core_structures.update(data['fucosylated type'].unique())
        for _, row in data.iterrows():
            core = row['fucosylated type']
            branches = eval(row['Branches'])  
            for branch in branches:
                pair = (core, branch)
                # print(pair)
                if pair not in co_occurrence_count:
                    co_occurrence_count[pair] = 0
                co_occurrence_count[pair] += 1
        df = pd.DataFrame([
            {
                'Branches': key[1],  
                'fucosylated type': key[0],  
                'Count': value   
            }
            for key, value in co_occurrence_count.items()
        ])
        df = df.sort_values(by=['Branches', 'Count'], ascending=[False, False])
        df = df[df['fucosylated type']!=' ']
        self.branches_structure_fucosylated_type = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['Branches']),pd.DataFrame(self.data['Ac/Gc'])], axis=1)
        co_occurrence_count = {}
        unique_branches = set()
        unique_core_structures = set()
        data = temp_data[['Branches', 'Ac/Gc']]          
        unique_branches.update(data['Branches'].apply(eval).explode().unique())
        unique_core_structures.update(data['Ac/Gc'].unique())
        for _, row in data.iterrows():
            core = row['Ac/Gc']
            branches = eval(row['Branches'])  
            for branch in branches:
                pair = (core, branch)
                # print(pair)
                if pair not in co_occurrence_count:
                    co_occurrence_count[pair] = 0
                co_occurrence_count[pair] += 1
        df = pd.DataFrame([
            {
                'Branches': key[1],  
                'Ac/Gc': key[0],  
                'Count': value   
            }
            for key, value in co_occurrence_count.items()
        ])
        df = df.sort_values(by=['Branches', 'Count'], ascending=[False, False])
        df = df[df['Ac/Gc']!=' ']
        self.branches_structure_acgc = df
        #
        type_counts = pd.DataFrame(self.data['Branches'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}),pd.DataFrame(self.data['Branches'])], axis=1)
        unique_branches = set()
        for branches_str in temp_data['Branches']:
            branches_list = literal_eval(branches_str)
            for branch in branches_list:
                if branch:  
                    unique_branches.add(branch)
        unique_branches_list = list(unique_branches)
        for branch in unique_branches_list:
            filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Branches', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['Branches'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.branches_structure_sialicacid_count = df
        #
        type_counts = pd.DataFrame(self.data['Branches'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')}),pd.DataFrame(self.data['Branches'])], axis=1)
        unique_branches = set()
        for branches_str in temp_data['Branches']:
            branches_list = literal_eval(branches_str)
            for branch in branches_list:
                if branch:  
                    unique_branches.add(branch)
        unique_branches_list = list(unique_branches)
        for branch in unique_branches_list:
            filtered_data = temp_data[temp_data['Branches'].apply(lambda x: branch in literal_eval(x))]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Branches', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['Branches'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.branches_structure_fucose_count = df
        
        #
        type_counts = pd.DataFrame(self.data['BranchNumber'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}),pd.DataFrame(self.data['BranchNumber'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['BranchNumber'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['BranchNumber', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['BranchNumber'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.branches_count_sialicacid_count = df
        #
        type_counts = pd.DataFrame(self.data['BranchNumber'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')}),pd.DataFrame(self.data['BranchNumber'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['BranchNumber'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['BranchNumber', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['BranchNumber'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.branches_count_fucose_count = df
        
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['lacdinac']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Lacdinac', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Lacdinac', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Lacdinac']!=' ']
        self.fucosylated_type_lacdinac = df
        self.fucosylated_type_lacdinac = self.fucosylated_type_lacdinac[self.fucosylated_type_lacdinac['Lacdinac']!='']
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame(self.identify_core_structure())[~pd.DataFrame(self.identify_core_structure()).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['Core_structure']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Core_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Core_structure', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.fucosylated_type_core_structure = df
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame(self.data['Glycan_type'])[~pd.DataFrame(self.data['Glycan_type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['Glycan_type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Glycan_type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Glycan_type', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.fucosylated_type_glycan_type = df
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame(self.data['Branches'])[~pd.DataFrame(self.data['Branches']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            count_dict = {}
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['Branches']
            item = filtered_data.value_counts().reset_index()
            for index, row in item.iterrows():
                branch_list = literal_eval(row['index'])
                count = row['Branches']
                for i in branch_list:
                    if i:
                        if i in count_dict:
                            count_dict[i] += count
                        else:
                            count_dict[i] = count
            item = pd.Series(count_dict)
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Branches', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branches', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.fucosylated_type_branches_structure = df
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame(self.data['BranchNumber'])[~pd.DataFrame(self.data['BranchNumber']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['BranchNumber']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Branch_Number', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branch_Number', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.fucosylated_type_branches_count = df
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')})[~pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.fucosylated_type_sialicacid_count = df
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')})[~pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('5')}).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.fucosylated_type_fucose_count = df
        #
        type_counts = self.fucosylated_type['FucosylatedType']
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['fucosylated type']), 
                             pd.DataFrame(self.data['Ac/Gc'])[~pd.DataFrame(self.data['Ac/Gc']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts:
            filtered_data = temp_data[temp_data['fucosylated type'] == branch]['Ac/Gc']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['FucosylatedType', 'Ac/Gc', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Ac/Gc', 'Count']
            temp_df['FucosylatedType'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.fucosylated_type_acgc = df
        
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             self.data['lacdinac'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).to_frame('lacdinac'), 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['lacdinac']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Lacdinac', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Lacdinac', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        df = df[df['Lacdinac']!=' ']
        self.acgc_lacdinac = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame(self.identify_core_structure())[~pd.DataFrame(self.identify_core_structure()).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['Core_structure']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Core_structure', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Core_structure', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.acgc_core_structure = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame(self.data['Glycan_type'])[~pd.DataFrame(self.data['Glycan_type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['Glycan_type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Glycan_type', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Glycan_type', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.acgc_glycan_type = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame(self.data['Branches'])[~pd.DataFrame(self.data['Branches']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            count_dict = {}
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['Branches']
            item = filtered_data.value_counts().reset_index()
            for index, row in item.iterrows():
                branch_list = literal_eval(row['index'])
                count = row['Branches']
                for i in branch_list:
                    if i:
                        if i in count_dict:
                            count_dict[i] += count
                        else:
                            count_dict[i] = count
            item = pd.Series(count_dict)
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Branches', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branches', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.acgc_branches_structure = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame(self.data['BranchNumber'])[~pd.DataFrame(self.data['BranchNumber']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['BranchNumber']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Branch_Number', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Branch_Number', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.acgc_branches_count = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')})[~pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('3')}).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['sialicacid_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Sialicacid_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Sialicacid_count', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.acgc_sialicacid_count = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame({'fucose_count': self.data['structure_coding'].str.count('5')})[~pd.DataFrame({'sialicacid_count': self.data['structure_coding'].str.count('5')}).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['fucose_count']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'Fucose_count', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['Fucose_count', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        self.acgc_fucose_count = df
        #
        type_counts = pd.DataFrame(self.data['Ac/Gc'].value_counts())
        result = {}
        temp_data = pd.merge(pd.DataFrame(self.data['Ac/Gc']), 
                             pd.DataFrame(self.data['fucosylated type'])[~pd.DataFrame(self.data['fucosylated type']).index.duplicated()], 
                             left_index=True, right_index=True, how='left')
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['Ac/Gc'] == branch]['fucosylated type']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['Ac/Gc', 'FucosylatedType', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['FucosylatedType', 'Count']
            temp_df['Ac/Gc'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        df = df[df['Ac/Gc']!=' ']
        df = df[df['FucosylatedType']!=' ']
        self.acgc_fucosylated_type = df
        
        self.data_manager.log_params('StrucGAP_GlycanStructure', 'cor', {})
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_lacdinac', self.struc_cod_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_core_structure', self.struc_cod_core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_glycan_type', self.struc_cod_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_branches_structure', self.struc_cod_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_branches_count', self.struc_cod_branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_sialicacid_count', self.struc_cod_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_fucose_count', self.struc_cod_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_fucosylated_type', self.struc_cod_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'struc_cod_acgc', self.struc_cod_acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_core_structure', self.lacdinac_core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_glycan_type', self.lacdinac_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_branches_structure', self.lacdinac_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_branches_count', self.lacdinac_branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_sialicacid_count', self.lacdinac_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_fucose_count', self.lacdinac_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_fucosylated_type', self.lacdinac_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'lacdinac_acgc', self.lacdinac_acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_lacdinac', self.core_structure_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_glycan_type', self.core_structure_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_branches_structure', self.core_structure_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_branches_count', self.core_structure_branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_sialicacid_count', self.core_structure_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_fucose_count', self.core_structure_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_fucosylated_type', self.core_structure_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'core_structure_acgc', self.core_structure_acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_lacdinac', self.glycan_type_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_core_structure', self.glycan_type_core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_branches_structure', self.glycan_type_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_branches_count', self.glycan_type_branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_sialicacid_count', self.glycan_type_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_fucose_count', self.glycan_type_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_fucosylated_type', self.glycan_type_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_type_acgc', self.glycan_type_acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure_core_structure', self.branches_structure_core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure_glycan_type', self.branches_structure_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure_fucosylated_type', self.branches_structure_fucosylated_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure_acgc', self.branches_structure_acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure_sialicacid_count', self.branches_structure_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_structure_fucose_count', self.branches_structure_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_count_sialicacid_count', self.branches_count_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'branches_count_fucose_count', self.branches_count_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_lacdinac', self.fucosylated_type_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_core_structure', self.fucosylated_type_core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_glycan_type', self.fucosylated_type_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_branches_structure', self.fucosylated_type_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_branches_count', self.fucosylated_type_branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_sialicacid_count', self.fucosylated_type_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_fucose_count', self.fucosylated_type_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'fucosylated_type_acgc', self.fucosylated_type_acgc)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_lacdinac', self.acgc_lacdinac)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_core_structure', self.acgc_core_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_glycan_type', self.acgc_glycan_type)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_branches_structure', self.acgc_branches_structure)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_branches_count', self.acgc_branches_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_sialicacid_count', self.acgc_sialicacid_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_fucose_count', self.acgc_fucose_count)
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'acgc_fucosylated_type', self.acgc_fucosylated_type)        
        return self
        
    def isoforms(self):
        """
        Investigates isomeric variants within identical glycan compositions.
        
        Parameters:
            None.
        
        Returns:
            self.glycan_composition_isoforms. 
            
        Return type:
            dataframe
        
        """
        #
        type_counts = pd.DataFrame(self.data['GlycanComposition'].value_counts())
        result = {}
        temp_data = pd.concat([pd.DataFrame(self.data['structure_coding']),pd.DataFrame(self.data['GlycanComposition'])], axis=1)
        for branch in type_counts.index:
            filtered_data = temp_data[temp_data['GlycanComposition'] == branch]['structure_coding']
            item = filtered_data.value_counts()
            result[branch] = item
        df = pd.DataFrame(columns=['GlycanComposition', 'structure_coding', 'Count'])
        for core_structure, series in result.items():
            temp_df = pd.DataFrame(series).reset_index()
            temp_df.columns = ['structure_coding', 'Count']
            temp_df['GlycanComposition'] = core_structure
            df = pd.concat([df, temp_df], ignore_index=True)
        self.glycan_composition_isoforms = df
        #
        self.data_manager.log_params('StrucGAP_GlycanStructure', 'isoforms', {})
        self.data_manager.log_output('StrucGAP_GlycanStructure', 'glycan_composition_isoforms', self.glycan_composition_isoforms)
        
        return self
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        
        with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_GlycanStructure.xlsx'), engine='xlsxwriter') as writer:
            export_data = self.GlycanComposition_rank.copy()
            export_data['GlycanComposition_ratio'] = export_data['GlycanComposition_count'] / export_data['GlycanComposition_count'].sum()
            export_data.to_excel(writer, sheet_name='GlycanComposition_rank')
            
            export_data = self.structure_coding_rank.copy()
            export_data['Structure_coding_ratio'] = export_data['Structure_coding_count'] / export_data['Structure_coding_count'].sum()
            export_data.to_excel(writer, sheet_name='structure_coding_rank')
            
            export_data = self.glycan_composition_isoforms.copy()
            export_data['Ratio'] = export_data.groupby('GlycanComposition')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_composition_isoforms'[:31])
            
            export_data = self.core_structure.copy()
            export_data['Core_structure_ratio'] = export_data['Core_structure_count'] / export_data['Core_structure_count'].sum()
            export_data.to_excel(writer, sheet_name='core_structure')
            export_data = self.glycan_type.copy()
            export_data['Glycan_type_ratio'] = export_data['Glycan_type_count'] / export_data['Glycan_type_count'].sum()
            export_data.to_excel(writer, sheet_name='glycan_type')
            export_data = self.branches_structure.copy()
            export_data['Branches_ratio'] = export_data['Branches_count'] / export_data['Branches_count'].sum()
            export_data.to_excel(writer, sheet_name='branches_structure')
            export_data = self.branches_count.copy()
            export_data['BranchNumber_ratio'] = export_data['BranchNumber_count'] / export_data['BranchNumber_count'].sum()
            export_data.to_excel(writer, sheet_name='branches_count')
            export_data = self.sialicacid_count.copy()
            export_data['Ratio'] = export_data['Number'] / export_data['Number'].sum()
            export_data.to_excel(writer, sheet_name='sialicacid_count')
            export_data = self.fucose_count.copy()
            export_data['Ratio'] = export_data['Number'] / export_data['Number'].sum()
            export_data.to_excel(writer, sheet_name='fucose_count')
            export_data = self.sialicacid_structure.copy()
            export_data['Ratio'] = export_data['Number'] / export_data['Number'].sum()
            export_data.to_excel(writer, sheet_name='sialicacid_structure')
            export_data = self.fucose_structure.copy()
            export_data['Ratio'] = export_data['Number'] / export_data['Number'].sum()
            export_data.to_excel(writer, sheet_name='fucose_structure')
            export_data = self.acgc.copy()
            export_data['Ratio'] = export_data['Ac/Gc_count'] / export_data['Ac/Gc_count'].sum()
            export_data['Global ratio'] = export_data['Ac/Gc_count'] / self.data.shape[0]
            export_data.to_excel(writer, sheet_name='acgc')
            export_data = self.lacdinac_count.copy()
            export_data['Ratio'] = export_data['lacdinac_count'] / export_data['lacdinac_count'].sum()
            export_data['Global ratio'] = export_data['lacdinac_count'] / self.data.shape[0]
            export_data.to_excel(writer, sheet_name='lacdinac_count')
            export_data = self.fucosylated_type.copy()
            export_data['Ratio'] = export_data['FucosylatedType_count'] / export_data['FucosylatedType_count'].sum()
            export_data['Global ratio'] = export_data['FucosylatedType_count'] / self.data.shape[0]
            export_data.to_excel(writer, sheet_name='fucosylated_type')
            export_data = self.fsg.copy()
            export_data['Ratio'] = export_data['FSG_count'] / export_data['FSG_count'].sum()
            export_data.to_excel(writer, sheet_name='fsg')
            
            export_data = self.struc_cod_core_structure.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_core_structure'[:31])
            export_data = self.struc_cod_glycan_type.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_glycan_type'[:31])
            export_data = self.struc_cod_branches_structure.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_branches_structure'[:31])
            export_data = self.struc_cod_branches_count.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_branches_count'[:31])
            export_data = self.struc_cod_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_sialicacid_count'[:31])
            export_data = self.struc_cod_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_fucose_count'[:31])
            export_data = self.struc_cod_fucosylated_type.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_fucosylated_type'[:31])
            export_data = self.struc_cod_acgc.copy()
            export_data['Ratio'] = export_data.groupby('structure_coding')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='struc_cod_acgc'[:31])
            
            export_data = self.lacdinac_core_structure.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_core_structure'[:31])
            export_data = self.lacdinac_glycan_type.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_glycan_type'[:31])
            export_data = self.lacdinac_branches_structure.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_branches_structure'[:31])
            export_data = self.lacdinac_branches_count.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_branches_count'[:31])
            export_data = self.lacdinac_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_sialicacid_count'[:31])
            export_data = self.lacdinac_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_fucose_count'[:31])
            export_data = self.lacdinac_fucosylated_type.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_fucosylated_type'[:31])
            export_data = self.lacdinac_acgc.copy()
            export_data['Ratio'] = export_data.groupby('Lacdinac')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='lacdinac_acgc'[:31])
            
            export_data = self.core_structure_lacdinac.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_lacdinac'[:31])
            export_data = self.core_structure_glycan_type.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_glycan_type'[:31])
            export_data = self.core_structure_branches_structure.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_branches_structure'[:31])
            export_data = self.core_structure_branches_count.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_branches_count'[:31])
            export_data = self.core_structure_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_sialicacid_count'[:31])
            export_data = self.core_structure_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_fucose_count'[:31])
            export_data = self.core_structure_fucosylated_type.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_fucosylated_type'[:31])
            export_data = self.core_structure_acgc.copy()
            export_data['Ratio'] = export_data.groupby('Core_structure')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='core_structure_acgc'[:31])
            
            export_data = self.glycan_type_lacdinac.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_lacdinac'[:31])
            export_data = self.glycan_type_core_structure.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_core_structure'[:31])
            export_data = self.glycan_type_branches_structure.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_branches_structure'[:31])
            export_data = self.glycan_type_branches_count.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_branches_count'[:31])
            export_data = self.glycan_type_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_sialicacid_count'[:31])
            export_data = self.glycan_type_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_fucose_count'[:31])
            export_data = self.glycan_type_fucosylated_type.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_fucosylated_type'[:31])
            export_data = self.glycan_type_acgc.copy()
            export_data['Ratio'] = export_data.groupby('Glycan_type')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='glycan_type_acgc'[:31])
            
            export_data = self.branches_structure_core_structure.copy()
            export_data['Ratio'] = export_data.groupby('Branches')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_structure_core_structure'[:31])
            export_data = self.branches_structure_glycan_type.copy()
            export_data['Ratio'] = export_data.groupby('Branches')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_structure_glycan_type'[:31])
            export_data = self.branches_structure_fucosylated_type.copy()
            export_data['Ratio'] = export_data.groupby('Branches')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_structure_fucosylated_type'[:31])
            export_data = self.branches_structure_acgc.copy()
            export_data['Ratio'] = export_data.groupby('Branches')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_structure_acgc'[:31])
            export_data = self.branches_structure_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('Branches')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_structure_sialicacid_count'[:31])
            export_data = self.branches_structure_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('Branches')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_structure_fucose_count'[:31])
                        
            export_data = self.branches_count_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('BranchNumber')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_count_sialicacid_count'[:31])
            export_data = self.branches_count_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('BranchNumber')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='branches_count_fucose_count'[:31])
            
            export_data = self.fucosylated_type_lacdinac.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_lacdinac'[:31])
            export_data = self.fucosylated_type_core_structure.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_core_structure'[:31])
            export_data = self.fucosylated_type_glycan_type.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_glycan_type'[:31])
            export_data = self.fucosylated_type_branches_structure.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_branches_structure'[:31])
            export_data = self.fucosylated_type_branches_count.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_branches_count'[:31])
            export_data = self.fucosylated_type_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_sialicacid_count'[:31])
            export_data = self.fucosylated_type_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_fucose_count'[:31])
            export_data = self.fucosylated_type_acgc.copy()
            export_data['Ratio'] = export_data.groupby('FucosylatedType')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='fucosylated_type_acgc'[:31])
            
            export_data = self.acgc_lacdinac.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_lacdinac'[:31])
            export_data = self.acgc_core_structure.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_core_structure'[:31])
            export_data = self.acgc_glycan_type.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_glycan_type'[:31])
            export_data = self.acgc_branches_structure.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_branches_structure'[:31])
            export_data = self.acgc_branches_count.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_branches_count'[:31])
            export_data = self.acgc_sialicacid_count.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_sialicacid_count'[:31])
            export_data = self.acgc_fucose_count.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_fucose_count'[:31])
            export_data = self.acgc_fucosylated_type.copy()
            export_data['Ratio'] = export_data.groupby('Ac/Gc')['Count'].transform(lambda x: x / x.sum())
            export_data.to_excel(writer, sheet_name='acgc_fucosylated_type'[:31])
            
            
    
