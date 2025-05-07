# -*- coding: utf-8 -*-
## 功能注释与关联性分析模块--81
class StrucGAP_FunctionAnnotation:
    def __init__(self, gs_data, data_manager):
        if hasattr(gs_data, 'sample_group'):
            self.sample_group = gs_data.sample_group
        else:
            self.sample_group = module4.sample_group
        
        if isinstance(gs_data, pd.DataFrame):
            self.gs_data = gs_data
            # self.gs_data = self.gs_data[~self.gs_data['PeptideSequence+structure_coding+ProteinID'].duplicated()]
            self.fc_data = self.gs_data
            self.data_type = ''
        
        elif (gs_data.__class__.__name__ == 'StrucGAP_Preprocess')|(gs_data.__class__.__name__ == 'StrucGAP_GlycoPeptideQuant'):
        # if hasattr(gs_data, 'data'): # module6.__class__.__name__
            self.gs_data = gs_data.data
            self.gs_data = self.gs_data[~self.gs_data['PeptideSequence+structure_coding+ProteinID'].duplicated()]
            if hasattr(gs_data, 'fc_result'):
                self.fc_data = pd.merge(self.gs_data, gs_data.fc_result, left_index=True, right_index=True, how='left')
                self.fc_data = self.fc_data[~self.fc_data['fc'].isnull()]
            elif not hasattr(gs_data, 'fc_result'):
                self.fc_data = self.gs_data
            self.data_type = ''
        elif gs_data.__class__.__name__ == 'StrucGAP_GlycoNetwork':
            data_type = input("Please enter your expected data in the current module (such as: 'protein_no_glyco_up'): ")
            expected_options = ['protein_up_glyco_up', 'protein_up_glyco_no', 'protein_up_glyco_down', 
                                'protein_no_glyco_up', 'protein_no_glyco_no', 'protein_no_glyco_down',
                                'protein_down_glyco_up', 'protein_down_glyco_no', 'protein_down_glyco_down',]
            matches = get_close_matches(data_type, expected_options, n=1, cutoff=0.5)
            if matches:
                data_type = matches[0]
                print(f"Using '{data_type}' as the input.")
            else:
                print("No close match found. Using 'protein_no_glyco_up' as the input.")
                data_type = 'protein_no_glyco_up'

            self.gs_data = getattr(gs_data, data_type, None)
            self.fc_data = self.gs_data
            if self.gs_data is not None:
                self.fc_data = self.fc_data.rename(columns={'p_g': 'pvalue_ttest_mannwhitneyu', 'fc_g': 'fc'})
                self.fc_data = self.fc_data[~self.fc_data['fc'].isnull()]
                self.data_type = data_type
            else:
                raise ValueError(f"No data found for {data_type} in gs_data.")
        
        # if GO_data_dir is None:
        #     GO_data_dir = self.download_go_file()
        
        # self.GO_data = self.parse_go_obo(GO_data_dir)
        
        self.data_manager = data_manager    
        self.data_manager.register_module('StrucGAP_FunctionAnnotation', self, {})
        # self.data_manager.log_params('StrucGAP_FunctionAnnotation', 'GO_data', {'GO_data': GO_data_dir})
        
    # def download_go_file(self):
    #     """An auxiliary function."""
    #     url = "http://current.geneontology.org/ontology/go-basic.obo"
    #     file_path = os.path.join(os.getcwd(), "go-basic.obo")
    #     if not os.path.exists(file_path):
    #         print("Can not find GO file, is downloading go-basic.obo ...")
    #         response = requests.get(url, stream=True)
    #         if response.status_code == 200:
    #             with open(file_path, "wb") as f:
    #                 for chunk in response.iter_content(chunk_size=8192):
    #                     f.write(chunk)
    #             print("Downloading finished. ")
    #         else:
    #             raise Exception("Can not download GO file, please try again. ")
    #     else:
    #         print("GO file already exists, skip downloading. ")
    #     return file_path
        
    # def parse_go_obo(self, file_path):
    #     """An auxiliary function."""
    #     nodes1 = []
    #     names1 = []
    #     namespaces = []
    #     relationships = []
    #     nodes2 = []
    #     names2 = []
    #     def append_all_lists():
    #         nodes1.append(node1)
    #         names1.append(name1)
    #         namespaces.append(namespace)
    #         relationships.append(relationship)
    #         nodes2.append(node2)
    #         names2.append(name2)
    #     with open(file_path, 'r', encoding='utf-8') as f:
    #         lines = f.readlines()
    #     for line in tqdm(lines):
    #         line = line.strip()
    #         if line.startswith('id:'):
    #             node1 = line[line.find('GO:'):].strip()
    #         elif line.startswith('name:'):
    #             name1 = line[5:].strip()
    #         elif line.startswith('namespace:'):
    #             namespace = line[10:].strip()
    #         elif line.startswith('is_a:'):
    #             relationship = 'is a'
    #             node2 = line.split('!')[0][5:].strip()
    #             name2 = line.split('!')[1].strip()
    #             append_all_lists()
    #         elif line.startswith('relationship:'):
    #             relationship = line.split('GO:')[0][13:].strip()
    #             node2 = line[line.find('GO:'):].split('!')[0].strip()
    #             name2 = line.split('!')[1].strip()
    #             append_all_lists()
    #     GO_data = pd.DataFrame({
    #         'GO_id1': nodes1,
    #         'GO_name1': names1,
    #         'type': namespaces,
    #         'relationship': relationships,
    #         'GO_id2': nodes2,
    #         'GO_name2': names2
    #     })
    #     return GO_data
        
    def compute_weighted_fc(self, df, up_down_fc_threshold):
        """An auxiliary function called by other functions to compute weight FC."""
        protein_fc = {}  
        upregulated_proteins = []    
        downregulated_proteins = []  
    
        protein_groups = df.groupby('GeneName')
        for protein_id, group in protein_groups:
            fc_values = group['fc'].values
            fc_up = group[group['fc'] > up_down_fc_threshold]
            fc_down = group[group['fc'] < 1/up_down_fc_threshold]
            if not fc_up.empty and not fc_down.empty:
                intensities = []
                for index, row in group.iterrows():
                    fc = row['fc']
                    matched_reporter_ions_str = row['Matched_Reporter_Ions']
                    if pd.isnull(matched_reporter_ions_str):
                        continue
                    try:
                        matched_reporter_ions = ast.literal_eval(matched_reporter_ions_str)
                    except (ValueError, SyntaxError):
                        continue
                    total_intensity = sum([ion_info[2] for ion_info in matched_reporter_ions.values()])
                    intensities.append((fc, total_intensity))
                if not intensities:
                    continue
                total_weight = sum([intensity for fc, intensity in intensities])
                weighted_fc = sum([fc * intensity for fc, intensity in intensities]) / total_weight
                protein_fc[protein_id] = weighted_fc
            else:
                fc_mean = group['fc'].mean()
                protein_fc[protein_id] = fc_mean
    
        protein_fc_unique = {}
    
        for protein_id_combined, fc_value in protein_fc.items():
            protein_ids = protein_id_combined.split(';')
            for pid in protein_ids:
                pid = pid.strip()
                if pid == 'Unknown':
                    continue
                if pid in protein_fc_unique:
                    protein_fc_unique[pid]['fc_values'].append(fc_value)
                    protein_fc_unique[pid]['weights'].append(1)
                else:
                    protein_fc_unique[pid] = {'fc_values': [fc_value], 'weights': [1]}
    
        final_protein_fc = {}
        for pid, data in protein_fc_unique.items():
            fc_values = data['fc_values']
            weights = data['weights']
            total_weight = sum(weights)
            weighted_fc = sum([fc * w for fc, w in zip(fc_values, weights)]) / total_weight
            final_protein_fc[pid] = weighted_fc
            if weighted_fc > up_down_fc_threshold:
                upregulated_proteins.append(pid)
            elif weighted_fc < 1/up_down_fc_threshold:
                downregulated_proteins.append(pid)
    
        upregulated_proteins = list(set(upregulated_proteins))
        downregulated_proteins = list(set(downregulated_proteins))
    
        return final_protein_fc, upregulated_proteins, downregulated_proteins
        
    def ora(self, organism=None, up_down_fc_threshold=None, background_input=False, pvalue_type='pvalue_ttest_mannwhitneyu',
            selected_terms=['GO:MF', 'GO:CC', 'GO:BP']):
        """
        Performs over-representation analysis using the g:Profiler API.
        
        Parameters:
            organism: database from ['mmusculus','hsapiens', 'rnorvegicus'].
            up_down_fc_threshold: FC threshold used to differentiate up and down regulated features.
            background_input: use both proteins as background or not.
            pvalue_type: 'pvalue_ttest', 'pvalue_mannwhitneyu' or 'pvalue_ttest_mannwhitneyu'.
            selected_terms: enrichment database.
        
        Returns:
            self.final_protein_fc
            self.upregulated_proteins
            self.downregulated_proteins
            self.both_proteins
            self.ora_no_background_up_result
            self.ora_no_background_down_result
            self.ora_background_up_result
            self.ora_background_down_result
            self.ora_no_background_both_proteins_result
            self.ora_background_both_proteins_result 
                
        Return type:
            dataframe
        
        """
        if up_down_fc_threshold == None:
            up_down_fc_threshold = float(input("Please enter the fc threshold (such as 1.5) to differentiate up and down regulated features: "))
 
        if organism == None:
            organism = input("Please enter database from ['mmusculus','hsapiens', 'rnorvegicus']: ")
            expected_options = ['mmusculus','hsapiens', 'rnorvegicus']
            matches = get_close_matches(organism, expected_options, n=1, cutoff=0.5)
            if matches:
                organism = matches[0]
                print(f"Using '{organism}' as the input.")
            else:
                print("No close match found. Using 'mmusculus' as the input.")
                organism = 'mmusculus'

            if organism not in ['mmusculus','hsapiens', 'rnorvegicus']:
                raise ValueError("Please enter the correct database name")
        # enrichr
        terms = ["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP"]
        print("Please select a term by entering a number, and separate multiple terms with a comma: ")
        for i, term in enumerate(terms, start=1):
            print(f"{i}. {term}")
        selected_indices = input("Please enter the term number you want to select (e.g. 1,3,5): ")
        indices = selected_indices.split(',')
        selected_terms = []
        for idx in indices:
            idx = idx.strip()  
            if idx.isdigit():
                idx_int = int(idx)
                if 1 <= idx_int <= len(terms):
                    selected_terms.append(terms[idx_int - 1])  
                else:
                    print(f"number {idx_int} out of range. ")
            else:
                print(f"input number '{idx}' not a valid number. ")
        # gene_list
        self.final_protein_fc = None  
        self.upregulated_proteins = None 
        self.downregulated_proteins = None 
        self.both_proteins = None 
        if pvalue_type in self.fc_data.columns:       
            fc_data_p = self.fc_data[self.fc_data[pvalue_type]<0.05]
            self.final_protein_fc, self.upregulated_proteins, self.downregulated_proteins = self.compute_weighted_fc(fc_data_p, up_down_fc_threshold)
        elif pvalue_type not in self.fc_data.columns:
            fc_data_p = self.fc_data
            self.both_proteins = list(set(fc_data_p['GeneName']))
        
        background = pd.DataFrame(self.gs_data[~self.gs_data['GeneName'].duplicated()]['GeneName'])
        background = background['GeneName'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).to_frame('GeneName')
        background = background[background['GeneName']!='Unknown']
        
        # 初始化结果属性
        self.ora_no_background_up_result = pd.DataFrame()
        self.ora_no_background_down_result = pd.DataFrame()
        self.ora_background_up_result = pd.DataFrame()
        self.ora_background_down_result = pd.DataFrame()
        self.ora_no_background_both_proteins_result = pd.DataFrame()
        self.ora_background_both_proteins_result = pd.DataFrame()
        self.gsea_result = pd.DataFrame()
            
        if background_input == False:
            if self.upregulated_proteins:  
                gpf = GProfiler(return_dataframe=True)
                enr = gpf.profile(organism=organism,
                       no_evidences=False,
                       sources = selected_terms,
                       user_threshold = 1,
                       all_results = True,
                       background = None,
                       query=self.upregulated_proteins)
                enr = enr[['source', 'native', 'name', 'p_value', 'term_size', 'intersection_size', 'parents', 'intersections']]
                enr['Term'] = enr.apply(lambda row: f"{row['name']} ({row['native']})", axis=1)
                enr['Overlap'] = enr.apply(lambda row: f"{row['intersection_size']}/{row['term_size']}", axis=1)
                enr = enr[['source', 'Term', 'name', 'p_value', 'Overlap', 'parents', 'intersections']]
                enr['intersections'] = enr['intersections'].apply(lambda x: ';'.join([item.upper() for item in x]))
                enr.columns = ['Gene_set', 'Term', 'GO_name', 'P-value', 'Overlap', 'parents', 'Genes']
                enr['GO_name'] = enr['GO_name'].str.upper()
                
                self.ora_no_background_up_result = enr
            if self.downregulated_proteins: 
                gpf = GProfiler(return_dataframe=True)
                enr = gpf.profile(organism=organism,
                       no_evidences=False,
                       sources = selected_terms,
                       user_threshold = 1,
                       all_results = True,
                       background = None,
                       query=self.downregulated_proteins)
                enr = enr[['source', 'native', 'name', 'p_value', 'term_size', 'intersection_size', 'parents', 'intersections']]
                enr['Term'] = enr.apply(lambda row: f"{row['name']} ({row['native']})", axis=1)
                enr['Overlap'] = enr.apply(lambda row: f"{row['intersection_size']}/{row['term_size']}", axis=1)
                enr = enr[['source', 'Term', 'name', 'p_value', 'Overlap', 'parents', 'intersections']]
                enr['intersections'] = enr['intersections'].apply(lambda x: ';'.join([item.upper() for item in x]))
                enr.columns = ['Gene_set', 'Term', 'GO_name', 'P-value', 'Overlap', 'parents', 'Genes']
                enr['GO_name'] = enr['GO_name'].str.upper()
                
                self.ora_no_background_down_result = enr
            if self.both_proteins:  
                gpf = GProfiler(return_dataframe=True)
                enr = gpf.profile(organism=organism,
                       no_evidences=False,
                       sources = selected_terms,
                       user_threshold = 1,
                       all_results = True,
                       background = None,
                       query=self.both_proteins)
                enr = enr[['source', 'native', 'name', 'p_value', 'term_size', 'intersection_size', 'parents', 'intersections']]
                enr['Term'] = enr.apply(lambda row: f"{row['name']} ({row['native']})", axis=1)
                enr['Overlap'] = enr.apply(lambda row: f"{row['intersection_size']}/{row['term_size']}", axis=1)
                enr = enr[['source', 'Term', 'name', 'p_value', 'Overlap', 'parents', 'intersections']]
                enr['intersections'] = enr['intersections'].apply(lambda x: ';'.join([item.upper() for item in x]))
                enr.columns = ['Gene_set', 'Term', 'GO_name', 'P-value', 'Overlap', 'parents', 'Genes']
                enr['GO_name'] = enr['GO_name'].str.upper()
                
                self.ora_no_background_both_proteins_result = enr
        
        elif background_input == True:
            if self.upregulated_proteins:  
                gpf = GProfiler(return_dataframe=True)
                enr_bg = gpf.profile(organism=organism,
                       no_evidences=False,
                       sources = selected_terms,
                       user_threshold = 1,
                       all_results = True,
                       background = list(background['GeneName']),
                       query=self.upregulated_proteins)
                enr_bg = enr_bg[['source', 'native', 'name', 'p_value', 'term_size', 'intersection_size', 'parents', 'intersections']]
                enr_bg['Term'] = enr_bg.apply(lambda row: f"{row['name']} ({row['native']})", axis=1)
                enr_bg['Overlap'] = enr_bg.apply(lambda row: f"{row['intersection_size']}/{row['term_size']}", axis=1)
                enr_bg = enr_bg[['source', 'Term', 'name', 'p_value', 'Overlap', 'parents', 'intersections']]
                enr_bg['intersections'] = enr_bg['intersections'].apply(lambda x: ';'.join([item.upper() for item in x]))
                enr_bg.columns = ['Gene_set', 'Term', 'GO_name', 'P-value', 'Overlap', 'parents', 'Genes']
                enr_bg['GO_name'] = enr_bg['GO_name'].str.upper()
                
                self.ora_background_up_result = enr_bg
            if self.downregulated_proteins:  
                gpf = GProfiler(return_dataframe=True)
                enr_bg = gpf.profile(organism=organism,
                       no_evidences=False,
                       sources = selected_terms,
                       user_threshold = 1,
                       all_results = True,
                       background = list(background['GeneName']),
                       query=self.downregulated_proteins)
                enr_bg = enr_bg[['source', 'native', 'name', 'p_value', 'term_size', 'intersection_size', 'parents', 'intersections']]
                enr_bg['Term'] = enr_bg.apply(lambda row: f"{row['name']} ({row['native']})", axis=1)
                enr_bg['Overlap'] = enr_bg.apply(lambda row: f"{row['intersection_size']}/{row['term_size']}", axis=1)
                enr_bg = enr_bg[['source', 'Term', 'name', 'p_value', 'Overlap', 'parents', 'intersections']]
                enr_bg['intersections'] = enr_bg['intersections'].apply(lambda x: ';'.join([item.upper() for item in x]))
                enr_bg.columns = ['Gene_set', 'Term', 'GO_name', 'P-value', 'Overlap', 'parents', 'Genes']
                enr_bg['GO_name'] = enr_bg['GO_name'].str.upper()
                
                self.ora_background_down_result = enr_bg
            if self.both_proteins:  
                gpf = GProfiler(return_dataframe=True)
                enr_bg = gpf.profile(organism=organism,
                       no_evidences=False,
                       sources = selected_terms,
                       user_threshold = 1,
                       all_results = True,
                       background = list(background['GeneName']),
                       query=self.both_proteins)
                enr_bg = enr_bg[['source', 'native', 'name', 'p_value', 'term_size', 'intersection_size', 'parents', 'intersections']]
                enr_bg['Term'] = enr_bg.apply(lambda row: f"{row['name']} ({row['native']})", axis=1)
                enr_bg['Overlap'] = enr_bg.apply(lambda row: f"{row['intersection_size']}/{row['term_size']}", axis=1)
                enr_bg = enr_bg[['source', 'Term', 'name', 'p_value', 'Overlap', 'parents', 'intersections']]
                enr_bg['intersections'] = enr_bg['intersections'].apply(lambda x: ';'.join([item.upper() for item in x]))
                enr_bg.columns = ['Gene_set', 'Term', 'GO_name', 'P-value', 'Overlap', 'parents', 'Genes']
                enr_bg['GO_name'] = enr_bg['GO_name'].str.upper()
                
                self.ora_background_both_proteins_result = enr_bg

        self.selected_terms = selected_terms
        self.organism = organism
        #
        self.data_manager.log_params('StrucGAP_FunctionAnnotation', 'ora', {'organism':organism, 
                                                                       'background_input':background_input,
                                                                       'terms':self.selected_terms, 
                                                                       'up_down_fc_threshold':up_down_fc_threshold})

        output_data = {
            'ora_no_background_up_result': self.ora_no_background_up_result,
            'ora_no_background_down_result': self.ora_no_background_down_result,
            'ora_background_up_result': self.ora_background_up_result,
            'ora_background_down_result': self.ora_background_down_result,
            'ora_no_background_both_proteins_result': self.ora_no_background_both_proteins_result,
            'ora_background_both_proteins_result': self.ora_background_both_proteins_result,
        }
        
        for key, value in output_data.items():
            if value is not None and not value.empty:
                self.data_manager.log_output('StrucGAP_FunctionAnnotation', key, value)

        return self
    
    def gsea(self):
        """
        Conducts GSEA using the GSEApy.
        
        Parameters:
            None.
        
        Returns:
            self.gsea_result
                
        Return type:
            dataframe
        
        """
        if self.organism == 'mmusculus' or 'rnorvegicus':
            organism = 'Mouse'
        elif self.organism == 'hsapiens':
            organism = 'Human'

        # enrichr
        terms = gp.get_library_name(organism=organism)
        print("Please select a term by entering a number, and separate multiple terms with a comma: ")
        for i, term in enumerate(terms, start=1):
            print(f"{i}. {term}")
        selected_indices = input("Please enter the term number you want to select (e.g. 1,3,5): ")
        indices = selected_indices.split(',')
        selected_terms = []
        for idx in indices:
            idx = idx.strip()  
            if idx.isdigit():
                idx_int = int(idx)
                if 1 <= idx_int <= len(terms):
                    selected_terms.append(terms[idx_int - 1])  
                else:
                    print(f"number {idx_int} out of range. ")
            else:
                print(f"input number '{idx}' not a valid number. ")
        
        protein_fc = pd.DataFrame.from_dict(self.final_protein_fc, orient='index')
        protein_fc.columns = ['fc']  
        protein_fc['fc'] = np.log2(protein_fc['fc'])
        protein_fc = protein_fc.sort_values(by='fc', ascending=False)
        protein_fc.index = protein_fc.index.str.upper()

        upregulated_protein_fc = pd.concat([pd.DataFrame(index=self.upregulated_proteins), protein_fc], axis=1, join='inner')
        upregulated_protein_fc['fc'] = np.log2(upregulated_protein_fc['fc'])
        upregulated_protein_fc = upregulated_protein_fc.sort_values(by='fc', ascending=False)
        downregulated_protein_fc = pd.concat([pd.DataFrame(index=self.downregulated_proteins), protein_fc], axis=1, join='inner')
        downregulated_protein_fc['fc'] = np.log2(downregulated_protein_fc['fc'])
        downregulated_protein_fc = downregulated_protein_fc.sort_values(by='fc', ascending=False)
        
        # # run prerank
        pre_res = gp.prerank(rnk=protein_fc, 
                             gene_sets=selected_terms,
                             threads=6,
                             min_size=3,
                             max_size=2000,
                             permutation_num=1000, 
                             outdir=None, 
                             seed=6,
                             verbose=True, 
                            )
        
        self.gsea_result = pre_res.res2d
        #    
        self.data_manager.log_params('StrucGAP_FunctionAnnotation', 'gsea', {})
        self.data_manager.log_output('StrucGAP_FunctionAnnotation', 'gsea_result', self.gsea_result)

        return self
    
    def identify_core_structure(self, data):  
        """An auxiliary function called by other functions to identify core structure."""
        temp_data = data[(data['GlycanComposition']!='N2H2')&(data['GlycanComposition']!='N2H2F1')]
        temp_data = temp_data[['structure_coding','MS2Scan', 'Corefucose', 'Bisection']] # data1 = self.data
        temp_data = temp_data[~temp_data['MS2Scan'].duplicated()]
        temp_data = temp_data.drop(columns='MS2Scan')
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
    
    def go_function_structure(self, function_data = None, p_value='P-value', cutoff=0.05): 
        """
        Performs integrated analysis between enriched GO terms and glycan substructural features.
        
        Parameters:
            function_data: enrichment result in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result'].
            p_value: ['P-value', 'Adjusted P-value'] used to screen enrichment terms.
            cutoff: p_value threshold.
                
        Returns:
            self.bp_core_structure
            self.mf_core_structure
            self.cc_core_structure
            self.bp_glycan_type
            self.mf_glycan_type
            self.cc_glycan_type
            self.bp_branches_structure
            self.mf_branches_structure
            self.cc_branches_structure
            self.bp_branches_count
            self.mf_branches_count
            self.cc_branches_count
            self.bp_sialicacid_count
            self.mf_sialicacid_count
            self.cc_sialicacid_count
            self.bp_fucose_count
            self.mf_fucose_count
            self.cc_fucose_count
            self.bp_sialicacid_structure
            self.mf_sialicacid_structure
            self.cc_sialicacid_structure
            self.bp_fucose_structure
            self.mf_fucose_structure
            self.cc_fucose_structure
            self.bp_lacdinac
            self.mf_lacdinac
            self.cc_lacdinac
            self.bp_structurecoding
            self.mf_structurecoding
            self.cc_structurecoding
            self.bp_fucosylated_type
            self.mf_fucosylated_type
            self.cc_fucosylated_type
            self.bp_acgc
            self.mf_acgc
            self.cc_acgc
            self.bp_core_structure_count
            self.mf_core_structure_count
            self.cc_core_structure_count
            self.bp_glycan_type_count
            self.mf_glycan_type_count
            self.cc_glycan_type_count
            self.bp_branches_structure_count
            self.mf_branches_structure_count
            self.cc_branches_structure_count
            self.bp_branches_count_count
            self.mf_branches_count_count
            self.cc_branches_count_count
            self.bp_sialicacid_count_count
            self.mf_sialicacid_count_count
            self.cc_sialicacid_count_count
            self.bp_fucose_count_count
            self.mf_fucose_count_count
            self.cc_fucose_count_count
            self.bp_sialicacid_structure_count
            self.mf_sialicacid_structure_count
            self.cc_sialicacid_structure_count
            self.bp_fucose_structure_count
            self.mf_fucose_structure_count
            self.cc_fucose_structure_count
            self.bp_lacdinac_count
            self.mf_lacdinac_count
            self.cc_lacdinac_count
            self.bp_structurecoding_count
            self.mf_structurecoding_count
            self.cc_structurecoding_count
            self.bp_fucosylated_type_count
            self.mf_fucosylated_type_count
            self.cc_fucosylated_type_count
            self.bp_acgc_count
            self.mf_acgc_count
            self.cc_acgc_count
                
        Return type:
            dataframe
        
        """
        self.bp_core_structure = pd.DataFrame()
        self.mf_core_structure = pd.DataFrame()
        self.cc_core_structure = pd.DataFrame()
        self.bp_glycan_type = pd.DataFrame()
        self.mf_glycan_type = pd.DataFrame()
        self.cc_glycan_type = pd.DataFrame()
        self.bp_branches_structure = pd.DataFrame()
        self.mf_branches_structure = pd.DataFrame()
        self.cc_branches_structure = pd.DataFrame()
        self.bp_branches_count = pd.DataFrame()
        self.mf_branches_count = pd.DataFrame()
        self.cc_branches_count = pd.DataFrame()
        self.bp_sialicacid_count = pd.DataFrame()
        self.mf_sialicacid_count = pd.DataFrame()
        self.cc_sialicacid_count = pd.DataFrame()
        self.bp_fucose_count = pd.DataFrame()
        self.mf_fucose_count = pd.DataFrame()
        self.cc_fucose_count = pd.DataFrame()
        self.bp_sialicacid_structure = pd.DataFrame()
        self.mf_sialicacid_structure = pd.DataFrame()
        self.cc_sialicacid_structure = pd.DataFrame()
        self.bp_fucose_structure = pd.DataFrame()
        self.mf_fucose_structure = pd.DataFrame()
        self.cc_fucose_structure = pd.DataFrame()
        self.bp_lacdinac = pd.DataFrame()
        self.mf_lacdinac = pd.DataFrame()
        self.cc_lacdinac = pd.DataFrame()
        self.bp_structurecoding = pd.DataFrame()
        self.mf_structurecoding = pd.DataFrame()
        self.cc_structurecoding = pd.DataFrame()
        self.bp_fucosylated_type = pd.DataFrame()
        self.mf_fucosylated_type = pd.DataFrame()
        self.cc_fucosylated_type = pd.DataFrame()
        self.bp_acgc = pd.DataFrame()
        self.mf_acgc = pd.DataFrame()
        self.cc_acgc = pd.DataFrame()
        self.bp_core_structure_count = pd.DataFrame()
        self.mf_core_structure_count = pd.DataFrame()
        self.cc_core_structure_count = pd.DataFrame()
        self.bp_glycan_type_count = pd.DataFrame()
        self.mf_glycan_type_count = pd.DataFrame()
        self.cc_glycan_type_count = pd.DataFrame()
        self.bp_branches_structure_count = pd.DataFrame()
        self.mf_branches_structure_count = pd.DataFrame()
        self.cc_branches_structure_count = pd.DataFrame()
        self.bp_branches_count_count = pd.DataFrame()
        self.mf_branches_count_count = pd.DataFrame()
        self.cc_branches_count_count = pd.DataFrame()
        self.bp_sialicacid_count_count = pd.DataFrame()
        self.mf_sialicacid_count_count = pd.DataFrame()
        self.cc_sialicacid_count_count = pd.DataFrame()
        self.bp_fucose_count_count = pd.DataFrame()
        self.mf_fucose_count_count = pd.DataFrame()
        self.cc_fucose_count_count = pd.DataFrame()
        self.bp_sialicacid_structure_count = pd.DataFrame()
        self.mf_sialicacid_structure_count = pd.DataFrame()
        self.cc_sialicacid_structure_count = pd.DataFrame()
        self.bp_fucose_structure_count = pd.DataFrame()
        self.mf_fucose_structure_count = pd.DataFrame()
        self.cc_fucose_structure_count = pd.DataFrame()
        self.bp_lacdinac_count = pd.DataFrame()
        self.mf_lacdinac_count = pd.DataFrame()
        self.cc_lacdinac_count = pd.DataFrame()
        self.bp_structurecoding_count = pd.DataFrame()
        self.mf_structurecoding_count = pd.DataFrame()
        self.cc_structurecoding_count = pd.DataFrame()
        self.bp_fucosylated_type_count = pd.DataFrame()
        self.mf_fucosylated_type_count = pd.DataFrame()
        self.cc_fucosylated_type_count = pd.DataFrame()
        self.bp_acgc_count = pd.DataFrame()
        self.mf_acgc_count = pd.DataFrame()
        self.cc_acgc_count = pd.DataFrame()
        
        if function_data == None:
            function_data = input("Please enter the data you want to execute function and structure correlation analysis (select from ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']): ")
            expected_options = ['ora_no_background_up_result', 'ora_no_background_down_result',
                                'ora_background_up_result', 'ora_background_up_result', 'gsea_result',
                                'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']
            matches = get_close_matches(function_data, expected_options, n=1, cutoff=0.5)
            if matches:
                function_data = matches[0]
                print(f"Using '{function_data}' as the input.")
            else:
                print("No close match found. Using 'ora_no_background_up_result' as the input.")
                function_data = 'ora_no_background_up_result'
            
            if function_data not in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']:
                raise ValueError('Please enter the correct data name')
        
        # GO_data = self.GO_data
        self.function_data_name = function_data
        
        #
        if function_data in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']:
            if function_data in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']:
                function_data = getattr(self, function_data, None)
                # function_data['GO_name'] = function_data['Term'].str.split(' \(', n=1, expand=True)[0].str.upper()
                function_data = function_data.set_index('GO_name',drop=False)
                if p_value is None:
                    p_value = input('Please enter the type of p_value (select from: P-value, Adjusted P-value) you would like to use to select the differential pathway: ')
                    expected_options = ['P-value', 'Adjusted P-value']
                    matches = get_close_matches(p_value, expected_options, n=1, cutoff=0.5)
                    if matches:
                        p_value = matches[0]
                        print(f"Using '{p_value}' as the input.")
                    else:
                        print("No close match found. Using 'P-value' as the input.")
                        p_value = 'P-value'
                if cutoff is None:
                    cutoff = input('Please enter the cutoff value (such as: 0.05, 0.01) you would like to use to select the differential pathway: ')
                function_data = function_data[function_data[p_value]<float(cutoff)]
            elif function_data == 'gsea_result':
                function_data = getattr(self, function_data, None)
                function_data['Gene_set'] = [x[0] for x in function_data['Term'].str.split('__')]
                function_data['GO_name'] = [x[1] for x in function_data['Term'].str.split('__')]
                function_data['GO_name'] = [x[0].upper() for x in function_data['GO_name'].str.split(' \(')]
                function_data = function_data.set_index('GO_name',drop=False)
                if p_value is None:
                    p_value = input('Please enter the type of p_value (select from: NOM p-val, FDR q-val) you would like to use to select the differential pathway: ')
                    expected_options = ['NOM p-val', 'FDR q-val']
                    matches = get_close_matches(p_value, expected_options, n=1, cutoff=0.5)
                    if matches:
                        p_value = matches[0]
                        print(f"Using '{p_value}' as the input.")
                    else:
                        print("No close match found. Using 'P-value' as the input.")
                        p_value = 'P-value'
                if cutoff is None:
                    cutoff = input('Please enter the cutoff value (such as: 0.05, 0.01) you would like to use to select the differential pathway: ')
                function_data = function_data[function_data[p_value]<float(cutoff)]
                function_data = function_data.rename(columns={'Lead_genes': 'Genes'})
            # GO_data = GO_data[GO_data['relationship']=='part_of']
            # GO_data = GO_data.set_index('GO_name1',drop=False)
            # GO_data['GO_name1'] = GO_data['GO_name1'].str.upper()
            # GO_data['GO_name2'] = GO_data['GO_name2'].str.upper()
            self.function_data = function_data
            
            #
            for i in ['biological_process', 'molecular_function', 'cellular_component']:
                if i == 'biological_process':
                    # bp_data = GO_data[GO_data['type']==i][['GO_name1', 'relationship', 'GO_name2']]
                    # bp_both_data = {}  # function_data[function_data['Gene_set']=='GO_Biological_Process_2023']['GO_name'].unique()     
                    # for j in tqdm(list(function_data[function_data['Gene_set']=='GO_Biological_Process_2023']['GO_name'].unique())):
                    #     bp_son_data = bp_data[bp_data['GO_name1']==j]
                    #     bp_son_data['chain'] = bp_son_data['relationship'] + '+' + bp_son_data['GO_name2']
                    #     bp_son_data.set_index('GO_name2', inplace=True, drop=True)
                    #     bp_son_data = bp_son_data.drop(columns='relationship')
                    #     bp_merge_data = bp_data[bp_data['GO_name1'].isin(list(bp_son_data.index))]
                    #     bp_merge_data.set_index('GO_name1', inplace=True, drop=True)
                    #     num = 0
                    #     while bp_data[bp_data['GO_name1'].isin(list(bp_merge_data.index))].shape[0] != 0:
                    #         bp_merge_data = bp_data[bp_data['GO_name1'].isin(list(bp_merge_data.index))]
                    #         bp_merge_data['chain'] = bp_merge_data['relationship'] + '+' + bp_merge_data['GO_name2']
                    #         bp_merge_data.set_index('GO_name1', inplace=True, drop=True)
                    #         num = num + 1
                    #         bp_merge_data = bp_merge_data.rename(columns={'chain':num})
                    #         bp_merge_data = bp_merge_data.drop(columns='relationship')
                    #         bp_son_data = pd.merge(bp_son_data, bp_merge_data, left_index=True, right_index=True, how='left')
                    #         bp_merge_data.set_index('GO_name2', inplace=True, drop=True)
                    #         bp_son_data.set_index('GO_name2', inplace=True, drop=True)
                    #         bp_son_data_1 = bp_son_data.reset_index() 
                    #         parent = []
                    #         for k in list(bp_son_data_1.index):
                    #             lst = list(bp_son_data_1.loc[k])
                    #             lst = [x for x in lst if pd.notnull(x)]
                    #             last_value = lst[-1]
                    #             target_value = None  
                    #             if last_value == 'is a+BIOLOGICAL_PROCESS':
                    #                 for item in reversed(lst[:-1]):
                    #                     if 'is a+' in item:
                    #                         target_value = item.split('+')[1]
                    #                         break
                    #             parent.append(target_value)
                    #         parent = [x for x in parent if x is not None]
                    #         parent = list(set(parent))
                    #     bp_both_data[j] = parent
                    # bp_both_data_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in bp_both_data.items() ]))
                    # bp_both_data_df = bp_both_data_df.T
                    bp_both_data_df = pd.DataFrame(function_data[function_data['Gene_set'].str.contains('GO:BP', na=False)]['GO_name'])
                    
                elif i == 'molecular_function':
                    # mf_data = GO_data[GO_data['type']==i][['GO_name1', 'relationship', 'GO_name2']]
                    # mf_both_data = {}
                    # for j in tqdm(list(function_data[function_data['Gene_set']=='GO_Molecular_Function_2023']['GO_name'].unique())):
                    #     mf_son_data = mf_data[mf_data['GO_name1']==j]
                    #     mf_son_data['chain'] = mf_son_data['relationship'] + '+' + mf_son_data['GO_name2']
                    #     mf_son_data.set_index('GO_name2', inplace=True, drop=True)
                    #     mf_son_data = mf_son_data.drop(columns='relationship')
                    #     mf_merge_data = mf_data[mf_data['GO_name1'].isin(list(mf_son_data.index))]
                    #     mf_merge_data.set_index('GO_name1', inplace=True, drop=True)
                    #     num = 0
                    #     while mf_data[mf_data['GO_name1'].isin(list(mf_merge_data.index))].shape[0] != 0:
                    #         mf_merge_data = mf_data[mf_data['GO_name1'].isin(list(mf_merge_data.index))]
                    #         mf_merge_data['chain'] = mf_merge_data['relationship'] + '+' + mf_merge_data['GO_name2']
                    #         mf_merge_data.set_index('GO_name1', inplace=True, drop=True)
                    #         num = num + 1
                    #         mf_merge_data = mf_merge_data.rename(columns={'chain':num})
                    #         mf_merge_data = mf_merge_data.drop(columns='relationship')
                    #         mf_son_data = pd.merge(mf_son_data, mf_merge_data, left_index=True, right_index=True, how='left')
                    #         mf_merge_data.set_index('GO_name2', inplace=True, drop=True)
                    #         mf_son_data.set_index('GO_name2', inplace=True, drop=True)
                    #         mf_son_data_1 = mf_son_data.reset_index() 
                    #         parent = []
                    #         for k in list(mf_son_data_1.index): 
                    #             lst = list(mf_son_data_1.loc[k])
                    #             lst = [x for x in lst if pd.notnull(x)]
                    #             last_value = lst[-1]
                    #             target_value = [] 
                    #             if last_value == 'is a+MOLECULAR_FUNCTION':
                    #                 for item in reversed(lst[:-1]):
                    #                     if 'is a+' in item:
                    #                         term = item.split('+')[1]
                    #                         target_value.append(term)
                    #                         if len(target_value) == 3:
                    #                             break
                    #             parent.extend(target_value)
                    #         parent = [x for x in parent if x is not None]
                    #         parent = list(set(parent))
                    #     mf_both_data[j] = parent
                    # mf_both_data_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in mf_both_data.items() ]))
                    # mf_both_data_df = mf_both_data_df.T
                    mf_both_data_df = pd.DataFrame(function_data[function_data['Gene_set'].str.contains('GO:MF', na=False)]['GO_name'])
                    
                elif i == 'cellular_component':
                    # cc_data = GO_data[GO_data['type']==i][['GO_name1', 'relationship', 'GO_name2']]
                    # cc_both_data = {}
                    # for j in tqdm(list(function_data[function_data['Gene_set']=='GO_Cellular_Component_2023']['GO_name'].unique())):
                    #     cc_son_data = cc_data[cc_data['GO_name1']==j]
                    #     cc_son_data['chain'] = cc_son_data['relationship'] + '+' + cc_son_data['GO_name2']
                    #     cc_son_data.set_index('GO_name2', inplace=True, drop=True)
                    #     cc_son_data = cc_son_data.drop(columns='relationship')
                    #     cc_merge_data = cc_data[cc_data['GO_name1'].isin(list(cc_son_data.index))]
                    #     cc_merge_data.set_index('GO_name1', inplace=True, drop=True)
                    #     num = 0
                    #     while cc_data[cc_data['GO_name1'].isin(list(cc_merge_data.index))].shape[0] != 0:
                    #         cc_merge_data = cc_data[cc_data['GO_name1'].isin(list(cc_merge_data.index))]
                    #         cc_merge_data['chain'] = cc_merge_data['relationship'] + '+' + cc_merge_data['GO_name2']
                    #         cc_merge_data.set_index('GO_name1', inplace=True, drop=True)
                    #         num = num + 1
                    #         cc_merge_data = cc_merge_data.rename(columns={'chain':num})
                    #         cc_merge_data = cc_merge_data.drop(columns='relationship')
                    #         cc_son_data = pd.merge(cc_son_data, cc_merge_data, left_index=True, right_index=True, how='left')
                    #         cc_merge_data.set_index('GO_name2', inplace=True, drop=True)
                    #         cc_son_data.set_index('GO_name2', inplace=True, drop=True)
                    #         cc_son_data_1 = cc_son_data.reset_index() 
                    #         parent = []
                    #         for k in list(cc_son_data_1.index): #k=1
                    #             lst = list(cc_son_data_1.loc[k])
                    #             lst = [x for x in lst if pd.notnull(x)]
                    #             last_value = lst[-1]
                    #             target_value = [] 
                    #             if last_value == 'is a+CELLULAR_COMPONENT':
                    #                 for item in lst:
                    #                     if 'part_of+' in item:
                    #                         term = item.split('+')[1]
                    #                         target_value.append(term)
                    #             parent.extend(target_value)
                    #         parent = [x for x in parent if x is not None]
                    #         parent = list(set(parent))
                    #     cc_both_data[j] = parent
                    # cc_both_data_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in cc_both_data.items() ]))
                    # cc_both_data_df = cc_both_data_df.T
                    cc_both_data_df = pd.DataFrame(function_data[function_data['Gene_set'].str.contains('GO:CC', na=False)]['GO_name'])
                    
            # data processing
            dataframes = [bp_both_data_df, mf_both_data_df, cc_both_data_df]  
            prefixes = ['bp', 'mf', 'cc']
            
            for both_data_df, prefix in zip(dataframes, prefixes): 
                if both_data_df is None or both_data_df.empty:
                    print(f"{prefix}_both_data_df is empty, skipping...")
                    continue
                #
                # both_data_df = both_data_df[~both_data_df[0].isnull()]
                # all_values = both_data_df.stack().unique()
                # both_data = pd.DataFrame(columns=all_values)
                # df_long = both_data_df.stack().reset_index()
                # df_long.columns = ['pathway', 'attribute_num', 'attribute']
                # attribute_to_pathways = df_long.groupby('attribute')['pathway'].apply(list)
                # attribute_to_pathways_df = attribute_to_pathways.reset_index()
                # attribute_to_pathways_df.columns = ['attribute', 'pathways']
                # max_len = attribute_to_pathways_df['pathways'].apply(len).max()
                # function_structure = pd.DataFrame(index=range(max_len))
                # for index, row in attribute_to_pathways_df.iterrows():
                #     attribute = row['attribute']
                #     pathways = row['pathways']
                #     pathways += [pd.NA] * (max_len - len(pathways))
                #     function_structure[attribute] = pathways
                # function_structure = function_structure.fillna(np.nan)
                function_structure = both_data_df.T
                
                gene_dict = {}
                for column, series in function_structure.items():
                    genes_list = function_data[function_data.index.isin(list(series.dropna()))]['Genes'].tolist()
                    all_genes = []
                    for item in genes_list:
                        if ';' in item:
                            genes = [gene.strip() for gene in item.split(';')]
                            all_genes.extend(genes)
                        else:
                            all_genes.append(item.strip())
                    unique_genes = list(set(all_genes))
                    gene_dict[column] = unique_genes
                    
                if gene_dict:
                    max_length = max(len(genes) for genes in gene_dict.values())
                else:
                    max_length = 0  # 或者根据需要设置一个默认值
                    
                gene_list = pd.DataFrame(index=range(max_length))
                for column, genes in gene_dict.items():
                    genes_extended = genes + [pd.NA] * (max_length - len(genes))
                    gene_list[column] = genes_extended
                # structure_statistics  
                structure_data = copy.deepcopy(self.fc_data)
                structure_data['GeneName'] = structure_data['GeneName'].str.upper()
                structure_data = structure_data.drop(columns='PeptideSequence+structure_coding+ProteinID')
                structure_data['GeneName'] = structure_data['GeneName'].str.split(';')
                structure_data = structure_data.explode('GeneName')
                
                #
                core_structure_result = pd.DataFrame(index=['A2B2C1D1dD1','A2B2C1D1dD2dD1','A2B2C1D1dD1dcbB5','A2B2C1D1dD2dD1dcbB5'])
                core_structure_count = pd.DataFrame(index=['A2B2C1D1dD1','A2B2C1D1dD2dD1','A2B2C1D1dD1dcbB5','A2B2C1D1dD2dD1dcbB5'])
                glycan_type_result = pd.DataFrame(index=structure_data['Glycan_type'].unique())
                glycan_type_count = pd.DataFrame(index=structure_data['Glycan_type'].unique())
                
                unique_branches = set()
                for branches_str in structure_data['Branches'].unique():
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
                branches_structure_result = pd.DataFrame(index=unique_branches_list) 
                branches_structure_count = pd.DataFrame(index=unique_branches_list)
                
                branches_count_result = pd.DataFrame(index=structure_data['BranchNumber'].unique())
                branches_count_count = pd.DataFrame(index=structure_data['BranchNumber'].unique())
                sialicacid_count_result = pd.DataFrame(index=structure_data['structure_coding'].str.count('3').value_counts().index)
                sialicacid_count_count = pd.DataFrame(index=structure_data['structure_coding'].str.count('3').value_counts().index)
                fucose_count_result = pd.DataFrame(index=structure_data['structure_coding'].str.count('5').value_counts().index)
                fucose_count_count = pd.DataFrame(index=structure_data['structure_coding'].str.count('5').value_counts().index)
                sialicacid_structure_result = pd.DataFrame(index=structure_data[structure_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().index)
                sialicacid_structure_count = pd.DataFrame(index=structure_data[structure_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().index)
                fucose_structure_result = pd.DataFrame(index=structure_data[structure_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().index)
                fucose_structure_count = pd.DataFrame(index=structure_data[structure_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().index)
                lacdinac_result = pd.DataFrame(index=[item for item in list(structure_data['lacdinac'].str.split(', ').explode().unique()) if item != ' '])
                lacdinac_count = pd.DataFrame(index=[item for item in list(structure_data['lacdinac'].str.split(', ').explode().unique()) if item != ' '])
                structurecoding_result = pd.DataFrame(index=structure_data['structure_coding'].unique())
                structurecoding_count = pd.DataFrame(index=structure_data['structure_coding'].unique())
                fucosylated_type_result = pd.DataFrame(index=structure_data['fucosylated type'].unique())
                fucosylated_type_result = fucosylated_type_result[fucosylated_type_result.index!=' ']
                fucosylated_type_count = pd.DataFrame(index=structure_data['fucosylated type'].unique())
                fucosylated_type_count = fucosylated_type_count[fucosylated_type_count.index!=' ']
                acgc_result = pd.DataFrame(index=structure_data['Ac/Gc'].unique())
                acgc_result = acgc_result[acgc_result.index!=' ']
                acgc_count = pd.DataFrame(index=structure_data['Ac/Gc'].unique())
                acgc_count = acgc_count[acgc_count.index!=' ']
                
                self.structure_data = structure_data
                self.gene_list = gene_list
                for column, series in gene_list.items():
                    pathway_structure_data = structure_data[structure_data['GeneName'].isin(list(series.dropna()))]
                    # core_structure
                    core_structure = self.identify_core_structure(pathway_structure_data)['Core_structure'].value_counts().reset_index()
                    core_structure.columns = ['Core_structure', column]
                    core_structure = core_structure[(core_structure['Core_structure']=='A2B2C1D1dD1dcbB5')|
                                                    (core_structure['Core_structure']=='A2B2C1D1dD2dD1dcbB5')|
                                                    (core_structure['Core_structure']=='A2B2C1D1dD1')|
                                                    (core_structure['Core_structure']=='A2B2C1D1dD2dD1')]
                    core_structure_c = core_structure.set_index('Core_structure', drop=False)
                    core_structure_count = pd.concat([core_structure_count, core_structure_c[column]], axis=1)
                    setattr(self, f'{prefix}_core_structure_count', core_structure_count.reset_index().rename(columns={'index':'Core_structure'}))
                    core_structure[column] = core_structure[column] / core_structure[column].sum()
                    core_structure.set_index('Core_structure', inplace=True, drop=False)
                    core_structure_result = pd.concat([core_structure_result, core_structure[column]], axis=1)
                    setattr(self, f'{prefix}_core_structure', core_structure_result.reset_index().rename(columns={'index':'Core_structure'}))
                    # glycan_type
                    glycan_type = pd.DataFrame(pathway_structure_data['Glycan_type'].value_counts().reset_index())
                    glycan_type.columns = ['Glycan_type', column]
                    glycan_type_c = glycan_type.set_index('Glycan_type', drop=False)
                    glycan_type_count = pd.concat([glycan_type_count, glycan_type_c[column]], axis=1)
                    setattr(self, f'{prefix}_glycan_type_count', glycan_type_count.reset_index().rename(columns={'index':'Glycan_type'}))
                    glycan_type[column] = glycan_type[column] / glycan_type[column].sum()
                    glycan_type.set_index('Glycan_type', inplace=True, drop=False)
                    glycan_type_result = pd.concat([glycan_type_result, glycan_type[column]], axis=1)
                    setattr(self, f'{prefix}_glycan_type', glycan_type_result.reset_index().rename(columns={'index':'Glycan_type'}))
                    # branches_structure
                    unique_branches = set()
                    for branches_str in pathway_structure_data['Branches'].value_counts().index:
                        branches_list = literal_eval(branches_str)
                        for branch in branches_list:
                            if branch:  
                                unique_branches.add(branch)
                    unique_branches_list = list(unique_branches)
                    resultlist = {}
                    for branch in unique_branches_list:
                        filtered_data = pathway_structure_data[pathway_structure_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                        item = filtered_data.shape[0]
                        resultlist[branch] = item
                    if not resultlist:
                        branches_structure = pd.DataFrame(columns = ['Branches', column])
                    if resultlist:
                        branches_structure = pd.DataFrame(index=resultlist.keys(),data=resultlist.values()).reset_index()
                        branches_structure.columns = ['Branches', column]
                    branches_structure_c = branches_structure.set_index('Branches', drop=False)
                    branches_structure_count = pd.concat([branches_structure_count, branches_structure_c[column]], axis=1)
                    setattr(self, f'{prefix}_branches_structure_count', branches_structure_count.reset_index().rename(columns={'index':'Branches'}))
                    branches_structure[column] = branches_structure[column] / branches_structure[column].sum()
                    branches_structure.set_index('Branches', inplace=True, drop=False)
                    branches_structure_result = pd.concat([branches_structure_result, branches_structure[column]], axis=1)
                    setattr(self, f'{prefix}_branches_structure', branches_structure_result.reset_index().rename(columns={'index':'Branches'}))
                    # branches_count
                    branches_count = pd.DataFrame(pathway_structure_data['BranchNumber'].value_counts().reset_index())
                    branches_count.columns = ['BranchNumber', column]
                    branches_count_c = branches_count.set_index('BranchNumber', drop=False)
                    branches_count_count = pd.concat([branches_count_count, branches_count_c[column]], axis=1)
                    setattr(self, f'{prefix}_branches_count_count', branches_count_count.reset_index().rename(columns={'index':'BranchNumber'}))
                    branches_count[column] = branches_count[column] / branches_count[column].sum()
                    branches_count.set_index('BranchNumber', inplace=True, drop=False)
                    branches_count_result = pd.concat([branches_count_result, branches_count[column]], axis=1)
                    setattr(self, f'{prefix}_branches_count', branches_count_result.reset_index().rename(columns={'index':'BranchNumber'}))
                    # sialicacid_count
                    sialicacid_count = pd.DataFrame(pathway_structure_data['structure_coding'].str.count('3').value_counts().reset_index())
                    sialicacid_count.columns = ['Sialicacid_count', column]
                    sialicacid_count_c = sialicacid_count.set_index('Sialicacid_count', drop=False)
                    sialicacid_count_count = pd.concat([sialicacid_count_count, sialicacid_count_c[column]], axis=1)
                    setattr(self, f'{prefix}_sialicacid_count_count', sialicacid_count_count.reset_index().rename(columns={'index':'Sialicacid_count'}))
                    sialicacid_count[column] = sialicacid_count[column] / sialicacid_count[column].sum()
                    sialicacid_count.set_index('Sialicacid_count', inplace=True, drop=False)
                    sialicacid_count_result = pd.concat([sialicacid_count_result, sialicacid_count[column]], axis=1)
                    setattr(self, f'{prefix}_sialicacid_count', sialicacid_count_result.reset_index().rename(columns={'index':'Sialicacid_count'}))
                    # fucose_count
                    fucose_count = pd.DataFrame(pathway_structure_data['structure_coding'].str.count('5').value_counts().reset_index())
                    fucose_count.columns = ['Fucose_count', column]
                    fucose_count_c = fucose_count.set_index('Fucose_count', drop=False)
                    fucose_count_count = pd.concat([fucose_count_count, fucose_count_c[column]], axis=1)
                    setattr(self, f'{prefix}_fucose_count_count', fucose_count_count.reset_index().rename(columns={'index':'Fucose_count'}))
                    fucose_count[column] = fucose_count[column] / fucose_count[column].sum()
                    fucose_count.set_index('Fucose_count', inplace=True, drop=False)
                    fucose_count_result = pd.concat([fucose_count_result, fucose_count[column]], axis=1)
                    setattr(self, f'{prefix}_fucose_count', fucose_count_result.reset_index().rename(columns={'index':'Fucose_count'}))
                    # sialicacid_structure
                    sialicacid_structure = pd.DataFrame(pathway_structure_data[pathway_structure_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().reset_index())
                    sialicacid_structure.columns = ['Sialicacid_structure', column]
                    sialicacid_structure_c = sialicacid_structure.set_index('Sialicacid_structure', drop=False)
                    sialicacid_structure_count = pd.concat([sialicacid_structure_count, sialicacid_structure_c[column]], axis=1)
                    setattr(self, f'{prefix}_sialicacid_structure_count', sialicacid_structure_count.reset_index().rename(columns={'index':'Sialicacid_structure'}))
                    sialicacid_structure[column] = sialicacid_structure[column] / sialicacid_structure[column].sum()
                    sialicacid_structure.set_index('Sialicacid_structure', inplace=True, drop=False)
                    sialicacid_structure_result = pd.concat([sialicacid_structure_result, sialicacid_structure[column]], axis=1)
                    setattr(self, f'{prefix}_sialicacid_structure', sialicacid_structure_result.reset_index().rename(columns={'index':'Sialicacid_structure'}))
                    # fucose_structure
                    fucose_structure = pd.DataFrame(pathway_structure_data[pathway_structure_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().reset_index())
                    fucose_structure.columns = ['Fucose_structure', column]
                    fucose_structure_c = fucose_structure.set_index('Fucose_structure', drop=False)
                    fucose_structure_count = pd.concat([fucose_structure_count, fucose_structure_c[column]], axis=1)
                    setattr(self, f'{prefix}_fucose_structure_count', fucose_structure_count.reset_index().rename(columns={'index':'Fucose_structure'}))
                    fucose_structure[column] = fucose_structure[column] / fucose_structure[column].sum()
                    fucose_structure.set_index('Fucose_structure', inplace=True, drop=False)
                    fucose_structure_result = pd.concat([fucose_structure_result, fucose_structure[column]], axis=1)
                    setattr(self, f'{prefix}_fucose_structure', fucose_structure_result.reset_index().rename(columns={'index':'Fucose_structure'}))
                    # lacdinac
                    lacdinac = pathway_structure_data.assign(lacdinac=pathway_structure_data['lacdinac'].str.split(', ')).explode('lacdinac')
                    lacdinac = lacdinac.drop_duplicates()
                    lacdinac = lacdinac['lacdinac'].value_counts().reset_index()
                    lacdinac = lacdinac[lacdinac['index']!=' ']
                    lacdinac.columns = ['Lacdinac', column]
                    lacdinac_c = lacdinac.set_index('Lacdinac', drop=False)
                    lacdinac_count = pd.concat([lacdinac_count, lacdinac_c[column]], axis=1)
                    setattr(self, f'{prefix}_lacdinac_count', lacdinac_count.reset_index().rename(columns={'index':'Lacdinac'}))
                    lacdinac[column] = lacdinac[column] / lacdinac[column].sum()
                    lacdinac.set_index('Lacdinac', inplace=True, drop=False)
                    lacdinac_result = pd.concat([lacdinac_result, lacdinac[column]], axis=1)
                    setattr(self, f'{prefix}_lacdinac', lacdinac_result.reset_index().rename(columns={'index':'Lacdinac'}))
                    # structure_coding
                    structurecoding = pd.DataFrame(pathway_structure_data['structure_coding'].value_counts().reset_index())
                    structurecoding.columns = ['Structure_coding', column]
                    structurecoding_c = structurecoding.set_index('Structure_coding', drop=False)
                    structurecoding_count = pd.concat([structurecoding_count, structurecoding_c[column]], axis=1)
                    setattr(self, f'{prefix}_structurecoding_count', structurecoding_count.reset_index().rename(columns={'index':'Structure_coding'}))
                    structurecoding[column] = structurecoding[column] / structurecoding[column].sum()
                    structurecoding.set_index('Structure_coding', inplace=True, drop=False)
                    structurecoding_result = pd.concat([structurecoding_result, structurecoding[column]], axis=1)
                    setattr(self, f'{prefix}_structurecoding', structurecoding_result.reset_index().rename(columns={'index':'Structure_coding'}))
                    # fucosylated_type
                    fucosylated_type = pd.DataFrame(pathway_structure_data['fucosylated type'].value_counts().reset_index())
                    fucosylated_type.columns = ['Fucosylated_type', column]
                    fucosylated_type = fucosylated_type[fucosylated_type['Fucosylated_type']!=' ']
                    fucosylated_type_c = fucosylated_type.set_index('Fucosylated_type', drop=False)
                    fucosylated_type_count = pd.concat([fucosylated_type_count, fucosylated_type_c[column]], axis=1)
                    setattr(self, f'{prefix}_fucosylated_type_count', fucosylated_type_count.reset_index().rename(columns={'index':'Fucosylated_type'}))
                    fucosylated_type[column] = fucosylated_type[column] / fucosylated_type[column].sum()
                    fucosylated_type.set_index('Fucosylated_type', inplace=True, drop=False)
                    fucosylated_type_result = pd.concat([fucosylated_type_result, fucosylated_type[column]], axis=1)
                    setattr(self, f'{prefix}_fucosylated_type', fucosylated_type_result.reset_index().rename(columns={'index':'Fucosylated_type'}))
                    # acgc
                    acgc = pd.DataFrame(pathway_structure_data['Ac/Gc'].value_counts().reset_index())
                    acgc.columns = ['Ac/Gc', column]
                    acgc = acgc[acgc['Ac/Gc']!=' ']
                    acgc_c = acgc.set_index('Ac/Gc', drop=False)
                    acgc_count = pd.concat([acgc_count, acgc_c[column]], axis=1)
                    setattr(self, f'{prefix}_acgc_count', acgc_count.reset_index().rename(columns={'index':'Ac/Gc'}))
                    acgc[column] = acgc[column] / acgc[column].sum()
                    acgc.set_index('Ac/Gc', inplace=True, drop=False)
                    acgc_result = pd.concat([acgc_result, acgc[column]], axis=1)
                    setattr(self, f'{prefix}_acgc', acgc_result.reset_index().rename(columns={'index':'Ac/Gc'}))

        else:
            raise ValueError('Please enter the correct data name')
        #
        self.data_manager.log_params('StrucGAP_FunctionAnnotation', 'go_function_structure', {'function_data':self.function_data_name})

        # 将所有数据存储在字典中
        output_data = {
            'go_function_data': self.function_data,
            'go_structure_data': self.structure_data,
            'bp_core_structure': self.bp_core_structure,
            'mf_core_structure': self.mf_core_structure,
            'cc_core_structure': self.cc_core_structure,
            'bp_glycan_type': self.bp_glycan_type,
            'mf_glycan_type': self.mf_glycan_type,
            'cc_glycan_type': self.cc_glycan_type,
            'bp_branches_structure': self.bp_branches_structure,
            'mf_branches_structure': self.mf_branches_structure,
            'cc_branches_structure': self.cc_branches_structure,
            'bp_branches_count': self.bp_branches_count,
            'mf_branches_count': self.mf_branches_count,
            'cc_branches_count': self.cc_branches_count,
            'bp_sialicacid_count': self.bp_sialicacid_count,
            'mf_sialicacid_count': self.mf_sialicacid_count,
            'cc_sialicacid_count': self.cc_sialicacid_count,
            'bp_fucose_count': self.bp_fucose_count,
            'mf_fucose_count': self.mf_fucose_count,
            'cc_fucose_count': self.cc_fucose_count,
            'bp_sialicacid_structure': self.bp_sialicacid_structure,
            'mf_sialicacid_structure': self.mf_sialicacid_structure,
            'cc_sialicacid_structure': self.cc_sialicacid_structure,
            'bp_fucose_structure': self.bp_fucose_structure,
            'mf_fucose_structure': self.mf_fucose_structure,
            'cc_fucose_structure': self.cc_fucose_structure,
            'bp_lacdinac': self.bp_lacdinac,
            'mf_lacdinac': self.mf_lacdinac,
            'cc_lacdinac': self.cc_lacdinac,
            'bp_fucosylated_type': self.bp_fucosylated_type,
            'mf_fucosylated_type': self.mf_fucosylated_type,
            'cc_fucosylated_type': self.cc_fucosylated_type,
            'bp_acgc': self.bp_acgc,
            'mf_acgc': self.mf_acgc,
            'cc_acgc': self.cc_acgc,
        }
        
        # 循环遍历字典中的每个项，并检查其值是否为空
        for key, value in output_data.items():
            if value is not None and not value.empty:  # 确保数据不为空且非空DataFrame
                self.data_manager.log_output('StrucGAP_FunctionAnnotation', key, value)

        return self
    
    def kegg_function_structure(self, function_data = None, p_value='P-value', cutoff=0.05):
        """
        Performs integrated analysis between enriched KEGG terms and glycan substructural features.
        
        Parameters:
            function_data: enrichment result in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result'].
            p_value: ['P-value', 'Adjusted P-value'] used to screen enrichment terms.
            cutoff: p_value threshold.
        
        Returns:
            self.kegg_core_structure
            self.kegg_glycan_type
            self.kegg_branches_structure
            self.kegg_branches_count
            self.kegg_sialicacid_count
            self.kegg_fucose_count
            self.kegg_sialicacid_structure
            self.kegg_fucose_structure
            self.kegg_lacdinac
            self.kegg_structurecoding
            self.kegg_fucosylated_type
            self.kegg_acgc
                
        Return type:
            dataframe
        
        """
        self.kegg_core_structure = pd.DataFrame()
        self.kegg_glycan_type = pd.DataFrame()
        self.kegg_branches_structure = pd.DataFrame()
        self.kegg_branches_count = pd.DataFrame()
        self.kegg_sialicacid_count = pd.DataFrame()
        self.kegg_fucose_count = pd.DataFrame()
        self.kegg_sialicacid_structure = pd.DataFrame()
        self.kegg_fucose_structure = pd.DataFrame()
        self.kegg_lacdinac = pd.DataFrame()
        self.kegg_structurecoding = pd.DataFrame()
        self.kegg_fucosylated_type = pd.DataFrame()
        self.kegg_acgc = pd.DataFrame()
        
        if function_data == None:
            function_data = input("Please enter the data you want to execute function and structure correlation analysis (select from ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']): ")
            expected_options = ['ora_no_background_up_result', 'ora_no_background_down_result',
                                'ora_background_up_result', 'ora_background_up_result', 'gsea_result',
                                'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']
            matches = get_close_matches(function_data, expected_options, n=1, cutoff=0.5)
            if matches:
                function_data = matches[0]
                print(f"Using '{function_data}' as the input.")
            else:
                print("No close match found. Using 'ora_no_background_up_result' as the input.")
                function_data = 'ora_no_background_up_result'
            
            if function_data not in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']:
                raise ValueError('Please enter the correct data name')
                
        self.function_data_name = function_data

        # function_data = data1
        if function_data in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'gsea_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']:
            if function_data in ['ora_no_background_up_result', 'ora_no_background_down_result', 'ora_background_up_result', 'ora_background_up_result', 'ora_no_background_both_proteins_result', 'ora_background_both_proteins_result']:
                function_data = getattr(self, function_data, None)
                function_data['Term'] = function_data['Term'].str.upper()
                function_data = function_data.set_index('Term',drop=False)
                if p_value is None:
                    p_value = input('Please enter the type of p_value (select from: P-value, Adjusted P-value) you would like to use to select the differential pathway: ')
                    expected_options = ['P-value', 'Adjusted P-value']
                    matches = get_close_matches(p_value, expected_options, n=1, cutoff=0.5)
                    if matches:
                        p_value = matches[0]
                        print(f"Using '{p_value}' as the input.")
                    else:
                        print("No close match found. Using 'P-value' as the input.")
                        p_value = 'P-value'
                if cutoff is None:
                    cutoff = input('Please enter the cutoff value (such as: 0.05, 0.01) you would like to use to select the differential pathway: ')
                function_data = function_data[function_data[p_value]<float(cutoff)]
            elif function_data == 'gsea_result':
                function_data = getattr(self, function_data, None)
                function_data['Term'] = [x[1].upper() for x in function_data['Term'].str.split('__')]
                function_data = function_data.set_index('Term',drop=False)
                if p_value is None:
                    p_value = input('Please enter the type of p_value (select from: NOM p-val, FDR q-val) you would like to use to select the differential pathway: ')
                    expected_options = ['NOM p-val', 'FDR q-val']
                    matches = get_close_matches(p_value, expected_options, n=1, cutoff=0.5)
                    if matches:
                        p_value = matches[0]
                        print(f"Using '{p_value}' as the input.")
                    else:
                        print("No close match found. Using 'P-value' as the input.")
                        p_value = 'P-value'
                if cutoff is None:
                    cutoff = input('Please enter the cutoff value (such as: 0.05, 0.01) you would like to use to select the differential pathway: ')
                function_data = function_data[function_data[p_value]<float(cutoff)]
                function_data = function_data.rename(columns={'Lead_genes': 'Genes'})
                
            self.function_data = function_data
            # data processing
            gene_dict = {}
            for i in function_data.index:
                genes_list = pd.DataFrame(function_data.loc[i]).T['Genes'].tolist()
                all_genes = []    
                for item in genes_list:
                    if ';' in item:
                        genes = [gene.strip() for gene in item.split(';')]
                        all_genes.extend(genes)
                    else:
                        all_genes.append(item.strip())
                unique_genes = list(set(all_genes))
                gene_dict[i] = unique_genes
                
            if gene_dict:
                max_length = max(len(genes) for genes in gene_dict.values())
            else:
                max_length = 0  # 或者根据需要设置一个默认值

            gene_list = pd.DataFrame(index=range(max_length))
            for column, genes in gene_dict.items():
                genes_extended = genes + [pd.NA] * (max_length - len(genes))
                gene_list[column] = genes_extended
            # structure_statistics
            structure_data = copy.deepcopy(module5.fc_data)
            structure_data['GeneName'] = structure_data['GeneName'].str.upper()
            structure_data = structure_data.drop(columns='PeptideSequence+structure_coding+ProteinID')
            structure_data['GeneName'] = structure_data['GeneName'].str.split(';')
            structure_data = structure_data.explode('GeneName')
            
            #
            core_structure_result = pd.DataFrame(index=['A2B2C1D1dD1','A2B2C1D1dD2dD1','A2B2C1D1dD1dcbB5','A2B2C1D1dD2dD1dcbB5'])
            glycan_type_result = pd.DataFrame(index=structure_data['Glycan_type'].unique())
            
            unique_branches = set()
            for branches_str in structure_data['Branches'].unique():
                branches_list = literal_eval(branches_str)
                for branch in branches_list:
                    if branch:  
                        unique_branches.add(branch)
            unique_branches_list = list(unique_branches)
            branches_structure_result = pd.DataFrame(index=unique_branches_list)   

            branches_count_result = pd.DataFrame(index=structure_data['BranchNumber'].unique())
            sialicacid_count_result = pd.DataFrame(index=structure_data['structure_coding'].str.count('3').value_counts().index)
            fucose_count_result = pd.DataFrame(index=structure_data['structure_coding'].str.count('5').value_counts().index)
            sialicacid_structure_result = pd.DataFrame(index=structure_data[structure_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().index)
            fucose_structure_result = pd.DataFrame(index=structure_data[structure_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().index)
            lacdinac_result = pd.DataFrame(index=[item for item in list(structure_data['lacdinac'].str.split(', ').explode().unique()) if item != ' '])
            structurecoding_result = pd.DataFrame(index=structure_data['structure_coding'].unique())
            fucosylated_type_result = pd.DataFrame(index=structure_data['fucosylated type'].unique())
            fucosylated_type_result = fucosylated_type_result[fucosylated_type_result.index!=' ']
            acgc_result = pd.DataFrame(index=structure_data['Ac/Gc'].unique())
            acgc_result = acgc_result[acgc_result.index!=' ']

            self.structure_data = structure_data
            for column, series in gene_list.items():
                pathway_structure_data = structure_data[structure_data['GeneName'].isin(list(series.dropna()))]
                # core_structure
                core_structure = self.identify_core_structure(pathway_structure_data)['Core_structure'].value_counts().reset_index()
                core_structure.columns = ['Core_structure', column]
                core_structure = core_structure[(core_structure['Core_structure']=='A2B2C1D1dD1dcbB5')|
                                                          (core_structure['Core_structure']=='A2B2C1D1dD2dD1dcbB5')|
                                                          (core_structure['Core_structure']=='A2B2C1D1dD1')|
                                                          (core_structure['Core_structure']=='A2B2C1D1dD2dD1')]
                core_structure[column] = core_structure[column] / core_structure[column].sum()
                core_structure.set_index('Core_structure', inplace=True, drop=False)
                core_structure_result = pd.concat([core_structure_result, core_structure[column]], axis=1)
                setattr(self, 'kegg_core_structure', core_structure_result.reset_index().rename(columns={'index':'Core_structure'}))
                # glycan_type
                glycan_type = pd.DataFrame(pathway_structure_data['Glycan_type'].value_counts().reset_index())
                glycan_type.columns = ['Glycan_type', column]
                glycan_type[column] = glycan_type[column] / glycan_type[column].sum()
                glycan_type.set_index('Glycan_type', inplace=True, drop=False)
                glycan_type_result = pd.concat([glycan_type_result, glycan_type[column]], axis=1)
                setattr(self, 'kegg_glycan_type', glycan_type_result.reset_index().rename(columns={'index':'Glycan_type'}))
                # branches_structure
                unique_branches = set()
                for branches_str in pathway_structure_data['Branches'].value_counts().index:
                    branches_list = literal_eval(branches_str)
                    for branch in branches_list:
                        if branch:  
                            unique_branches.add(branch)
                unique_branches_list = list(unique_branches)
            
                resultlist = {}
                for branch in unique_branches_list:
                    filtered_data = pathway_structure_data[pathway_structure_data['Branches'].apply(lambda x: branch in literal_eval(x))]
                    item = filtered_data.shape[0]
                    resultlist[branch] = item
                if not resultlist:
                    branches_structure = pd.DataFrame(columns = ['Branches', column])
                if resultlist:
                    branches_structure = pd.DataFrame(index=resultlist.keys(),data=resultlist.values()).reset_index()
                    branches_structure.columns = ['Branches', column]
                branches_structure[column] = branches_structure[column] / branches_structure[column].sum()
                branches_structure.set_index('Branches', inplace=True, drop=False)
                branches_structure_result = pd.concat([branches_structure_result, branches_structure[column]], axis=1)
                setattr(self, 'kegg_branches_structure', branches_structure_result.reset_index().rename(columns={'index':'Branches'}))
                # branches_count
                branches_count = pd.DataFrame(pathway_structure_data['BranchNumber'].value_counts().reset_index())
                branches_count.columns = ['BranchNumber', column]
                branches_count[column] = branches_count[column] / branches_count[column].sum()
                branches_count.set_index('BranchNumber', inplace=True, drop=False)
                branches_count_result = pd.concat([branches_count_result, branches_count[column]], axis=1)
                setattr(self, 'kegg_branches_count', branches_count_result.reset_index().rename(columns={'index':'BranchNumber'}))
                # sialicacid_count
                sialicacid_count = pd.DataFrame(pathway_structure_data['structure_coding'].str.count('3').value_counts().reset_index())
                sialicacid_count.columns = ['Sialicacid_count', column]
                sialicacid_count[column] = sialicacid_count[column] / sialicacid_count[column].sum()
                sialicacid_count.set_index('Sialicacid_count', inplace=True, drop=False)
                sialicacid_count_result = pd.concat([sialicacid_count_result, sialicacid_count[column]], axis=1)
                setattr(self, 'kegg_sialicacid_count', sialicacid_count_result.reset_index().rename(columns={'index':'Sialicacid_count'}))
                # fucose_count
                fucose_count = pd.DataFrame(pathway_structure_data['structure_coding'].str.count('5').value_counts().reset_index())
                fucose_count.columns = ['Fucose_count', column]
                fucose_count[column] = fucose_count[column] / fucose_count[column].sum()
                fucose_count.set_index('Fucose_count', inplace=True, drop=False)
                fucose_count_result = pd.concat([fucose_count_result, fucose_count[column]], axis=1)
                setattr(self, 'kegg_fucose_count', fucose_count_result.reset_index().rename(columns={'index':'Fucose_count'}))
                # sialicacid_structure
                sialicacid_structure = pd.DataFrame(pathway_structure_data[pathway_structure_data['structure_coding'].str.contains('3')]['structure_coding'].value_counts().reset_index())
                sialicacid_structure.columns = ['Sialicacid_structure', column]
                sialicacid_structure[column] = sialicacid_structure[column] / sialicacid_structure[column].sum()
                sialicacid_structure.set_index('Sialicacid_structure', inplace=True, drop=False)
                sialicacid_structure_result = pd.concat([sialicacid_structure_result, sialicacid_structure[column]], axis=1)
                setattr(self, 'kegg_sialicacid_structure', sialicacid_structure_result.reset_index().rename(columns={'index':'Sialicacid_structure'}))
                # fucose_structure
                fucose_structure = pd.DataFrame(pathway_structure_data[pathway_structure_data['structure_coding'].str.contains('5')]['structure_coding'].value_counts().reset_index())
                fucose_structure.columns = ['Fucose_structure', column]
                fucose_structure[column] = fucose_structure[column] / fucose_structure[column].sum()
                fucose_structure.set_index('Fucose_structure', inplace=True, drop=False)
                fucose_structure_result = pd.concat([fucose_structure_result, fucose_structure[column]], axis=1)
                setattr(self, 'kegg_fucose_structure', fucose_structure_result.reset_index().rename(columns={'index':'Fucose_structure'}))
                # lacdinac
                lacdinac = pathway_structure_data.assign(lacdinac=pathway_structure_data['lacdinac'].str.split(', ')).explode('lacdinac')
                lacdinac = lacdinac.drop_duplicates()
                lacdinac = lacdinac['lacdinac'].value_counts().reset_index()
                lacdinac = lacdinac[lacdinac['index']!='']
                lacdinac.columns = ['Lacdinac', column]
                lacdinac[column] = lacdinac[column] / lacdinac[column].sum()
                lacdinac.set_index('Lacdinac', inplace=True, drop=False)
                lacdinac_result = pd.concat([lacdinac_result, lacdinac[column]], axis=1)
                setattr(self, f'kegg_lacdinac', lacdinac_result.reset_index().rename(columns={'index':'Lacdinac'}))
                # structure_coding
                structurecoding = pd.DataFrame(pathway_structure_data['structure_coding'].value_counts().reset_index())
                structurecoding.columns = ['Structure_coding', column]
                structurecoding[column] = structurecoding[column] / structurecoding[column].sum()
                structurecoding.set_index('Structure_coding', inplace=True, drop=False)
                structurecoding_result = pd.concat([structurecoding_result, structurecoding[column]], axis=1)
                setattr(self, f'kegg_structurecoding', structurecoding_result.reset_index().rename(columns={'index':'Structure_coding'}))
                # fucosylated_type
                fucosylated_type = pd.DataFrame(pathway_structure_data['fucosylated type'].value_counts().reset_index())
                fucosylated_type.columns = ['Fucosylated_type', column]
                fucosylated_type = fucosylated_type[fucosylated_type['Fucosylated_type']!=' ']
                fucosylated_type[column] = fucosylated_type[column] / fucosylated_type[column].sum()
                fucosylated_type.set_index('Fucosylated_type', inplace=True, drop=False)
                fucosylated_type_result = pd.concat([fucosylated_type_result, fucosylated_type[column]], axis=1)
                setattr(self, f'kegg_fucosylated_type', fucosylated_type_result.reset_index().rename(columns={'index':'Fucosylated_type'}))
                # acgc
                acgc = pd.DataFrame(pathway_structure_data['Ac/Gc'].value_counts().reset_index())
                acgc.columns = ['Ac/Gc', column]
                acgc = acgc[acgc['Ac/Gc']!=' ']
                acgc[column] = acgc[column] / acgc[column].sum()
                acgc.set_index('Ac/Gc', inplace=True, drop=False)
                acgc_result = pd.concat([acgc_result, acgc[column]], axis=1)
                setattr(self, f'kegg_acgc', acgc_result.reset_index().rename(columns={'index':'Ac/Gc'}))

        else:
            raise ValueError('Please enter the correct data name')
        #
        self.data_manager.log_params('StrucGAP_FunctionAnnotation', 'kegg_function_structure', {'function_data':function_data})

        # 将所有数据存储在字典中
        output_data = {
            'kegg_function_data': self.function_data,
            'kegg_structure_data': self.structure_data,
            'kegg_core_structure': self.kegg_core_structure,
            'kegg_glycan_type': self.kegg_glycan_type,
            'kegg_branches_structure': self.kegg_branches_structure,
            'kegg_branches_count': self.kegg_branches_count,
            'kegg_sialicacid_count': self.kegg_sialicacid_count,
            'kegg_fucose_count': self.kegg_fucose_count,
            'kegg_sialicacid_structure': self.kegg_sialicacid_structure,
            'kegg_fucose_structure': self.kegg_fucose_structure,
            'kegg_lacdinac': self.kegg_lacdinac,
            'kegg_fucosylated_type': self.kegg_fucosylated_type,
            'kegg_acgc': self.kegg_acgc,
        }
        
        # 循环遍历字典中的每个项，并检查其值是否为空
        for key, value in output_data.items():
            if value is not None and not value.empty:  # 确保数据不为空且非空DataFrame
                self.data_manager.log_output('StrucGAP_FunctionAnnotation', key, value)
        
        return self
    
    def output(self):
        """Outputs both analysis results."""
        output_dir = './analysis_result'
        os.makedirs(output_dir, exist_ok=True)
        
        term_name = 'Unknown'
        
        if self.selected_terms[0] in ['GO:BP','GO:CC','GO:MF']:
            term_name = 'GO'
        elif self.selected_terms[0] == 'KEGG':
            term_name = 'KEGG'
        
        with pd.ExcelWriter(os.path.join(output_dir, f'StrucGAP_FunctionAnnotation_{term_name}_{self.data_type}_{self.function_data_name}.xlsx'), engine='xlsxwriter') as writer:

            if self.ora_no_background_up_result is not None and not self.ora_no_background_up_result.empty:
                self.ora_no_background_up_result.to_excel(writer, sheet_name=f'{term_name}_ora_no_bg_up_result'[:31])
            if self.ora_no_background_down_result is not None and not self.ora_no_background_down_result.empty:
                self.ora_no_background_down_result.to_excel(writer, sheet_name=f'{term_name}_ora_no_bg_down_result'[:31])
            if self.ora_background_up_result is not None and not self.ora_background_up_result.empty:
                self.ora_background_up_result.to_excel(writer, sheet_name=f'{term_name}_ora_bg_up_result'[:31])
            if self.ora_background_down_result is not None and not self.ora_background_down_result.empty:
                self.ora_background_down_result.to_excel(writer, sheet_name=f'{term_name}_ora_bg_down_result'[:31])
            if self.ora_no_background_both_proteins_result is not None and not self.ora_no_background_both_proteins_result.empty:
                self.ora_no_background_both_proteins_result.to_excel(writer, sheet_name=f'{term_name}_ora_no_bg_both_proteins_result'[:31])
            if self.ora_background_both_proteins_result is not None and not self.ora_background_both_proteins_result.empty:
                self.ora_background_both_proteins_result.to_excel(writer, sheet_name=f'{term_name}_ora_bg_both_proteins_result'[:31])
            if self.gsea_result is not None and not self.gsea_result.empty:
                self.gsea_result.to_excel(writer, sheet_name=f'{term_name}_gsea_result'[:31])
                
            output_data = {
                'go_function_data': self.function_data,      
                'go_structure_data': self.structure_data,      
                'bp_core_structure_ratio': self.bp_core_structure if hasattr(self, 'bp_core_structure') else None,
                'bp_core_structure_count'[:31]: self.bp_core_structure_count if hasattr(self, 'bp_core_structure_count') else None,
                'mf_core_structure_ratio': self.mf_core_structure if hasattr(self, 'mf_core_structure') else None,
                'mf_core_structure_count'[:31]: self.mf_core_structure_count if hasattr(self, 'mf_core_structure_count') else None,
                'cc_core_structure_ratio': self.cc_core_structure if hasattr(self, 'cc_core_structure') else None,
                'cc_core_structure_count'[:31]: self.cc_core_structure_count if hasattr(self, 'cc_core_structure_count') else None,
                'bp_glycan_type_ratio': self.bp_glycan_type if hasattr(self, 'bp_glycan_type') else None,
                'bp_glycan_type_count'[:31]: self.bp_glycan_type_count if hasattr(self, 'bp_glycan_type_count') else None,
                'mf_glycan_type_ratio': self.mf_glycan_type if hasattr(self, 'mf_glycan_type') else None,
                'mf_glycan_type_count'[:31]: self.mf_glycan_type_count if hasattr(self, 'mf_glycan_type_count') else None,
                'cc_glycan_type_ratio': self.cc_glycan_type if hasattr(self, 'cc_glycan_type') else None,
                'cc_glycan_type_count'[:31]: self.cc_glycan_type_count if hasattr(self, 'cc_glycan_type_count') else None,
                'bp_branches_structure_ratio': self.bp_branches_structure if hasattr(self, 'bp_branches_structure') else None,
                'bp_branches_structure_count'[:31]: self.bp_branches_structure_count if hasattr(self, 'bp_branches_structure_count') else None,
                'mf_branches_structure_ratio': self.mf_branches_structure if hasattr(self, 'mf_branches_structure') else None,
                'mf_branches_structure_count'[:31]: self.mf_branches_structure_count if hasattr(self, 'mf_branches_structure_count') else None,
                'cc_branches_structure_ratio': self.cc_branches_structure if hasattr(self, 'cc_branches_structure') else None,
                'cc_branches_structure_count'[:31]: self.cc_branches_structure_count if hasattr(self, 'cc_branches_structure_count') else None,
                'bp_branches_count_ratio': self.bp_branches_count if hasattr(self, 'bp_branches_count') else None,
                'bp_branches_count_count'[:31]: self.bp_branches_count_count if hasattr(self, 'bp_branches_count_count') else None,
                'mf_branches_count_ratio': self.mf_branches_count if hasattr(self, 'mf_branches_count') else None,
                'mf_branches_count_count'[:31]: self.mf_branches_count_count if hasattr(self, 'mf_branches_count_count') else None,
                'cc_branches_count_ratio': self.cc_branches_count if hasattr(self, 'cc_branches_count') else None,
                'cc_branches_count_count'[:31]: self.cc_branches_count_count if hasattr(self, 'cc_branches_count_count') else None,
                'bp_sialicacid_count_ratio': self.bp_sialicacid_count if hasattr(self, 'bp_sialicacid_count') else None,
                'bp_sialicacid_count_count'[:31]: self.bp_sialicacid_count_count if hasattr(self, 'bp_sialicacid_count_count') else None,
                'mf_sialicacid_count_ratio': self.mf_sialicacid_count if hasattr(self, 'mf_sialicacid_count') else None,
                'mf_sialicacid_count_count'[:31]: self.mf_sialicacid_count_count if hasattr(self, 'mf_sialicacid_count_count') else None,
                'cc_sialicacid_count_ratio': self.cc_sialicacid_count if hasattr(self, 'cc_sialicacid_count') else None,
                'cc_sialicacid_count_count'[:31]: self.cc_sialicacid_count_count if hasattr(self, 'cc_sialicacid_count_count') else None,
                'bp_fucose_count_ratio': self.bp_fucose_count if hasattr(self, 'bp_fucose_count') else None,
                'bp_fucose_count_count'[:31]: self.bp_fucose_count_count if hasattr(self, 'bp_fucose_count_count') else None,
                'mf_fucose_count_ratio': self.mf_fucose_count if hasattr(self, 'mf_fucose_count') else None,
                'mf_fucose_count_count'[:31]: self.mf_fucose_count_count if hasattr(self, 'mf_fucose_count_count') else None,
                'cc_fucose_count_ratio': self.cc_fucose_count if hasattr(self, 'cc_fucose_count') else None,
                'cc_fucose_count_count'[:31]: self.cc_fucose_count_count if hasattr(self, 'cc_fucose_count_count') else None,
                'bp_sialicacid_structure_ratio': self.bp_sialicacid_structure if hasattr(self, 'bp_sialicacid_structure') else None,
                'bp_sialicacid_structure_count'[:31]: self.bp_sialicacid_structure_count if hasattr(self, 'bp_sialicacid_structure_count') else None,
                'mf_sialicacid_structure_ratio': self.mf_sialicacid_structure if hasattr(self, 'mf_sialicacid_structure') else None,
                'mf_sialicacid_structure_count'[:31]: self.mf_sialicacid_structure_count if hasattr(self, 'mf_sialicacid_structure_count') else None,
                'cc_sialicacid_structure_ratio': self.cc_sialicacid_structure if hasattr(self, 'cc_sialicacid_structure') else None,
                'cc_sialicacid_structure_count'[:31]: self.cc_sialicacid_structure_count if hasattr(self, 'cc_sialicacid_structure_count') else None,
                'bp_fucose_structure_ratio': self.bp_fucose_structure if hasattr(self, 'bp_fucose_structure') else None,
                'bp_fucose_structure_count'[:31]: self.bp_fucose_structure_count if hasattr(self, 'bp_fucose_structure_count') else None,
                'mf_fucose_structure_ratio': self.mf_fucose_structure if hasattr(self, 'mf_fucose_structure') else None,
                'mf_fucose_structure_count'[:31]: self.mf_fucose_structure_count if hasattr(self, 'mf_fucose_structure_count') else None,
                'cc_fucose_structure_ratio': self.cc_fucose_structure if hasattr(self, 'cc_fucose_structure') else None,
                'cc_fucose_structure_count'[:31]: self.cc_fucose_structure_count if hasattr(self, 'cc_fucose_structure_count') else None,
                'bp_lacdinac_ratio': self.bp_lacdinac if hasattr(self, 'bp_lacdinac') else None,
                'bp_lacdinac_count'[:31]: self.bp_lacdinac_count if hasattr(self, 'bp_lacdinac_count') else None,
                'mf_lacdinac_ratio': self.mf_lacdinac if hasattr(self, 'mf_lacdinac') else None,
                'mf_lacdinac_count'[:31]: self.mf_lacdinac_count if hasattr(self, 'mf_lacdinac_count') else None,
                'cc_lacdinac_ratio': self.cc_lacdinac if hasattr(self, 'cc_lacdinac') else None,
                'cc_lacdinac_count'[:31]: self.cc_lacdinac_count if hasattr(self, 'cc_lacdinac_count') else None,
                'bp_structurecoding_ratio': self.bp_structurecoding if hasattr(self, 'bp_structurecoding') else None,
                'bp_structurecoding_count'[:31]: self.bp_structurecoding_count if hasattr(self, 'bp_structurecoding_count') else None,
                'mf_structurecoding_ratio': self.mf_structurecoding if hasattr(self, 'mf_structurecoding') else None,
                'mf_structurecoding_count'[:31]: self.mf_structurecoding_count if hasattr(self, 'mf_structurecoding_count') else None,
                'cc_structurecoding_ratio': self.cc_structurecoding if hasattr(self, 'cc_structurecoding') else None,
                'cc_structurecoding_count'[:31]: self.cc_structurecoding_count if hasattr(self, 'cc_structurecoding_count') else None,
                'bp_fucosylated_type_ratio': self.bp_fucosylated_type if hasattr(self, 'bp_fucosylated_type') else None,
                'bp_fucosylated_type_count'[:31]: self.bp_fucosylated_type_count if hasattr(self, 'bp_fucosylated_type_count') else None,
                'mf_fucosylated_type_ratio': self.mf_fucosylated_type if hasattr(self, 'mf_fucosylated_type') else None,
                'mf_fucosylated_type_count'[:31]: self.mf_fucosylated_type_count if hasattr(self, 'mf_fucosylated_type_count') else None,
                'cc_fucosylated_type_ratio': self.cc_fucosylated_type if hasattr(self, 'cc_fucosylated_type') else None,
                'cc_fucosylated_type_count'[:31]: self.cc_fucosylated_type_count if hasattr(self, 'cc_fucosylated_type_count') else None,
                'bp_acgc_ratio': self.bp_acgc if hasattr(self, 'bp_acgc') else None,
                'bp_acgc_count'[:31]: self.bp_acgc_count if hasattr(self, 'bp_acgc_count') else None,
                'mf_acgc_ratio': self.mf_acgc if hasattr(self, 'mf_acgc') else None,
                'mf_acgc_count'[:31]: self.mf_acgc_count if hasattr(self, 'mf_acgc_count') else None,
                'cc_acgc_ratio': self.cc_acgc if hasattr(self, 'cc_acgc') else None,
                'cc_acgc_count'[:31]: self.cc_acgc_count if hasattr(self, 'cc_acgc_count') else None,
            }
            
            # 循环遍历字典中的每个项，并检查其值是否为空
            for key, value in output_data.items():
                try:
                    if value is not None and not value.empty: 
                        value.to_excel(writer, sheet_name=key)
                except AttributeError:
                    # 如果属性未定义，则跳过该项
                    continue
                    
            output_data = {
                'kegg_function_data': self.function_data,
                'kegg_structure_data': self.structure_data,
                'kegg_core_structure': self.kegg_core_structure if hasattr(self, 'kegg_core_structure') else None,
                'kegg_glycan_type': self.kegg_glycan_type if hasattr(self, 'kegg_glycan_type') else None,
                'kegg_branches_structure': self.kegg_branches_structure if hasattr(self, 'kegg_branches_structure') else None,
                'kegg_branches_count': self.kegg_branches_count if hasattr(self, 'kegg_branches_count') else None,
                'kegg_sialicacid_count': self.kegg_sialicacid_count if hasattr(self, 'kegg_sialicacid_count') else None,
                'kegg_fucose_count': self.kegg_fucose_count if hasattr(self, 'kegg_fucose_count') else None,
                'kegg_sialicacid_structure': self.kegg_sialicacid_structure if hasattr(self, 'kegg_sialicacid_structure') else None,
                'kegg_fucose_structure': self.kegg_fucose_structure if hasattr(self, 'kegg_fucose_structure') else None,
                'kegg_lacdinac': self.kegg_lacdinac if hasattr(self, 'kegg_lacdinac') else None,
                'kegg_structurecoding': self.kegg_structurecoding if hasattr(self, 'kegg_structurecoding') else None,
                'kegg_fucosylated_type': self.kegg_fucosylated_type if hasattr(self, 'kegg_fucosylated_type') else None,
                'kegg_acgc': self.kegg_acgc if hasattr(self, 'kegg_acgc') else None,
            }
            
            # 循环遍历字典中的每个项，并检查其值是否为空
            for key, value in output_data.items():
                try:
                    if value is not None and not value.empty: 
                        value.to_excel(writer, sheet_name=key)
                except AttributeError:
                    # 如果属性未定义，则跳过该项
                    continue
            

                      
