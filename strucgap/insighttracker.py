# -*- coding: utf-8 -*-
## 数据分析管理模块
class StrucGAP_InsightTracker:
    def __init__(self):
        self.module_records = {}
        self.analysis_params = {}
        self.outputs = {}
        
    def register_module(self, module_name, module_instance, params):
        """
        Registers a module instance, records the initial parameters, and output locations.

        Parameters:
            module_name (str): The name of the module being registered.
            module_instance (object): The instance of the module to be registered.
            params (dict): The initial parameters used to configure the module.
    
        Returns:
            None
            
        """
        self.module_records[module_name] = {'instance': module_instance, 'params': params, 'outputs': {}}
        self.analysis_params[module_name] = {}

    def log_params(self, module_name, function_name, params):
        """
        Logs the parameters for each function call.
    
        Parameters:
            module_name (str): The name of the module.
            function_name (str): The name of the function being logged.
            params (dict): The parameters passed to the function.
    
        Returns:
            None
            
        """
        if module_name not in self.analysis_params:
            self.analysis_params[module_name] = {}
        self.analysis_params[module_name][function_name] = params

    def log_output(self, module_name, output_name, output_data):
        """
        Logs the output data of a module.
    
        Parameters:
            module_name (str): The name of the module.
            output_name (str): The name of the output being logged.
            output_data (object): The output data produced by the module.
    
        Returns:
            None
            
        """
        if module_name in self.module_records:
            self.module_records[module_name]['outputs'][output_name] = output_data
            self.outputs[output_name] = output_data

    def retrieve_data(self, module_name, output_name):
        """
        Retrieves output data from a specified module for use in other modules.
    
        Parameters:
            module_name (str): The name of the module.
            output_name (str): The name of the output to retrieve.
    
        Returns:
            object or None: The output data if found, otherwise None.
            
        """
        return self.module_records[module_name]['outputs'].get(output_name, None)

    def show_params(self, module_name):
        """
        Displays the analysis parameters for a specified module.

        Parameters:
            module_name (str): The name of the module whose parameters are to be displayed.
    
        Returns:
            dict or None: The analysis parameters for the specified module, or None if not found.
        """
        return self.analysis_params.get(module_name, None)

    def get_all_data(self):
        """
        Returns all registered modules' parameters and output data.
    
        Parameters:
            None.
    
        Returns:
            dict: A dictionary containing all registered module parameters and output data.
            
        """
        return self.module_records
    
    def output_analysis_params(self, output_dir='./analysis_result', output_file='GAP_analysis_params.xlsx'):
        """
        Outputs the analysis parameters to an Excel file.
    
        Parameters:
            output_dir (str): The directory where the output file will be saved (default is './analysis_result').
            output_file (str): The name of the output Excel file (default is 'GAP_analysis_params.xlsx').
    
        Returns:
            None
            
        """
        os.makedirs(output_dir, exist_ok=True)  # 创建输出目录
        output_path = os.path.join(output_dir, output_file)

        # 使用ExcelWriter管理多个Sheet
        with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
            for module_name, module_params in self.analysis_params.items():
                rows = []
                for function_name, params in module_params.items():
                    if isinstance(params, dict) and len(params) != 0:  
                        for sub_key, sub_value in params.items():
                            rows.append([function_name, sub_key, sub_value])
                    else:  # 如果不是字典，直接存储
                        rows.append([function_name, np.nan, np.nan])
                # 创建一个DataFrame用于存储结果
                df = pd.DataFrame(rows, columns=["Function", "Parameter", "Value"])
                df.to_excel(writer, sheet_name=module_name, index=False)
                
    def key_information_extraction(self, module):
        """
        Mines outputs from both StrucGAP_GlycoPeptideQuant and StrucGAP_FunctionAnnotation to identify structurally and functionally relevant glycan substructural features.
        
        Parameters:
            module:  StrucGAP_GlycoPeptideQuant or StrucGAP_FunctionAnnotation.
        
        Returns:
            None (output key information directly as a table).
        
        """
        if module not in ['StrucGAP_GlycoPeptideQuant', 'StrucGAP_FunctionAnnotation']:
            print("Please select module in ['StrucGAP_GlycoPeptideQuant' or 'StrucGAP_FunctionAnnotation']!")
            return
        
        if module == 'StrucGAP_GlycoPeptideQuant' and 'StrucGAP_GlycoPeptideQuant' in data_manager.module_records.keys():
            StrucGAP_GlycoPeptideQuant_key_set = {}
            print('StrucGAP_GlycoPeptideQuant module has been executed, key information extracting ...')
            for i in ['core_structure','branches_structure','glycan_type','branches_count',
                      'glycan_composition','lacdinac','fucosylated_type','acgc']: # i='glycan_type'
                ratio_sheetname = f"result_{i}_ratio"
                up_sheetname = f"result_{i}_up"
                up_ratio_sheetname = f"result_{i}_up_ratio"
                down_sheetname = f"result_{i}_down"
                down_ratio_sheetname = f"result_{i}_down_ratio"
                #
                data_ratio = getattr(module4, ratio_sheetname)
                data_ratio = data_ratio.replace(0, np.nan)
                data_up = getattr(module4, up_sheetname)
                data_up = data_up.replace(0, np.nan)
                data_up_ratio = getattr(module4, up_ratio_sheetname)
                data_up_ratio = data_up_ratio.replace(0, np.nan)
                data_down = getattr(module4, down_sheetname)
                data_down = data_down.replace(0, np.nan)
                data_down_ratio = getattr(module4, down_ratio_sheetname)
                data_down_ratio = data_down_ratio.replace(0, np.nan)
                
                def is_monotonic(row):
                    values = row[1:].astype(float)  
                    increasing = all(values[i] <= values[i+1] for i in range(len(values)-1))
                    decreasing = all(values[i] >= values[i+1] for i in range(len(values)-1))
                    return increasing or decreasing
                
                output = pd.DataFrame()
                
                def stack_with_empty_row(final_df, new_df):
                    if not new_df.empty:
                        column_names = new_df.columns.tolist()
                        new_df.columns = [0,1,2,3,4,5]
                        new_df.loc[-1] = column_names  
                        new_df.index = new_df.index + 1 
                        new_df = new_df.sort_index()
                        # 堆叠数据
                        final_df = pd.concat([final_df, new_df], ignore_index=True)
                        # 插入空行
                        empty_row = pd.DataFrame({col: np.nan for col in new_df.columns}, index=[0])
                        final_df = pd.concat([final_df, empty_row], ignore_index=True)
                    return final_df
                
                # ratio
                df = pd.DataFrame(data_ratio)
                for col in df.columns[1:]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df_filtered = df.dropna()
                result = df_filtered[df_filtered.apply(is_monotonic, axis=1)]
                if not result.empty:
                    output = stack_with_empty_row(output, result)
                    
                # up
                df = pd.DataFrame(data_up_ratio)
                for col in df.columns[1:]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df_filtered = df.dropna()
                result = df_filtered[df_filtered.apply(is_monotonic, axis=1)]
                if not result.empty:
                    result1 = data_up[data_up[data_up.columns[0]].isin(result[data_up.columns[0]])]
                    result1['min_value'] = result1.iloc[:, 2:].min(axis=1)
                    result1 = result1[result1['min_value'] >= 3]
                    result1 = result1.drop(columns=['min_value'])
                    output = stack_with_empty_row(output, result1)
                    if not result1.empty:
                        result = result[result[result.columns[0]].isin(result1[result1.columns[0]])]
                        output = stack_with_empty_row(output, result)
                        
                # down
                df = pd.DataFrame(data_down_ratio)
                for col in df.columns[1:]:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                df_filtered = df.dropna()
                result = df_filtered[df_filtered.apply(is_monotonic, axis=1)]
                if not result.empty:
                    result1 = data_down[data_down[data_down.columns[0]].isin(result[data_down.columns[0]])]
                    result1['min_value'] = result1.iloc[:, 2:].min(axis=1)
                    result1 = result1[result1['min_value'] >= 3]
                    result1 = result1.drop(columns=['min_value'])
                    output = stack_with_empty_row(output, result1)
                    if not result1.empty:
                        result = result[result[result.columns[0]].isin(result1[result1.columns[0]])]
                        output = stack_with_empty_row(output, result)
                        
                StrucGAP_GlycoPeptideQuant_key_set[i] = output
                output_dir = './analysis_result'
                os.makedirs(output_dir, exist_ok=True)
                with pd.ExcelWriter(os.path.join(output_dir, 'StrucGAP_GlycoPeptideQuant_key_information.xlsx'), engine='xlsxwriter') as writer:
                    for sheet_name, df in StrucGAP_GlycoPeptideQuant_key_set.items():
                        df.iloc[:-1,:].to_excel(writer, sheet_name=sheet_name, index=True)
                        
        if module == 'StrucGAP_FunctionAnnotation' and 'StrucGAP_FunctionAnnotation' in data_manager.module_records.keys():
            StrucGAP_FunctionAnnotation_key_set = {}
            function_data_type = data_manager.analysis_params['StrucGAP_FunctionAnnotation']['go_function_structure']['function_data']
            database = data_manager.analysis_params['StrucGAP_FunctionAnnotation']['ora']['terms']
            if 'GO' in database[0]:
                if len(database) == 3:
                    database = 'GO'
                else:
                    database = database
            if 'KEGG' in database[0]:
                database = 'KEGG'
            
            if database == 'GO':
                term_list = ['bp_core_structure','mf_core_structure','cc_core_structure',
                      'bp_glycan_type','mf_glycan_type','cc_glycan_type',
                      'bp_branches_structure','mf_branches_structure','cc_branches_structure',
                      'bp_branches_count','mf_branches_count','cc_branches_count',
                      'bp_sialicacid_count','mf_sialicacid_count','cc_sialicacid_count',
                      'bp_fucose_count','mf_fucose_count','cc_fucose_count',
                      'bp_sialicacid_structure','mf_sialicacid_structure','cc_sialicacid_structure',
                      'bp_fucose_structure','mf_fucose_structure','cc_fucose_structure',
                      'bp_lacdinac','mf_lacdinac','cc_lacdinac',
                      'bp_structurecoding','mf_structurecoding','cc_structurecoding',
                      'bp_fucosylated_type','mf_fucosylated_type','cc_fucosylated_type',
                      'bp_acgc','mf_acgc','cc_acgc',
                      ]
            elif database == 'KEGG':
                term_list = ['kegg_core_structure',
                      'kegg_glycan_type',
                      'kegg_branches_structure',
                      'kegg_branches_count',
                      'kegg_sialicacid_count',
                      'kegg_fucose_count',
                      'kegg_sialicacid_structure',
                      'kegg_fucose_structure',
                      'kegg_lacdinac',
                      'kegg_structurecoding',
                      'kegg_fucosylated_type',
                      'kegg_acgc',
                      ]
            
            for i in term_list: # i='mf_branches_count'
                
                print(i)
                count_sheetname = f"{i}_count"
                ratio_sheetname = f"{i}"
                
                data_count = getattr(module5, count_sheetname)
                data_count = data_count.replace(0, np.nan)
                data_ratio = getattr(module5, ratio_sheetname)
                data_ratio = data_ratio.replace(0, np.nan)
                
                def extract_columns_with_dynamic_threshold(df, target_percentage=0.05):
                    num_columns = len(df.columns)-1
                    if num_columns <= 10:
                        target_percentage=0.5
                    target_columns = int(num_columns * target_percentage)
                    first_column = df.columns[0]
                    thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                    
                    best_threshold = None
                    best_column_count = float('inf')
                    
                    for threshold in thresholds:
                        selected_columns = [col for col in df.columns[1:] if df[col].max() > threshold]
                        if abs(len(selected_columns) - target_columns) < abs(best_column_count - target_columns):
                            best_threshold = threshold
                            best_column_count = len(selected_columns)
                    
                    selected_columns = [first_column] + [col for col in df.columns[1:] if df[col].max() > best_threshold]
                    
                    print(f"Using threshold {best_threshold}, selected {len(selected_columns)} columns out of {num_columns}")
                    return df[selected_columns]
                
                def filter_columns_by_min_value(df, min_value=3):
                    first_column = df.columns[0]
                    selected_columns = [first_column] + [col for col in df.columns[1:] if df[col].min() > min_value]
                    return df[selected_columns]
                
                output = pd.DataFrame()
                
                def stack_with_empty_row(final_df, new_df):
                    if not new_df.empty:
                        column_names = new_df.columns.tolist()
                        new_df.columns = list(range(0,len(new_df.columns)))
                        new_df.loc[-1] = column_names  
                        new_df.index = new_df.index + 1 
                        new_df = new_df.sort_index()
                        # 堆叠数据
                        final_df = pd.concat([final_df, new_df], ignore_index=True)
                        # 插入空行
                        empty_row = pd.DataFrame({col: np.nan for col in new_df.columns}, index=[0])
                        final_df = pd.concat([final_df, empty_row], ignore_index=True)
                    return final_df
                
                df = pd.DataFrame(data_ratio)
                columns_to_keep = df.iloc[:, 1:].columns[~df.iloc[:, 1:].isin([1]).any()]
                df = df[[df.columns[0]] + list(columns_to_keep)]
                result = extract_columns_with_dynamic_threshold(df)
                if not result.iloc[:,1:].empty:
                    result1 = data_count[result.columns]
                    result1 = filter_columns_by_min_value(result1)
                    output = stack_with_empty_row(output, result1)
                    if not result1.empty:
                        result = result[list(result1.loc[0])]
                        output = stack_with_empty_row(output, result)
                
                if not output.iloc[:,1:].empty:
                    StrucGAP_FunctionAnnotation_key_set[i] = output
                output_dir = './analysis_result'
                os.makedirs(output_dir, exist_ok=True)
                with pd.ExcelWriter(os.path.join(output_dir, f'StrucGAP_FunctionAnnotation_{database}_{function_data_type}_key_information.xlsx'), engine='xlsxwriter') as writer:
                    for sheet_name, df in StrucGAP_FunctionAnnotation_key_set.items():
                        df.iloc[:-1,:].to_excel(writer, sheet_name=sheet_name, index=True)
    
    def output_pickle(self):
        """
        Serializes and saves all picklable variables (excluding Flask-related objects) to a file.
    
        Returns:
            None
        """
        # 定义一个函数来检查对象是否可以被序列化
        def is_picklable(obj):
            try:
                pickle.dumps(obj)
            except (pickle.PicklingError, TypeError):
                return False
            return True
        # 过滤出可序列化且非Flask相关的变量
        with open('all_variables.pkl', 'wb') as f:
            variables_to_save = {}
            for k, v in globals().items():
                # 先过滤掉 Flask 相关的上下文代理对象和模块类型
                if isinstance(v, werkzeug.local.LocalProxy) or isinstance(v, types.ModuleType):
                    continue
                # 只保存可以序列化的对象
                if is_picklable(v):
                    variables_to_save[k] = v
            pickle.dump(variables_to_save, f)
            
    def read_pickle(self):
        """
        Loads and restores previously serialized variables from a pickle file.
    
        Returns:
            None
        """
        # 确保在加载前定义所需的自定义函数和类
        def is_picklable(obj):
            try:
                pickle.dumps(obj)
            except (pickle.PicklingError, TypeError):
                return False
            return True
        # 从保存的文件中加载变量
        with open('all_variables.pkl', 'rb') as f:
            loaded_variables = pickle.load(f)
        # 将变量重新加载到当前全局命名空间中
        globals().update(loaded_variables)
    
    
 