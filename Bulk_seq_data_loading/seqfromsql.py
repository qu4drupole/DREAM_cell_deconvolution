class RNAseq_data:
    # import modules
    pd = __import__('pandas')
    np = __import__('numpy')
    re = __import__('re')
    #from sqlalchemy import create_engine 
    import sqlalchemy as sqal
    import mysql.connector
    
    engine = sqal.create_engine('mysql+mysqlconnector://dream_user:dream_sql_pw@192.168.144.21/test_dream')
    #engine = sqal.create_engine('mysql+mysqlconnector://Simon:Bane@localhost/test_dream')
    
    # basic attributes
    def __init__(self, ct=None, norm='FPKM',  scope='fine'):
        
        course_bcells = ['b_cells', 'naive_b_cells', 'activated_b_cells', 'memory_b_cells']
        course_cd4 = ['']

        cell_types = ['fibroblast', 'Naive CD4 T cells', 'PBMC', 'endothelial',
       'Monocytes', 'Macrophage', 'B cells', 'GC B cells', 'CD4 T cells',
       'T cells', 'Granulocytes', 'Memory CD4 T cells', 'NK cells',
       'Basophils', 'Central Memory CD8 T cells', 'Effector CD4 T cells',
       'Effector Memory CD8 T cells', 'Follicular helper T cells',
       'Memory B cells', 'Myeloid Dendritic cells', 'Naive B cells',
       'Naive CD8 T cells', 'Neutrophils', 'gamma delta T cells',
       'Th1 cells', 'Th17 cells', 'Th2 cells', 'Tregs', 'Plasmablast',
       'CD8 T cells', 'Plasmacytoid Dendritic Cells', 'Dendritic cells',
       'Activated T cells', 'White blood cells', 'Eosinophils',
       'Naive T cells', 'Central Memory', 'Effector Memory',
       'Central Memory T cells', 'Memory CD8 T cells', 'Plasma cells',
       'Memory T cells', 'NKT cells', 'Central Memory CD4 T cells',
       'Effector Memory T cells', 'Activated B cells',
       'Naive T effector cells']
        cell_types = list(map(lambda x: self.re.sub(' ','_',x.lower()), cell_types))
        norm_types = ['FPKM','RPKM','TPM','RAW','unknown']
        
        if any(param == None for param in [ct, norm]):
            print('Instantiate class with: RNAseq_data(CELL TYPE, NORMALIZATION, SCOPE) \n \n'
                  'where CELL TYPE is a list of cell types. Must be one or more of:'+str(cell_types)+'\n\n'
                  'where NORMALIZATION is how the counts are normalized. Must be one of: '+str(norm_types)+'\n\n'
                  'where SCOPE is the cell type specificity. Must be either \'fine\'(default) or \'course\'')
            return
        
        if isinstance(ct,(list,)):
            self.ct = ct
        else:
            print('Usage Error: cell type must be a list with one or more of: '+str(cell_types))
            return
        
        if norm in norm_types:
            self.norm = norm
        else:
            print('Usage Error: normalization method must be one of: '+str(norm_types))
            return
        
        if scope in ['course','fine']:
            self.scope = scope
        else:
            print('Usage Error: scope must be of type \'course\' or \'fine\' (default)')
            return
        
    def allData(self):
        df_dict = {}
        for celltype in self.ct:
            df = self.pd.read_sql_table(celltype, con=self.engine, )
            df.drop('index', 1, inplace=True)
            df_dict[celltype] = df
        return df_dict
    
    def normData(self):
        df_dict = {}
        for celltype in self.ct:
            df = self.pd.read_sql_table(celltype, con=self.engine)
            df.drop('index', 1, inplace=True)
            sampleNames = df.select_dtypes(exclude=['object']).columns.to_numpy()
            sampleNames_idx = list(map(lambda x: self.re.search('norm_(.*)',x).group(1) == self.norm, sampleNames))
            sampleNames_norm = sampleNames[sampleNames_idx]
            sampleNames_norm = self.np.insert(sampleNames_norm,0,'gene_symbol_sql')
            if len(sampleNames_norm) == 1:
                df = None
            else: 
                df = df[sampleNames_norm]
            df_dict[celltype] = df
        return df_dict
        
    def mergeCellTypes(self):
        # holds the merged df and cell type list
        merge_dict = {}
        merge_dict['cellTypes'] = [] 
        initial = 0
        ct_dfs = self.normData()
        for celltype in ct_dfs.keys():
            if initial == 1:
                merge_df = ct_dfs[celltype]
                merge_df.set_index(merge_df['gene_symbol_sql'], inplace=True)
                merge_df.drop(['gene_symbol_sql'], 1, inplace=True)
                df = df.join(merge_df, how = 'inner')
                merge_dict['cellTypes'].extend([celltype] * len(merge_df.columns))
            else:
                df = ct_dfs[celltype]
                df.set_index(df['gene_symbol_sql'], inplace=True)
                df.drop(['gene_symbol_sql'], 1, inplace=True)
                merge_dict['cellTypes'].extend([celltype] * len(df.columns))
                initial = 1
        df.dropna(inplace=True)
        merge_dict['merged_df'] = df
        return merge_dict