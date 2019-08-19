class RNAseq_data:
    # import modules
    pd = __import__('pandas')
    np = __import__('numpy')
    re = __import__('re')
    pk = __import__('pickle')
    #from sqlalchemy import create_engine 
    import sqlalchemy as sqal
    import mysql.connector
    
    #engine = sqal.create_engine('mysql+mysqlconnector://dream_user:dream_sql_pw@192.168.144.21/test_dream')
    engine = sqal.create_engine('mysql+mysqlconnector://Simon:Bane@localhost/test_dream')
    
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
                  'where SCOPE is the cell type specificity. Must be either \'fine\'(default) or \'coarse\'')
            return
        
        with open('picklejar/CellType.pkl', 'rb') as ctDictFile:
            cellDict = self.pk.load(ctDictFile)
            coarse_ctDict = cellDict['Coarse']['Main']
            coarseCells = coarse_ctDict.keys()
        
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
        
        self.ctDict = {}
        if scope in ['coarse','fine']:
            self.scope = scope
            if scope == 'coarse':
                if all(cell in coarseCells for cell in ct):
                    for cell in ct:
                        self.ctDict[cell] = coarse_ctDict[cell]
                else:
                    print('Usage Error: Specified cell type not a coarse cell type \n')
                    print('Must be one or more of: \n'+str(coarseCells))
                    return
        else:
            print('Usage Error: scope must be of type \'coarse\' or \'fine\' (default)')
            return
    
    def primaryDf(self,df,gseSet):
        minNa = 500000
        gse_name = None
        for gse in gseSet:
            gse_col = df.columns[df.columns.str.contains(gse)]
            gse_df = df[gse_col]
            gseNAnum = gse_df.isna().any(axis=1).sum()
            if gseNAnum < minNa:
                minNa = gseNAnum
                gse_name = gse
        return gse_name
            
    def mergeTest(self,df1,df2):
        #check too make sure not too many NAs exist...
        minRow = df1.shape[0]
        df_merge = df1.join(df2, how='inner')
        if df_merge.shape[0] < minRow*0.2:
            return 0
        else:
            return 1
        
    def condenseGSE(self,df):
        gseSet = set(map(lambda x: self.re.search('gse_(GSE\d+)',x).group(1), df.columns))
        gse1_name = self.primaryDf(df,gseSet)
        gseSet.remove(gse1_name)
        df1 = df[df.columns[df.columns.str.contains(gse1_name)]].copy()
        df1.dropna(inplace=True)
        for gse in gseSet:
            df2 = df[df.columns[df.columns.str.contains(gse)]].copy()
            df2.dropna(inplace=True)
            if self.mergeTest(df1, df2) == 1:
                df1 = df1.join(df2, how='inner')
        df = df1
        return df
    
    def allData(self):
        df_dict = {}
        if len(self.ctDict) > 0:
            for celltype in self.ctDict.keys():
                token = 0
                for subcell in self.ctDict[celltype]:
                    try:
                        df = self.pd.read_sql_table(subcell, con=self.engine)
                        if sum(df.duplicated(subset='gene_symbol_sql')) > 0:
                            df.drop_duplicates(subset='gene_symbol_sql', inplace=True)
                        df.drop('index', 1, inplace=True)
                        df.set_index('gene_symbol_sql', inplace=True)
                        df = self.condenseGSE(df)
                        if token == 0:
                            token = 1
                            df_dict[celltype] = df
                        else:
                            df_dict[celltype] = df_dict[celltype].join(df, how='inner')
                    except:
                        pass
        else: 
            for celltype in self.ct:
                df = self.pd.read_sql_table(celltype, con=self.engine)
                if sum(df.duplicated(subset='gene_symbol_sql')) > 0:
                    df.drop_duplicates(subset='gene_symbol_sql', inplace=True)
                df.drop('index', 1, inplace=True)
                df.set_index('gene_symbol_sql', inplace=True)
                df = self.condenseGSE(df)
                df_dict[celltype] = df
        return df_dict
    
    def normData(self):
        df_dict = {}
        if len(self.ctDict) > 0:
            for celltype in self.ctDict.keys():
                token = 0
                for subcell in self.ctDict[celltype]:
                    try:
                        df = self.pd.read_sql_table(subcell, con=self.engine)
                        if sum(df.duplicated(subset='gene_symbol_sql')) > 0:
                            df.drop_duplicates(subset='gene_symbol_sql', inplace=True)
                        df.drop('index', 1, inplace=True)
                        df.set_index('gene_symbol_sql', inplace=True)
                        sampleNames_norm = df.columns[df.columns.str.contains(self.norm)]
                        if len(sampleNames_norm) == 1:
                            df = None
                            continue
                        else:
                            df = df[sampleNames_norm]
                            df = self.condenseGSE(df)
                            #insert merge check to dropnas
                        if token == 0:
                            token = 1
                            df_dict[celltype] = df
                        else:
                            df_dict[celltype] = df_dict[celltype].join(df, how='inner')
                    except:
                        pass
        else:
            for celltype in self.ct:
                try:
                    df = self.pd.read_sql_table(celltype, con=self.engine)
                    if sum(df.duplicated(subset='gene_symbol_sql')) > 0:
                        df.drop_duplicates(subset='gene_symbol_sql', inplace=True)
                    df.drop('index', 1, inplace=True)
                    df.set_index('gene_symbol_sql', inplace=True)
                    sampleNames_norm = df.columns[df.columns.str.contains(self.norm)]
                    if len(sampleNames_norm) == 1:
                        df = None
                        continue
                    else: 
                        df = df[sampleNames_norm]
                        df = self.condenseGSE(df)
                    df_dict[celltype] = df
                except:
                    pass
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
                df = df.join(merge_df, how = 'inner')
                merge_dict['cellTypes'].extend([celltype] * len(merge_df.columns))
            else:
                df = ct_dfs[celltype]
                merge_dict['cellTypes'].extend([celltype] * len(df.columns))
                initial = 1
        df.dropna(inplace=True)
        merge_dict['merged_df'] = df
        return merge_dict