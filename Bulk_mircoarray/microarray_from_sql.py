class microarray_data:
    # import modules
    pd = __import__('pandas')
    np = __import__('numpy')
    re = __import__('re')
    pk = __import__('pickle')
    #from sqlalchemy import create_engine 
    import sqlalchemy as sqal
    import GEOparse
    import mysql.connector
    
    #engine = sqal.create_engine('mysql+mysqlconnector://dream_user:dream_sql_pw@192.168.144.21/test_dream')
    engine = sqal.create_engine('mysql+mysqlconnector://Simon:Bane@localhost/microarray_data')
    
    # basic attributes
    def __init__(self, ct=None, norm='RMA',  scope='coarse'):

        cell_types = ['NK_cells', 'fibroblast', 'Monocytes', 'CD8_T_cells', 'B_cells',
       'Dendritic_cells', 'endothelial', 'PBMC', 'Memory_CD4_T_cells',
       'Th17_cells', 'Naive_CD4_T_cells', 'Th1_cells', 'Th2_cells',
       'CD4_T_cells', 'Macrophage', 'Neutrophils', 'Granulocytes',
       'T_cells', 'Central_Memory', 'Plasmacytoid_Dendritic_Cells',
       'Plasma_cells', 'Tregs', 'Basophils', 'Mast_cells',
       'Naive_B_cells', 'GC_B_cells', 'Memory_B_cells', 'Plasmablast',
       'Immature_B_cells', 'White_blood_cells', 'Activated_B_cells',
       'Immature_Dendritic_cells', 'Myeloid_Dendritic_cells',
       'Eosinophils', 'Effector_Memory_T_cells', 'Central_Memory_T_cells',
       'Effector_Memory', 'Plasma_B_cells', 'Naive_CD8_T_cells',
       'Central_Memory_CD4_T_cells', 'Effector_Memory_CD4_T_cells',
       'gamma_delta_T_cells', 'Resting_B_cells', 'Activated_T_cells',
       'Memory_T_cells', 'Naive_T_cells', 'NKT_cells',
       'Effector_Memory_CD8_T_cells', 'Memory_CD8_T_cells',
       'Follicular_helper_T_cells', 'Th1Th17_cells',
       'Central_Memory_CD8_T_cells', 'Activated_Macrophages',
       'Effector_CD4_T_cells', 'pre-GC_B_cells']
        cell_types = list(map(lambda x: x.lower(), cell_types))
        
        norm_types = ['RMA','gcRMA','RMA-quantile','MAS5','unknown']
        
        if any(param == None for param in [ct, norm]):
            print('Instantiate class with: RNAseq_data(CELL TYPE, NORMALIZATION, SCOPE) \n \n'
                  'where CELL TYPE is a list of cell types. Must be one or more of:'+str(cell_types)+'\n\n'
                  'where NORMALIZATION is how the counts are normalized. Must be one of: '+str(norm_types)+'\n\n'
                  'where SCOPE is the cell type specificity. Must be either \'fine\'(default) or \'coarse\'')
            return
        
        with open('picklejar/ctDict.pckl', 'rb') as ctDictFile:
            cellDict = self.pk.load(ctDictFile)
            self.allCtDict = cellDict
            coarse_ctDict = cellDict['Coarse']
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
    
    # INTERNAL USE ONLY
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
    
    # INTERNAL USE ONLY
    def mergeTest(self,df1,df2):
        #check too make sure not too many NAs exist...
        minRow = df1.shape[0]
        df_merge = df1.join(df2, how='inner')
        if df_merge.shape[0] < minRow*0.2:
            return 0
        else:
            return 1
    
    # INTERNAL USE ONLY  
    def condenseGSE(self,df):
        gseSet = set(map(lambda x: self.re.search('GSE_(GSE\d+)',x).group(1), df.columns))
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
    
    def allData(self, ct):
        df_dict = {}
        # for working with 'Coarse' cells
        if len(self.ctDict) > 0:
            if ct != self.ct:
                ct_query = ct
            else:
                ct_query = self.ctDict.keys()
            for celltype in ct_query:
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
        # for working with 'Fine' cells
        else: 
            for celltype in ct:
                df = self.pd.read_sql_table(celltype, con=self.engine)
                if sum(df.duplicated(subset='gene_symbol_sql')) > 0:
                    df.drop_duplicates(subset='gene_symbol_sql', inplace=True)
                df.drop('index', 1, inplace=True)
                df.set_index('gene_symbol_sql', inplace=True)
                df = self.condenseGSE(df)
                df_dict[celltype] = df
        return df_dict
    
    def normData(self, ct):
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
                        if len(sampleNames_norm) < 2:
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
                        continue
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
    
    # takes argument of cell type--but should know whether coarse of fine
    def convertUnknown(self, ct_df, forMixing = 0):
        unknown_gse = self.getNormalizationPlus(ct_df)
        unknown_gsms = []
        for gse, norm in unknown_gse.items():
            if norm == 'RAW':
                gsmCols = ct_df.columns[ct_df.columns.str.contains(gse)]
                unknown_gsms.extend(gsmCols)
        if forMixing == 1:
            if any(ct_df.columns.str.contains('RAW')):
                raw_gsms = ct_df.columns[ct_df.columns.str.contains('RAW')]
                unknown_gsms.extend(raw_gsms)
        df_unknown = ct_df[unknown_gsms]
        df_unknown.dropna(inplace=True)
        if forMixing == 0:
            df_unknown = df_unknown.loc[df_unknown.index.intersection(self.genelens.index)]
            df_unknown = self.convertRawtoFPKM(df_unknown, self.genelens)
        return df_unknown
        
    # INTERNAL USE ONLY
    # this function is for making random mixtures with only 1 target ct
    def createMixRatio(self, ctList, mix):
        # I'm specifically writing this to only work with an object instantiated with
        # all coarse cell types--and ctList is length 1 with only the ct of interest
        random_total = 1-mix[0]
        # need to randomly mix the remaining 7 CTs
        other_mix = self.np.random.uniform(0, 1, len(ctList) - 1)
        other_mix = random_total * other_mix / other_mix.sum()
        new_mix = self.np.insert(other_mix, 0, mix[0])
        return new_mix
        
    def createMixDF(self, ctList, n_samp=10, mix = None):
        if mix is None:
            mix = [1/len(ctList)] * len(ctList)
        
        # usage is where you only specify 1 cell type
        if len(mix) == 1:
            other_ct = list(self.ctDict.keys())
            other_ct.remove(ctList[0])
            ctList.extend(other_ct)
            train_ctData = self.normData(ct=ctList[0])
            otherCts = [ct for ct in self.ct if ct != ctList[0]]
            train_mixData = self.normData(ct=otherCts)
            for i in range(n_samp):
#                 print('working on '+str(i)+'\n\n')
                ct_df = train_ctData[ctList[0]].sample(1, axis=1)
                colname = ct_df.columns.values[0]
                ct_df.rename(columns={colname:ctList[0]}, inplace = True)
                for ct in otherCts:
                    merge_df = train_mixData[ct].sample(1, axis = 1)
                    colname = merge_df.columns.values[0]
                    merge_df.rename(columns={colname:ct}, inplace = True)
                    ct_df = ct_df.join(merge_df, how='inner')
                # mix
                if len(mix) == len(ctList):
                    ct_mix = dict(zip(ctList,mix))
                else:
                    rand_mix = self.createMixRatio(ctList, mix)
                    ct_mix = dict(zip(ctList,rand_mix))
#                 print(str(ct_mix))
                for ct, ratio in ct_mix.items():
                    ct_df.loc[:,ct] = ct_df[ct]*ratio
#                 print(ct_df.head())
                if i == 0:
                    mixed_counts = ct_df.apply(sum, 1)
                    mix_sample_df = self.pd.DataFrame(mixed_counts, columns=['sample'+str(i)])
                else:
                    mixed_counts = ct_df.apply(sum, 1)
                    merge_sample_df = self.pd.DataFrame(mixed_counts, columns=['sample'+str(i)])
                    mix_sample_df = mix_sample_df.join(merge_sample_df, how='inner')
#                 print(mix_sample_df.head()) 
        return mix_sample_df
        
        
        