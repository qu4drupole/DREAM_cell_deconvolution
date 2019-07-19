<<<<<<< HEAD
# basic class to read data from SQL database

class RNAseq_data:
    # import modules
    pd = __import__('pandas')
    np = __import__('numpy')
    re = __import__('re')
    #from sqlalchemy import create_engine 
    import sqlalchemy as sqal
    import mysql.connector
    
    # basic attributes
    def __init__(self, ct=None, norm=None,  scope='fine'):

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
        
    def data(self):
        engine = self.sqal.create_engine('mysql+mysqlconnector://Simon:Bane@localhost/test_dream')
        df = self.pd.read_sql_table(self.ct[0], con=engine)
        return df
=======
# basic class to read data from SQL database

class RNAseq_data:
    # import modules
    pd = __import__('pandas')
    np = __import__('numpy')
    re = __import__('re')
    #from sqlalchemy import create_engine 
    import sqlalchemy as sqal
    import mysql.connector
    
    # basic attributes
    def __init__(self, ct=None, norm=None,  scope='fine'):

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
        
    def data(self):
        engine = self.sqal.create_engine('mysql+mysqlconnector://Simon:Bane@localhost/test_dream')
        df = self.pd.read_sql_table(self.ct[0], con=engine)
        return df
>>>>>>> 1b230e48b06b5b780618589572042b11524ce0f0
        