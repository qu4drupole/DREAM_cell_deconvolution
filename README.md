# DREAM_cell_deconvolution 
Scripts, tools, and notes about DREAM challenge for cell type deconvolution

==========================================================================

!!!!    YOU CAN NOT ACCESS THE DATABASE UNLESS YOU GIVE ME YOUR IP ADDRESS    !!!!

==========================================================================

files in 'Bulk_seq_data_loading'....
  - There are three jupyter notebooks that outline the data processing procedure:
    - 'data_import_preprocess' does the heavy lifting
    - 'RNAseq_data_class' describes a python data object that can read data from the SQL database--_not completed!_
    - 'seqfromsql.py' is what you want to import to use the RNAseq_data class
    - 'EDA_LOGFILE-SQL' is exploration of bulk data loading process--_not completed!_
  - 'gseList.csv' is a csv file of all trainging GSEs and relevant attributes
  - 'test_log.txt' is a log file from the data import process
  - 'success_gse.pckl' is a pickled set of GSE records successfully uploaded to the SQL database
