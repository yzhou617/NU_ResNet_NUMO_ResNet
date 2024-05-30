
from PDB_test_data_ENTRNA_function_V3 import ENTRNA_testing_function_train_tRNA_rRNA
import pandas as pd

# 1. GroupIintron
print('GroupIintron: ') 
RNA_testing_data_1 = 'PDB_data_GroupIintron_with_info.csv'
RNA_testing_data_classification_results_1 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_1)

# 2. GroupIIintron
print('GroupIIintron: ')
RNA_testing_data_2 = 'PDB_data_GroupIIintron_with_info.csv'
RNA_testing_data_classification_results_2 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_2)

# 3. HairpinRibozyme
print('HairpinRibozyme: ')
RNA_testing_data_3 = 'PDB_data_HairpinRibozyme_with_info.csv'
RNA_testing_data_classification_results_3 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_3)

# 4. HammerheadRibozyme
print('HammerheadRibozyme: ')
RNA_testing_data_4 = 'PDB_data_HammerheadRibozyme_with_info.csv'
RNA_testing_data_classification_results_4 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_4)

# 5. otherRibozyme
print('otherRibozyme: ')
RNA_testing_data_5 = 'PDB_data_otherRibozyme_with_info.csv'
RNA_testing_data_classification_results_5 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_5)

# 6. RNasePRNA
print('RNasePRNA: ')
RNA_testing_data_6 = 'PDB_data_RNasePRNA_with_info.csv'
RNA_testing_data_classification_results_6 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_6)

# 7. SRPRNA
print('SRPRNA: ')
RNA_testing_data_7 = 'PDB_data_SRPRNA_with_info.csv'
RNA_testing_data_classification_results_7 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_7)

# 8. ViralandPhage
print('ViralandPhage: ')
RNA_testing_data_8 = 'PDB_data_ViralandPhage_with_info.csv'
RNA_testing_data_classification_results_8 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_8)

# 9. SmallnuclearRNA
print('SmallnuclearRNA: ')
RNA_testing_data_9 = 'PDB_data_SmallnuclearRNA_with_info.csv'
RNA_testing_data_classification_results_9 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_9)

# 10. InternalRibosomeEntrySite
print('InternalRibosomeEntrySite: ')
RNA_testing_data_10 = 'PDB_data_InternalRibosomeEntrySite_with_info.csv'
RNA_testing_data_classification_results_10 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_10)

# 11. Combined data from different RNA families
print('Combined data from different RNA families: ')
RNA_testing_data_11 = 'PDB_testing_data_different_RNA_families_concatenation.csv'
RNA_testing_data_classification_results_11 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_11)

# 12. otherRNA
print('otherRNA: ')
RNA_testing_data_12 = 'PDB_data_otherRNA_with_info.csv'
RNA_testing_data_classification_results_12 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_12)

# 13. Combined data from different RNA families and other RNA
print('Combined data from different RNA families and other RNA: ')
RNA_testing_data_13 = 'PDB_testing_data_different_RNA_families_otherRNA_concatenation.csv'
RNA_testing_data_classification_results_13 = ENTRNA_testing_function_train_tRNA_rRNA(testing_data_file_name=RNA_testing_data_13)


