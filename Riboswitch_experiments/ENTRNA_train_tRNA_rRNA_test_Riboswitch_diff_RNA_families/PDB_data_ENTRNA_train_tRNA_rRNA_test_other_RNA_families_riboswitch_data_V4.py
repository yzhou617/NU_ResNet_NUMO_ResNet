
from PDB_test_data_ENTRNA_function_riboswitch_data_V3 import ENTRNA_testing_function_train_tRNA_rRNA_test_riboswitch_data
from PDB_test_data_ENTRNA_function_riboswitch_data_other_RNA_families_V3 import ENTRNA_testing_function_train_tRNA_rRNA_test_other_RNA_families_riboswitch_data
import pandas as pd
import numpy as np

'''
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
'''


# 14. Riboswitch
# obtain the riboswitch data

Riboswitch_RNA_Data = pd.read_csv('NU_ResNet_NUMO_ResNet_trained_tRNA_rRNA_data_test_Riboswitch_data_analysis_with_GC_percent.csv')

Riboswitch_RNA_Data['RNA_seq']=np.nan

Riboswitch_RNA_Data['RNA_seq'] = Riboswitch_RNA_Data['RNA_seq'].astype("string")

num_row_Riboswitch_RNA_Data = Riboswitch_RNA_Data.shape[0]

for i in range(num_row_Riboswitch_RNA_Data):
        Riboswitch_RNA_Data.at[i, 'RNA_seq']=Riboswitch_RNA_Data.at[i, 'sequence'].upper()

Riboswitch_RNA_Data_obtained = Riboswitch_RNA_Data.rename(columns={"name":"RNA_name", "GT_struct":"RNA_secondary_structure"})

Riboswitch_RNA_data_obtained = Riboswitch_RNA_Data_obtained[['RNA_name', 'RNA_seq', 'RNA_secondary_structure']]

print('Riboswitch: ')
RNA_testing_data_14 = Riboswitch_RNA_data_obtained
RNA_testing_data_classification_results_14 = ENTRNA_testing_function_train_tRNA_rRNA_test_riboswitch_data(testing_data=RNA_testing_data_14)

# 11. Combined data from different RNA families
print('Combined data from different RNA families: ')
RNA_testing_data_11 = 'PDB_testing_data_different_RNA_families_concatenation.csv'
RNA_testing_data_classification_results_11 = ENTRNA_testing_function_train_tRNA_rRNA_test_other_RNA_families_riboswitch_data(testing_data_file_name=RNA_testing_data_11, riboswitch_data=RNA_testing_data_14)

# 13. Combined data from different RNA families and other RNA
print('Combined data from different RNA families and other RNA: ')
RNA_testing_data_13 = 'PDB_testing_data_different_RNA_families_otherRNA_concatenation.csv'
RNA_testing_data_classification_results_13 = ENTRNA_testing_function_train_tRNA_rRNA_test_other_RNA_families_riboswitch_data(testing_data_file_name=RNA_testing_data_13, riboswitch_data=RNA_testing_data_14)









