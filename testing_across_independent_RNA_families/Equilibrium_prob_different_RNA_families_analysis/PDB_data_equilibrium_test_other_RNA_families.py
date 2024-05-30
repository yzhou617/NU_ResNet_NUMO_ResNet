
from PDB_test_data_ensemble_analysis_function_V2 import equilibrium_prob_testing_function
import pandas as pd

# 1. GroupIintron
print('GroupIintron: ') 
RNA_testing_data_1 = 'PDB_data_GroupIintron_with_info.csv'
RNA_testing_data_classification_results_1 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_1, testing_data_name='GroupIintron')

# 2. GroupIIintron
print('GroupIIintron: ')
RNA_testing_data_2 = 'PDB_data_GroupIIintron_with_info.csv'
RNA_testing_data_classification_results_2 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_2, testing_data_name='GroupIIintron')

# 3. HairpinRibozyme
print('HairpinRibozyme: ')
RNA_testing_data_3 = 'PDB_data_HairpinRibozyme_with_info.csv'
RNA_testing_data_classification_results_3 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_3, testing_data_name='HairpinRibozyme')

# 4. HammerheadRibozyme
print('HammerheadRibozyme: ')
RNA_testing_data_4 = 'PDB_data_HammerheadRibozyme_with_info.csv'
RNA_testing_data_classification_results_4 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_4, testing_data_name='HammerheadRibozyme')

# 5. otherRibozyme
print('otherRibozyme: ')
RNA_testing_data_5 = 'PDB_data_otherRibozyme_with_info.csv'
RNA_testing_data_classification_results_5 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_5, testing_data_name='otherRibozyme')

# 6. RNasePRNA
print('RNasePRNA: ')
RNA_testing_data_6 = 'PDB_data_RNasePRNA_with_info.csv'
RNA_testing_data_classification_results_6 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_6, testing_data_name='RNasePRNA')

# 7. SRPRNA
print('SRPRNA: ')
RNA_testing_data_7 = 'PDB_data_SRPRNA_with_info.csv'
RNA_testing_data_classification_results_7 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_7, testing_data_name='SRPRNA')

# 8. ViralandPhage
print('ViralandPhage: ')
RNA_testing_data_8 = 'PDB_data_ViralandPhage_with_info.csv'
RNA_testing_data_classification_results_8 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_8, testing_data_name='ViralandPhage')

# 9. SmallnuclearRNA
print('SmallnuclearRNA: ')
RNA_testing_data_9 = 'PDB_data_SmallnuclearRNA_with_info.csv'
RNA_testing_data_classification_results_9 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_9, testing_data_name='SmallnuclearRNA')

# 10. InternalRibosomeEntrySite
print('InternalRibosomeEntrySite: ')
RNA_testing_data_10 = 'PDB_data_InternalRibosomeEntrySite_with_info.csv'
RNA_testing_data_classification_results_10 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_10, testing_data_name='InternalRibosomeEntrySite')

# 11. Combined data from different RNA families
print('Combined data from different RNA families: ')
RNA_testing_data_11 = 'PDB_testing_data_different_RNA_families_concatenation.csv'
RNA_testing_data_classification_results_11 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_11, testing_data_name='combined_data_from_different_RNA_families')

# 12. otherRNA
print('otherRNA: ')
RNA_testing_data_12 = 'PDB_data_otherRNA_with_info.csv'
RNA_testing_data_classification_results_12 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_12, testing_data_name='otherRNA')

# 13. Combined data from different RNA families and other RNA
print('Combined data from different RNA families and other RNA: ')
RNA_testing_data_13 = 'PDB_testing_data_different_RNA_families_otherRNA_concatenation.csv'
RNA_testing_data_classification_results_13 = equilibrium_prob_testing_function(testing_data_file_name=RNA_testing_data_13, testing_data_name='combined_data_from_different_RNA_families_and_otherRNA')


