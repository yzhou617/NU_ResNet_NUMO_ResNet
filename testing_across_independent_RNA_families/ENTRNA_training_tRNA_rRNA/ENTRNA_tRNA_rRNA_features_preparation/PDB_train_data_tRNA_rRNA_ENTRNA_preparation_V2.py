
import numpy as np
import pandas as pd
from extract_features import extract_features_pseudoknot_free


# 1. Load the data
Original_RNA_Data_train = pd.read_csv('PDB_data_tRNA_rRNA_training_data.csv')

# 2. Select columns
Original_RNA_Data_train = Original_RNA_Data_train[['RNA_name', 'RNA_seq', 'RNA_secondary_structure', 'RNA_generated_seq']]

#Original_RNA_Data_train = Original_RNA_Data_train.rename(columns={'MCC_RNAfold' : 'MCC'})

# 3. Normalize the MCC
#Original_RNA_Data_train['Normalized_MCC']=np.nan

#Original_RNA_Data_train['Normalized_MCC'] = (Original_RNA_Data_train['MCC']+1)/2

# 4. Obtain the pair of RNA sequence and correct RNA secondary structure  
Original_RNA_Data_actual_RNA_structure = Original_RNA_Data_train[['RNA_name', 'RNA_seq', 'RNA_secondary_structure']].copy(deep=True)

#Original_RNA_Data_actual_RNA_structure['Normalized_MCC']=np.nan

#Original_RNA_Data_actual_RNA_structure['Normalized_MCC']=1

# 5. Obtain the pair of RNA sequence and predicted RNA secondary structure whose normalized MCC less than 1.
Original_RNA_Data_predicted_RNA_structure = Original_RNA_Data_train[['RNA_name', 'RNA_generated_seq', 'RNA_secondary_structure']].copy(deep=True)

#Original_RNA_Data_predicted_RNA_structure = Original_RNA_Data_predicted_RNA_structure[Original_RNA_Data_predicted_RNA_structure['Normalized_MCC'] < 1]

# 6. Obtain the label and rename the corresponding column name

Original_RNA_Data_actual_RNA_structure['RNA_Label'] = 1

Original_RNA_Data_predicted_RNA_structure['RNA_Label'] = 0

Original_RNA_Data_actual_RNA_structure = Original_RNA_Data_actual_RNA_structure.rename(columns={"RNA_seq":"RNA_sequence", "RNA_secondary_structure":"RNA_sec_structure"})

Original_RNA_Data_predicted_RNA_structure = Original_RNA_Data_predicted_RNA_structure.rename(columns={"RNA_generated_seq":"RNA_sequence", "RNA_secondary_structure":"RNA_sec_structure"})

# 7. Concatenate two data sets
Original_RNA_Data_combined = pd.concat([Original_RNA_Data_actual_RNA_structure, Original_RNA_Data_predicted_RNA_structure], axis=0, ignore_index=True) 

# 8. Obtain the label for each sample  
#Original_RNA_Data_combined['RNA_Label'] = np.where(Original_RNA_Data_combined['Normalized_MCC'] < 1, 0, 1)

# 9. Obtain the features for ENTRNA
Original_RNA_Data_combined['ent_3']=np.nan
Original_RNA_Data_combined['gc_perentage']=np.nan
Original_RNA_Data_combined['ensemble_diversity']=np.nan
Original_RNA_Data_combined['expected_accuracy']=np.nan
Original_RNA_Data_combined['fe_per']=np.nan

num_row_Original_RNA_Data_combined = Original_RNA_Data_combined.shape[0]

for i in range(num_row_Original_RNA_Data_combined):
    ENTRNA_features = extract_features_pseudoknot_free(seq=Original_RNA_Data_combined.at[i, 'RNA_sequence'].upper(), sec_str=Original_RNA_Data_combined.at[i, 'RNA_sec_structure'])
    Original_RNA_Data_combined.at[i, 'ent_3'] = ENTRNA_features.at[0, 'ent_3']
    Original_RNA_Data_combined.at[i, 'gc_perentage'] = ENTRNA_features.at[0, 'gc_perentage']
    Original_RNA_Data_combined.at[i, 'ensemble_diversity'] = ENTRNA_features.at[0, 'ensemble_diversity']
    Original_RNA_Data_combined.at[i, 'expected_accuracy'] = ENTRNA_features.at[0, 'expected_accuracy']
    Original_RNA_Data_combined.at[i, 'fe_per'] = ENTRNA_features.at[0, 'fe_per']
    
# 10. Shuffle the data
Original_RNA_Data_combined = Original_RNA_Data_combined.sample(frac=1, random_state=182).reset_index(drop=True)

# 11. Save the results in .csv file
Original_RNA_Data_combined.to_csv('PDB_training_data_tRNA_rRNA_ENTRNA_features_20240510.csv', index=False)










