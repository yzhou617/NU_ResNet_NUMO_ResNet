
import pandas as pd
import numpy as np

#%%

def analysis_results_equilibrium_prob_function(testing_data_file_name):
    
    Original_RNA_Data = pd.read_csv(testing_data_file_name)
    Original_RNA_Data['indicator_FE_p_larger_MFE_p']=np.nan
    Original_RNA_Data['indicator_FE_p_equal_to_MFE_p']=np.nan
    Original_RNA_Data['indicator_FE_p_less_FE_n']=np.nan
    Original_RNA_Data['indicator_prob_p_larger_prob_n']=np.nan
    
    num_row_Original_RNA_Data = Original_RNA_Data.shape[0] 
    
    for i in range(num_row_Original_RNA_Data):
        if Original_RNA_Data.at[i, 'FE_p_RNA_str'] > Original_RNA_Data.at[i, 'FE_p_MFE_RNA_str']:
            Original_RNA_Data.at[i, 'indicator_FE_p_larger_MFE_p'] = 1
        else:
            Original_RNA_Data.at[i, 'indicator_FE_p_larger_MFE_p'] = 0

        if Original_RNA_Data.at[i, 'FE_p_RNA_str'] == Original_RNA_Data.at[i, 'FE_p_MFE_RNA_str']:
            Original_RNA_Data.at[i, 'indicator_FE_p_equal_to_MFE_p'] = 1
        else:
            Original_RNA_Data.at[i, 'indicator_FE_p_equal_to_MFE_p'] = 0

        if Original_RNA_Data.at[i, 'FE_p_RNA_str'] < Original_RNA_Data.at[i, 'FE_n_RNA_str']:
            Original_RNA_Data.at[i, 'indicator_FE_p_less_FE_n'] = 1
        else:
            Original_RNA_Data.at[i, 'indicator_FE_p_less_FE_n'] = 0
            
        if Original_RNA_Data.at[i, 'positive_sample_ensemble_str_prob'] > Original_RNA_Data.at[i, 'negative_sample_ensemble_str_prob']:
            Original_RNA_Data.at[i, 'indicator_prob_p_larger_prob_n'] = 1
        else:
            Original_RNA_Data.at[i, 'indicator_prob_p_larger_prob_n'] = 0
            
        percent_FE_p_larger_MFE_p = Original_RNA_Data['indicator_FE_p_larger_MFE_p'].sum()/num_row_Original_RNA_Data
        percent_FE_p_equal_to_MFE_p = Original_RNA_Data['indicator_FE_p_equal_to_MFE_p'].sum()/num_row_Original_RNA_Data
        percent_FE_p_less_FE_n = Original_RNA_Data['indicator_FE_p_less_FE_n'].sum()/num_row_Original_RNA_Data
        percent_prob_p_larger_prob_n = Original_RNA_Data['indicator_prob_p_larger_prob_n'].sum()/num_row_Original_RNA_Data
        
    print('The percentage that the FE of positive sample greater than the MFE structure is ', percent_FE_p_larger_MFE_p)
    print('The percentage that the FE of positive sample equal to the MFE structure is ', percent_FE_p_equal_to_MFE_p)
    print('The percentage that the FE of positive sample less than the FE of negative sample is ', percent_FE_p_less_FE_n)    
    print('The percentage that the equilibrium prob of positive samples greater than the equilibrium prob of negative sample is ', percent_prob_p_larger_prob_n)
   
    return percent_FE_p_larger_MFE_p, percent_prob_p_larger_prob_n
    
    




