
import pandas as pd
import numpy as np 
from model_files.NUMO_ResNet_model_function_V2 import NUMO_ResNet_model
from model_files_NU_ResNet.NU_ResNet_model_function_V2 import NU_ResNet_model


padding_len = 410.0
NU_ResNet_model_path_utilized = "./model_files_NU_ResNet/Best_validation_loss_model_NU_ResNet_batch_20_lr_0.0001_weightd_0.15_epochs_100_scheduler_gamma_0.95_PDB_data_tRNA_rRNA.pth"
NUMO_ResNet_model_path_utilized = "./model_files/Best_validation_loss_model_NUMO_ResNet_batch_20_lr_0.0001_weightd_0.1_epochs_100_scheduler_gamma_0.95_PDB_data_tRNA_rRNA.pth"

Original_RNA_Data_obtained = pd.read_excel('41592_2022_1605_MOESM2_ESM.xlsx', sheet_name = "12.SHAPEdirectedFolding")
#Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam.head(10)

# Obtain the Riboswitch data

Original_RNA_Data_test_rfam_obtained = Original_RNA_Data_obtained[['name', 'sequence', 'GT_struct']]
#Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam.rename(columns={'MCC_RNAfold' : 'MCC'})

Original_RNA_Data_test_rfam_selected = Original_RNA_Data_test_rfam_obtained.iloc[[1, 3, 6, 16]]

Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam_selected.reset_index(drop=True)


Original_RNA_Data_test_rfam['NU_ResNet_score_correctss']=np.nan
Original_RNA_Data_test_rfam['NUMO_ResNet_score_correctss']=np.nan
Original_RNA_Data_test_rfam['RNA_length']=np.nan
Original_RNA_Data_test_rfam['GC_percent']=np.nan

num_row_Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam.shape[0]
    
for i in range(num_row_Original_RNA_Data_test_rfam):
    print(i)
    Original_RNA_Data_test_rfam.at[i, "RNA_length"] = len(Original_RNA_Data_test_rfam.at[i, "sequence"])
    
    NU_ResNet_score_correct_ss = NU_ResNet_model(RNA_sequence=Original_RNA_Data_test_rfam.at[i, "sequence"].upper(), 
                                       RNA_sec_structure=Original_RNA_Data_test_rfam.at[i, "GT_struct"], 
                                       padding_length=padding_len, 
                                       model_name=NU_ResNet_model_path_utilized)
    NU_ResNet_score_correct_ss_np = NU_ResNet_score_correct_ss.numpy()
    
    Original_RNA_Data_test_rfam.at[i, "NU_ResNet_score_correctss"] = NU_ResNet_score_correct_ss_np[0, 1]

    NUMO_ResNet_score_correct_ss = NUMO_ResNet_model(RNA_sequence=Original_RNA_Data_test_rfam.at[i, "sequence"].upper(), 
                                       RNA_sec_structure=Original_RNA_Data_test_rfam.at[i, "GT_struct"], 
                                       padding_length=padding_len, 
                                       model_name=NUMO_ResNet_model_path_utilized,
                                       embed=False)
    NUMO_ResNet_score_correct_ss_np = NUMO_ResNet_score_correct_ss.numpy()
    
    Original_RNA_Data_test_rfam.at[i, "NUMO_ResNet_score_correctss"] = NUMO_ResNet_score_correct_ss_np[0, 1]
    RNA_seq = Original_RNA_Data_test_rfam.at[i, "sequence"].upper()
    Original_RNA_Data_test_rfam.at[i, 'GC_percent'] = (RNA_seq.count('G') + RNA_seq.count('C')) / float(len(RNA_seq))
    
Original_RNA_Data_test_rfam.to_csv("NU_ResNet_NUMO_ResNet_trained_tRNA_rRNA_data_test_Riboswitch_data_analysis_with_GC_percent.csv", index=False)



 