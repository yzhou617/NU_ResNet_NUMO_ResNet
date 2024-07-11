
import pandas as pd
from data_preparation_RNA_mat_3d_nt_localized_info_mat_function_PDB_data import data_prep_RNA_mat_3d_nt_localized_info_mat_PDB_data
import torch
import torch.nn as nn
from train_val_test_plot_functions import create_weighted_sampler, dataloader_prep_with_sampler, dataloader_prep, train_val_ResNet_expert, test_ResNet_expert, train_val_figure_plot_function
from ResNet_architecture_grayscale_mat_nt_localized_info_mat import ResNet_18_pair_grayscale_mat_nt_localized_info_mat
import numpy as np
from PDB_data_RNA_family_riboswitch_data_testing_NUMO_ResNet import RNA_family_testing_function_riboswitch_data_NUMO_ResNet

#%%
###### First Part: Data Preparation
# 1.1 Read the data, select the utilized columns, and rename the column name.

# Training data
#Original_RNA_Data_train = pd.read_csv('PDB_data_tRNA_rRNA_training_data.csv')
#Original_RNA_Data_train = Original_RNA_Data_train.head(10)
'''
Original_RNA_Data_train['RNA_seq_upper']=np.nan
Original_RNA_Data_train['RNA_seq_upper'] = Original_RNA_Data_train['RNA_seq_upper'].astype("string")

num_row_Original_RNA_Data_train = Original_RNA_Data_train.shape[0]

for i in range(num_row_Original_RNA_Data_train):
    Original_RNA_Data_train.at[i, "RNA_seq_upper"] = Original_RNA_Data_train.at[i, "RNA_seq"].upper() 
'''

#Original_RNA_Data_train = Original_RNA_Data_train[['name', 'sequence', 'correct_structure', 'RNAfold_structure', 'MCC_RNAfold']]
#Original_RNA_Data_train = Original_RNA_Data_train.rename(columns={'MCC_RNAfold' : 'MCC'})


# Validation data
#Original_RNA_Data_val = pd.read_csv('PDB_data_tRNA_rRNA_validation_data.csv')
#Original_RNA_Data_val = Original_RNA_Data_val.head(10)
'''
Original_RNA_Data_val['RNA_seq_upper']=np.nan
Original_RNA_Data_val['RNA_seq_upper'] = Original_RNA_Data_val['RNA_seq_upper'].astype("string")

num_row_Original_RNA_Data_val = Original_RNA_Data_val.shape[0]

for i in range(num_row_Original_RNA_Data_val):
    Original_RNA_Data_val.at[i, "RNA_seq_upper"] = Original_RNA_Data_val.at[i, "RNA_seq"].upper()
'''

#Original_RNA_Data_val = Original_RNA_Data_val[['name', 'sequence', 'correct_structure', 'RNAfold_structure', 'MCC_RNAfold']]
#Original_RNA_Data_val = Original_RNA_Data_val.rename(columns={'MCC_RNAfold' : 'MCC'})

'''
# Testing data
Original_RNA_Data_test = pd.read_csv('PDB_data_reliable_negative_samples_testing_data.csv')
#Original_RNA_Data_test = Original_RNA_Data_test.head(10)

Original_RNA_Data_test['RNA_seq_upper']=np.nan
Original_RNA_Data_test['RNA_seq_upper'] = Original_RNA_Data_test['RNA_seq_upper'].astype("string")

num_row_Original_RNA_Data_test = Original_RNA_Data_test.shape[0]

for i in range(num_row_Original_RNA_Data_test):
    Original_RNA_Data_test.at[i, "RNA_seq_upper"] = Original_RNA_Data_test.at[i, "RNA_seq"].upper()

Original_RNA_Data_test = Original_RNA_Data_test[~Original_RNA_Data_test.RNA_seq_upper.str.contains('|'.join(["P"]))].reset_index(drop=True)
#Original_RNA_Data_test = Original_RNA_Data_test[['name', 'sequence', 'correct_structure', 'RNAfold_structure', 'MCC_RNAfold']]

#Original_RNA_Data_test = Original_RNA_Data_test.rename(columns={'MCC_RNAfold' : 'MCC'})
'''
'''
# Rfam testing data
Original_RNA_Data_test_rfam = pd.read_csv('rfam_mcc_summary.csv')
#Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam.head(10)

Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam[['name', 'sequence', 'correct_structure', 'RNAfold_structure', 'MCC_RNAfold']]

Original_RNA_Data_test_rfam = Original_RNA_Data_test_rfam.rename(columns={'MCC_RNAfold' : 'MCC'})
'''

#%%

'''
# Training data
train_x_color_mat_np, train_x_nt_localized_info_mat_np, train_y_np = data_prep_RNA_mat_3d_nt_localized_info_mat_PDB_data(Original_RNA_Data = Original_RNA_Data_train, padding_length=410.0)

train_x_color_mat_4D_torch = torch.from_numpy(train_x_color_mat_np)
train_x_color_mat_4D_torch = train_x_color_mat_4D_torch.type(torch.float)

train_x_nt_localized_info_mat_4D_torch = torch.from_numpy(train_x_nt_localized_info_mat_np)
train_x_nt_localized_info_mat_4D_torch = train_x_nt_localized_info_mat_4D_torch.type(torch.float)

train_y_torch = torch.from_numpy(train_y_np)
train_y_torch = train_y_torch.type(torch.long)

# Validation data
val_x_color_mat_np, val_x_nt_localized_info_mat_np, val_y_np = data_prep_RNA_mat_3d_nt_localized_info_mat_PDB_data(Original_RNA_Data = Original_RNA_Data_val, padding_length=410.0)

val_x_color_mat_4D_torch = torch.from_numpy(val_x_color_mat_np)
val_x_color_mat_4D_torch = val_x_color_mat_4D_torch.type(torch.float)

val_x_nt_localized_info_mat_4D_torch = torch.from_numpy(val_x_nt_localized_info_mat_np)
val_x_nt_localized_info_mat_4D_torch = val_x_nt_localized_info_mat_4D_torch.type(torch.float)

val_y_torch = torch.from_numpy(val_y_np)
val_y_torch = val_y_torch.type(torch.long)
'''
'''
# Testing data
test_x_color_mat_np, test_x_nt_localized_info_mat_np, test_y_np = data_prep_RNA_mat_3d_nt_localized_info_mat_PDB_data(Original_RNA_Data = Original_RNA_Data_test, padding_length=410.0)

test_x_color_mat_4D_torch = torch.from_numpy(test_x_color_mat_np)
test_x_color_mat_4D_torch = test_x_color_mat_4D_torch.type(torch.float)

test_x_nt_localized_info_mat_4D_torch = torch.from_numpy(test_x_nt_localized_info_mat_np)
test_x_nt_localized_info_mat_4D_torch = test_x_nt_localized_info_mat_4D_torch.type(torch.float)

test_y_torch = torch.from_numpy(test_y_np)
test_y_torch = test_y_torch.type(torch.long)
'''

'''
# Rfam testing data
rfam_test_x_color_mat_np, rfam_test_x_nt_localized_info_mat_np, rfam_test_y_np = data_prep_grayscale_mat_nt_localized_info_mat(Original_RNA_Data = Original_RNA_Data_test_rfam)

rfam_test_x_color_mat_4D_torch = torch.from_numpy(rfam_test_x_color_mat_np)
rfam_test_x_color_mat_4D_torch = rfam_test_x_color_mat_4D_torch.type(torch.float)

rfam_test_x_nt_localized_info_mat_4D_torch = torch.from_numpy(rfam_test_x_nt_localized_info_mat_np)
rfam_test_x_nt_localized_info_mat_4D_torch = rfam_test_x_nt_localized_info_mat_4D_torch.type(torch.float)

rfam_test_y_torch = torch.from_numpy(rfam_test_y_np)
rfam_test_y_torch = rfam_test_y_torch.type(torch.long)
'''

#%%

#sampler_weight_data_loader_utilized = create_weighted_sampler(y_numpy = train_y_np)

device_utilized = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

NN_architecture_utilized = ResNet_18_pair_grayscale_mat_nt_localized_info_mat(nt_info_mat_input_num_channels=1, grayscale_mat_input_num_channels=4)

NN_model_utilized = NN_architecture_utilized.to(device_utilized)

# Training Part 2.3: Set hyperparameters, loss function, optimizer, and scheduler.
learning_rate_value_utilized = 0.0001
weight_decay_value_utilized = 0.1
num_epochs_utilized = 100
batch_size_value_utilized = 20
gamma_value_utilized = 0.95

optimizer_utilized = torch.optim.Adam(NN_model_utilized.parameters(), lr=learning_rate_value_utilized, weight_decay=weight_decay_value_utilized)
scheduler_for_optimizer_utilized = torch.optim.lr_scheduler.ExponentialLR(optimizer_utilized, gamma=gamma_value_utilized, last_epoch=-1)

loss_function_utilized = nn.CrossEntropyLoss()

model_name_utilized = 'NUMO_ResNet_batch_'+str(batch_size_value_utilized)+'_lr_'+str(learning_rate_value_utilized)+\
    '_weightd_'+str(weight_decay_value_utilized)+'_epochs_'+str(num_epochs_utilized)+'_scheduler_gamma_'+\
        str(gamma_value_utilized)+'_PDB_data_tRNA_rRNA.pth'



# Training Part 2.4: Prepare the Dataloader
'''
training_data_utilized = dataloader_prep_with_sampler(x_color_mat_stack_4D_tensor = train_x_color_mat_4D_torch,
                                                      x_nt_localized_info_mat_stack_4D_tensor = train_x_nt_localized_info_mat_4D_torch,
                                                      y_tensor = train_y_torch, 
                                                      batch_size_value = batch_size_value_utilized, 
                                                      sampler_applied = sampler_weight_data_loader_utilized)


val_data_utilized = dataloader_prep(x_color_mat_stack_4D_tensor = val_x_color_mat_4D_torch,
                                    x_nt_localized_info_mat_stack_4D_tensor = val_x_nt_localized_info_mat_4D_torch,
                                    y_tensor = val_y_torch, 
                                    batch_size_value = batch_size_value_utilized)
'''
'''
testing_data_utilized = dataloader_prep(x_color_mat_stack_4D_tensor = test_x_color_mat_4D_torch,
                                        x_nt_localized_info_mat_stack_4D_tensor = test_x_nt_localized_info_mat_4D_torch,
                                        y_tensor = test_y_torch, 
                                        batch_size_value = batch_size_value_utilized)
'''
'''
rfam_testing_data_utilized = dataloader_prep(x_color_mat_stack_4D_tensor = rfam_test_x_color_mat_4D_torch, 
                                             x_nt_localized_info_mat_stack_4D_tensor = rfam_test_x_nt_localized_info_mat_4D_torch, 
                                             y_tensor = rfam_test_y_torch, 
                                             batch_size_value = batch_size_value_utilized)
'''

#%%

'''
train_val_record_obtained, saved_models_metrics_obtained = train_val_ResNet_expert(device = device_utilized, 
                                                                                   NN_model = NN_model_utilized, 
                                                                                   num_epochs = num_epochs_utilized, 
                                                                                   optimizer = optimizer_utilized, 
                                                                                   scheduler_for_optimizer = scheduler_for_optimizer_utilized, 
                                                                                   loss_function = loss_function_utilized, 
                                                                                   model_name = model_name_utilized, 
                                                                                   training_data_used = training_data_utilized, 
                                                                                   val_data_used = val_data_utilized)

print('The metrics of best validation accuracy model is as follows.')
print(saved_models_metrics_obtained['best_val_acc_model'])

print('The metrics of best validation loss model is as follows.')
print(saved_models_metrics_obtained['best_val_loss_model'])

print('The metrics of the model from last epoch is as follows.')
print(saved_models_metrics_obtained['last_epoch_model'])
'''

# Save the results from training into csv files.
'''
train_val_record_file_name = 'NUMOResNet_train_val_record_batch_'+str(batch_size_value_utilized)+'_lr_'+str(learning_rate_value_utilized)+\
    '_weightd_'+str(weight_decay_value_utilized)+'_epochs_'+str(num_epochs_utilized)+'_scheduler_gamma_'+\
        str(gamma_value_utilized)+'PDB_data_tRNA_rRNA.csv'

df_train_val_record = pd.DataFrame.from_dict(train_val_record_obtained)

df_train_val_record.to_csv(train_val_record_file_name, index=False)

saved_models_metrics_file_name = 'NUMOResNet_saved_models_metrics_batch_'+str(batch_size_value_utilized)+'_lr_'+str(learning_rate_value_utilized)+\
    '_weightd_'+str(weight_decay_value_utilized)+'_epochs_'+str(num_epochs_utilized)+'_scheduler_gamma_'+\
        str(gamma_value_utilized)+'PDB_data_tRNA_rRNA.csv'

df_saved_models_metrics = pd.DataFrame.from_dict(saved_models_metrics_obtained)

df_saved_models_metrics.to_csv(saved_models_metrics_file_name, index=False)
'''

#%%
# Load the trained models
# Best val accuracy model

best_val_acc_model_name_utilized = 'Best_validation_accuracy_model_'+model_name_utilized

NN_archit_test_best_acc = ResNet_18_pair_grayscale_mat_nt_localized_info_mat(nt_info_mat_input_num_channels=1, grayscale_mat_input_num_channels=4)

NN_test_model_best_acc = NN_archit_test_best_acc.to(device_utilized)

NN_test_model_best_acc.load_state_dict(torch.load(best_val_acc_model_name_utilized, map_location=device_utilized)['saved_model'])
print('The best validation accuracy model is saved from epoch ', torch.load(best_val_acc_model_name_utilized, map_location=device_utilized)['epoch_num'])

# Best val loss model

best_val_loss_model_name_utilized = 'Best_validation_loss_model_'+model_name_utilized

NN_archit_test_best_loss = ResNet_18_pair_grayscale_mat_nt_localized_info_mat(nt_info_mat_input_num_channels=1, grayscale_mat_input_num_channels=4)

NN_test_model_best_loss = NN_archit_test_best_loss.to(device_utilized)

NN_test_model_best_loss.load_state_dict(torch.load(best_val_loss_model_name_utilized, map_location=device_utilized)['saved_model'])
print('The best validation loss model is saved from epoch ', torch.load(best_val_loss_model_name_utilized, map_location=device_utilized)['epoch_num'])

# Model from last epoch

model_from_last_epoch_name_utilized = 'Model_from_last_epoch_'+model_name_utilized

NN_archit_test_model_last_epoch = ResNet_18_pair_grayscale_mat_nt_localized_info_mat(nt_info_mat_input_num_channels=1, grayscale_mat_input_num_channels=4)

NN_test_model_last_epoch = NN_archit_test_model_last_epoch.to(device_utilized)

NN_test_model_last_epoch.load_state_dict(torch.load(model_from_last_epoch_name_utilized, map_location=device_utilized)['saved_model'])
print('The model from last epoch is saved from epoch ', torch.load(model_from_last_epoch_name_utilized, map_location=device_utilized)['epoch_num'])

'''
#%%
# 1. GroupIintron
print('GroupIintron: ')

GroupI_loss_b_v_a, GroupI_metric_b_v_a, GroupI_loss_b_v_l, GroupI_metric_b_v_l, GroupI_loss_model_l_e, GroupI_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_GroupIintron_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='GroupIintron',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 2. GroupIIintron
print('GroupIIintron: ')

GroupII_loss_b_v_a, GroupII_metric_b_v_a, GroupII_loss_b_v_l, GroupII_metric_b_v_l, GroupII_loss_model_l_e, GroupII_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_GroupIIintron_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='GroupIIintron',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 3. SRPRNA
print('SRPRNA: ')

SRP_loss_b_v_a, SRP_metric_b_v_a, SRP_loss_b_v_l, SRP_metric_b_v_l, SRP_loss_model_l_e, SRP_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_SRPRNA_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='SRPRNA',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 4. HairpinRibozyme
print('HairpinRibozyme: ')

HR_loss_b_v_a, HR_metric_b_v_a, HR_loss_b_v_l, HR_metric_b_v_l, HR_loss_model_l_e, HR_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_HairpinRibozyme_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='HairpinRibozyme',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 5. HammerheadRibozyme
print('HammerheadRibozyme: ')

HHR_loss_b_v_a, HHR_metric_b_v_a, HHR_loss_b_v_l, HHR_metric_b_v_l, HHR_loss_model_l_e, HHR_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_HammerheadRibozyme_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='HammerheadRibozyme',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)





#%%
# 6. otherRibozyme
print('otherRibozyme: ')

OR_loss_b_v_a, OR_metric_b_v_a, OR_loss_b_v_l, OR_metric_b_v_l, OR_loss_model_l_e, OR_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_otherRibozyme_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='otherRibozyme',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 7. ViralandPhage
print('ViralandPhage: ')

VP_loss_b_v_a, VP_metric_b_v_a, VP_loss_b_v_l, VP_metric_b_v_l, VP_loss_model_l_e, VP_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_ViralandPhage_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='ViralandPhage',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 8. SmallnuclearRNA
print('SmallnuclearRNA: ')

SN_loss_b_v_a, SN_metric_b_v_a, SN_loss_b_v_l, SN_metric_b_v_l, SN_loss_model_l_e, SN_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_SmallnuclearRNA_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='SmallnuclearRNA',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

#%%
# 9. InternalRibosomeEntrySite
print('InternalRibosomeEntrySite: ')

IRES_loss_b_v_a, IRES_metric_b_v_a, IRES_loss_b_v_l, IRES_metric_b_v_l, IRES_loss_model_l_e, IRES_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_InternalRibosomeEntrySite_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='InternalRibosomeEntrySite',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)


#%%
# 10. RNasePRNA
print('RNasePRNA: ')

RNaseP_loss_b_v_a, RNaseP_metric_b_v_a, RNaseP_loss_b_v_l, RNaseP_metric_b_v_l, RNaseP_loss_model_l_e, RNaseP_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_RNasePRNA_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='RNasePRNA',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)


#%%
# 11. otherRNA
print('otherRNA: ')

O_loss_b_v_a, O_metric_b_v_a, O_loss_b_v_l, O_metric_b_v_l, O_loss_model_l_e, O_metric_model_l_e = RNA_family_testing_function(test_data_file='PDB_data_otherRNA_with_info.csv',
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='otherRNA',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)


'''
'''
#%%
# concatenate the data from different RNA families
print('concatenate the data from different RNA families: ')

GroupIintron_RNA_Data = pd.read_csv('PDB_data_GroupIintron_with_info.csv')

GroupIIintron_RNA_Data = pd.read_csv('PDB_data_GroupIIintron_with_info.csv')

SRP_RNA_Data = pd.read_csv('PDB_data_SRPRNA_with_info.csv')

HairpinRibozyme_RNA_Data = pd.read_csv('PDB_data_HairpinRibozyme_with_info.csv')

HammerheadRibozyme_RNA_Data = pd.read_csv('PDB_data_HammerheadRibozyme_with_info.csv')

otherRibozyme_RNA_Data = pd.read_csv('PDB_data_otherRibozyme_with_info.csv')

ViralandPhage_RNA_Data = pd.read_csv('PDB_data_ViralandPhage_with_info.csv')

SmallnuclearRNA_RNA_Data = pd.read_csv('PDB_data_SmallnuclearRNA_with_info.csv')

InternalRibosomeEntrySite_RNA_Data = pd.read_csv('PDB_data_InternalRibosomeEntrySite_with_info.csv')

RNasePRNA_RNA_Data = pd.read_csv('PDB_data_RNasePRNA_with_info.csv')

concatenated_RNA_Data = pd.concat([GroupIintron_RNA_Data,
                                   GroupIIintron_RNA_Data,
                                   SRP_RNA_Data,
                                   HairpinRibozyme_RNA_Data,
                                   HammerheadRibozyme_RNA_Data,
                                   otherRibozyme_RNA_Data,
                                   ViralandPhage_RNA_Data,
                                   SmallnuclearRNA_RNA_Data,
                                   InternalRibosomeEntrySite_RNA_Data,
                                   RNasePRNA_RNA_Data], axis=0, ignore_index=True)

concatenated_RNA_Data.to_csv('concatenated_different_RNA_families_testing_data.csv', index=False)

'''


#%%
# obtain the riboswitch data

Riboswitch_RNA_Data = pd.read_csv('NU_ResNet_NUMO_ResNet_trained_tRNA_rRNA_data_test_Riboswitch_data_analysis_with_GC_percent.csv')

Riboswitch_RNA_Data['RNA_seq_upper']=np.nan

Riboswitch_RNA_Data['RNA_seq_upper'] = Riboswitch_RNA_Data['RNA_seq_upper'].astype("string")

num_row_Riboswitch_RNA_Data = Riboswitch_RNA_Data.shape[0]

for i in range(num_row_Riboswitch_RNA_Data):
        Riboswitch_RNA_Data.at[i, 'RNA_seq_upper']=Riboswitch_RNA_Data.at[i, 'sequence'].upper()

Riboswitch_RNA_Data_obtained = Riboswitch_RNA_Data.rename(columns={"name":"RNA_name", "GT_struct":"RNA_secondary_structure"})

Riboswitch_RNA_Data_obtained = Riboswitch_RNA_Data_obtained[['RNA_name', 'RNA_seq_upper', 'RNA_secondary_structure']]


concat_loss_b_v_a, concat_metric_b_v_a, concat_loss_b_v_l, concat_metric_b_v_l, concat_loss_model_l_e, concat_metric_model_l_e = RNA_family_testing_function_riboswitch_data_NUMO_ResNet(test_data_file='concatenated_different_RNA_families_testing_data.csv',
                                                                                                                                                                                         riboswitch_data_set = Riboswitch_RNA_Data_obtained,
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='concatenated_RNA_data',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

'''
#%%
# concatenate the data from different RNA families and other RNA
print('concatenate the data from different RNA families and other RNA: ')

otherRNA_RNA_Data = pd.read_csv('PDB_data_otherRNA_with_info.csv')

concatenated_different_RNA_families_otherRNA = pd.concat([concatenated_RNA_Data,
                                                          otherRNA_RNA_Data], axis=0, ignore_index=True)

concatenated_different_RNA_families_otherRNA.to_csv('concatenated_different_RNA_families_otherRNA_testing_data.csv', index=False)

'''

concat_2_loss_b_v_a, concat_2_metric_b_v_a, concat_2_loss_b_v_l, concat_2_metric_b_v_l, concat_2_loss_model_l_e, concat_2_metric_model_l_e = RNA_family_testing_function_riboswitch_data_NUMO_ResNet(test_data_file='concatenated_different_RNA_families_otherRNA_testing_data.csv',
                                                                                                                                                                                                     riboswitch_data_set = Riboswitch_RNA_Data_obtained,
                                                                                                                                                             device_utilized_para=device_utilized,
                                                                                                                                                             loss_function_utilized_para=loss_function_utilized,
                                                                                                                                                             testing_data_name_para='concatenated_data_from_different_RNA_families_otherRNA',
                                                                                                                                                             NN_test_model_best_acc_para=NN_test_model_best_acc,
                                                                                                                                                             best_val_acc_model_name_para=best_val_acc_model_name_utilized,
                                                                                                                                                             NN_test_model_best_loss_para=NN_test_model_best_loss,
                                                                                                                                                             best_val_loss_model_name_para=best_val_loss_model_name_utilized,
                                                                                                                                                             NN_test_model_last_epoch_para=NN_test_model_last_epoch,
                                                                                                                                                             model_from_last_epoch_name_para=model_from_last_epoch_name_utilized)

'''
#%% 

figure_name_utilized = 'NUMO_ResNet_batch_'+str(batch_size_value_utilized)+'_lr_'+str(learning_rate_value_utilized)+\
    '_weightd_'+str(weight_decay_value_utilized)+'_epochs_'+str(num_epochs_utilized)+'_scheduler_gamma_'+\
        str(gamma_value_utilized)+'_PDB_data_tRNA_rRNA'


train_val_figure_plot_function(figure_name = figure_name_utilized, 
                               loss_value_train = train_val_record_obtained['loss_train'], 
                               loss_value_val = train_val_record_obtained['loss_val'], 
                               accuracy_train = train_val_record_obtained['acc_train'], 
                               accuracy_val = train_val_record_obtained['acc_val'], 
                               auc_roc_train = train_val_record_obtained['aucroc_train'], 
                               auc_roc_val = train_val_record_obtained['aucroc_val'], 
                               auc_roc_train_check = train_val_record_obtained['aucroc_train_check'], 
                               auc_roc_val_check = train_val_record_obtained['aucroc_val_check'])

                               


'''
























