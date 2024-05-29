
import pandas as pd
from data_preparation_RNA_mat_3d_nt_localized_info_mat_function_PDB_data import data_prep_RNA_mat_3d_nt_localized_info_mat_PDB_data
import torch
import torch.nn as nn
from train_val_test_plot_functions_color_mat_ResNet_18 import create_weighted_sampler, dataloader_prep_with_sampler, dataloader_prep, train_val_ResNet_expert_color_mat, test_ResNet_expert, train_val_figure_plot_function
import numpy as np

def RNA_family_testing_function(test_data_file, device_utilized_para, loss_function_utilized_para, testing_data_name_para, NN_test_model_best_acc_para, best_val_acc_model_name_para, NN_test_model_best_loss_para, best_val_loss_model_name_para, NN_test_model_last_epoch_para, model_from_last_epoch_name_para):

    # 1.1 Read the data, select the utilized columns, and rename the column name.
    # Testing data
    Original_RNA_Data_test = pd.read_csv(test_data_file)
    #Original_RNA_Data_test = Original_RNA_Data_test.head(10)
    '''
    Original_RNA_Data_test['RNA_seq_upper']=np.nan
    Original_RNA_Data_test['RNA_seq_upper'] = Original_RNA_Data_test['RNA_seq_upper'].astype("string")

    num_row_Original_RNA_Data_test = Original_RNA_Data_test.shape[0]

    for i in range(num_row_Original_RNA_Data_test):
        Original_RNA_Data_test.at[i, "RNA_seq_upper"] = Original_RNA_Data_test.at[i, "RNA_seq"].upper()
    '''
    Original_RNA_Data_test = Original_RNA_Data_test[~Original_RNA_Data_test.RNA_seq_upper.str.contains('|'.join(["P"]))].reset_index(drop=True)

    # Testing data
    test_x_color_mat_np, test_x_nt_localized_info_mat_np, test_y_np = data_prep_RNA_mat_3d_nt_localized_info_mat_PDB_data(Original_RNA_Data = Original_RNA_Data_test, padding_length=410.0)

    test_x_color_mat_4D_torch = torch.from_numpy(test_x_color_mat_np)
    test_x_color_mat_4D_torch = test_x_color_mat_4D_torch.type(torch.float)

    test_x_nt_localized_info_mat_4D_torch = torch.from_numpy(test_x_nt_localized_info_mat_np)
    test_x_nt_localized_info_mat_4D_torch = test_x_nt_localized_info_mat_4D_torch.type(torch.float)

    test_y_torch = torch.from_numpy(test_y_np)
    test_y_torch = test_y_torch.type(torch.long)

    # Training Part 2.4: Prepare the Dataloader

    testing_data_utilized = dataloader_prep(x_color_mat_stack_4D_tensor = test_x_color_mat_4D_torch,
                                            x_nt_localized_info_mat_stack_4D_tensor = test_x_nt_localized_info_mat_4D_torch,
                                            y_tensor = test_y_torch, 
                                            batch_size_value = 1)

    # Test trained models on testing data

    loss_best_val_acc, metric_best_val_acc = test_ResNet_expert(device = device_utilized_para, NN_model = NN_test_model_best_acc_para, 
                                                                                            loss_function = loss_function_utilized_para, model_name = best_val_acc_model_name_para, 
                                                                                            testinging_data_used = testing_data_utilized, testing_data_name = testing_data_name_para)

    print('The metrics of best validation accuracy model is as follows.')
    print(metric_best_val_acc)

    loss_best_val_loss, metric_best_val_loss = test_ResNet_expert(device = device_utilized_para, NN_model = NN_test_model_best_loss_para, 
                                                                                            loss_function = loss_function_utilized_para, model_name = best_val_loss_model_name_para, 
                                                                                            testinging_data_used = testing_data_utilized, testing_data_name = testing_data_name_para)

    print('The metrics of best validation loss model is as follows.')
    print(metric_best_val_loss)

    loss_model_last_epoch, metric_model_last_epoch = test_ResNet_expert(device = device_utilized_para, NN_model = NN_test_model_last_epoch_para, 
                                                                                            loss_function = loss_function_utilized_para, model_name = model_from_last_epoch_name_para, 
                                                                                            testinging_data_used = testing_data_utilized, testing_data_name = testing_data_name_para)

    print('The metrics of the model from last epoch is as follows.')
    print(metric_model_last_epoch)

    return loss_best_val_acc, metric_best_val_acc, loss_best_val_loss, metric_best_val_loss, loss_model_last_epoch, metric_model_last_epoch 







