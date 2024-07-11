
import numpy as np
import torch
from .nt_localized_info_matrix_generation_function import nt_localized_info_mat_generator
from .RNA_matrix_3d import RNA_matrix_3d_generator_canonical_bp
from .ResNet_architecture_grayscale_mat_nt_localized_info_mat import ResNet_18_pair_grayscale_mat_nt_localized_info_mat



def NUMO_ResNet_model(RNA_sequence, RNA_sec_structure, padding_length, model_name, embed):
    # 1.1 Obtain the color matrix and Local info
    color_mat = RNA_matrix_3d_generator_canonical_bp(RNA_sequence, RNA_sec_structure)
    nt_info_mat = nt_localized_info_mat_generator(RNA_sequence, RNA_sec_structure)

    # 1.2 Obtain the size of color matrix
    color_mat_size = color_mat.shape[2]

    # 1.3 Set the size which the color matrix needs to be padded to
    size_maximum = padding_length

    # 1.4 Do the padding for the color matrix

    size_diff = size_maximum - color_mat_size
    if (size_diff % 2) == 0:
        num_padding_left_above = int(size_diff/2)
        num_padding_right_below = int(size_diff/2)
    else:
        num_padding_left_above = int(size_diff/2)
        num_padding_right_below = int(size_diff/2)+1
    
    # 1.5 Do the padding for the color matrix
        
    RNA_mat_3d_padding = np.pad(color_mat, ((0, 0), (num_padding_left_above, num_padding_right_below), (num_padding_left_above, num_padding_right_below)), 'constant')
    
    RNA_mat_3d_padding_4d = RNA_mat_3d_padding[np.newaxis,:,:,:]
    
    RNA_mat_3d_padding_torch = torch.from_numpy(RNA_mat_3d_padding_4d)
    
    RNA_mat_3d_padding_torch = RNA_mat_3d_padding_torch.type(torch.float)
    
    # 1.6 Do the padding for nt info

    RNA_nt_info_padding = np.pad(nt_info_mat, ((0, int(size_diff)), (0, 0)), 'constant')

    RNA_nt_info_padding_3d = RNA_nt_info_padding[np.newaxis, np.newaxis,:,:]

    RNA_nt_info_padding_torch = torch.from_numpy(RNA_nt_info_padding_3d)
    
    RNA_nt_info_padding_torch = RNA_nt_info_padding_torch.type(torch.float)

    # 1.7 Run Model on Sequence
    device_utilized = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    NN_architecture_utilized = ResNet_18_pair_grayscale_mat_nt_localized_info_mat(nt_info_mat_input_num_channels=1,grayscale_mat_input_num_channels=4)
    
    NN_model_utilized = NN_architecture_utilized.to(device_utilized)
    
    NN_model_utilized.load_state_dict(torch.load(model_name, map_location=device_utilized)['saved_model'])
    
    NN_model_utilized.eval()
    
    with torch.no_grad():
        RNA_mat_3d_padding_torch_utilized = RNA_mat_3d_padding_torch.to(device_utilized)

        RNA_nt_info_padding_torch_utilized = RNA_nt_info_padding_torch.to(device_utilized)

        if embed:
            result_1 = NN_model_utilized.resnet_18_color_mat_modified(RNA_mat_3d_padding_torch_utilized)
            result_2 = NN_model_utilized.resnet_18_nt_localized_mat_modified(RNA_nt_info_padding_torch_utilized)
            embedded_output = torch.cat((result_1, result_2), dim=1)

            return embedded_output
        else:
            pred_test = NN_model_utilized(RNA_mat_3d_padding_torch_utilized, RNA_nt_info_padding_torch_utilized)
            
            pred_test_score = torch.nn.functional.softmax(pred_test.cpu().detach(), dim=1)
        
            return pred_test_score

# Test DEPRICATED, CHANGE IMPORTS TO TEST OR GO OUTSIDE FOLDER
'''
RNA_seq = "CUCGUCUAGUCAUUUCUGGCCCCACUGGAGGUCGAG"
RNA_ss = "((((.((((......))))((((...)).)).))))"
padding_len = 410.0
model_path_utilized = "NUMO_ResNet_PDB_data/Best_validation_loss_model_NUMO_ResNet_batch_20_lr_0.0001_weightd_0.1_epochs_100_scheduler_gamma_0.95_PDB_data.pth"
print(NUMO_ResNet_model(RNA_seq, RNA_ss, padding_len, model_path_utilized, embed=True).size())
'''
