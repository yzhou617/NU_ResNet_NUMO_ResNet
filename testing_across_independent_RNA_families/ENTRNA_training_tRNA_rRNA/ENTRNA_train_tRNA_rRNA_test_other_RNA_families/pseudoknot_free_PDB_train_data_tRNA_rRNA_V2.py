import numpy as np
import pickle
import pandas as pd
from sklearn import preprocessing
from scipy.spatial import distance
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.model_selection import KFold
import itertools
from sklearn.metrics import accuracy_score
from scipy.spatial import distance
import extract_features
from sklearn.preprocessing import Normalizer


def bulid_model_PDB_train_data(real_RNA_loc):
    df_real = pd.read_csv(real_RNA_loc, header = 0)
    df_v1 = df_real[df_real['fe_per'] < 2]
    #file_name_list = df_v1['file_name'].tolist()
    #df_v1['length'] = df_v1['seq'].apply(lambda x:len(x))
    result_training_accuracy_list = []
    result_test_accuracy_list = []
    feature_set = ['ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy', 'fe_per']

    mat_all = df_v1[feature_set].values

    #train = list(range(len(file_name_list)))
    #mat_train_real = mat_all[train]
    #file_name_train_list = [file_name_list[i] for i in train]
    #x_negative_list = []
    #for ii in range(len(file_name_train_list)):
    #    x_true = mat_train_real[ii]
    #    df_temp = pd.read_csv(folder_simulation_result+file_name_train_list[ii]+'.csv', header = 0)
    #    mat_temp_v1 = df_temp[feature_set].as_matrix()

    #    distance_list = []
    #    for row_i in mat_temp_v1:
    #        if np.isnan(row_i).any() == True:
    #            distance_list.append(0)
    #        else:
    #            dst = distance.euclidean(row_i,x_true)
    #            distance_list.append(dst)
    #    x_negative_list.append(mat_temp_v1[np.argmax(distance_list)])
    x_train_raw = mat_all
    scaler = Normalizer().fit(x_train_raw)
    x_train = scaler.transform(x_train_raw)
    y_train = df_v1['RNA_Label'].values
    clf = LogisticRegression()        
    clf.fit(x_train,y_train)    
    train_acc = accuracy_score(y_train,clf.predict(x_train))

    return scaler, clf, train_acc

    
def pred_foldability(df_pred,clf, scaler):
    mat_x = df_pred[['ent_3', 'gc_perentage', 'ensemble_diversity', 'expected_accuracy', 'fe_per']].values
    for row_i in mat_x:
        if np.isnan(row_i).any() == True:
            return float('nan')
            
        else:
            return clf.predict_proba(scaler.transform([row_i,]))[0][1]


def entrna_main_PDB_train_data_tRNA_rRNA(seq,sec_str,real_RNA_loc = "PDB_training_data_tRNA_rRNA_ENTRNA_features_20240510.csv"):
    df_pred = extract_features.extract_features_pseudoknot_free(seq,sec_str)
    scaler, clf, train_acc = bulid_model_PDB_train_data(real_RNA_loc)
    foldability = pred_foldability(df_pred,clf, scaler)
    return foldability, train_acc

'''
def entrna_train(real_RNA_loc = "./util/RNASTRAND_pseudoknot_free_feature.csv", folder_simulation_result = "./util/RNASTRAND_extract_feature_pseudoknot_free/"):
    scaler, clf, train_acc = bulid_model(real_RNA_loc,folder_simulation_result)
    return train_acc
'''

