
from Equilibrium_prob_analysis_V2 import analysis_results_equilibrium_prob_function  
import pandas as pd



#%%

file_name_list=['GroupIintronequilibrium_prob_analysis.csv',
                'GroupIIintronequilibrium_prob_analysis.csv',
                'HairpinRibozymeequilibrium_prob_analysis.csv',
                'HammerheadRibozymeequilibrium_prob_analysis.csv',
                'otherRibozymeequilibrium_prob_analysis.csv',
                'RNasePRNAequilibrium_prob_analysis.csv',
                'SRPRNAequilibrium_prob_analysis.csv',
                'ViralandPhageequilibrium_prob_analysis.csv',
                'SmallnuclearRNAequilibrium_prob_analysis.csv',
                'InternalRibosomeEntrySiteequilibrium_prob_analysis.csv',
                'combined_data_from_different_RNA_familiesequilibrium_prob_analysis.csv',
                'otherRNAequilibrium_prob_analysis.csv',
                'combined_data_from_different_RNA_families_and_otherRNAequilibrium_prob_analysis.csv']

for file_name in file_name_list:
    
    print(file_name)
    result_1, result_2 = analysis_results_equilibrium_prob_function(testing_data_file_name=file_name)
    










