Link-prediction gene list classification 
===================================  


Classify gene fusion pairs into different groups according to frequency and probability
-----------------------------------



### 1. Import necessary packages  
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from matplotlib.colors import ListedColormap
    import datetime
    from sklearn import svm
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import make_moons, make_circles, make_classification
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
    from sklearn.naive_bayes import GaussianNB
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
    from sklearn.metrics import classification_report, confusion_matrix  
    import unittest
    from sklearn.svm import SVC
  
### 2. Set working directory to path where Link-prediction file is saved
    file_destination = "genefusion.predictions.2019.03.08.xlsx"
    df = pd.read_excel(file_destination)
    df.reset_index(inplace=True)
    df1 = {}


### 3. Make a list for gene and delete the duplicate
    genes = []
    count = 0
    id_value ="g1_id"

    def remove(duplicate): 
        final_list = [] 
        for num in duplicate: 
            if num not in final_list: 
                final_list.append(num) 
        return final_list 

    for i in range(0, 2960):
	    genes.append(df.loc[i][id_value])
    genes_wo_repeats = remove(genes)
    for i in range(0, len(genes_wo_repeats) - 1):
	    if type(genes_wo_repeats[i]) is not str:
		    del genes_wo_repeats[i]
          
### 4. Count the gene frequency    
    frequency = []
    for i in range(0, len(genes_wo_repeats)):
	    count = genes.count(genes_wo_repeats[i])
	
	    frequency.append(count)
  
### 5. Calculate gene frequency and divide it into top gene and bottom gene 
    def averagenum(frequency):
        nsum = 0
        for i in range(len(frequency)):
            nsum += frequency[i]
        return nsum / len(frequency)
    #Average of freqency is 6.7, so bigger than 7 is top gene, otherwise, is bottom  
    
    def create_sequence_dictionary(frequency, wo_repeats, lessthan, limit):
	    dictionary = {}
	    for i in range(0, len(frequency) - 1):
		    dictionary.update({wo_repeats[i] : frequency[i]})
	    for item in list(dictionary.keys()):
		    if lessthan == True:
			    if dictionary[item] > limit:
				    del dictionary[item]	
		    else:
			    if dictionary[item] < limit:
				    del dictionary[item]
	    return dictionary
    top_genes = create_sequence_dictionary(frequency, genes_wo_repeats, False, 7)
    bottom_genes = create_sequence_dictionary(frequency, genes_wo_repeats, True, 6)
    bot_gene_values = list(bottom_genes.keys())
    top_gene_values = list(top_genes.keys())
          
### 6. Calculate the fusion probability of genes
    top_gene_values_max = max(list(top_genes.values()))
    max_key = [key for key, value in top_genes.items() if value == top_gene_values_max]
    print(max_key[0])

    max_key_index = top_gene_values.index(max_key[0])
    print(max_key_index)

    all_gene = list(bottom_genes.keys()) + list(top_genes.keys())
    def find_probabilities(dictionary, index):
        dictionary_keys = list(dictionary.keys())   
        probabilities = []
        for item in dictionary_keys:
            indexes = df.loc[df[id_value] == item].index
            probabilities.append(df.loc[indexes[index]]["prediction"])
        return probabilities
    low_freq_prob = find_probabilities(bottom_genes, 0)
    high_freq_prob = find_probabilities(top_genes, 0)
    print(len(low_freq_prob))
    print(len(high_freq_prob))

    probabilities = low_freq_prob + high_freq_prob
    print(len(probabilities))
    frequencies = list(bottom_genes.values()) + list(top_genes.values())
    print(len(frequencies))

    probabilities_sorted = sorted(probabilities, reverse=True)
    indeces = []
    for i in range(0, 15):
        indeces.append(probabilities.index(probabilities_sorted[i]))
    frequencies_sorted = sorted(frequencies, reverse=True)
    indece_s = []
    for i in range(0, 20):
        indece_s.append(frequencies.index(frequencies_sorted[i]))
  
 ### 7. Genes are divided into high frequency groups and high probability groups
    desired_freq_values = []
    desired_prob_values = []

    for i in range(0, len(frequencies)):
	if(frequencies[i] >= 11):
		desired_freq_values.append(i)
	if(probabilities[i] >= .55):
		desired_prob_values.append(i)
    print(desired_prob_values)
    print(desired_freq_values)
 
 ### 8. Genes are divided into four groups according to the median of probability and frequency
 
 
 
 
 
 
 
 
 
 
 
