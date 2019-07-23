Link-prediction gene list classification 
===================================  


Classify gene fusion pairs into different groups according to frequency and proabality.
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
          
### 多行文本框    
    这是一个有多行的文本框  
    你可以写入代码等,每行文字只要输入两个Tab再输入文字即可  
    这里你可以输入一段代码  
  
### 比如我们可以在多行文本框里输入一段代码,来一个Java版本的HelloWorld吧  
 
