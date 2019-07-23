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
    median = np.median(probabilities_sorted)

    y = []
    for i in range(0, len(frequencies)):
        if probabilities[i] >= median and frequencies[i] <= 6: 
            y.append(0) # This is high prob low freq group
        elif probabilities[i] >= median and frequencies[i] > 6:
            y.append(1) # This is high prob high freq group
        elif probabilities[i] < median and frequencies[i] <= 6: 
            y.append(2) # This is low prob low freq group
        elif probabilities[i] < median and frequencies[i] > 6: 
            y.append(3) # This is low prob high freq group

    X = np.array([probabilities, frequencies])
    X = np.transpose(X)
    print(X.shape)
    print(y.shape)
 
### 9. Train the model for the classifier algorithm
    X_train = [list(np.random.rand(len(X)))]
    X_train_freq = list(np.random.randint(1, 3, size=len(list(bottom_genes.values())))) + list(np.random.randint(12, 30,     size=len(list(top_genes.values()))))

    y_results = list((np.zeros(len(low_freq_prob), dtype=int))) + list(np.random.randint(1,2, size=len(high_freq_prob)))


    h = 0.2
    names = ["Nearest Neighbors", "Linear SVM", "RBF SVM", "Decision Tree", "Random Forest", "AdaBoost",
    "Naive Bayes", "LDA", "QDA"]
    classifiers = [
    KNeighborsClassifier(2), svm.SVC(kernel="linear"), svm.SVC(gamma=2, C=1), DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1), AdaBoostClassifier(), GaussianNB(),
    LDA(), QDA()]

    clf = classifiers[1]
    clf.fit(X, y)
 
### 10. Make SVM classification scatter plot
    def make_meshgrid(x, y, h=0.2):
        x_min, x_max = x.min() -1, x.max() +1
        y_min, y_max = y.min() - 1, y.max() + 1
        xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
        return xx, yy
    def plot_contours(ax, clf, xx, yy, **params):
        Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        out = ax.contourf(xx, yy, Z, **params)
    title = ('SVM')
    fig, ax = plt.subplots()
    X0 = X[:,0]
    X1 = X[:,1]
    xx, yy = make_meshgrid(X0, X1)
    plot_contours(ax, clf, xx, yy, cmap=plt.cm.coolwarm, alpha=1)
    ax.scatter(X0, X1, c=y, cmap=plt.cm.coolwarm, s=20, edgecolors='k')
    ax.set_ylabel('frequency')
    ax.set_xlabel('probability')
    ax.set_xticks(())
    ax.set_yticks(())
    ax.set_title(title)
    ax.legend()
    for i in range(0, len(desired_prob_values)):
        ax.annotate(all_gene[desired_prob_values[i]], xy=(X0[desired_prob_values[i]], X1[desired_prob_values[i]]),
            xytext=(X0[desired_prob_values[i]], X1[desired_prob_values[i]]), fontsize=12)
    for i in range(0, len(desired_freq_values)):
        ax.annotate(all_gene[desired_freq_values[i]], xy=(X0[desired_freq_values[i]], X1[desired_freq_values[i]]),
            xytext=(X0[desired_freq_values[i]], X1[desired_freq_values[i]]), fontsize=12)
    plt.show()
 
 
 
 
 
 
