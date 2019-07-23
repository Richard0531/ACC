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
### 单行文本框  
    这是一个单行的文本框,只要两个Tab再输入文字即可  
          
### 多行文本框    
    这是一个有多行的文本框  
    你可以写入代码等,每行文字只要输入两个Tab再输入文字即可  
    这里你可以输入一段代码  
  
### 比如我们可以在多行文本框里输入一段代码,来一个Java版本的HelloWorld吧  
 
