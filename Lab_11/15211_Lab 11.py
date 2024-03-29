import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Q1) I) i
# Loading iris dataset
iris = pd.read_csv('iris.csv')
iris_df = pd.DataFrame(iris)
print(iris_df)

# Extracting X and y variables from iris dataset
X = iris_df.iloc[:,:4]
y = iris_df['variety']
print(f"X values: \n{X}")

# How do you know which data to select?

# Q1) I) ii
# Conveting the categorical data to numerical
from sklearn.preprocessing import LabelEncoder

print(f"Y values: \n{y}")
le = LabelEncoder()
y_en = le.fit_transform(y) # Label encoded species
print(y_en)

iris_df['variety encoded'] = y_en
print(iris_df)
Y = iris_df['variety encoded']
print(Y)

# Splitting data into test and train
from sklearn.model_selection import train_test_split

x_train, x_test, y_train, y_test = train_test_split(X, Y,
                                                    test_size = 0.2,
                                                    random_state = 42)

print(f"Shape of x_train: {x_train.shape}\n"
      f"Shape of x_test: {x_test.shape}\n"
      f"Shape of y_train: {y_train.shape}\n"
      f"Shape of y_test: {y_test.shape}\n")

# Standardizing the iris feature values
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
x_train_sc = scaler.fit_transform(x_train)
x_test_sc = scaler.fit_transform(x_test)

# Training the model
from sklearn.neural_network import MLPClassifier

clf = MLPClassifier(hidden_layer_sizes = [10,10,10],
                    learning_rate_init = 0.001, max_iter = 1000,
                    verbose = True)
clf.fit(x_train_sc, y_train)

# Visualizing the validation loss
plt.figure()
plt.plot(clf.loss_curve_)
plt.title('Validation Loss')
plt.show()

# Predicting the species labels from the test data
y_pred = clf.predict(x_test_sc)
print(y_pred)

# Evaluating the model on the test data using
# the classification report and confusion matrix
from sklearn.metrics import classification_report, confusion_matrix
report = classification_report(y_test, y_pred)
print(report)

c_matrix = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:\n")
print(c_matrix)
print('\n')

from sklearn.metrics import multilabel_confusion_matrix
mlc_matrix = multilabel_confusion_matrix(y_test, y_pred)
print(mlc_matrix)
print('\n')

# Predict the species of new iris plants using
# trained model
new_data = [[5.9,3.0,7.0,5.0],[4.6, 3.0, 1.5, 0.2],[6.2, 3.0,4.1,1.2]]
new_data_sc = scaler.fit_transform(new_data)
new_pred = clf.predict(new_data_sc)
print(new_pred)

############################################################################

# Q1) V
# Training the model
from sklearn.neural_network import MLPClassifier

clf = MLPClassifier(hidden_layer_sizes = [2,2,2],
                    learning_rate_init = 0.001, max_iter = 1000,
                    verbose = True)
clf.fit(x_train_sc, y_train)

# Visualizing the validation loss
plt.figure()
plt.plot(clf.loss_curve_)
plt.title('Validation Loss')
plt.show()

# Predicting the species labels from the test data
y_pred = clf.predict(x_test_sc)
print(y_pred)

# Evaluating the model on the test data using
# precision and accuracy scores
report = classification_report(y_test, y_pred)
print(report)

c_matrix = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:")
print(c_matrix)
print('\n')

mlc_matrix = multilabel_confusion_matrix(y_test, y_pred)
print(mlc_matrix)
print('\n')




