import numpy as np
import pandas as pd
from sklearn.datasets import load_iris

pd.set_option('display.max_columns', None)

# Q2) I
# importing iris dataset
iris = load_iris()
x_train = iris.data
y_train = iris.target
print(x_train)
print(y_train)

# test dataset
x_test = np.array([[4.6,3.0,1.5,0.2], [6.2, 3.0, 4.1, 1.2]])

# Standardizing the training and testing dataset
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
x_train_sc = scaler.fit_transform(x_train)
x_test_sc = scaler.fit_transform(x_test)

print(f"x_train_sc: \n"
      f"{x_train_sc}\n"
      f"x_test_sc: \n"
      f"{x_test_sc}\n")

# Q2) II) i
# KNN with number of nearest neighbours = 2
# Extracting the sepal features
print("Training KNN model with sepal data\n\n")
sepal_x_train_sc = x_train_sc[:,:2]

# Finding the 2 nearest neighbors for plant 1
from sklearn.neighbors import NearestNeighbors
neighbours = NearestNeighbors(n_neighbors = 2)
neighbours.fit(sepal_x_train_sc)

sepal_x_test_sc = x_test_sc[:1,:2]

# Indices of the 2 nearest neighbours
indices = neighbours.kneighbors(sepal_x_test_sc,
                                return_distance = False)[:,:][0]

columns = ['sepal.length','sepal.width','petal.length','petal.width']
neigh_data = pd.DataFrame([x_train[indices[0],:],x_train[indices[1],:]],
                          columns=columns)
neigh_data['species_label'] = [y_train[indices[0]],y_train[indices[1]]]
neigh_data['indices'] = indices

# Function to find the species based on the species label
# that is given in the iris documentation
def find_species(pred):
    species = []
    for i in pred:
        if i == 0:
            species.append('Setosa')
        elif i == 1:
            species.append('Versicolor')
        else:
            species.append('Virginica')

    return species
neigh_data['species'] = find_species(neigh_data['species_label'])

print(neigh_data)

# Q2) II) ii
# Training the KNN model using only sepal data
from sklearn.neighbors import KNeighborsClassifier

knn = KNeighborsClassifier(n_neighbors = 5)
knn.fit(sepal_x_train_sc, y_train)

probability = []
for plant in sepal_x_train_sc:
    #print(plant)
    probability.append(knn.predict_proba([plant])[0])

probability = pd.DataFrame(probability,
                           columns = ['Setosa','Versicolor','Virginica'])
print(probability)

# Predicting the values of the 2 unknown plants using sepal
# values only
pred = knn.predict(x_test_sc[:,:2]) # Taking only the sepal values
x_test_df = pd.DataFrame(x_test[:,:2],
                         columns = ['sepal.length','sepal.width'])
x_test_df['predicted_label'] = pred
x_test_df['Species'] = find_species(pred)

print(x_test_df)

######################################################################################################################

# Q2) III) i
# KNN with number of nearest neighbours = 2
# Extracting the petal features
print("\n\nTraining KNN model with petal data\n")
petal_x_train_sc = x_train_sc[:,2:4]

# Finding the 2 nearest neighbors for plant 1
from sklearn.neighbors import NearestNeighbors
neighbours = NearestNeighbors(n_neighbors = 2)
neighbours.fit(petal_x_train_sc)

petal_x_test_sc = x_test_sc[:1,2:4]

indices = neighbours.kneighbors(petal_x_test_sc,
                                return_distance = False)[0]


columns = ['sepal.length','sepal.width','petal.length','petal.width']
neigh_data = pd.DataFrame([x_train[indices[0],:],x_train[indices[1],:]],
                          columns=columns)
neigh_data['species_label'] = [y_train[indices[0]],y_train[indices[1]]]
neigh_data['indices'] = indices
neigh_data['species'] = find_species(neigh_data['species_label'])

print(neigh_data)

# Q2) III) ii
# Training the KNN model using only petal data
from sklearn.neighbors import KNeighborsClassifier

knn = KNeighborsClassifier(n_neighbors = 5)
knn.fit(petal_x_train_sc, y_train)

# Finding the probability that each plant will
# categorize into each species
probability = []
for plant in petal_x_train_sc:
    #print(plant)
    probability.append(knn.predict_proba([plant])[0])

probability = pd.DataFrame(probability,
                             columns = ['Setosa','Versicolor','Virginica'])
print('\n')
print(probability)


# Predicting the values of the 2 unknown plants using petal
# values only
pred = knn.predict(x_test_sc[:,2:4]) # Taking only the petal values
x_test_df = pd.DataFrame(x_test[:,2:4],
                         columns = ['petal.length','petal.width'])
x_test_df['predicted_label'] = pred
x_test_df['Species'] = find_species(pred)

print(x_test_df)

######################################################################################################################

# Q2) IV) i
# KNN with number of nearest neighbours = 2
print("\n\nTraining KNN model with both sepal and petal data\n")


# Finding the 2 nearest neighbors for plant 1
from sklearn.neighbors import NearestNeighbors
neighbours = NearestNeighbors(n_neighbors = 2)
neighbours.fit(x_train_sc)

indices = neighbours.kneighbors(x_test_sc,
                                return_distance = False)[:,:][0]

columns = ['sepal.length','sepal.width','petal.length','petal.width']
neigh_data = pd.DataFrame([x_train[indices[0],:],x_train[indices[1],:]],
                          columns=columns)
neigh_data['species_label'] = [y_train[indices[0]],y_train[indices[1]]]
neigh_data['indices'] = indices
neigh_data['species'] = find_species(neigh_data['species_label'])

print(neigh_data)

# Q2) IV) ii
# Training the KNN model using both sepal and petal data
from sklearn.neighbors import KNeighborsClassifier

knn = KNeighborsClassifier(n_neighbors = 5)
knn.fit(x_train_sc, y_train)

# Finding the probability that each plant will
# categorize into each species
probability = []
for plant in x_train_sc:
    #print(plant)
    probability.append(knn.predict_proba([plant])[0])

probability = pd.DataFrame(probability,
                             columns = ['Setosa','Versicolor','Virginica'])
print('\n')
print(probability)


# Predicting the values of the 2 unknown plants using petal
# values only
pred = knn.predict(x_test_sc) # Taking only the petal values
x_test_df = pd.DataFrame(x_test,
                         columns = ['sepal.length','sepal.width',
                                    'petal.length','petal.width'])
x_test_df['predicted_label'] = pred
x_test_df['Species'] = find_species(pred)

print(x_test_df)








