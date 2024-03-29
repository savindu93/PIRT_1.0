import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Q1) I
data = pd.read_csv('iris.csv')
iris = pd.DataFrame(data)

print(iris)

# Q1) II
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Drop any Nan entries in df columns
iris_1 = iris.dropna().drop(['petal.length','petal.width','variety','v_short'],
                            axis = 1)
print(iris_1)

# Scaling the data
# scaler = StandardScaler()
# iris_sc = scaler.fit_transform(iris_1)
# print(iris_sc)

# Finding the optimal number of clusters
# Parameters for KMeans
km = {"init":"random", "n_init": 10, "random_state":1}

sse = []
for k in range (1, 11):
    kmeans = KMeans(n_clusters = k, **km)
    kmeans.fit(iris_1)
    sse.append(kmeans.inertia_)

plt.plot(range(1,11), sse)
plt.xticks(range(1,11))
plt.xlabel("No. of Clusters")
plt.ylabel("SSE")
plt.show(block = False)

# Performing k-means with the found optimal no. of clusters
kmeans = KMeans(init = 'random', n_clusters = 3, n_init = 10,
                random_state = 1)

kmeans = kmeans.fit(iris_1)

iris['Species'] = kmeans.labels_
print(iris)

print(kmeans.cluster_centers_)

centers = kmeans.cluster_centers_
# centers = scaler.inverse_transform(centers)
# print(centers)

center_df = pd.DataFrame(centers, columns = ['sepal.length','sepal.width'])
print(center_df)


# Q1) III
# Predicting with new values
test_data = [[4.6,3.0], [6.2, 3.0]]
#test_sc = scaler.fit_transform(test_data)

pred = kmeans.predict(test_data)
print(f"Predited labels: {pred}")

# Function to find the species based on the species
# labels determined by Kmeans clustering
def find_species(predictions):
    species = []
    for i in pred:
        if i == 2:
            species.append('Setosa')
        elif i == 0:
            species.append('Versicolor')
        elif i == 1:
            species.append('Virginica')

    return species

pred_df = pd.DataFrame(test_data, columns = ['sepal.length','sepal.width'])
pred_df['predicted species'] = find_species(pred)
print(pred_df)


import seaborn as sns
# plt.figure()
# sns.scatterplot(data = iris, x = iris['sepal.length'], y = iris['sepal.width'],
#                 hue = iris['variety'])
# sns.scatterplot(data = center_df, x = center_df['sepal.length'],
#                 y = center_df['sepal.width'], c = "red", marker = "^", s = 50)
# sns.scatterplot(data = pred_df, x = pred_df['sepal.length'],
#                 y = pred_df['sepal.width'], c = "purple", marker = "o",
#                 s = 50)
# plt.show()

# Q1) IV
# PLotting the scatterplot with the relevant data
fig = plt.figure()
ax = fig.add_subplot()
x = iris['sepal.length']
y = iris['sepal.width']

# The color for the dataplots of each species is
# based on the cluster label determined by the
# Kmeans clustering algorithm
colors = iris['Species']
# Scatterplot with the sepal length and sepal width
ax.scatter(x, y, s = 20, c = colors)


# Annotating the scatterplot with labels determined
# by the kmeans model
annotations = iris['Species']
for xi, yi, text in zip(x, y, annotations):
    ax.annotate(text,
                xy = (xi,yi), xycoords = 'data',
                xytext = (-5,-5), textcoords = "offset points",
                fontsize = 7)

# Annotating the scatterplot with the species short form
# as labels from the original dataset
v_short = iris['v_short']
for xi, yi, text in zip(x, y, v_short):
    ax.annotate(text, xy = (xi, yi), xycoords = 'data',
                xytext = (5,17), textcoords = 'offset pixels',
                fontsize = 7,
                arrowprops = dict(facecolor = 'black',arrowstyle = '->'),
                horizontalalignment = 'left',
                verticalalignment = 'bottom')

# Plotting cluster centers
x = center_df['sepal.length']
y = center_df['sepal.width']
ax.scatter(x, y, s = 40, marker = '^', c = 'red')

# Plotting predicted values
x = pred_df['sepal.length']
y = pred_df['sepal.width']
ax.scatter(x, y, s = 40, c = 'green')

plt.show()

####################################################################################################################

# Carrying out Kmeans clustering with petal length and width data

# Q1) V
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Drop any Nan entries in df columns
iris_1 = iris.dropna().drop(['sepal.length','sepal.width','variety','v_short'],
                            axis = 1)
print(iris_1)

# Scaling the data
# scaler = StandardScaler()
# iris_sc = scaler.fit_transform(iris_1)
# print(iris_sc)

# Finding the optimal number of clusters
# Parameters for KMeans
km = {"init":"random", "n_init": 10, "random_state":1}

sse = []
for k in range (1, 11):
    kmeans = KMeans(n_clusters = k, **km)
    kmeans.fit(iris_1)
    sse.append(kmeans.inertia_)

plt.plot(range(1,11), sse)
plt.xticks(range(1,11))
plt.xlabel("No. of Clusters")
plt.ylabel("SSE")
plt.show(block = False)

# Performing k-means with the found optimal no. of clusters
kmeans = KMeans(init = 'random', n_clusters = 3, n_init = 10,
                random_state = 1)

kmeans = kmeans.fit(iris_1)

iris['Species'] = kmeans.labels_
print(iris)

print(kmeans.cluster_centers_)

centers = kmeans.cluster_centers_
# centers = scaler.inverse_transform(centers)
# print(centers)

center_df = pd.DataFrame(centers, columns = ['petal.length','petal.width'])
print(center_df)

# Predicting with new values
test_data = [[1.5,0.2], [4.1, 1.2]]
#test_sc = scaler.fit_transform(test_data)

pred = kmeans.predict(test_data)
print(f"Predited labels: {pred}")

pred_df = pd.DataFrame(test_data, columns = ['petal.length','petal.width'])

def find_species(predictions):
    species = []
    for i in pred:
        if i == 2:
            species.append('Setosa')
        elif i == 0:
            species.append('Versicolor')
        elif i == 1:
            species.append('Virginica')

    return species

pred_df['predicted species'] = find_species(pred)
print(pred_df)


import seaborn as sns
# plt.figure()
# sns.scatterplot(data = iris, x = iris['sepal.length'], y = iris['sepal.width'],
#                 hue = iris['variety'])
# sns.scatterplot(data = center_df, x = center_df['sepal.length'],
#                 y = center_df['sepal.width'], c = "red", marker = "^", s = 50)
# sns.scatterplot(data = pred_df, x = pred_df['sepal.length'],
#                 y = pred_df['sepal.width'], c = "purple", marker = "o",
#                 s = 50)
# plt.show()

# PLotting the scatterplot with the relevant data
fig = plt.figure()
ax = fig.add_subplot()
x = iris['petal.length']
y = iris['petal.width']
colors = iris['Species']
# Scatterplot with the sepal length and sepal width
ax.scatter(x, y, s = 20, c = colors)


# Annotating the scatterplot with labels determined
# by the kmeans model
annotations = iris['Species']
for xi, yi, text in zip(x, y, annotations):
    ax.annotate(text,
                xy = (xi,yi), xycoords = 'data',
                xytext = (-5,-5), textcoords = "offset points",
                fontsize = 7)

# Annotating the scatterplot with the species short form
# as labels from the original dataset
v_short = iris['v_short']
for xi, yi, text in zip(x, y, v_short):
    ax.annotate(text, xy = (xi, yi), xycoords = 'data',
                xytext = (5,17), textcoords = 'offset pixels',
                fontsize = 7,
                arrowprops = dict(facecolor = 'black',arrowstyle = '->'),
                horizontalalignment = 'left',
                verticalalignment = 'bottom')

# Plotting cluster centers
x = center_df['petal.length']
y = center_df['petal.width']
ax.scatter(x, y, s = 40, marker = '^', c = 'red')

# Plotting predicted values
x = pred_df['petal.length']
y = pred_df['petal.width']
ax.scatter(x, y, s = 40, c = 'green')

plt.show()




