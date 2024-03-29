import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Q1) I
# Import the iris data into a dataframe
data = pd.read_csv('iris.csv')
iris = pd.DataFrame(data)

print(iris)

# Partition test and training datasets
setosa_train = iris.loc[(iris['variety'] == 'Setosa')][0:int(0.9*49)]
versi_train = iris.loc[(iris['variety'] == 'Versicolor')][0:int(0.9*50)]
virgi_train = iris.loc[(iris['variety'] == 'Virginica')][0:int(0.9*50)]

setosa_test = iris.loc[(iris['variety'] == 'Setosa')][int(0.9*49)+1:-1]
versi_test = iris.loc[(iris['variety'] == 'Versicolor')][int(0.9*50)+1:-1]
virgi_test = iris.loc[(iris['variety'] == 'Virginica')][int(0.9*50)+1:-1]

train_data = pd.concat([setosa_train,versi_train,virgi_train],
                       ignore_index = True)
test_data = pd.concat([setosa_test,versi_test,virgi_test],
                      ignore_index = True)
print(train_data)
print(test_data)

# Q1) II
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder

# Drop any Nan entries in df columns
train = train_data.dropna().drop(['petal.length','petal.width','variety','v_short'],
                            axis = 1)
test = test_data.dropna().drop(['petal.length','petal.width','v_short'],
                            axis = 1)
print(train)

# Label encode categorical species data into numerical
# data in the test dataset
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder()
label = le.fit_transform(test['variety'])
test['Species'] = label
print(test)

# Scaling the data
scaler = StandardScaler()
train_sc = scaler.fit_transform(train)
test_sc = scaler.fit_transform(test[['sepal.length','sepal.width']])
print(train_sc)
print(test_sc)

# # Finding the optimal number of clusters
# # Parameters for KMeans
# km = {"init":"random", "n_init": 10, "random_state":1}
#
# sse = []
# for k in range (1, 11):
#
#     kmeans = KMeans(n_clusters = k, **km)
#     kmeans.fit(iris_sc)
#     sse.append(kmeans.inertia_)
#
# plt.plot(range(1,11), sse)
# plt.xticks(range(1,11))
# plt.xlabel("No. of Clusters")
# plt.ylabel("SSE")
# plt.show(block = False)
#
# Performing k-means with the found optimal no. of clusters
kmeans = KMeans(init = 'random', n_clusters = 3, n_init = 10,
                random_state = 1)
kmeans_1 = KMeans(init = 'random', n_clusters = 3, n_init = 10,
                random_state = 1)

# With training data
kmeans_scaled = kmeans.fit(train_sc)
kmeans = kmeans_1.fit(train)

train_data['Species_sc'] = kmeans_scaled.labels_
train_data['Species'] = kmeans.labels_
print(train_data)

print('KMeans cluster centers with scaled training data\n')
print(kmeans_scaled.cluster_centers_)
print('\nKMeans cluster centers with training data\n')
print(kmeans.cluster_centers_)

# With test data
pred_sc = kmeans_scaled.predict(test_sc)
pred = kmeans.predict(test[['sepal.length','sepal.width']])

test['Pred_sc'] = pred_sc
test['Pred'] = pred

print(test)

# Testing the accuracy and precision of the model
from sklearn.metrics import precision_score,accuracy_score

print(f"\n\nAccuracy (scaled data): {accuracy_score(test['Species'],test['Pred_sc'])}\n"
      f"Accuracy : {accuracy_score(test['Species'],test['Pred'])}\n"
      f"Precision (scaled data): {precision_score(test['Species'],test['Pred_sc'], average = 'weighted')}\n"
      f"Precision : {precision_score(test['Species'],test['Pred'],  average = 'weighted')}\n")



# Plot the cluster centers with the training data
# centers = kmeans.cluster_centers_
# centers = scaler.inverse_transform(centers)
# print(centers)
#
# center_df = pd.DataFrame(centers, columns = ['sepal.length','sepal.width'])
# print(center_df)


# Q1) III
# Predicting with new values
new_data = [[4.6,3.0], [6.2, 3.0]]

pred_1 = kmeans.predict(new_data)
print(pred_1)

#
# pred_df = pd.DataFrame(test_data, columns = ['sepal.length','sepal.width'])
# print(pred_df.iloc[:2])
#
#
#
# import seaborn as sns
#
# plt.figure()
# sns.scatterplot(data = iris, x = iris['sepal.length'], y = iris['sepal.width'],
#                 hue = iris['variety'], style = iris['variety'])
# sns.scatterplot(data = center_df, x = center_df['sepal.length'],
#                 y = center_df['sepal.width'], c = "red", marker = "^", s = 100)
# sns.scatterplot(data = pred_df, x = pred_df['sepal.length'],
#                 y = pred_df['sepal.width'], c = "purple", marker = "o",
#                 s = 50)
# plt.show()





