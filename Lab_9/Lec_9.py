import numpy as np
#
# arr_1 = np.array([1,2,3,4])
# print(arr_1)
#
# arr_2 = np.array([[1,2,3],[4,5,6],[0.7,0.8,0.9]])
# print(arr_2)
# print(arr_2.shape)
# print(arr_2[0:,])
# print(arr_2.dtype)
#
# arr_3 = np.array([1,2,3], dtype = int )
# print(arr_3)
# print(np.array(arr_3, dtype = float))
#
# r1 = np.random.rand(1000)
# print(r1)
#
# r2 = np.random.randn(5,2)
# print(r2)

import matplotlib.pyplot as plt

r3 = np.random.randn(100000)

# plt.hist(r3)
# plt.title('Histogram of Normal Population')
# plt.show()

# fig, axes = plt.subplots(1, 2, figsize = (12,3))
# axes[0].hist(r3)
# axes[0].set_title('Normal Pop')
# axes[1].hist(r1)
# axes[1].set_title('Random Pop')
# plt.show()

# import seaborn as sns
#
# fig, axes = plt.subplots(1, 2, figsize = (12,3))
# sns.histplot(r3, ax = axes[0], bins = 50, kde = True)
# axes[0].set_title('Normal Pop')
# sns.histplot(r1, ax = axes[1], bins = 50, kde = True)
# axes[1].set_title('Random Pop')
# plt.show()

import seaborn as sns
# import pandas as pd
#
# marks = pd.DataFrame([['Bhagya',97,78],
#                       ['Maheshi',95,85],
#                       ['Ruchila',88,98]])
# column_headers = ['Name','Chemistry','Biology']
#
# marks.columns = column_headers
#
# print(marks)
# print(marks.info())
# print(marks.describe())

# fig, axes = plt.subplots(figsize = (12,9))
# sns.boxplot(data = marks[['Chemistry','Biology']], ax = axes)
# plt.show()

#statistics
npop = np.random.normal(30, 5, 2000)

# fig, axes = plt.subplots(figsize = (9,4))
# sns.histplot(data = npop, ax=axes, kde =True)
# plt.show()

# from statsmodels.graphics.gofplots import qqplot
#
# fig, axes = plt.subplots(1, 2, figsize = (9,4))
# qqplot(npop, line = 's', ax= axes[0])
# axes[0].set_title('QQ plot for normal population')
# qqplot(r1, line = 's', ax = axes[1])
# axes[1].set_title('QQ plot for radnom population')
# plt.show()

from scipy import stats
#
# stat, p = stats.shapiro(npop)
# print('Stat val = %.3f, p = %.3f' %(stat,p))
#
# stat, p = stats.shapiro(r1)
# print('Stat val = %.3f, p = %.3f' %(stat,p))
#
mew = 20
t, p = stats.ttest_1samp(npop,mew)
print(f"Stat val = {t: .3f}, p = {p: .3f}")





