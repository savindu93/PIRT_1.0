import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from statsmodels.graphics.gofplots import qqplot

# Q1) III
temp = pd.read_csv('Temperature.csv')
temp_df = pd.DataFrame(temp)
print(temp_df)

temp_info = temp_df.describe()
print(temp_info)

# Q1) IV
fig, axes = plt.subplots(1, 2, figsize = (12,9))
sns.histplot(data = temp_df[['temperature']], ax = axes[0], kde = True)
axes[0].set_title('Histogram of Body Temperature')
qqplot(temp_df['temperature'], line = 's', ax = axes[1])
axes[1].set_title('QQplot for Body Temperature')
plt.show()

# Normality test
stat, p = stats.shapiro(temp_df[['temperature']])
print(f"Stat val: {stat: .3f}, p value: {p: .3f}")

# Q1) V
# 1 sample t-test
pop_temp = 98.6  # Putative body temperature as taught in schools
t, p = stats.ttest_1samp(temp_df[['temperature']], pop_temp)
print(f"t-value: {t[0]: .3f}, p-value: {p[0]: .3f}")

#------------------------------------------------------------------------------

# Q2) III
length = pd.read_csv('HornedLizards.csv')
length_df = pd.DataFrame(length).dropna()

print(length_df)

# Horn lengths of surviving horned lizards
surv_length_df = length_df.loc[(length_df['Survive'] == 'survived')]
print(surv_length_df)

# Horn lengths of dead horned lizards
dead_length_df = length_df.loc[(length_df['Survive'] == 'dead')].reset_index(drop = True)
print(dead_length_df)

print(surv_length_df.describe())
print(dead_length_df.describe())

# Histogram and QQplots for the 2 samples
fig, axes = plt.subplots(2, 2, figsize = (12,9))
fig.tight_layout(pad = 5.0)
sns.histplot(length_df, x = 'Squamosal horn length', hue = 'Survive',
             ax = axes[0][0], kde = True)
axes[0][0].set_title('Histogram for Horn Lengths of Surviving/ Dead Horned Lizards')
qqplot(surv_length_df['Squamosal horn length'], line = 's', ax = axes[1][0])
axes[1][0].set_title('QQplot for Horned Lengths of Surviving Horned Lizards')
qqplot(dead_length_df['Squamosal horn length'], line = 's', ax = axes[1][1])
axes[1][1].set_title('QQplot for Horned Lengths of Dead Horned Lizards')
plt.show()

# Normality test
stat, p = stats.shapiro(surv_length_df[['Squamosal horn length']])
print(f"Surviving horned izards; Stat val : {stat: .3f}, p-value : {p: .3f}")

stat, p = stats.shapiro(dead_length_df[['Squamosal horn length']])
print(f"Dead horned lizards; Stat val : {stat: .3f}, p-value : {p: .3f}")

# Testing for homogeneity of variance
stat, p = stats.levene(surv_length_df[['Squamosal horn length']],
                       dead_length_df[['Squamosal horn length']])
print(f"Homogeneity of variance; Stat val : {stat[0]: .3f}, p-value : {p[0]: .3f}")

# 2 sample t-test
t, p = stats.ttest_ind(surv_length_df[['Squamosal horn length']],
                          dead_length_df[['Squamosal horn length']],
                          equal_var = True)
print(f"2 sample t-test; t : {t[0]: .3f}, p-value : {p[0]: .3f}")

# Non-parametric: Mann-Whitney test
stat, p = stats.mannwhitneyu(surv_length_df[['Squamosal horn length']],
                             dead_length_df[['Squamosal horn length']])
print(f"Mann-Whitney test; Stat val: {stat[0]: .3f}, p-value: {p[0]: .3f}")

# Boxplots and violin plots of 2 samples
fig = plt.plot(figsize = (12,9))
sns.boxplot(data = length_df,
            x = 'Survive',
            y = 'Squamosal horn length',
            hue = 'Survive')
plt.show()

fig = plt.plot(figsize = (12,9))
sns.violinplot(data = length_df,
               x = 'Survive',
               y = 'Squamosal horn length',
               hue = 'Survive')
plt.title('violinplot of Horn Lengths')
plt.show()

#------------------------------------------------------------------------------

# Q3) II
roa = pd.read_csv('BlackbirdTestosterone.csv')
roa_df = pd.DataFrame(roa)

roa_df = roa_df.drop(['Before','After','dif'], axis = 1)

print(roa_df)

print(roa_df['log before'].describe())
print(roa_df['log after'].describe())
print(roa_df['dif in logs'].describe())

# Q3) III
# Constructing the histogram and QQplot of log difference
fig, axes = plt.subplots(1, 2, figsize = (12,9))
sns.histplot(data = roa_df['dif in logs'], ax = axes[0], kde = True)
axes[0].set_title('Histogram of dif in logs')
qqplot(data = roa_df['dif in logs'], line = 's', ax = axes[1])
axes[1].set_title('QQplot of dif in logs')
plt.show(block = False)

# Normality test
stat, p = stats.shapiro(roa_df['dif in logs'])
print(f"Stat val: {stat: .3f}, p-value: {p: .3f}")

# Q3) IV
# Constructing boxplots and violoin plots of the 2 samples
roa_df = roa_df.drop(['dif in logs'], axis = 1)

roa_melt_df = pd.melt(roa_df, var_name = 'Implant State',
                      value_name = 'Log Antibody Production Rate')
print(roa_melt_df)

fig, axes = plt.subplots(1, 2, figsize = (12, 9))
sns.boxplot(x = 'Implant State', y = 'Log Antibody Production Rate',
            data = roa_melt_df, ax = axes[0])
axes[0].set_title('Boxplot of Log Antibody Production Rates')
sns.violinplot(x = 'Implant State', y = 'Log Antibody Production Rate',
            data = roa_melt_df, ax = axes[1])
axes[1].set_title('Violin plot of Log Antibody Production Rates ')
plt.show()

#Performing paired sample t-test
stat, p = stats.ttest_rel(roa_df['log before'],
                            roa_df['log after'],
                          alternative = 'less')
print(f"Stat val: {stat: .3f}, p-value: {p: .3f}")

#------------------------------------------------------------------------------

# Q4) III
# Creating the contingency table as a dataframe

data = [[1, 10, 37, 48],[49,35,9, 93],[50, 45, 46, 141]]
index = ['Eaten by birds','Not eaten by birds','Column total']
columns = ['Uninfected','Lightly infected','Highly infected','Row total']

df = pd.DataFrame(data, index = index, columns = columns)
print(df)

# Q4) iV
# Creating a mosaic plot with given values
df_1 = df.drop(['Row total'], axis = 1).drop(['Column total'], axis = 0)
print(df_1)

from statsmodels.graphics.mosaicplot import mosaic

mosaic(df_1.stack())
plt.show()

# Q4) V
# Performing chi-square contingency test
from scipy.stats import chi2_contingency

observed_val = np.array([df_1.iloc[0][0:3],
                         df_1.iloc[1][0:3]])

stat, p, dof, exp_freq = chi2_contingency(observed_val)
print(f"\nStat val: {stat: .3f}, p-value: {p: 0.3f}, Degree of freedom: {dof}\n")

# Q4) VI
# Obtaining expected values
exp_freq_df = pd.DataFrame(exp_freq,
                           index = index[0:2], columns = columns[0:3])

print(exp_freq_df)






