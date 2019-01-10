import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from sklearn.mixture import GaussianMixture

# read data, print first 5 to check
data_path = '/Users/durand.dc/Documents/works/Misc.Projects/wen/rawdata.csv'
rawDf = pd.read_csv(data_path, header=0)
print(rawDf.head())

# Find best fit for most obvious distribution

# data preparation
oneDataSet = rawDf.Qc.dropna()

# figure preparation
fig, ax = plt.subplots(1, 1)
# plot the origional distribution
sns.distplot(oneDataSet)

# fit multiple normal (gaussian) distribution (mixture)
oneSeries = oneDataSet.values[:, np.newaxis]
'''
# Fit function needs two axis:
>>> print(oneDataSet.values)
[28.434 31.43  31.426 ... 32.749 30.148 32.601]
>>> print(oneSeries)
[[28.434]
 [31.43]
 [31.426]
 ...
 [32.749]
 [30.148]
 [32.601]]
'''
max = oneSeries.max()
min = oneSeries.min()
N = 6
models = []  # model container
for i in range(1, N + 1):  # get all fits with 1 to N components
    models.append(GaussianMixture(n_components=i,
                                  init_params='kmeans').fit(oneSeries))
# get quality value of all fits:
# AIC, the lower the better
AIC = [m.aic(oneSeries) for m in models]
bestN = np.argmin(AIC)
bestGMM = models[bestN]

# prepare data for plot the fitted lines
# x axis data
x = np.linspace(min, max, 1000)
# weighted log probabilities:
logprob = bestGMM.score_samples(x[:, np.newaxis])
# pdf = Probability Density Functions
pdf = np.exp(logprob)
# Predict posterior probability of each component given the data.
responsibilities = bestGMM.predict_proba(x[:, np.newaxis])
pdf_individual = responsibilities * pdf[:, np.newaxis]

ax.plot(x, pdf, linestyle='-', color='k')
ax.plot(x, pdf_individual, linestyle='--', color='k')
ax.set_title('Fit of core alignment 4Q')


plt.show()

# -------------------------------------------------------------------------------
# Plotting all data

fig, axs = plt.subplots(3, 2, figsize=(8, 9))
titles = {
    'Eh': 'Hexagonal TODs alignment result 4E',
    'Qh': 'Hexagonal TODs alignment result 4Q',
    'Ec': 'Core alignment result 4E',
    'Qc': 'Core alignment result 4Q',
    'Et': 'Trimer alignment result 4E',
    'Qt': 'Trimer alignment result 4Q'
}
for col, ax in zip(rawDf.columns, axs.ravel()):
    ax.set_title(titles[col])
    oneDataSet = rawDf.loc[:, col].dropna()
    sns.distplot(oneDataSet, ax=ax)
    oneSeries = oneDataSet.values[:, np.newaxis]
    max = oneSeries.max()
    min = oneSeries.min()
    if col.endswith('c'):  # do multiple fit
        # fit multiple normal (gaussian) distribution (mixture)
        xlims = [24, 40]  # remove unwanted data in the fit
        oneDataSet = oneDataSet[(oneDataSet > xlims[0]) &
                                (oneDataSet < xlims[1])].dropna()
        oneSeries = oneDataSet.values[:, np.newaxis]
        n = (2 if col == 'Ec' else 3)
        model = GaussianMixture(n_components=n,
                                init_params='kmeans').fit(oneSeries)
        x = np.linspace(min, max, 1000)
        # weighted log probabilities:
        logprob = model.score_samples(x[:, np.newaxis])
        # pdf = Probability Density Functions
        pdf = np.exp(logprob)
        ax.plot(x, pdf, linestyle='-', color='k')
        # Predict posterior probability of each component given the data.
        responsibilities = model.predict_proba(x[:, np.newaxis])
        pdf_individual = responsibilities * pdf[:, np.newaxis]
        ax.plot(x, pdf_individual, linestyle='--', color='k')

    elif col.endswith('h'):
        # fit single normal (gaussian) distribution (mixture)
        xlims = [24, 45]  # remove unwanted data in the fit
        oneDataSet = oneDataSet[(oneDataSet > xlims[0]) &
                                (oneDataSet < xlims[1])].dropna()
        oneSeries = oneDataSet.values[:, np.newaxis]
        n = 1
        model = GaussianMixture(n_components=n,
                                init_params='kmeans').fit(oneSeries)
        x = np.linspace(min, max, 1000)
        # weighted log probabilities:
        logprob = model.score_samples(x[:, np.newaxis])
        # pdf = Probability Density Functions
        pdf = np.exp(logprob)
        ax.plot(x, pdf, linestyle='-', color='k')
        responsibilities = model.predict_proba(x[:, np.newaxis])
        pdf_individual = responsibilities * pdf[:, np.newaxis]
        ax.plot(x, pdf_individual, linestyle='--', color='k')
    elif col.endswith('t'):
        # fit single normal (gaussian) distribution (mixture)
        xlims = [12.5, 27.5]  # remove unwanted data in the fit
        oneDataSet = oneDataSet[(oneDataSet > xlims[0]) &
                                (oneDataSet < xlims[1])].dropna()
        oneSeries = oneDataSet.values[:, np.newaxis]
        n = 1
        model = GaussianMixture(n_components=n,
                                init_params='kmeans').fit(oneSeries)
        x = np.linspace(min, max, 1000)
        # weighted log probabilities:
        logprob = model.score_samples(x[:, np.newaxis])
        # pdf = Probability Density Functions
        pdf = np.exp(logprob)
        ax.plot(x, pdf, linestyle='-', color='k')
        responsibilities = model.predict_proba(x[:, np.newaxis])
        pdf_individual = responsibilities * pdf[:, np.newaxis]
        ax.plot(x, pdf_individual, linestyle='--', color='k')

    print(f'\n{"-"*100}')
    print(f'{titles[col]} fit {n} population result:')
    print('weights_            \n', model.weights_)
    print('means_              \n', model.means_.flatten())
    print('covariances_        \n', model.covariances_.flatten())
    print('precisions_         \n', model.precisions_.flatten())


plt.tight_layout()
plt.show()
