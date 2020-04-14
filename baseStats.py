#  minimal example for a first inspection of our variants


#       scaffolds
#       variant density



#####################################

# import random     # to set the seed
#random.seed(42)
import sys
import allel as al
import numpy as np
#np.random.seed(42)
import pandas as pd
import zarr
#import scipy.parse
#import scipy.spatial
%matplotlib
import matplotlib.pyplot as plt
#import matplotlib as mpl
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')

######################################
variants = al.VariantChunkedTable(subs['variants'], index= 'CHROM')


## count (and plot) the number of variants per scaffold
scafs, scaf_counts = np.unique(variants['CHROM'], return_counts= True)
scafdf = pd.DataFrame({'scaffold': scafs, 'nVariables': scaf_counts})

#plt.figure(figsize= (16,4))
#sns.countplot(scafdf['nVariables'])
## barplot of number of variants per scaffold
fig, axes = plt.subplots(figsize= (14,4))
axes.plot(np.arange(len(scafdf['nVariables'])), scafdf['nVariables'])
axes.set_xlabel('scaffolds')
axes.set_ylabel('count')
axes.set_title('Number of variants per scaffold')


## histogram of number of variables per scaffold
fig, axes = plt.subplots()
sns.distplot(scafdf['nVariables'], kde= True)   # kde= False for counts on y-axis, and not the ratio
axes.set_xlabel('number of variants')
axes.set_ylabel('proportion')
axes.set_title('Distribution of variants per scaffold')

## summary of number of variants per scaffold

scafdf.describe()






###
gtsub.is_het()



# plot windoed variant density

def plot_windowed_variant_density(pos, window_size, title=None):

    # setup windows
    bins = np.arange(0, pos.max(), window_size)

    # use window midpoints as x coordinate
    x = (bins[1:] + bins[:-1])/2

    # compute variant density in each window
    h, _ = np.histogram(pos, bins=bins)
    y = h / window_size

    # plot
    fig, ax = plt.subplots(figsize=(12, 3))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y)
    ax.set_xlabel('Chromosome position (bp)')
    ax.set_ylabel('Variant density (bp$^{-1}$)')
    if title:
        ax.set_title(title)


pos = variants['POS'][:]

plot_windowed_variant_density(pos, window_size=100000, title='Raw variant density')


