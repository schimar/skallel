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
sns.set_style('whitegrid')
sns.set_style('ticks')

######################################

## create VariantChunkedTable object
variants = al.VariantChunkedTable(subs['variants']) #, index= 'CHROM')


## count the number of variants per scaffold
scafs, scaf_counts = np.unique(variants['CHROM'], return_counts= True)
scafdf = pd.DataFrame({'scaffold': scafs, 'nVariables': scaf_counts})


## barplot of number of variants per scaffold
fig, ax = plt.subplots(figsize= (14,4))
ax.plot(np.arange(len(scafdf['nVariables'])), scafdf['nVariables'])
ax.set_xlabel('scaffolds')
ax.set_ylabel('count')
ax.set_title('Number of variants per scaffold')


## histogram of number of variables per scaffold
fig, ax = plt.subplots()
sns.distplot(scafdf['nVariables'], kde= True)   # kde= False for counts on y-axis, and not the ratio
ax.set_xlabel('number of variants')
ax.set_ylabel('proportion')
ax.set_title('Distribution of variants per scaffold')

## summary of number of variants per scaffold   (this is simiar to R's summary() function, except that it doesn't give you a median; use np.median(scafdf) for that)

scafdf.describe()


## count the number (and proportion) of heterozygous calls

gtsub.count_het(axis=0)     # axis (0 = across loci, 1 = across samples)

propHets = pd.Series(gtsub.count_het(axis= 0)/len(gtsub))



def plotPropHets(propHets, ids):
    # plot the proportion of heterozygous genotypes
    plt.subplots(figsize= (20,5))
    ax = sns.barplot(np.arange(len(propHets)), propHets, hue= ids.pops, dodge= False)
    ax.set_ylim([0,1])
    ax.set_xlabel('samples')
    ax.set_ylabel('proportion heterozygous')
    ax.set_title('proportion of heterozygous genotype calls')
    #ax.set_xticks(np.arange(len(propHets)))
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True)
    ax.set_xticklabels(ids.id, rotation= 40, ha= 'right', fontsize= 8)


plotPropHets(propHets, ids)


##################

vars0 = variants['CHROM' == 'scaffold_0'][:]

pos0 = vars0['POS']
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





###
gtsub.is_het()




