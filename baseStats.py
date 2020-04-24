#  minimal example for a first inspection of our variants

#       Depth of Coverage
#       number of variants per scaffold
#       variant density (not yet working)



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
#sns.set_style('ticks')

######################################
## explore depth of coverage (DP) in your data

# if you wanted to get overall DP, you could access subs['variants/DP']
pd.Series(subs['variants/DP'][:]).describe()
# this gives you DP across all individuals

# if, however, you want to get DP for each individual, then access 'calldata' instead of the 'variants' tree
dpInds = pd.DataFrame(subs['calldata/DP'][:], columns= ids.id)

# get the basic stats for DP per individual
dpInds.describe()

# NOTE: you could use sns to plot (e.g.) the mean for each individual. For a relatively large DataFrame, however, this will be very slow
#sns.violinplot(dpInds, estimator= np.mean)
#sns.barplot(dpInds, estimator= np.mean)

# instead, we can first use the apply method with our favorite function, and then plot these
dpMu = np.array(dpInds.apply(np.mean))
dpSd = np.array(dpInds.apply(np.std))


def plotDP(dpMu, dpSd, ids):
    plt.subplots(figsize= (20,5))
    ax = sns.barplot(np.arange(len(dpMu)), dpMu, hue= ids.pops, dodge= False)
    ax.set_xlabel('samples')
    ax.set_ylabel('mean DP')
    ax.set_title('Depth of Coverage per individual')
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=True)
    ax.set_xticklabels(ids.id, rotation= 40, ha= 'right', fontsize= 8)
    ax.errorbar(np.arange(len(dpMu)), y= dpMu, yerr=dpSd, fmt= 'none', ecolor= 'grey')
    plt.tight_layout()

plotDP(dpMu, dpSd, ids)


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
    plt.tight_layout()

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




