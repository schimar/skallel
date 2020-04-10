#  minimal example for the estimation of diversity statistics from a vcf file

# In this script, I calculate different types of Fst estimators:



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


## first, we subset all variants with segregating alleles

subpops = {
    'all': list(range(len(ids))),
    'A': ids[ids.pops == 'A'].index.tolist(),
    'N': ids[ids.pops == 'N'].index.tolist(),
    'S': ids[ids.pops == 'S'].index.tolist(),
}

ac_subpops = gtsub.count_alleles_subpops(subpops, max_allele= 1)

segAll = ac_subpops['all'].is_segregating()[:]

gtseg = gtsub.compress(segAll, axis= 0)


########

## Weir & Cockerham's Fst pfor each locus

a, b, c, = al.weir_cockerham_fst(gtseg, list(subpops.values())[1:])

# estimate theta (a.k.a. Fst) for each variant & allele directly:
fst = a / (a + b + c)




# compare Hudson's and Weir & Cockerham's per locus Fst:

# only take variants that are segregating between the two pops
acu = al.AlleleCountsArray(ac_subpops['S'][:] + ac_subpops['N'][:])
flt = acu.is_segregating() & (acu.max_allele() == 1)
print('retaining', np.count_nonzero(flt), 'SNPs')

ac1 = al.AlleleCountsArray(ac_subpops['S'].compress(flt, axis=0)[:, :2])
ac2 = al.AlleleCountsArray(ac_subpops['N'].compress(flt, axis=0)[:, :2])

genotype = gtsub.compress(flt, axis=0)
#genotype


pop1_idx = subpops['S']
pop2_idx = subpops['N']
a, b, c = al.weir_cockerham_fst(genotype, subpops=[pop1_idx, pop2_idx], max_allele=1)
snp_fst_wc = (a / (a + b + c))[:, 0]

num, den = al.stats.fst.hudson_fst(ac1, ac2)
snp_fst_hudson = num / den
snp_fst_hudson

# plot the two Fst estimators
fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(ax=ax, offset=5)
ax.plot(snp_fst_hudson, snp_fst_wc, color='k', marker='.', linestyle=' ')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('Hudson $F_{ST}$')
ax.set_ylabel('Weir & Cockerham $F_{ST}$')
ax.set_title('%s (%s) vs %s (%s), SNP $F_{ST}$' % ('S', len(subpops['S']), 'N', len(subpops['N'])));

## block-jacknife Fst

# W & C's Fst with standard error obtained with block-jackknife
fst_wc, se_wc, vb_wc, _ = al.stats.fst.blockwise_weir_cockerham_fst(genotype, subpops=[subpops['S'], subpops['N']], blen=10000, max_allele=1)
print('%.04f +/- %.04f (Weir & Cockerham)' % (fst_wc, se_wc))

# A-S  = 0.1590 +/- 0.0053
# S-N  = 0.1816 +/- 0.0060
# N-A  = 0.1075 +/- 0.0052


# Hudson's Fst
fst, fst_se, _, _ = al.stats.fst.blockwise_hudson_fst(ac_segS, ac_segN, blen=100000)
print("Hudson's Fst: %.3f +/- %.3f (Hudson)" % (fst, fst_se))

# A-S  = 0.155 +/- 0.015
# S-N  = 0.178 +/- 0.018
# N-A  = 0.103 +/- 0.016

def plot_fst(ac1, ac2, pos, blen=2000):

    fst, se, vb, _ = al.stats.fst.blockwise_hudson_fst(ac1, ac2, blen=blen)

    # use the per-block average Fst as the Y coordinate
    y = vb

    # use the block centres as the X coordinate
    x = al.stats.window.moving_statistic(pos, statistic=lambda v: (v[0] + v[-1]) / 2, size=blen)

    # plot
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y, 'k-', lw=.5)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())


plot_fst(ac_segS, ac_segN, variants['POS'][:])
# NOTE: lengths between x & y differ!!!


