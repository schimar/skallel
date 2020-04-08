# minimal example for PCA

# In this script, we do the following:
#       - load in the variants from a zarr path (and optionally create this from a vcf first),
#       - have a very quick look at the zarr data structure
#       - create a GenotypeArray from the genotypes in the zarr path
#       - count alleles for all variants and for each population, respectively (here, we only need the allele counts for all, but to see how this can be done for each pop in one step (which you would need to calculate e.g. Fst)
#       - subset only the alleles which are segregating across all individuals (since PCA doesn't like rows with no variance)   (this gives us another GenotypeArray)
#       - create a matrix (or better: a 2-dimensional numpy array) with the number of alt alleles for each locus and individual
#       - have a quick look at LD, subsample loci, and perform LD pruning
#       - perform several PCAs:
#                               Singular Value Decomposition (SVD) with Patterson's variance scaling for the LD-pruned allele counts
#                               SVD + LD pruning & Patterson's scaling for a random subset of 100000 loci
#                               SVD + Patterson's scaling without LD pruning (i.e. for all loci)
#                               SVD + LD pruning, without Patterson's scaling
#                               randomized PCA + LD pruning + Patterson's scaling
#                               (and link to PCA with even numbers of samples and evaluation of lower-level PCs)


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

# get the data
## (I'd suggest converting to zarr)

# subs = al.read_vcf('subs85f1byPop.vcf', numbers= {'GT': 2, 'ALT': 1}, fields= '*')

subsvcfPath = 'subs85f1byPop.vcf'
subszarrPath = 'subs85f1byPop.zarr'

# convert from vcf to zarr
al.vcf_to_zarr(subsvcfPath, subszarrPath, fields='*', log=sys.stdout, overwrite=True)

# read data
# NOTE: in the future, you can simply do this, which is much faster and memory-efficient
subs = zarr.open_group(subszarrPath, mode='r')

# see what's in there:
    subs.tree(expand= True)
    # or, to show the main keys of the dict:
    sorted(subs.keys())




## read id file and perform some basic checks   (NOTE: this is a space-delimited file with columns [ids, nest, pops])
ids = pd.read_csv('tapopsOrd.csv', delimiter= ' ')


    # count number of inds per pop and nest
    ids.groupby(by= ['nest', 'pops']).count()

    # is the order of samples in ids and subs the same?
    np.all(list(subs['samples']) == ids['id'].values)

    # if not, then we can add a column to ids with the index in which they appear in subs:
    samples = list(subs['samples'])
    subsIndex = [samples.index(s) for s in ids['id']]
    ids['subsIndex'] = subsIndex



###############
# create genotype array

gtsub = al.GenotypeArray(subs['calldata/GT'])

    # NOTE: you can simply get a boolean 2d-array of inds by heterozygous loci:
    gtsub.is_het()


# get id indices (for all, and the three respective pops [A, N, S])
# this is a dictionary, with the three keys [all, A, N, S] and values are lists of indices for each respective id
subpops = {
    'all': list(range(len(ids))),
    'A': ids[ids.pops == 'A'].index.tolist(),
    'N': ids[ids.pops == 'N'].index.tolist(),
    'S': ids[ids.pops == 'S'].index.tolist(),
}


# count alleles
# this is also a dict, with the same keys as subpops [all, A, N, S] and each of the four 'values' is an AlleleCountsArray
ac_subpops_subs = gtsub.count_alleles_subpops(subpops, max_allele=1)

    # how many biallelic singletons
    np.count_nonzero((ac_subpops_subs['all'].max_allele() == 1) & ac_subpops_subs['all'].is_singleton(1))


# locate segregating alleles and subset (which is a np 1d array - see segAll_subs.shape)
segAll_subs = ac_subpops_subs['all'].is_segregating()[:]

    # check segregating alleles per pop:
    for pop in ac_subpops_subs.keys():
        print(pop, ac_subpops_subs[pop].count_segregating())


# subset loci with segregating alleles
# NOTE: we here create a new GenotypeArray
gtseg_subs = gtsub.compress(segAll_subs, axis=0)


# create matrix (or rather a np 2d array) with loci x individuals containing icontaining the number of alt alleles (0, 1 or 2)
nAltSub = gtseg_subs.to_n_alt()

# check for LD (first 1000 loci)
def plot_ld(gn, title):
    m = al.rogers_huff_r(gn) ** 2
    ax = al.plot_pairwise_ld(m)
    ax.set_title(title)

plot_ld(nAltSub[:1000], 'Figure 1. Pairwise LD.')


# random subsampling of loci
n = 100000  # number of SNPs to choose randomly
vidx = np.random.choice(nAltSub.shape[0], n, replace=False)
vidx.sort()
gnr = nAltSub.take(vidx, axis=0)

plot_ld(gnr[:1000], 'Figure 2. Pairwise LD after random downsampling.')



# LD pruning    (i.e. sliding a window along the data, computing pairwise LD between all SNPs within each window, then removing one SNP from each correlated pair)
def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = al.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn

gnu = ld_prune(nAltSub, size=200, step=50, threshold=.1, n_iter=5)

plot_ld(gnu[:1000], 'Figure 3. Pairwise LD after LD pruning.')



###############

# PCA using Singular Value Decomposition (SVD) and Patterson's scaling on LD-pruned data (see gnu.shape for dimensions)
coords1, model1 = al.pca(gnu, n_components=10, scaler='patterson')

populations = ids.pops.unique()
pop_colours = {
    'A': sns.color_palette()[0],
    'N': sns.color_palette()[1],
    'S': sns.color_palette()[2],
}


def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop],
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))


def fig_pca(coords, model, title, sample_population=None):
    if sample_population is None:
        sample_population = ids.pops.values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()

fig_pca(coords1, model1, 'Figure 4. Conventional PCA.')




# pca without LD pruning for the random subset of 100k loci with Patterson's scaling
coords2, model2 = al.pca(gnr, n_components=10, scaler='patterson')
fig_pca(coords2, model2, 'Figure 5. Conventional PCA without LD pruning.')

# now for the full set (gtseg) with Patterson's scaling
#   NOTE: probably do not run this on your laptop
coords2all, model2all = al.pca(nAltSub, n_components=10, scaler='patterson')
fig_pca(coords2all, model2all, 'Conventional PCA without LD pruning.')

# pca + LD-pruning, without Patterson's scaling
coords3, model3 = al.pca(gnu, n_components=10, scaler=None)
fig_pca(coords3, model3, 'Figure 6. Conventional PCA without variance scaling.')

# randomized PCA with LD-pruning and Patterson's scaling
coords5, model5 = al.randomized_pca(gnu, n_components=10, scaler='patterson')
fig_pca(coords5, model5, 'Figure 8. Randomized PCA.')

# pca with even sample sizes NOTE: not really needed here, see alimanfoo's Fast PCA site
# (https://alimanfoo.github.io/2015/09/28/fast-pca.html)
# also: see alimanfoo's Fast-PCA post on an evaluation of the lower PCs in randomized PCA


