# skallel
This repo contains a set of scripts for specific tasks in scikit-allel. Please note that this is not intended to be a full documentation for scikit-allel, but simply a collection of scripts for specific purposes. For the full documentation of scikit-allel, please see the [official scikit-allel documentation](https://scikit-allel.readthedocs.io/en/stable/index.html). This code has been run and tested on Python version 3.7.5 with scikit-allel version 1.2.1. Some of this code comes from [Alistair Miles' github](alimanfoo.github.io), where part of it was changed to conform to scikit-allel version 1.2.1. Most of the scripts in here are kept really short, to keep the focus on a specific purpose. Please see the the relevant tutorials by Alistair Miles (which are mentioned in the scripts, where applicable), if you want more detailed tutorials. 
   
Most scripts assume that you have already loaded your vcf file.  
This means that you'll already have the following objects loaded:  
	- ***subs*** (zarr.hierarchy.Group object of your vcf)  
	- ***ids*** (a pandas DataFrame, with individual ids (corresponding to the headers in your vcf) and population as columns  
	- ***gtsub*** (a GenotypeArray based on 'subs')  


Some imports, needed for the below code to run:
``` 
import allel as al
import pandas as pd
import zarr
```
  
I recommend converting your vcf to zarr for easier reading, which can be done as follows:
```

vcfPath = 'your_file.vcf'
subs = al.read_vcf('your_file.vcf', numbers= {'GT': 2, 'ALT': 1}, fields= '*')		# with all (*) fields read

zarrPath = 'your_file.zarr'
al.vcf_to_zarr(subsvcfPath, subszarrPath, fields='*', log=sys.stdout, overwrite=True)
```
  
So next time, you can quickly load the data from zarr, without storing the full dataset in your memory:
```
subs = zarr.open_group(subszarrPath, mode='r')
```

Then, you create ***gtsub***:  
```
gtsub = al.GenotypeArray(subs['calldata/GT'])
```
  
And load ***ids***:
```
ids = pd.read_csv('your_ids.csv', delimiter= 'your delimiter')
```  


## Some useful things to do  

You can easily inspect the structure of ***subs*** with the following:
```
subs.tree(expand= True)
```

Further, check if the order of samples in ***ids*** and ***subs*** are the same:
```
np.all(list(subs['samples']) == ids['id'].values)
```
  
if not, then you can add a column to ids with the index in which they appear in subs:
```
samples = list(subs['samples'])
subsIndex = [samples.index(s) for s in ids['id']]
ids['subsIndex'] = subsIndex
```
