# skallel
Set of scripts for specific tasks in scikit-allel 

Most scripts assume that you have already loaded your vcf file.  
This means that you'll already have the following objects loaded:
	- subs (zarr.hierarchy.Group)
	- ids (a pandas DataFrame, with individual ids (corresponding to the headers in your vcf) and population as columns 
	- gtsub (a GenotypeArray based on 'subs' ```gtsub = al.GenotypeArray(subs['calldata/GT'])```) 


	<br/>
I recommend converting your vcf to zarr for easier reading, which can be done as follows:
```

vcfPath = 'your_file.vcf'
subs = al.read_vcf('your_file.vcf', numbers= {'GT': 2, 'ALT': 1}, fields= '*')		# with all (*) fields read

zarrPath = 'your_file.zarr'
al.vcf_to_zarr(subsvcfPath, subszarrPath, fields='*', log=sys.stdout, overwrite=True)
```
So next time, you can quickly load the data from zarr:
```
subs = zarr.open_group(subszarrPath, mode='r')
```


You can easily inspect your data structure with the following:
```
subs.tree(expand= True)
```
