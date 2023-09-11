# Nextflow Workflow for Estimating Mutation Age Using DMLE

The typical usage is
- Edit the `nextflow.config` file with your input parameters
- And run the workflow in the order
```
nextflow run getHaplotypeFrequencyFiles.nf
```
```
nextflow run createDmleInputFile.nf
```
```
nextflow run estimateMutationAge.nf
```

The results should be written to the output directory specified in the config file.
