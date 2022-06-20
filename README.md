# Building Containerized Workflows

## Execution 

### Requirement

- snakemake 
```
conda activate /cluster/customapps/nexus/sharedutils/anaconda_envs/snakemake_v6.15.2/
```
- singularity 
```
snakemake --use-singularity --singularity-args '-B /cluster/work/nexus/' --cores 3
```
