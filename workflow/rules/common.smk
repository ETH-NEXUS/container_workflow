import glob
import pandas as pd

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

sample_names = samples["sample"].tolist()

# Function to return input fastqs
# ToDo: Modify to make pipeline run when only R1 is specified
def get_fastq(wildcards):
    fastq_R1 = dict(zip(samples['sample'], samples['R1']))
    fastq_R2 = dict(zip(samples['sample'], samples['R2']))
    return {"fq1": f"{fastq_R1[wildcards.sample]}", "fq2": f"{fastq_R2[wildcards.sample]}"}

