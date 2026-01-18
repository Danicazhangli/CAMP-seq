import shutil
import subprocess
from random import random
from time import sleep
from pathlib import Path

# obtain parameters from snakemake
from snakemake.script import snakemake
sample_name = snakemake.wildcards.sample_name
R1_file = Path(snakemake.input.R1_file).absolute()
R2_file = Path(snakemake.input.R2_file).absolute()
chem = snakemake.params.chem
rt8codes_fasta = Path(snakemake.params.rt8codes_fasta).absolute()
R1R2_Cutadapt_Done = Path(snakemake.output.R1R2_Cutadapt)
result_path = R1R2_Cutadapt_Done.parent.absolute()
log_file = result_path / f'demultiplexing_{sample_name}.log'
dry_run = snakemake.params.get('dry_run', True)

shutil.rmtree(result_path, ignore_errors=True)
result_path.mkdir(parents=True, exist_ok=True)


with open(log_file, 'w') as f:
    f.write(f'chem: {chem}\n')
    f.write(f'{R1_file.name}\t{R2_file.name}\n')

if dry_run:
    sleep(random() * 2)
    R1R2_Cutadapt_Done.touch()
elif chem != 'CapitalbioSeq-CB4':
    R1R2_Cutadapt_Done.touch()
else:
    action = 'trim'  # action: none, trim
    command = f'cutadapt \
        --cores 0 \
        --error-rate 1 \
        --no-indels \
        --action={action} \
        --front ^file:{rt8codes_fasta} \
        --quality-cutoff 10 \
        --output {result_path}/{{name}}_R2.demultiplexing.fastq.gz \
        --paired-output {result_path}/{{name}}_R1.demultiplexing.fastq.gz \
        {R2_file} \
        {R1_file}'
    try:
        print(f"Running command: {command}")
        result = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        with open(log_file, 'a') as f:
            f.write(f'{command}\n')
            f.write(result.decode('utf-8'))
        R1R2_Cutadapt_Done.touch()
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output.decode('utf-8')}")
