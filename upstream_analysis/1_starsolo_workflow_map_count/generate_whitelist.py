import os, sys
import subprocess
import itertools
from random import random
from time import sleep
from snakemake.script import snakemake
from joblib import Parallel, delayed
from pathlib import Path

# obtain parameters from snakemake
from snakemake.script import snakemake
dry_run = snakemake.config.get('dry_run', True)
sample_name = snakemake.wildcards.sample_name
chem = snakemake.params.chem
readsFiles = Path(snakemake.input.readsFiles).absolute()
result_whitelist1_file = Path(snakemake.output.whitelist1).absolute()
result_whitelist2_file = Path(snakemake.output.whitelist2).absolute()
result_whitelist3_file = Path(snakemake.output.whitelist3).absolute()

result_whitelist1_path = result_whitelist1_file.parent.absolute()
result_whitelist2_path = result_whitelist2_file.parent.absolute()
result_whitelist3_path = result_whitelist3_file.parent.absolute()

subset_reads = 1000000000
bc_pattern1 = 'CCCCCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXNNNNNNNN'
bc_pattern2 = 'XXXXXXXXXXXXXXXXXXCCCCCCCCCCCCXXXXXXXXXXXXXXXXXXNNNNNNNN'
bc_pattern3 = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXCCCCCCCCCCCCNNNNNNNN'

def generate_whitelist(data_path, result_path, bc_pattern: str, subset_reads=1000000000):
    Path(result_path).mkdir(exist_ok=True)

    command = f'cd {result_path} && umi_tools whitelist \
        --stdin {data_path} \
        --bc-pattern={bc_pattern} \
        --method=reads \
        --plot-prefix=expect_whitelist \
        --subset-reads={subset_reads} \
        --log2stderr > whitelist.txt'

    try:
        result = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        with open(Path(result_path) / 'umitools_whitelist.log', 'w') as f:
            f.write(result.decode('utf-8'))

        input_file = result_path / "whitelist.txt"
        output_file = result_path / "whitelist_starsolo.txt"
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                first_word = line.split()[0]
                outfile.write(first_word + '\n')

    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output.decode('utf-8')}")

if dry_run:
    sleep(random() * 2)
    result_whitelist1_file.touch()
    result_whitelist2_file.touch()
    result_whitelist3_file.touch()
elif chem == "DropSeq":
    result_whitelist1_file.touch()
    result_whitelist2_file.touch()
    result_whitelist3_file.touch()
elif chem == "CapitalbioSeq":
    tasks = []
    readsFile_R1 = readsFiles / f"{sample_name}_R1.patched.fastq.gz"
    tasks.append(delayed(generate_whitelist)(readsFile_R1, result_whitelist1_path, bc_pattern1, subset_reads))
    tasks.append(delayed(generate_whitelist)(readsFile_R1, result_whitelist2_path, bc_pattern2, subset_reads))
    tasks.append(delayed(generate_whitelist)(readsFile_R1, result_whitelist3_path, bc_pattern3, subset_reads))
    multi_work = Parallel(n_jobs=-1, backend='loky')
    res = multi_work(tasks)
elif chem == "CapitalbioSeq-CB4":
    tasks = []
    readsFile_R1 = readsFiles / f"{sample_name}_R1.fastq.gz"
    tasks.append(delayed(generate_whitelist)(readsFile_R1, result_whitelist1_path, bc_pattern1, subset_reads))
    tasks.append(delayed(generate_whitelist)(readsFile_R1, result_whitelist2_path, bc_pattern2, subset_reads))
    tasks.append(delayed(generate_whitelist)(readsFile_R1, result_whitelist3_path, bc_pattern3, subset_reads))
    multi_work = Parallel(n_jobs=-1, backend='loky')
    res = multi_work(tasks)
