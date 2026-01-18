import re
import subprocess
import shutil
from random import random
from time import sleep
from joblib import Parallel, delayed
from pathlib import Path

# obtain parameters from snakemake
from snakemake.script import snakemake
sample_name = snakemake.wildcards.sample_name
R1_file = Path(snakemake.input.R1_file).absolute()
R2_file = Path(snakemake.input.R2_file).absolute()
pattern = snakemake.params.pattern
R1_merged = Path(snakemake.output.R1_merged).absolute()
R2_merged = Path(snakemake.output.R2_merged).absolute()
result_path = R1_merged.parent.absolute()
log_file = result_path / f'merge_fastq_{sample_name}.log'
dry_run = snakemake.params.get('dry_run', True)

if R1_file.is_dir():
    assert R1_file == R2_file
    R1_files = {}
    R2_files = {}
    for f in R1_file.iterdir():
        match = re.match(pattern, f.name)
        if match:
            read_order = match.group(1)
            read_type = match.group(2)
            if read_type == '1':
                R1_files[int(read_order)] = f
            elif read_type == '2':
                R2_files[int(read_order)] = f
    R1_files = sorted(R1_files.items(), key=lambda x: x[0])
    R2_files = sorted(R2_files.items(), key=lambda x: x[0])
    R1_files = [f[1] for f in R1_files]
    R2_files = [f[1] for f in R2_files]
else:
    R1_files = [R1_file]
    R2_files = [R2_file]

with open(log_file, 'w') as f:
    for f1p, f2p in zip(R1_files, R2_files):
        f.write(f'{f1p.name}\t{f2p.name}\n')

if dry_run:
    sleep(random() * 2)
    R1_merged.touch()
    R2_merged.touch()
else:
    if len(R1_files) == 1:
        para_task_jobs = []
        para_task_jobs.append(delayed(shutil.copyfile)(R1_files[0], R1_merged))
        para_task_jobs.append(delayed(shutil.copyfile)(R2_files[0], R2_merged))
        multi_work = Parallel(n_jobs=-1, backend='loky')
        multi_work(para_task_jobs)
    else:
        def task_job(fastq_paths, result_path):
            command = f"cat {' '.join([str(fastq_path.absolute()) for fastq_path in fastq_paths])} > {result_path}"
            try:
                print(f"Running command: {command}")
                subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)

            except subprocess.CalledProcessError as e:
                print(f"Error: {e.output.decode('utf-8')}")

        para_task_jobs = []
        para_task_jobs.append(delayed(task_job)(R1_files, R1_merged))
        para_task_jobs.append(delayed(task_job)(R2_files, R2_merged))
        multi_work = Parallel(n_jobs=-1, backend='loky')
        multi_work(para_task_jobs)
