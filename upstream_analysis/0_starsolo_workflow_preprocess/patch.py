import shutil
import subprocess
import fastq as fq
import miniFasta as fa
import itertools
from joblib import Parallel, delayed
from random import random
from time import sleep
from pathlib import Path


# obtain parameters from snakemake
from snakemake.script import snakemake
sample_name = snakemake.wildcards.sample_name
chem = snakemake.params.chem
rt8codes_fasta = Path(snakemake.params.rt8codes_fasta).absolute()
data_path = Path(snakemake.input.R1R2_Cutadapt).parent.absolute()
R1_patched = Path(snakemake.output.R1_patched).absolute()
R2_patched = Path(snakemake.output.R2_patched).absolute()
snakemake_output_path = R1_patched.parent.parent
result_path = R1_patched.parent
log_file = result_path / f'patch_{sample_name}.log'
cleanup = snakemake.params.get('cleanup', False)
dry_run = snakemake.params.get('dry_run', True)


R1_unknown = result_path / f"{R1_patched.name.replace('.patched.fastq.gz', '.unknown.fastq.gz')}"
R2_unknown = result_path / f"{R2_patched.name.replace('.patched.fastq.gz', '.unknown.fastq.gz')}"
shutil.rmtree(result_path, ignore_errors=True)
result_path.mkdir(parents=True, exist_ok=True)


with open(log_file, 'w') as f:
    f.write(f'chem: {chem}\n')
    f.write(f'{sample_name}\n')

if dry_run:
    sleep(random() * 2)
    R1_unknown.touch()
    R2_unknown.touch()
    R1_patched.touch()
    R2_patched.touch()
elif chem != 'CapitalbioSeq-CB4':
    para_task_jobs = []
    para_task_jobs.append(delayed(shutil.move)(snakemake_output_path / f"{sample_name}_R1.fastq.gz", R1_patched))
    para_task_jobs.append(delayed(shutil.move)(snakemake_output_path / f"{sample_name}_R2.fastq.gz", R2_patched))
    multi_work = Parallel(n_jobs=-1, backend='loky')
    multi_work(para_task_jobs)
else:
    para_task_jobs = []
    para_task_jobs.append(delayed(shutil.move)(snakemake_output_path / f"{sample_name}_R1.fastq.gz", R1_patched.parent))
    para_task_jobs.append(delayed(shutil.move)(snakemake_output_path / f"{sample_name}_R2.fastq.gz", R2_patched.parent))
    multi_work = Parallel(n_jobs=-1, backend='loky')
    multi_work(para_task_jobs)

    # unknown
    unknown_R1_demultiplexing = (data_path / f"unknown_R1.demultiplexing.fastq.gz").absolute()
    unknown_R2_demultiplexing = (data_path / f"unknown_R2.demultiplexing.fastq.gz").absolute()
    para_task_jobs = []
    para_task_jobs.append(delayed(shutil.move)(unknown_R1_demultiplexing, R1_unknown))
    para_task_jobs.append(delayed(shutil.move)(unknown_R2_demultiplexing, R2_unknown))
    multi_work = Parallel(n_jobs=-1, backend='loky')
    multi_work(para_task_jobs)

    # append rt8codes to R1
    R1_files = []
    R2_files = []
    for f in data_path.glob('*_R1.demultiplexing.fastq.gz'):
        if "unknown" in f.name:
            continue
        R1_files.append(f)
        R2_files.append(f.parent / f.name.replace('_R1', '_R2'))

    adapter_map = {}
    for val in itertools.islice(fa.read(rt8codes_fasta), None):
        adapter_map[val.getHead()[1:]] = val.getSeq()

    def append_rt8codes(R1_file, adapter, insert_behind=-1):
        adapter_qstr = ''.join(['I' for _ in range(len(adapter))])
        R1_file_patched = f"{R1_file.parent / R1_file.name.removesuffix('.fastq.gz')}.patched.fastq.gz"

        commands = []
        if insert_behind == -1:
            commands.append(f"zcat {R1_file} | sed '2~4 s/^/{adapter}/;4~4 s/^/{adapter_qstr}/' | gzip > {R1_file_patched}")
        else:
            commands.append(f"zcat {R1_file} | sed '2~4 s/^\(.\{{{insert_behind}\}}\)/\\1{adapter}/;4~4 s/^\(.\{{{insert_behind}\}}\)/\\1{adapter_qstr}/' | gzip > {R1_file_patched}")

        try:
            for command in commands:
                _ = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e.output.decode('utf-8')}")

    tasks = []
    for R1_file in R1_files:
        adapter_key = '_'.join(R1_file.name.split('_')[:2])
        tasks.append(delayed(append_rt8codes)(R1_file, f"{adapter_map[adapter_key]}", insert_behind=-1))

    multi_work = Parallel(n_jobs=-1, backend='loky')
    multi_work(tasks)

    # merge R1 files and R2 files
    def task_job(fastq_paths, result_path):
        command = f"cat {' '.join([str(fastq_path.absolute()) for fastq_path in fastq_paths])} > {result_path}"
        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e.output.decode('utf-8')}")

    R1_files = []
    R2_files = []
    for f in data_path.glob('*_R2.demultiplexing.fastq.gz'):
        if "unknown" in f.name:
            continue
        R1_files.append(f.parent / f"{f.name.removesuffix('.fastq.gz').replace('_R2', '_R1')}.patched.fastq.gz")
        R2_files.append(f)

    para_task_jobs = []
    para_task_jobs.append(delayed(task_job)(R1_files, R1_patched))
    para_task_jobs.append(delayed(task_job)(R2_files, R2_patched))
    multi_work = Parallel(n_jobs=-1, backend='loky')
    multi_work(para_task_jobs)

if cleanup:
    Path(snakemake_output_path / f"{sample_name}_R1.fastq.gz").unlink(missing_ok=True)
    Path(snakemake_output_path / f"{sample_name}_R2.fastq.gz").unlink(missing_ok=True)

    for f in data_path.glob('*.fastq.gz'):
        f.unlink(missing_ok=True)
