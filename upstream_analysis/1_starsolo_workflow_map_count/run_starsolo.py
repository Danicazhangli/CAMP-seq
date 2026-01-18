import os, sys
import subprocess
import itertools
import shutil
from random import random
from time import sleep
from joblib import Parallel, delayed
from pathlib import Path

# obtain parameters from snakemake
from snakemake.script import snakemake
threads = snakemake.threads
dry_run = snakemake.config.get('dry_run', True)
sample_name = snakemake.wildcards.sample_name
strand = snakemake.params.strand
chem = snakemake.params.chem
starindex = snakemake.params.starindex
extra_whitelist = snakemake.params.get('extra_whitelist', None)
extra_whitelist = Path(extra_whitelist).absolute() if extra_whitelist else None
result_whitelist1 = Path(snakemake.input.whitelist1).absolute()
result_whitelist2 = Path(snakemake.input.whitelist2).absolute()
result_whitelist3 = Path(snakemake.input.whitelist3).absolute()
readsFiles = Path(snakemake.input.readsFiles).absolute()
R1R2_StarSolo = Path(snakemake.output.R1R2_StarSolo).absolute()


readsFile_R1_patched = readsFiles / f"{sample_name}_R1.patched.fastq.gz"
readsFile_R2_patched = readsFiles / f"{sample_name}_R2.patched.fastq.gz"
readsFile_R1_unknown = readsFiles / f"{sample_name}_R1.unknown.fastq.gz"
readsFile_R2_unknown = readsFiles / f"{sample_name}_R2.unknown.fastq.gz"
R1R2_StarSolo_path = R1R2_StarSolo.parent.absolute()

solo_params = {
    'CPUS': threads,
    'RAM': 100000000000,
    # 'BAM': "--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN sS sQ sM"
    # 'BAM': "--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
    'CBLEN': 12,  # CR V3+: CBLEN=16 UMILEN=12; CR V2: CBLEN=16 UMILEN=10; CR V1: CBLEN=14 UMILEN=10
    'UMILEN': 8,
    'BAM': "--outSAMtype None",
    'soloMultiMappers': "EM",  # Unique, Uniform, Rescue, PropUnique, EM
    'UNMAPPED': "",  # "--outReadsUnmapped Fastx"
    'GZIP': "--readFilesCommand zcat",
    'STRAND': strand,
    'soloCBposition': '0_0_0_6 0_7_0_18 0_25_0_36 0_43_0_54',
    'soloUMIposition': '0_55_0_62',
    'clipAdapterType': "Hamming",  # Hamming, CellRanger4, None
    'soloType': 'CB_UMI_Complex',  # None, CB_UMI_Simple, CB_UMI_Complex, CB_samTagOut, SmartSeq
    'soloAdapterMismatchesNmax': 1,  # int>0: maximum number of mismatches allowed in adapter sequence
    'soloCBmatchWLtype': "1MM",  # Exact, 1MM, 1MM_multi, 1MM_multi_pseudocounts, 1MM_multi_Nbase_pseudocounts, EditDist_2
    'soloUMIdedup': "1MM_All",  # 1MM_All, 1MM_Directional_UMItools, 1MM_Directional, Exact, NoDedup, 1MM_CR
    'soloUMIfiltering': "-",  # -, MultiGeneUMI, MultiGeneUMI_All, MultiGeneUMI_CR
    'soloCellFilter': "CellRanger2.2 3000 0.99 10",
    'outFilterScoreMin': 0,
    'outFilterScoreMinOverLread': 0.33,  # 0.33 0.66
    'outFilterMatchNminOverLread': 0.33,  # 0.33 0.66
    'chimSegmentMin': 5  # 5 0
}

def starsolo_task(task, solo_params):
    if Path(task['OUTPUT_DIR']).exists():
        shutil.rmtree(Path(task['OUTPUT_DIR']), ignore_errors=True)

    Path(task['OUTPUT_DIR']).mkdir(exist_ok=True, parents=True)

    command0 = f"STAR --runThreadN {solo_params['CPUS']} \
        --outFileNamePrefix {str(task['OUTPUT_DIR']) + '/'} \
        --genomeDir {str(task['INDEX_DIR'])} \
        --readFilesIn {str(task['FASTQ_R2'])} {str(task['FASTQ_R1'])} \
        --runDirPerm \"All_RWX\" {solo_params['GZIP']} {solo_params['BAM']} \
        --soloCBwhitelist {str(task['WHITELIST_PATH'])} \
        --soloType {solo_params['soloType']} \
        --soloAdapterMismatchesNmax {solo_params['soloAdapterMismatchesNmax']} \
        --soloStrand {solo_params['STRAND']} \
        --soloUMIdedup {solo_params['soloUMIdedup']} \
        --soloCBmatchWLtype {solo_params['soloCBmatchWLtype']} \
        --soloUMIfiltering {solo_params['soloUMIfiltering']} \
        --soloCellFilter {solo_params['soloCellFilter']} \
        --clipAdapterType {solo_params['clipAdapterType']} \
        --outFilterScoreMin {solo_params['outFilterScoreMin']} \
        --outFilterScoreMinOverLread {solo_params['outFilterScoreMinOverLread']} \
        --outFilterMatchNminOverLread {solo_params['outFilterMatchNminOverLread']} \
        --chimSegmentMin {solo_params['chimSegmentMin']} \
        --soloFeatures Gene GeneFull \
        --soloMultiMappers {solo_params['soloMultiMappers']} {solo_params['UNMAPPED']} \
        --limitBAMsortRAM {solo_params['RAM']}"

    if solo_params['soloType'] == 'CB_UMI_Simple':
        command0 += f" --soloBarcodeReadLength 0 \
            --soloCBstart 1 \
            --soloCBlen {solo_params['CBLEN']} \
            --soloUMIstart {solo_params['CBLEN'] + 1} \
            --soloUMIlen {solo_params['UMILEN']}"
    elif solo_params['soloType'] == 'CB_UMI_Complex':
        command0 += f" --soloCBposition {solo_params['soloCBposition']} \
        --soloUMIposition {solo_params['soloUMIposition']}"

    command0 += f" > {task['OUTPUT_DIR']}/process_log.txt"

    try:
        _ = subprocess.check_output(command0, shell=True, stderr=subprocess.STDOUT)

        aligned_bam = Path(task['OUTPUT_DIR']) / 'Aligned.sortedByCoord.out.bam'
        if aligned_bam.exists():
            command1 = f"samtools index --threads {solo_params['CPUS']} {aligned_bam}"
            _ = subprocess.check_output(command1, shell=True, stderr=subprocess.STDOUT)
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error: {e.output.decode('utf-8')}")
        return False

if dry_run:
    sleep(random() * 2)
    R1R2_StarSolo.touch()
elif chem == "DropSeq":
    solo_params['soloType'] = 'CB_UMI_Simple'

    solo_params['CBLEN'] = 12
    solo_params['UMILEN'] = 8
    task = {
        'OUTPUT_DIR': R1R2_StarSolo_path / "starsolo_outputs",
        'FASTQ_R1': readsFile_R1_patched,
        'FASTQ_R2': readsFile_R2_patched,
        'INDEX_DIR': starindex,
        'WHITELIST_PATH': "None"
    }

    if starsolo_task(task, solo_params):
        R1R2_StarSolo.touch()

elif chem == "CapitalbioSeq":
    solo_params['soloType'] = 'CB_UMI_Complex'

    solo_params["soloCBposition"] = "0_0_0_11 0_18_0_29 0_36_0_47"
    solo_params["soloUMIposition"] = "0_48_0_55"
    task = {
        'OUTPUT_DIR': R1R2_StarSolo_path / "starsolo_outputs",
        'FASTQ_R1': readsFile_R1_patched,
        'FASTQ_R2': readsFile_R2_patched,
        'INDEX_DIR': starindex,
        'WHITELIST_PATH': " ".join([str(result_whitelist1), str(result_whitelist2), str(result_whitelist3)])
    }

    if starsolo_task(task, solo_params):
        R1R2_StarSolo.touch()
elif chem == "CapitalbioSeq-CB4":
    if readsFile_R1_unknown.exists() and readsFile_R2_unknown.exists():
        solo_params['soloType'] = 'CB_UMI_Complex'

        solo_params["soloCBposition"] = "0_0_0_11 0_18_0_29 0_36_0_47"
        solo_params["soloUMIposition"] = "0_48_0_55"
        task = {
            'OUTPUT_DIR': R1R2_StarSolo_path / "starsolo_outputs_UNKNOWN",
            'FASTQ_R1': readsFile_R1_unknown,
            'FASTQ_R2': readsFile_R2_unknown,
            'INDEX_DIR': starindex,
            'WHITELIST_PATH': " ".join([str(result_whitelist1), str(result_whitelist2), str(result_whitelist3)])
        }
        run_unknown_done = starsolo_task(task, solo_params)

        solo_params["soloCBposition"] = "0_0_0_6 0_7_0_18 0_25_0_36 0_43_0_54"
        solo_params["soloUMIposition"] = "0_55_0_62"
        task = {
            'OUTPUT_DIR': R1R2_StarSolo_path / "starsolo_outputs_CB4",
            'FASTQ_R1': readsFile_R1_patched,
            'FASTQ_R2': readsFile_R2_patched,
            'INDEX_DIR': starindex,
            'WHITELIST_PATH': " ".join([str(extra_whitelist), str(result_whitelist1), str(result_whitelist2), str(result_whitelist3)])
        }
        run_cb4_done = starsolo_task(task, solo_params)

        if run_unknown_done and run_cb4_done:
            R1R2_StarSolo.touch()
