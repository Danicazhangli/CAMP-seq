import os
import h5py
import pickle
import typing as t
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as ssp
import ahocorasick

import scanpy as sc
import anndata as ad
from anndata import AnnData

from tqdm.auto import tqdm
from joblib import Parallel, delayed

class CB4_UNKNOWN_Merger:
    def __init__(self, h5ad_CB4: AnnData, h5ad_UNKNOWN: AnnData, prefix: str = "NNNNNNN", n_jobs: int = 24):
        assert (h5ad_CB4.var == h5ad_UNKNOWN.var).all(axis=1).all()

        self.cb4_var = h5ad_CB4.var
        self.cb4_obs = h5ad_CB4.obs
        self.cb4_X = h5ad_CB4.X.tocsc()
        self.unknown_obs = h5ad_UNKNOWN.obs
        self.unknown_X = h5ad_UNKNOWN.X.tocsc()

        self.cb4_idx = self.cb4_obs.index.to_list()
        self.unknown_idx = self.unknown_obs.index.to_list()

        self.prefix = prefix
        self.n_jobs = n_jobs
        if self.n_jobs < 1:
            self.n_jobs = os.cpu_count() or 1
        self.hits_table = self._build_hits_map()

    def run(self):
        results = self._parallel_process()
        merged_X, merged_idx = self._merge_results(results)

        merged_X.data = np.round(merged_X.data)
        h5ad_merged = ad.AnnData(merged_X)
        h5ad_merged_obs = pd.DataFrame(index=merged_idx)
        h5ad_merged_obs.index.name = 'barcode'
        h5ad_merged.obs = h5ad_merged_obs
        h5ad_merged.var = self.cb4_var

        return h5ad_merged

    def _build_hits_map(self):
        """Build Hits Map

        Returns:
            list[np.array[int]]: Hits Table
        """
        # 1) Build an automaton with UNKNOWN as the "model"
        AC = ahocorasick.Automaton()
        for pat_id, pat in enumerate(self.unknown_idx):
            AC.add_word(pat, pat_id)
        AC.make_automaton()

        # 2) Scan CB4 once and write the reverse matches into hits.
        hits = [list() for _ in range(len(self.unknown_idx))]
        for cb4_row, name in enumerate(self.cb4_idx):
            for _, pat_id in AC.iter(name):
                hits[pat_id].append(cb4_row)

        return [np.asarray(lst, dtype=int) for lst in hits]

    def _process_one(self, i):
        """Core func parallel processing

        Args:
            i (int): unknown_idx index

        Returns:
            tuple[rows, cols, vals, prefix, csc_matrix_row]: _description_
        """
        # unknown_name = self.unknown_idx[i]
        # hits = np.where([unknown_name in s for s in self.cb4_idx])[0]
        hits = self.hits_table[i]
        if hits.size:
            row_coo = self.unknown_X.getrow(i).tocoo()
            nnz = row_coo.nnz
            rows = np.repeat(hits, nnz)
            cols = np.tile(row_coo.col, hits.size)
            vals = np.tile(row_coo.data / hits.size, hits.size)
            return rows, cols, vals, None, None
        else:
            return (
                np.empty(0, dtype=int),
                np.empty(0, dtype=int),
                np.empty(0, dtype=self.unknown_X.dtype),
                f"{self.prefix}_{self.unknown_idx[i]}",
                self.unknown_X.getrow(i),
            )

    def _process_chunk(self, chunk_i):
        """Process Chunk of unknown_X

        Args:
            list[int]: chunk indexes

        Returns:
            list[func]: List of results from _process_one
        """
        results = []
        for i in chunk_i:
            results.append(self._process_one(i))
        return results

    def _parallel_process(self):
        """Process in parallel

        Returns:
            tuple[rows, cols, vals, prefix, csc_matrix_row]: _description_
        """
        unknown_X_index = np.arange(self.unknown_X.shape[0])
        chunk_indexes = np.array_split(unknown_X_index, self.n_jobs)
        chunk_results = Parallel(n_jobs=self.n_jobs, backend="loky")(
            delayed(self._process_chunk)(chunk_i) for chunk_i in chunk_indexes
        )
        results = []
        for chunk in chunk_results:
            results.extend(chunk)
        return results

    def _merge_results(self, results):
        """Merge results from parallel processing

        Args:
            tuple[rows, cols, vals, prefix, csc_matrix_row]: _description_
        """
        upd_rows, upd_cols, upd_vals = [], [], []
        extras_rows = []
        extras_idx = []
        cb4_nrow, cb4_ncol = self.cb4_X.shape

        for rows, cols, vals, ex_idx, ex_row in results:
            if rows.size:
                upd_rows.append(rows)
                upd_cols.append(cols)
                upd_vals.append(vals)
            if ex_idx is not None:  # no hit
                extras_idx.append(ex_idx)
                extras_rows.append(ex_row)
        if upd_rows:
            upd_rows = np.concatenate(upd_rows)
            upd_cols = np.concatenate(upd_cols)
            upd_vals = np.concatenate(upd_vals)
            delta = ssp.coo_matrix(
                (upd_vals, (upd_rows, upd_cols)), shape=(cb4_nrow, cb4_ncol)
            ).tocsc()
            merged_X = self.cb4_X + delta
        else:
            merged_X = self.cb4_X.copy()

        merged_idx = list(self.cb4_idx)
        if extras_rows:
            extras_mat = ssp.vstack(extras_rows, format="csc")
            merged_X = ssp.vstack([merged_X, extras_mat], format="csc")
            merged_idx.extend(extras_idx)

        return merged_X, merged_idx

def convert_solomtx_to_h5ad(mtx_dir, gtf_pickle=None):
    mtx_dir = Path(mtx_dir)

    matrix_fn = mtx_dir / 'matrix.mtx'
    feature_fn = mtx_dir / 'features.tsv'
    barcodes_fn = mtx_dir / 'barcodes.tsv'
    h5ad_fn = mtx_dir / 'matrix.h5ad'

    feature_df = pd.read_csv(feature_fn, sep='\t', header=None)
    var = feature_df
    var.index = var[0]
    var = var.rename_axis(index=None, columns=None)
    var.columns = ['gene_ids', 'genome', 'feature_types']
    var = var[['gene_ids', 'feature_types', 'genome']]
    var['genome'] = var['genome'].apply(lambda x: x.split('_')[0])

    if gtf_pickle is not None:
        gtf_pickle = Path(gtf_pickle)

    if gtf_pickle is not None and gtf_pickle.exists():
        with open(gtf_pickle, 'rb') as f:
            gtf_data = pickle.load(f)

        select_pair = gtf_data.drop_duplicates(subset='gene_id', keep='first')[['gene_id', 'gene_name', 'gene_biotype']]
        merge_biotype = pd.merge(var, select_pair, left_on='gene_ids', right_on='gene_id', how='left').drop('gene_id', axis=1)
        merge_biotype.index = var.index
        var = merge_biotype

    barcodes_df = pd.read_csv(barcodes_fn, sep='\t', header=None, index_col=0)
    obs = barcodes_df
    obs.index.name = 'barcode'

    adata = sc.read_mtx(matrix_fn)
    # Without copy - downstream does not work correctly
    # adata.X (csr_matrix) -> adata.X.T (csc_matrix)
    adata = adata.T.copy()

    adata.obs = obs
    adata.var = var

    adata.write_h5ad(h5ad_fn, compression='gzip')
    return h5ad_fn.absolute()

def write_10X_h5(adata, file, chemistry_description="Dropseq CBLEN=12 UMILEN=8"):
    """Writes adata to a 10X-formatted h5 file.

    https://www.10xgenomics.com/support/software/cell-ranger/analysis/outputs/cr-outputs-h5-matrices

    Note that this function is not fully tested and may not work for all cases.
    It will not write the following keys to the h5 file compared to 10X:
    '_all_tag_keys', 'pattern', 'read', 'sequence'

    Args:
        adata (AnnData object): AnnData object to be written.
        file (str): File name to be written to. If no extension is given, '.h5' is appended.

    Raises:
        FileExistsError: If file already exists.

    Returns:
        None
    """

    if '.h5' not in file:
        file = f'{file}.h5'

    if Path(file).exists():
        print(f"There already is a file `{file}`, detele and regenerate it.")
        os.remove(file)

    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)

    def str_max(x):
        return max([len(i) for i in x])

    w = h5py.File(file, 'w')

    w.attrs['chemistry_description'] = chemistry_description
    w.attrs['filetype'] = 'matrix'
    w.attrs['library_ids'] = np.array([b'STAR-avx2'])
    w.attrs['original_gem_groups'] = np.array([1])
    w.attrs['software_version'] = 'STAR-2.7.11a'
    w.attrs['version'] = np.int64(2)

    adata_csr = adata.X.tocsr()
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata_csr.data, dtype=f'<i{int_max(adata_csr.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("_all_tag_keys", data=np.array([b'genome']))
    ftrs.create_dataset("feature_type", data=np.array(adata.var.feature_types, dtype=f'|S{str_max(adata.var.feature_types)}'))
    ftrs.create_dataset("genome", data=np.array(adata.var.genome, dtype=f'|S{str_max(adata.var.genome)}'))
    ftrs.create_dataset("id", data=np.array(adata.var.gene_ids, dtype=f'|S{str_max(adata.var.gene_ids)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata_csr.indices, dtype=f'<i{int_max(adata_csr.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata_csr.indptr, dtype=f'<i{int_max(adata_csr.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata_csr.shape)[::-1], dtype=f'<i{int_max(adata_csr.shape)}'))

def write_gtf(df: pl.DataFrame, export_path: Path, headers: t.List[str] = None, COMMONS_COL=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame']):
    def commons_cols(row) -> str:
        return "\t".join([str(row[field] or '.') for field in COMMONS_COL])

    def custom_fields(row) -> str:
        return "; ".join([f'{field} "{row[field]}"' for field in row.keys() if (field not in COMMONS_COL) and (row[field])])

    headers = headers or []
    with open(export_path, 'w') as f:
        for header in headers:
            f.write(f"{header}\n")
        for row in df.iter_rows(named=True):
            f.write(f"{commons_cols(row)}\t{custom_fields(row)}\n")

def create_ensembl_row(
        seqname="",
        source="",
        feature="",
        start=1,
        end=1,
        score=None,
        strand="",
        frame=0,
        gene_id="",
        gene_version="",
        gene_name="",
        gene_source="",
        gene_biotype="",
        transcript_id="",
        transcript_version="",
        transcript_name="",
        transcript_source="",
        transcript_biotype="",
        tag="",
        transcript_support_level="",
        exon_number="",
        exon_id="",
        exon_version="",
        protein_id="",
        protein_version="",
        ccds_id=""
):
    # seqname, source, feature, start, end, score, strand, frame, gene_id, gene_version, gene_name, gene_source, gene_biotype, transcript_id, transcript_version, transcript_name, transcript_source, transcript_biotype, tag, transcript_support_level, exon_number, exon_id, exon_version, protein_id, protein_version, ccds_id
    return {
        'seqname': seqname,
        'source': source,
        'feature': feature,
        'start': start,
        'end': end,
        'score': score,
        'strand': strand,
        'frame': frame,
        'gene_id': gene_id,
        'gene_version': gene_version,
        'gene_name': gene_name,
        'gene_source': gene_source,
        'gene_biotype': gene_biotype,
        'transcript_id': transcript_id,
        'transcript_version': transcript_version,
        'transcript_name': transcript_name,
        'transcript_source': transcript_source,
        'transcript_biotype': transcript_biotype,
        'tag': tag,
        'transcript_support_level': transcript_support_level,
        'exon_number': exon_number,
        'exon_id': exon_id,
        'exon_version': exon_version,
        'protein_id': protein_id,
        'protein_version': protein_version,
        'ccds_id': ccds_id
    }

def ensembl_dicts_to_pldf(ensembl_pddf):
    return pl.from_dicts(ensembl_pddf).with_columns([
        pl.col('seqname').cast(pl.Categorical).alias('seqname'),
        pl.col('source').cast(pl.Categorical).alias('source'),
        pl.col('feature').cast(pl.Categorical).alias('feature'),
        pl.col('start').cast(pl.Int64).alias('start'),
        pl.col('end').cast(pl.Int64).alias('end'),
        pl.when(pl.col('score') == 'nan').then(None).otherwise(pl.col('score')).cast(pl.Float32).alias('score'),
        pl.col('strand').cast(pl.Categorical).alias('strand'),
        pl.col('frame').cast(pl.Int64).alias('frame')
    ])

def sparse_equal(A: ssp.csc_matrix, B: ssp.csc_matrix) -> bool:
    if A.shape != B.shape:
        return False

    return (A != B).nnz == 0
