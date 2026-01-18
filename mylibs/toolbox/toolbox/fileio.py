import json
import csv
import hashlib
import h5py
import platform
import socket
import getpass
from datetime import datetime
from pathlib import Path
from collections import OrderedDict
from . import util


def ensure_dir(dirname):
    dirname = Path(dirname)
    if not dirname.is_dir():
        dirname.mkdir(parents=True, exist_ok=False)

def get_file_md5(file_path):
    with open(file_path, 'rb') as f:
        md5_obj = hashlib.md5()
        while True:
            data = f.read(8192)
            if not data:
                break
            md5_obj.update(data)
    return md5_obj.hexdigest()

def generate_cache_file(filename, caches_dir='caches', child_dir=None, filetype='', suffix='', prefix=''):
    if child_dir is not None:
        cache = Path(caches_dir) / Path(child_dir)
    else:
        cache = Path(caches_dir)
    cache.mkdir(parents=True, exist_ok=True)
    return cache / Path(
        f"{prefix}{'_' if prefix != '' else ''}{filename}{'_' if suffix != '' else ''}{suffix}{'.' if filetype != '' else ''}{filetype}")

def generate_random_cache_file(caches_dir='caches', child_dir=None, filetype='', suffix='', prefix=''):
    return generate_cache_file(util.generate_uuid(), caches_dir, child_dir, filetype, suffix, prefix)

def generate_timedate_cache_file(caches_dir='caches', child_dir=None, filetype='', suffix='', prefix=''):
    return generate_cache_file(util.generate_datetime(), caches_dir, child_dir, filetype, suffix, prefix)

def read_json(fname):
    fname = Path(fname)
    with fname.open('rt') as handle:
        return json.load(handle, object_hook=OrderedDict)

def write_json(content, fname):
    fname = Path(fname)
    with fname.open('wt', encoding='utf-8') as handle:
        json.dump(content, handle, indent=4, sort_keys=False, ensure_ascii=False)

def get_stem0(file_path: Path):
    return file_path.stem.split('.')[0]

def read_csv_to_string(file_path):
    with open(file_path, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        result = '\n'.join([','.join(row) for row in reader])
    return result

def read_h5py(fname, mode='r', timedate_prefix=True, owner=None, include_sysinfo=True):
    fname = Path(fname)
    if not fname.exists() and timedate_prefix:
        fname = fname.parent / Path(f"{util.generate_datetime()}_{fname.stem}.h5")
    else:
        fname = fname.parent / Path(f"{fname.stem}.h5")
    h5f = h5py.File(fname, mode)
    curr_dt = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if owner is None:
        owner = getpass.getuser()
    if include_sysinfo and 'r' not in mode:
        if 'system' not in h5f.attrs:
            h5f.attrs['system'] = platform.uname().system
            h5f.attrs['release'] = platform.uname().release
            h5f.attrs['os'] = platform.uname().node
            h5f.attrs['os_version'] = platform.uname().version
            h5f.attrs['arch'] = platform.uname().machine
            h5f.attrs['machine'] = socket.gethostname()
            h5f.attrs['created_time'] = curr_dt
            h5f.attrs['created_by'] = owner
        h5f.attrs['last_modified_time'] = curr_dt
    h5f.flush()
    return h5f
