import os, glob, argparse, warnings

import numpy as np
import pandas as pd
from nibabel.freesurfer.io import read_annot

from freesurfer import fsaveragetools

MNI_COLS = ['MNI_x', 'MNI_y', 'MNI_z']
T1_COLS = ['T1_x', 'T1_y', 'T1_z']
FS_COLS = ['fs_x', 'fs_y', 'fs_z']
OUT_COLS = ['labels', 'T1_R'] + MNI_COLS + T1_COLS + \
           ['T1_AnatomicalRegion', 'hemi', 'hemi_method'] + FS_COLS + \
           ['fs_vertex', 'HCP_region', 'dist_to_pial_mm', 'subj_vertex']


def load_subject(subj_dir):
    '''
    loads coordinates.csv and last_valid_electrode.txt from subj_dir
    and keeps the first last_valid_electrode rows.
    rows without coordinates are kept and flagged in 'valid_coords'.
    '''
    subj = os.path.basename(os.path.normpath(subj_dir))
    df = pd.read_csv(os.path.join(subj_dir, 'coordinates.csv'))
    with open(os.path.join(subj_dir, 'last_valid_electrode.txt'), 'r') as f:
        last_valid = int(f.read().strip())
    if last_valid > len(df):
        raise ValueError(f'{subj}: last_valid_electrode={last_valid} '
                         f'but coordinates.csv has {len(df)} rows')
    df = df.iloc[:last_valid].reset_index(drop=True)
    # depth electrodes are not expected in the kept rows
    is_depth = df['labels'].astype(str).str.match(r'[dD]') | (df['T1_R'] == 'D')
    if is_depth.any():
        warnings.warn(f'{subj}: possible depth electrodes in kept rows: '
                      f'{list(df.loc[is_depth, "labels"])}')
    df['valid_coords'] = df[MNI_COLS + T1_COLS].notna().all(axis=1)
    return df


def _majority(hemis):
    ''' returns 'lh' or 'rh' by majority vote, None if empty or tie '''
    n_lh = np.sum(hemis == 'lh')
    n_rh = np.sum(hemis == 'rh')
    if n_lh == n_rh:
        return None
    return 'lh' if n_lh > n_rh else 'rh'


def assign_hemisphere(df, tol=5.0, subj=''):
    '''
    assigns hemisphere from the sign of MNI_x (negative -> lh).
    electrodes with |MNI_x| <= tol are ambiguous and get the majority
    hemisphere of their device (label without trailing digits), then
    the subject majority. rows without coordinates get NaN.
    '''
    x = df['MNI_x'].to_numpy(dtype=float)
    valid = df['valid_coords'].to_numpy()
    hemi = np.full(len(df), np.nan, dtype=object)
    method = np.full(len(df), np.nan, dtype=object)
    # confident electrodes by sign
    confident = valid & (np.abs(x) > tol)
    hemi[confident] = np.where(x[confident] < 0, 'lh', 'rh')
    method[confident] = 'sign'
    # ambiguous electrodes near the midline
    device = df['labels'].astype(str).str.replace(r'\d+$', '', regex=True).to_numpy()
    subj_major = _majority(hemi[confident])
    for i in np.where(valid & ~confident)[0]:
        dev_major = _majority(hemi[confident & (device == device[i])])
        if dev_major is not None:
            hemi[i], method[i] = dev_major, 'device'
        elif subj_major is not None:
            hemi[i], method[i] = subj_major, 'subject'
        else:
            hemi[i], method[i] = ('lh' if x[i] < 0 else 'rh'), 'sign_deadband'
    if np.any(hemi == 'lh') and np.any(hemi == 'rh'):
        warnings.warn(f'{subj}: electrodes assigned to both hemispheres '
                      f'(lh={np.sum(hemi == "lh")}, rh={np.sum(hemi == "rh")})')
    df['hemi'] = hemi
    df['hemi_method'] = method
    return df


def map_hemi_to_fsaverage(elec_locs, subj_dir, fs_dir, hemi):
    '''
    maps T1 electrode locations of one hemisphere to fsaverage.
    closest subject pial vertex -> same vertex on subject sphere.reg
    -> closest fsaverage sphere.reg vertex -> fsaverage pial.
    returns fs coords, fs vertex index, subject vertex index and
    distance (mm) from the electrode to the subject pial.
    '''
    subj = os.path.basename(os.path.normpath(subj_dir))
    fst = fsaveragetools(subj, hemi,
                         fsaverage_pial=os.path.join(fs_dir, f'{hemi}.pial'),
                         fsaverage_sphere=os.path.join(fs_dir, f'{hemi}.sphere.reg'),
                         subj_pial=os.path.join(subj_dir, f'{hemi}.pial'),
                         subj_sphere=os.path.join(subj_dir, f'{hemi}.sphere.reg'))
    fst.load_fs_pial()
    fst.load_fs_sphere()
    fst.load_subj_pial()
    fst.load_subj_shpere()
    if fst.subj_pial_verts.shape[0] != fst.subj_sphr_verts.shape[0]:
        raise ValueError(f'{subj}: {hemi}.pial and {hemi}.sphere.reg '
                         'have different number of vertices')
    elec_locs = np.asarray(elec_locs, dtype=float)
    n = elec_locs.shape[0]
    subj_ind = np.zeros(n, dtype=int)
    fs_ind = np.zeros(n, dtype=int)
    dist = np.zeros(n)
    for i in range(n):
        d2 = np.sum((fst.subj_pial_verts - elec_locs[i, :])**2, axis=1)
        subj_ind[i] = np.argmin(d2)
        dist[i] = np.sqrt(d2[subj_ind[i]])
        fs_ind[i] = np.argmin(np.sum((fst.fs_sphr_verts -
                                      fst.subj_sphr_verts[subj_ind[i], :])**2, axis=1))
    return fst.fs_pial_verts[fs_ind, :], fs_ind, subj_ind, dist


def lookup_hcp(fs_ind, hemi, fs_dir, labels=None, subj=''):
    '''
    returns the HCP-MMP1 region name for each fsaverage vertex.
    unlabeled vertices ('???' or -1) are returned as 'unknown'
    with a warning.
    '''
    annot_ids, _, names = read_annot(os.path.join(fs_dir, f'{hemi}.HCP-MMP1.annot'))
    names = np.array([n.decode('utf-8') for n in names], dtype=object)
    ids = annot_ids[fs_ind]
    regions = np.where(ids >= 0, names[np.clip(ids, 0, None)], 'unknown')
    regions = np.where(regions == '???', 'unknown', regions).astype(object)
    unknown = regions == 'unknown'
    if unknown.any():
        who = list(np.asarray(labels)[unknown]) if labels is not None else int(unknown.sum())
        warnings.warn(f'{subj} {hemi}: electrodes mapped to unlabeled '
                      f'HCP-MMP1 vertices, set to unknown: {who}')
    return regions


def process_subject(subj_dir, fs_dir, tol=5.0, out_name='coordinates_hcp.csv'):
    '''
    runs the full mapping for one subject and saves the output
    csv next to the input.
    '''
    subj = os.path.basename(os.path.normpath(subj_dir))
    df = load_subject(subj_dir)
    df = assign_hemisphere(df, tol=tol, subj=subj)
    for c in FS_COLS + ['dist_to_pial_mm']:
        df[c] = np.nan
    df['fs_vertex'] = pd.Series(pd.NA, index=df.index, dtype='Int64')
    df['subj_vertex'] = pd.Series(pd.NA, index=df.index, dtype='Int64')
    df['HCP_region'] = pd.Series(np.nan, index=df.index, dtype=object)
    for hemi in ['lh', 'rh']:
        m = (df['hemi'] == hemi).to_numpy()
        if not m.any():
            continue
        fs_locs, fs_ind, subj_ind, dist = map_hemi_to_fsaverage(
            df.loc[m, T1_COLS].to_numpy(dtype=float), subj_dir, fs_dir, hemi)
        df.loc[m, FS_COLS] = fs_locs
        df.loc[m, 'fs_vertex'] = fs_ind
        df.loc[m, 'subj_vertex'] = subj_ind
        df.loc[m, 'dist_to_pial_mm'] = dist
        df.loc[m, 'HCP_region'] = lookup_hcp(fs_ind, hemi, fs_dir,
                                             labels=df.loc[m, 'labels'], subj=subj)
    df = df[OUT_COLS]
    fn_out = os.path.join(subj_dir, out_name)
    df.to_csv(fn_out, index=False)
    print(f'{subj}: n={len(df)} lh={np.sum(df["hemi"] == "lh")} '
          f'rh={np.sum(df["hemi"] == "rh")} no_coords={df["hemi"].isna().sum()} '
          f'unknown={np.sum(df["HCP_region"] == "unknown")} -> {fn_out}')
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='map electrodes to HCP-MMP1 regions')
    parser.add_argument('--data_dir', default='./Data')
    parser.add_argument('--fs_dir', default='./sample_data/fsaverage')
    parser.add_argument('--subjects', nargs='*', default=None)
    parser.add_argument('--tol', type=float, default=5.0,
                        help='|MNI_x| (mm) below which the hemisphere is ambiguous')
    args = parser.parse_args()
    warnings.simplefilter('always')
    if args.subjects:
        subj_dirs = [os.path.join(args.data_dir, s) for s in args.subjects]
    else:
        subj_dirs = sorted(os.path.dirname(p) for p in
                           glob.glob(os.path.join(args.data_dir, '*', 'coordinates.csv')))
    for subj_dir in subj_dirs:
        process_subject(subj_dir, args.fs_dir, tol=args.tol)
