import os, glob, argparse, warnings

import numpy as np
import pandas as pd
from nibabel.freesurfer.io import read_annot

from freesurfer import fsaveragetools

MNI_COLS = ['MNI_x', 'MNI_y', 'MNI_z']
T1_COLS = ['T1_x', 'T1_y', 'T1_z']
FS_COLS = ['fs_x', 'fs_y', 'fs_z']
# output column -> fsaverage annotation name ({hemi}.{name}.annot)
ATLASES = {'HCP_region': 'HCP-MMP1',
           'Schaefer400_17Net': 'Schaefer2018_400Parcels_17Networks_order'}
# annotation names for unlabeled vertices (medial wall etc.)
UNLABELED = {'???', 'Background+FreeSurfer_Defined_Medial_Wall'}
OUT_COLS = ['labels'] + MNI_COLS + T1_COLS + \
           ['T1_R', 'AR_Percentage', 'T1_AnatomicalRegion', 'hemi', 'hemi_method'] + \
           FS_COLS + ['fs_vertex', 'subj_vertex', 'dist_to_pial_mm',
                      'HCP_region', 'Schaefer400_17Net',
                      'Schaefer400_17Net_network', 'reason']
# electrode types by T1_R
GRID_STRIP = {'G', 'EG', 'S'}
DEPTH = 'D'
# min percentage of a depth anatomical label to be used
CTX_THRESH = 30.0
# reasons that are mapped to fsaverage
MAPPED = {'grid_strip', 'depth_cortex'}


def load_subject(subj_dir, use_last_valid=False):
    '''
    loads coordinates.csv from subj_dir. if use_last_valid, also loads
    last_valid_electrode.txt and keeps the first last_valid_electrode rows.
    rows without coordinates are kept and flagged in 'valid_coords'.
    '''
    subj = os.path.basename(os.path.normpath(subj_dir))
    df = pd.read_csv(os.path.join(subj_dir, 'coordinates.csv'))
    if use_last_valid:
        with open(os.path.join(subj_dir, 'last_valid_electrode.txt'), 'r') as f:
            last_valid = int(f.read().strip())
        if last_valid > len(df):
            raise ValueError(f'{subj}: last_valid_electrode={last_valid} '
                             f'but coordinates.csv has {len(df)} rows')
        df = df.iloc[:last_valid].reset_index(drop=True)
    df['valid_coords'] = df[MNI_COLS + T1_COLS].notna().all(axis=1)
    return df


def _ar_pairs(row, pair_cols):
    '''
    returns the (pct, pct_str, label) anatomical pairs of a row in file order.
    pair_cols is [pct_col, label_col, pct_col, label_col, ...].
    pairs without a percentage or with an empty label are skipped.
    '''
    pairs = []
    for pc, lc in zip(pair_cols[0::2], pair_cols[1::2]):
        p, l = row[pc], row[lc]
        if isinstance(p, str) and p.strip().endswith('%') and \
           isinstance(l, str) and l.strip():
            pairs.append((float(p.strip()[:-1]), p, l.strip()))
    return pairs


def classify_electrodes(df, subj=''):
    '''
    sets 'reason' for each row:
    no_coords, grid_strip, depth_cortex, depth_ctx_unknown,
    depth_not_cortex or unknown_type.
    for depth electrodes the label is the first one (by descending
    percentage) with pct >= CTX_THRESH that is not Unknown or White-Matter.
    AR_Percentage and T1_AnatomicalRegion are set to that pair; if no
    label qualifies the original first pair is kept.
    '''
    pair_cols = ['AR_Percentage', 'T1_AnatomicalRegion'] + \
                [c for c in df.columns if c.startswith('AR_') and c != 'AR_Percentage']
    reason = np.full(len(df), 'unknown_type', dtype=object)
    for i, row in df.iterrows():
        if not row['valid_coords']:
            reason[i] = 'no_coords'
        elif row['T1_R'] in GRID_STRIP:
            reason[i] = 'grid_strip'
        elif row['T1_R'] == DEPTH:
            pairs = sorted(_ar_pairs(row, pair_cols), key=lambda x: -x[0])
            pick = next((pr for pr in pairs if pr[0] >= CTX_THRESH and
                         pr[2] != 'Unknown' and not pr[2].endswith('White-Matter')), None)
            if pick is None:
                reason[i] = 'depth_not_cortex'
                continue
            df.at[i, 'AR_Percentage'] = pick[1]
            df.at[i, 'T1_AnatomicalRegion'] = pick[2]
            if pick[2].startswith('ctx-') and pick[2].endswith('-unknown'):
                reason[i] = 'depth_ctx_unknown'
            elif pick[2].startswith('ctx-'):
                reason[i] = 'depth_cortex'
            else:
                reason[i] = 'depth_not_cortex'
    df['reason'] = reason
    return df


def check_depth_hemisphere(df, subj=''):
    '''
    warns when the assigned hemisphere of a depth_cortex electrode
    differs from the hemisphere in its ctx-{hemi}-* label.
    '''
    m = (df['reason'] == 'depth_cortex').to_numpy()
    ctx_hemi = df.loc[m, 'T1_AnatomicalRegion'].str.split('-').str[1]
    bad = ctx_hemi != df.loc[m, 'hemi']
    if bad.any():
        who = [f'{l} ({h} vs {r})' for l, h, r in
               zip(df.loc[m, 'labels'][bad], df.loc[m, 'hemi'][bad],
                   df.loc[m, 'T1_AnatomicalRegion'][bad])]
        warnings.warn(f'{subj}: depth electrodes with hemisphere different '
                      f'from their ctx label: {who}')


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


def lookup_annot(fs_ind, hemi, fs_dir, annot_name, labels=None, subj=''):
    '''
    returns the region name in {hemi}.{annot_name}.annot for each
    fsaverage vertex. unlabeled vertices (-1 or a name in UNLABELED)
    are returned as 'unknown' with a warning.
    '''
    annot_ids, _, names = read_annot(os.path.join(fs_dir, f'{hemi}.{annot_name}.annot'))
    names = np.array([n.decode('utf-8') for n in names], dtype=object)
    ids = annot_ids[fs_ind]
    regions = np.where(ids >= 0, names[np.clip(ids, 0, None)], 'unknown')
    regions = np.where(np.isin(regions, list(UNLABELED)), 'unknown', regions).astype(object)
    unknown = regions == 'unknown'
    if unknown.any():
        who = list(np.asarray(labels)[unknown]) if labels is not None else int(unknown.sum())
        warnings.warn(f'{subj} {hemi}: electrodes mapped to unlabeled '
                      f'{annot_name} vertices, set to unknown: {who}')
    return regions


def schaefer_network(parcel):
    '''
    extracts the network from a Schaefer parcel name,
    e.g. 17Networks_LH_DefaultA_PFCm_1 -> DefaultA
    '''
    if not isinstance(parcel, str) or parcel == 'unknown':
        return parcel
    return parcel.split('_')[2]


def process_subject(subj_dir, fs_dir, tol=5.0, use_last_valid=False,
                    out_name='coordinates_fs_atlas.csv'):
    '''
    runs the full mapping for one subject and saves the output
    csv next to the input.
    '''
    subj = os.path.basename(os.path.normpath(subj_dir))
    df = load_subject(subj_dir, use_last_valid=use_last_valid)
    df = classify_electrodes(df, subj=subj)
    df = assign_hemisphere(df, tol=tol, subj=subj)
    check_depth_hemisphere(df, subj=subj)
    for c in FS_COLS + ['dist_to_pial_mm']:
        df[c] = np.nan
    df['fs_vertex'] = pd.Series(pd.NA, index=df.index, dtype='Int64')
    df['subj_vertex'] = pd.Series(pd.NA, index=df.index, dtype='Int64')
    for col in ATLASES:
        df[col] = pd.Series(np.nan, index=df.index, dtype=object)
    # cortical depth electrodes in ctx-*-unknown are not mapped
    df.loc[df['reason'] == 'depth_ctx_unknown', list(ATLASES)] = 'unknown'
    to_map = df['reason'].isin(MAPPED).to_numpy()
    for hemi in ['lh', 'rh']:
        m = to_map & (df['hemi'] == hemi).to_numpy()
        if not m.any():
            continue
        fs_locs, fs_ind, subj_ind, dist = map_hemi_to_fsaverage(
            df.loc[m, T1_COLS].to_numpy(dtype=float), subj_dir, fs_dir, hemi)
        df.loc[m, FS_COLS] = fs_locs
        df.loc[m, 'fs_vertex'] = fs_ind
        df.loc[m, 'subj_vertex'] = subj_ind
        df.loc[m, 'dist_to_pial_mm'] = dist
        for col, annot_name in ATLASES.items():
            df.loc[m, col] = lookup_annot(fs_ind, hemi, fs_dir, annot_name,
                                          labels=df.loc[m, 'labels'], subj=subj)
    df['Schaefer400_17Net_network'] = df['Schaefer400_17Net'].map(schaefer_network)
    df = df[OUT_COLS]
    fn_out = os.path.join(subj_dir, out_name)
    df.to_csv(fn_out, index=False)
    print(f'{subj}: n={len(df)} {df["reason"].value_counts().to_dict()} '
          f'mapped_lh={np.sum(to_map & (df["hemi"] == "lh"))} '
          f'mapped_rh={np.sum(to_map & (df["hemi"] == "rh"))} '
          f'unknown_hcp={np.sum(df["HCP_region"] == "unknown")} '
          f'unknown_schaefer={np.sum(df["Schaefer400_17Net"] == "unknown")} -> {fn_out}')
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='map electrodes to HCP-MMP1 regions')
    parser.add_argument('--data_dir', default='./Data')
    parser.add_argument('--fs_dir', default='./sample_data/fsaverage')
    parser.add_argument('--subjects', nargs='*', default=None)
    parser.add_argument('--tol', type=float, default=5.0,
                        help='|MNI_x| (mm) below which the hemisphere is ambiguous')
    parser.add_argument('--use_last_valid', action='store_true',
                        help='keep only the first last_valid_electrode.txt rows')
    args = parser.parse_args()
    warnings.simplefilter('always')
    if args.subjects:
        subj_dirs = [os.path.join(args.data_dir, s) for s in args.subjects]
    else:
        subj_dirs = sorted(os.path.dirname(p) for p in
                           glob.glob(os.path.join(args.data_dir, '*', 'coordinates.csv')))
    for subj_dir in subj_dirs:
        process_subject(subj_dir, args.fs_dir, tol=args.tol,
                        use_last_valid=args.use_last_valid)
