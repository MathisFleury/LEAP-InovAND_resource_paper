#!/usr/bin/env python3
"""
Build Custom Subcortical Atlas for Yabplot
==========================================

One-time script: extracts 3D VTK meshes from the 4S156 NIfTI atlas
for the 4 subcortical atlas components:
  - CIT168Subcortical  (striatum, basal ganglia, hypothalamus...)
  - SubcorticalHCP     (hippocampus, amygdala)
  - Cerebellum         (10 cerebellar regions)
  - ThalamusHCP        (thalamic nuclei + LH/RH_Thalamus)

Saves .vtk mesh files to:
  _resources/custom_subcortical_atlas_4S156/

After running this script once, the visualization scripts in
6_functional_analysis and 7_cluster_functional_analysis can use:
  plot_subcortical(data=..., custom_atlas_path=CUSTOM_ATLAS_DIR, ...)
"""

from pathlib import Path
import pandas as pd
from yabplot import build_subcortical_atlas

# =============================================================================
# PATHS
# =============================================================================
NII_PATH = (
    '/Volumes/Imaging5/EEG_MRI-MF/LEAP/_converted_BIDS/derivatives/xcp_d/'
    'atlases/atlas-4S156Parcels/'
    'atlas-4S156Parcels_space-MNI152NLin2009cAsym_res-2_dseg.nii.gz'
)
TSV_PATH = (
    '/Volumes/Imaging5/EEG_MRI-MF/LEAP/_converted_BIDS/derivatives/xcp_d/'
    'atlases/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv'
)

OUT_DIR = str(Path(__file__).parent / 'custom_subcortical_atlas_4S156')

SUBCORTICAL_ATLASES = ['CIT168Subcortical', 'SubcorticalHCP', 'Cerebellum', 'ThalamusHCP']


def main():
    print("=" * 60)
    print("Building Custom Subcortical Atlas — 4S156")
    print("=" * 60)

    # Load TSV and filter to the 4 atlas names
    df = pd.read_csv(TSV_PATH, sep='\t')
    df_sub = df[df['atlas_name'].isin(SUBCORTICAL_ATLASES)][['index', 'label', 'atlas_name']]

    # Build labels_dict: {integer_id: label_name}
    labels_dict = dict(zip(df_sub['index'].astype(int), df_sub['label']))

    print(f"\nRegions to extract ({len(labels_dict)}):")
    for idx, name in sorted(labels_dict.items()):
        print(f"  {idx:4d}  {name}")

    print(f"\nOutput directory: {OUT_DIR}")

    # The atlas NIfTI is 4D (x, y, z, 1) — squeeze to 3D for marching cubes
    import nibabel as nib
    import numpy as np
    import tempfile, os
    img = nib.load(NII_PATH)
    data3d = np.squeeze(img.get_fdata())
    print(f"NIfTI shape: {img.shape} → squeezed to {data3d.shape}")
    tmp_nii = os.path.join(OUT_DIR, '_atlas_3d.nii.gz')
    Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(data3d, img.affine, img.header), tmp_nii)

    print("\nExtracting meshes (this may take a few minutes)...")

    build_subcortical_atlas(
        nii_path=tmp_nii,
        labels_dict=labels_dict,
        out_dir=OUT_DIR,
    )

    # Remove the temporary 3D NIfTI
    if os.path.exists(tmp_nii):
        os.remove(tmp_nii)

    # Count saved VTKs
    from pathlib import Path as P
    vtks = list(P(OUT_DIR).glob('*.vtk'))
    print(f"\nDone: {len(vtks)} VTK meshes saved to {OUT_DIR}")


if __name__ == '__main__':
    main()
