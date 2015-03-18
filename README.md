# mri_mask_metrics
A tool to assist in creating a whole-brain group mask for very large MRI data sets. This is a work in progress.

Given a list of masks and a target mask to compare to, this script will calcualte the following metrics:
- Overlap between each mask and the target (jaccard index and dice coefficient).
- Number of voxels that are different between mask and target.
- Number of voxels in target that are missing from mask.
- For each parcel in the Harvard-Oxford cortical and subcortical atlases:
  - Percentage of parcel voxels that are also in mask.
 
The path to the HO atlases will be pulled automatically using `$FSLDIR`. I've been using the MNI brain mask as the target so far, and it seems to work well. The code is set up to run on images with 2mm voxels, but this can be modified by pointing to a different set of atlas images.

## Dependencies
- Pandas
- NumPy
- NiBabel
- NiPype

## Usage
```
python calc_overlaps.py mask_list target_mask output_path
```
Parameters:
- `mask_list` is a list of paths to binary mask NIfTI iamges, one on each line.
- `target_mask` is the path to a binary mask NIfTI image.
- `output_path` the path/name of the output `.csv` file.
