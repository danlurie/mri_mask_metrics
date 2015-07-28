import sys
import os
import argparse
import numpy as np
import pandas as pd
import nibabel as nb
import xml.etree.ElementTree as ET
from nipype.algorithms import metrics
from nipype.interfaces import afni as afni

def divbz(array1, array2):
    """
    A hacky way to avoid divide-by-zero errors.
    
    Adds epsilon to any zero values in the input arrays, does a 
    floating point divide on the arrays, then rounds back to zero.

    Parameters
    ----------
    array1, array2 : array_like
        Input arrays.

    Returns
    -------
    array_like
        Array of element-wise quotients, as floats.
    """

    array1[array1 == 0] = np.spacing(1)
    array2[array2 == 0] = np.spacing(1)
    return np.round(np.true_divide(array1, array2), decimals=15)

def roi_overlap(atlas_data, mask_data):
    """
    Calcualte percent overlap between a mask and atlas.
    
    For each ROI in atlas_data, calculate the percentage of voxels
    that are also present in mask_data. See Notes before use.
    
    Parameters
    ----------
    atlas_data : array_like
        Atlas data, where each ROI is a different integer value.
    mask_data : array_like
        A boolean mask of the same dimensions as atlas_data.

    Returns
    -------
    pct_coverage : array_like
        For each ROI from 0 to 59, the percentage of voxels in that ROI
        that are also present in mask_data.

    Notes
    -----
    Voxel counts for each ROI are currently calculated by doing a
    bin-count on atlas_data, with number of bins hard-coded to 60.
    This is inelegant, and requires masking the output array later
    to match values with ROI names.

    Future versions will hopefully avoid this issue.
    """

    if atlas_data.shape == mask_data.shape:
        coverage_data = atlas_data * mask_data
        coverage_data = coverage_data.astype('uint8')
        nvox_atlas = np.bincount(atlas_data.ravel(), minlength=60)
        nvox_coverage = np.bincount(coverage_data.ravel(), minlength=60)
        pct_coverage = divbz(nvox_coverage, nvox_atlas)
        return pct_coverage 
    else:
        print('Atlas and mask must be the same dimensions!')

def get_atlas_labels(path_to_xml):
    tree = ET.parse(path_to_xml)
    root = tree.getroot()
    label_list = []
    for row in root[1]:
        label_list.append(row.text)
    return np.array(label_list)

# Read in required command line arguments.
mask_list_file, target_mask_path, output_file_path = sys.argv[1:]

# Load the list of masks.
mask_list = np.loadtxt(mask_list_file, dtype='S')

# Set up function to calculate overlap between individual mask and target_mask.
overlap = metrics.Overlap()
overlap.inputs.volume1 = target_mask_path
overlap.inputs.out_file = 'overlap.nii'

# Set up a function that does the same as above, but only considers voxels in target_mask.
overlap_m = metrics.Overlap()
overlap_m.inputs.volume1 = target_mask_path
overlap_m.inputs.mask_volume = target_mask_path
overlap_m.inputs.out_file = 'overlap_m.nii'

# Load Harvard-Oxford atlases (default is 2mm, 25%).
fsl_dir = os.getenv('FSLDIR')
ho_cort = nb.load('/'.join([fsl_dir,
                            'data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz'])).get_data()
ho_sub = nb.load('/'.join([fsl_dir,
                           'data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz'])).get_data()

# Create an array to store calculated metric values.
tup_store =[]

for mask_path in mask_list:
    # Calculate whole-image overlap.
    overlap.inputs.volume2 = mask_path
    ol_res = overlap.run()
    ol_diff_data = nb.load('overlap.nii').get_data()
    ol_diff_nvox = np.count_nonzero(ol_diff_data)

    # Calculate overlap within target_mask.
    overlap_m.inputs.volume2 = mask_path
    ol_m_res = overlap_m.run()
    ol_m_diff_data = nb.load('overlap_m.nii').get_data()
    ol_m_diff_nvox = np.count_nonzero(ol_m_diff_data)

    # Store target_mask overlap metrics.
    val_array = [mask_path, ol_res.outputs.dice, ol_m_res.outputs.dice, ol_res.outputs.jaccard,
                 ol_m_res.outputs.jaccard, ol_diff_nvox, ol_m_diff_nvox]

    # Calculate atlas overlap percentages, and mask output arrays to exclude non-ROI bin-counts.
    mask_data = nb.load(mask_path).get_data()
    cort_pct_array = roi_overlap(ho_cort, mask_data)[1:49]
    sub_pct_array = roi_overlap(ho_sub, mask_data)[np.array([10, 11, 12, 13, 16, 17, 18, 26,
                                                             49, 50, 51, 52, 53, 54, 58])] 
    
    # Store atlas overlap percentages.
    val_array.extend(cort_pct_array)
    val_array.extend(sub_pct_array)

    # Archive overlap metrics for this subject.
    tup_store.append(tuple(val_array))

    print "Metrics calcualted for %s" % mask_path

# Define default column names for output .csv
column_names = ['mask', 'dice', 'dice_m', 'jaccard', 'jaccard_m', 'diff_nvox', 'missing_nvox']

# Get ROI names and add them as column names.
column_names.extend(get_atlas_labels('/'.join([fsl_dir,'data/atlases/HarvardOxford-Cortical.xml'])))
ho_sub_cols = np.array([3,4,5,6,7,8,9,10,14,15,16,17,18,19,20]) 
column_names.extend(get_atlas_labels('/'.join([fsl_dir,'data/atlases/HarvardOxford-Subcortical.xml']))[ho_sub_cols])

# Create a data frame from metrics archive and write as a .csv
overlap_df = pd.DataFrame.from_records(tup_store, columns=column_names)
overlap_df.to_csv(output_file_path)
