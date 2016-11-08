import os
import subprocess
import sys

i = sys.argv[1]

IN_DIR = '/home/despo/dlurie/Projects/megarest_lag/data'
OUT_DIR = '/home/despo/dlurie/Projects/megarest_lag/analysis'
MAX_LAG = '10'
FWHM = '5'
SCRIPT_PATH = '/home/despo/dlurie/Projects/megarest_lag/gen_lag_maps_only.py'

i = str(i)
print('Processing subject {}...'.format(i))
# Get the subject input folder.
subject_in_dir = os.path.join(IN_DIR, "sub"+i)
# Create a folder for this subject in the output directory.
print('...creating subject output directory.')
subject_out_dir = os.path.join(OUT_DIR, "sub"+i)
if not os.path.isdir(subject_out_dir):
    os.mkdir(subject_out_dir)

for scan_condition in ['Aifo', 'Dlpfc', 'Sham']:
    # Get the condition input folder.
    print('Condition: {}...'.format(scan_condition))
    scan_condition_in_dir = os.path.join(subject_in_dir, scan_condition)
    # create the condition output folder
    print('...creationg condition output directory.')
    scan_condition_out_dir = os.path.join(subject_out_dir, scan_condition)
    if not os.path.isdir(scan_condition_out_dir):
        os.mkdir(scan_condition_out_dir)
    for scan_run in ['rest_1', 'rest_2', 'rest_3']:
        print('Scan: {}...'.format(scan_run))
        # Get the GM mask path.
        print('...getting GM mask')
        gm_mask_path = os.path.join(scan_condition_out_dir, 'sub{0}_{1}_{2}_minimal_gm_30p_mask.nii.gz'.format(i, scan_condition, scan_run))
        # Get the brain mask path.
        print('...getting brain mask')
        brain_mask_path = os.path.join(scan_condition_out_dir, 'sub{0}_{1}_{2}_minimal_brain_mask.nii.gz'.format(i, scan_condition, scan_run))
        # Generate lag map for CompCor image..
        print('...generating lag map for CompCor image.')
        compcor_4D_path = os.path.join(scan_condition_in_dir, 'sub{0}_{1}_{2}_compcor.nii.gz'.format(i, scan_condition, scan_run))
        subprocess.check_call(['python', SCRIPT_PATH, compcor_4D_path, brain_mask_path, gm_mask_path, MAX_LAG, scan_condition_out_dir,
            'sub{0}_{1}_{2}_compcor'.format(i, scan_condition, scan_run), FWHM])
   
 
