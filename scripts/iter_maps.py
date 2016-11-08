import os
import subprocess

IN_DIR = '/home/despo/dlurie/Projects/megarest_lag/data'
OUT_DIR = '/home/despo/dlurie/Projects/megarest_lag/analysis'
MAX_LAG = '6'
FWHM = '5'
SCRIPT_PATH = '/home/despo/dlurie/Projects/lag_maps/gen_lag_map.py'

#for i in range(101,115)+[116]+range(118,130):
for i in range(101,105):
    i = str(i)
    print('Processing subject {}...'.format(i))
    # Get the subject input folder.
    subject_in_dir = os.path.join(IN_DIR, "sub"+i)
    # Create a folder for this subject in the output directory.
    print('...creating subject output directory.')
    subject_out_dir = os.path.join(OUT_DIR, "sub"+i)
    if not os.path.isdir(subject_out_dir):
        os.mkdir(subject_out_dir)
    # Get the brain mask path.
    print('...getting brain mask.')
    brain_mask_path = os.path.join(subject_in_dir, 'sub{}_brain_mask.nii'.format(i))
    # Get the GM mask path.
    print('...getting 50% GM mask.')
    gm_mask_path = os.path.join(subject_in_dir, 'sub{}_gm_mask_p50.nii'.format(i))
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
            # Generate lag map for minimally preprocessed.
            print('...generating lag map for minimally preprocessed image.')
            minimal_4D_path = os.path.join(scan_condition_in_dir, 'sub{0}_{1}_{2}_minimal.nii.gz'.format(i, scan_condition, scan_run))
            subprocess.check_call(['python', SCRIPT_PATH, minimal_4D_path, brain_mask_path, gm_mask_path, MAX_LAG, scan_condition_out_dir,
                'sub{0}_{1}_{2}_minimal_lag_map'.format(i, scan_condition, scan_run), FWHM])
            # Generate lag map for CompCor
            print('...generating lag map for CompCor image.')
            compcor_4D_path = os.path.join(scan_condition_in_dir, 'sub{0}_{1}_{2}_compcor.nii.gz'.format(i, scan_condition, scan_run))
            subprocess.check_call(['python', SCRIPT_PATH, compcor_4D_path, brain_mask_path, gm_mask_path, MAX_LAG, scan_condition_out_dir,
                'sub{0}_{1}_{2}_compcor_lag_map'.format(i, scan_condition, scan_run), FWHM])
            # Generate lag map for CompCor+GSR
            print('...generating lag map for CompCor+GSR image.')
            compcor_gsr_4D_path = os.path.join(scan_condition_in_dir, 'sub{0}_{1}_{2}_compcor_gsr.nii.gz'.format(i, scan_condition, scan_run))
            subprocess.check_call(['python', SCRIPT_PATH, compcor_gsr_4D_path, brain_mask_path, gm_mask_path, MAX_LAG, scan_condition_out_dir,
                'sub{0}_{1}_{2}_compcor_gsr_lag_map'.format(i, scan_condition, scan_run), FWHM])

 
