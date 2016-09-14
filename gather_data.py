import os
import shutil
import subprocess
import glob

out_dir = '/home/despo/dlurie/Projects/megarest_lag/data'

for i in range(101,130):
    i = str(i)
    try:
        print('Processing subject {}...'.format(i))
        # Create a folder for this subject in the output directory.
        print('...creating subject output directory')
        subject_out_dir = os.path.join(out_dir, "sub"+i)
        if not os.path.isdir(subject_out_dir):
            os.mkdir(subject_out_dir)
        # Copy GM segmentation.
        print('...copying GM segmentation.')
        gm_seg_path_old = '/home/despo/arielle/megarest_sc/data/anat/{}_gm.nii'.format(i)
        gm_seg_path_new = os.path.join(subject_out_dir, 'sub{}_gm_seg.nii'.format(i))
        if not os.path.exists(gm_seg_path_new):
            shutil.copyfile(gm_seg_path_old, gm_seg_path_new)
       
        for scan_condition in ['Aifo', 'Dlpfc', 'Sham']:
            print('Condition: {}'.format(scan_condition))
            scan_condition_dir = os.path.join(subject_out_dir, scan_condition)
            print('Creationg condition directory')
            if not os.path.isdir(scan_condition_dir):
                os.mkdir(scan_condition_dir)
            print('Copying ROI')
            if scan_condition is 'Sham':
                roi_path_old = '/home/despo/arielle/megarest_sc/data/ROIs/TMS/{}-S1_L_mask.nii.gz'.format(i)
                roi_path_new = os.path.join(scan_condition_dir, 'sub{0}_S1_L_roi.nii.gz'.format(i))
            elif scan_condition is 'Aifo' or 'Dlpfc':
                roi_path_old = '/home/despo/arielle/megarest_sc/data/ROIs/TMS/{0}-{1}.nii.gz'.format(i, scan_condition)
                roi_path_new = os.path.join(scan_condition_dir, 'sub{0}_{1}_roi.nii.gz'.format(i, scan_condition))
            print(roi_path_old)
            print(roi_path_new)
            if not os.path.exists(roi_path_new):
                shutil.copyfile(roi_path_old, roi_path_new)

            for scan_run in ['rest_1', 'rest_2', 'rest_3']:
                print('Scan: {}'.format(scan_run))
                # Copy motion parameters.
                print('...copying motion parameters.')
                motion_params_path_old = glob.glob('/home/despo/arielle/megarest_sc/data/BOLD/{0}/{1}/{2}/rp*.txt'.format(i, scan_condition, scan_run))[0]
                motion_params_path_new = os.path.join(scan_condition_dir, 'sub{0}_{1}_{2}_motion_params.txt'.format(i, scan_condition, scan_run))
                if not os.path.exists(motion_params_path_new):
                    shutil.copyfile(motion_params_path_old, motion_params_path_new)
                # Stack minimally preprocessed images.
                print('...stacking minimally preprocessed')
                minimal_3D_paths = glob.glob('/home/despo/arielle/megarest_sc/data/BOLD/{0}/{1}/{2}/sraf*.nii'.format(i, scan_condition, scan_run))
                minimal_4D_path = os.path.join(scan_condition_dir, 'sub{0}_{1}_{2}_minimal.nii.gz'.format(i, scan_condition, scan_run))
                if not os.path.exists(minimal_4D_path):
                    subprocess.check_call(['fslmerge', '-t', minimal_4D_path] + minimal_3D_paths)
                # Stack CompCor images.
                print('...stacking CompCor')
                compcor_3D_paths = glob.glob('/home/despo/arielle/megarest_sc/data/BOLD/{0}/{1}/{2}/vSCInf_vFILT1_vSIMULT0_vPCAsimN*_sraf*.nii.gz'.format(i, scan_condition, scan_run)) 
                compcor_4D_path = os.path.join(scan_condition_dir, 'sub{0}_{1}_{2}_compcor.nii.gz'.format(i, scan_condition, scan_run))
                if not os.path.exists(compcor_4D_path):
                    subprocess.check_call(['fslmerge', '-t', compcor_4D_path] + compcor_3D_paths)
                # Stack CompCor+GSR images.
                print('...stacking CompCor+GSR')
                compcor_gsr_3D_paths = glob.glob('/home/despo/arielle/megarest_sc/data/BOLD/{0}/{1}/{2}/vSCInf_vFILT1_vSIMULT0_vPCAsimN*_vGS1_sraf*.nii'.format(i, scan_condition, scan_run)) 
                compcor_gsr_4D_path = os.path.join(scan_condition_dir, 'sub{0}_{1}_{2}_compcor_gsr.nii.gz'.format(i, scan_condition, scan_run))
                if not os.path.exists(compcor_gsr_4D_path):
                    subprocess.check_call(['fslmerge', '-t', compcor_gsr_4D_path] + compcor_gsr_3D_paths)
    except:
        print("ERROR ENCOUNTERED FOR SUBJECT {}".format(i))
        continue





    
