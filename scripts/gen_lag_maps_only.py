from __future__ import division, print_function
import numpy as np
import nibabel as nib
from scipy import signal as sci_signal
from nilearn import image, masking
import sys
import os


def get_numerator(signal_a, signal_b, lag):
    """
    Calculates the numerator of the cross-correlation equation.
    
    Parameters
    ----------
    signal_a : array_like (1D)
        Reference signal.
    signal_b : array_like (1D)
        Test signal. Must be the same length as signal_a.
    lag : int
        Lag by which signal_b will be shifted relative to signal_a.
        
    Returns
    -------
    array_like (1D)
        Element-wise product of matching time points in the lagged signals.
    """
    if lag == 0:
        numerator = np.multiply(signal_a, signal_b)
    # If lag is positive, shift signal_b forwards relative to signal_a.
    if lag > 0:
        numerator = np.multiply(signal_a[lag:], signal_b[0:-lag])
    # If lag is negative, shift signal_b backward relative to signal_a.
    if lag < 0:
        numerator = np.multiply(signal_b[-lag:], signal_a[0:lag])
    return numerator

def get_denominator(signal_a, signal_b):
    """
    Calculates the denominator of the cross-correlation equation.
    
    Parameters
    ----------
    signal_a : array_like (1D)
        Reference signal.
    signal_b : array_like (1D)
        Test signal. Must be the same length as signal_a.
        
    Returns
    -------
    float
        Product of the standard deviations of the input signals.
    """ 
    return np.std(signal_a) * np.std(signal_b)

def calc_xcorr(signal_a, signal_b, lag):
    """
    Calculate the cross-correlation of two signals at a given lag.
    
    Parameters
    ----------
    signal_a : array_like (1D)
        Reference signal.
    signal_b : array_like (1D)
        Test signal. Must be the same length as signal_a.
    lag : int
        Lag by which signal_b will be shifted relative to signal_a.
        
    Returns
    -------
    float
        Normalized cross-correlation.
    """ 
    xcorr = np.true_divide(1., len(signal_a)-np.absolute(lag)) * np.sum(np.true_divide(get_numerator(signal_a, signal_b, lag),
              get_denominator(signal_a, signal_b)))
    return xcorr

def sliding_xcorr(signal_a, signal_b, lags):
    """
    Calculate the cross-correlation of two signals over a range of lags.
    
    Parameters
    ----------
    signal_a : array_like (1D)
        Reference signal.
    signal_b : array_like (1D)
        Test signal. Must be the same length as signal_a.
    lags : array_like (1D)
        Lags by which signal_b will be shifted relative to signal_a.
        
    Returns
    -------
    array_like (1D)
        Normalized cross-correlation at each lag.
    """ 
    xcorr_vals = []
    for lag in lags:
        xcorr = calc_xcorr(signal_a, signal_b, lag)
        xcorr_vals.append(xcorr)
    return np.array(xcorr_vals)

# Adapted from https://gist.github.com/endolith/255291
def parabolic(sample_array, peak_index):
    """
    Quadratic interpolation for estimating the true position of an
    inter-sample local maximum when nearby samples are known.
   
    Parameters
    ----------
    sample_array : array_like (1D)
        Array of samples.
    peak_index : int
        Index for the local maximum in sample_array for which to estimate the inter-sample maximum.
   
    Returns
    -------
    tuple
        The (x,y) coordinates of the vertex of a parabola through peak_index and its two neighbors.
    """
    vertex_x = 1/2. * (sample_array[peak_index-1] - sample_array[peak_index+1]) / (sample_array[peak_index-1] - 2 * sample_array[peak_index] + sample_array[peak_index+1]) + peak_index
    vertex_y = sample_array[peak_index] - 1/4. * (sample_array[peak_index-1] - sample_array[peak_index+1]) * (vertex_x - peak_index)
    return (vertex_x, vertex_y)

def gen_lag_map(epi_img, brain_mask_img, gm_mask_img, lags):
    epi_gm_masked_data = masking.apply_mask(epi_img, gm_mask_img)
    gm_mean_signal = np.mean(epi_gm_masked_data, axis=1)
    epi_brain_masked_data = masking.apply_mask(epi_img, brain_mask_img).T
    lag_index_correction = np.sum(np.array(lags) > 0)
    xcorr_array = []
    for voxel in epi_brain_masked_data:
        vox_signal = voxel
        vox_xcorr = sliding_xcorr(gm_mean_signal, vox_signal, lags)
        xcorr_maxima = sci_signal.argrelmax(np.array(vox_xcorr), order=1)[0]
        if len(xcorr_maxima) == 0:
            interp_max = np.argmax(vox_xcorr) - lag_index_correction
        elif len(xcorr_maxima) == 1:
            interp_max = parabolic(vox_xcorr, xcorr_maxima[0])[0]
            interp_max = interp_max - lag_index_correction
        elif len(xcorr_maxima) > 1:
            xpeak = xcorr_maxima[np.argmax(vox_xcorr[xcorr_maxima])]
            interp_max = parabolic(vox_xcorr, xpeak)[0]
            interp_max = interp_max - lag_index_correction
        xcorr_array.append(interp_max)
    return(np.array(xcorr_array))


epi_path = sys.argv[1]
brain_mask_path = sys.argv[2]
gm_mask_path = sys.argv[3]
max_lag = int(sys.argv[4])
out_dir = sys.argv[5]
out_prefix = sys.argv[6]
smooth_kernel = sys.argv[7]

TR = 2
lags = range(-max_lag, max_lag+1)

# Load the epi.
epi_img = nib.load(epi_path)
# Load the gray matter mask.
gm_mask_img = nib.load(gm_mask_path)
# Load the brain mask.
brain_mask_img = nib.load(brain_mask_path)
# Generate lag map.
lag_map_data = gen_lag_map(epi_img, brain_mask_img, gm_mask_img, lags)
lag_map_image = masking.unmask(lag_map_data, brain_mask_img)
# Create a smoothed version of the lag map.
lag_map_image_smoothed = image.smooth_img(lag_map_image, float(smooth_kernel))
# Save lag map images.
nib.save(lag_map_image, os.path.join(out_dir, out_prefix+'_lag{}.nii.gz'.format(str(max_lag))))
nib.save(lag_map_image_smoothed, os.path.join(out_dir, out_prefix+'_lag{}_smoothed{}.nii.gz'.format(str(max_lag),smooth_kernel)))



