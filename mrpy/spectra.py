import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import levene
from .utils import loadnii


def load_spect(fn_base, num_strings=None, ext='.nii.gz'):
    if num_strings == None:
        num_strings = [f'{i:03d}' for i in range(1, 81)]

    aff, sample = loadnii(fn_base + num_strings[0] + ext)
    y, x, w = np.squeeze(sample).shape

    data = np.zeros((y, x, len(num_strings), w))
    for i, num in enumerate(num_strings):
        data[:, :, i] = np.squeeze(loadnii(fn_base + num + ext)[1])

    return aff, data


def get_hw(data, level):
    if data.ndim != 2:
        raise ValueError('Input data must be already masked / two-dimensional')

    # Computing average spectrum
    mean = data.mean(0)[5:-5]  # empirically disregarding tails at edges

    # Setting the threshold to "level" percent of the total peak-floor
    # height above the floor
    thresh = level / 100 * (mean.max() - mean.min()) + mean.min()
    above_thresh = np.where(mean >= thresh)[0]

    v0 = above_thresh[0]
    vf = above_thresh[-1]
    vmax = np.argmax(mean)
    hw = np.max([abs(vmax - v0), abs(vmax - vf)])
    return hw


def asym(data, hw=20, method='sum', include_peak=True):
    n = data.shape[-1]
    if n % 2 == 0:
        n0 = n//2 - 1  # 95 for our data

    if include_peak:
        lo = data[..., n0-hw:n0+1]
        hi = data[..., n0:n0+1+hw]
        if method == 'sum':
            t = data[..., n0-hw:n0+1+hw].sum(-1)
        elif method == 'trapz':
            t = np.trapz(data[..., n0-hw:n0+1+hw])
    else:
        lo = data[..., n0-hw:n0]
        hi = data[..., n0+1:n0+1+hw]
        if method == 'sum':
            t = lo.sum(-1) + hi.sum(-1)
        elif method == 'trapz':
            t = np.trapz(lo) + np.trapz(hi)

    if method == 'sum':
        asym = (hi.sum(-1) - lo.sum(-1)) / t
    elif method == 'trapz':
        asym = (np.trapz(hi) - np.trapz(lo)) / t

    return asym


def default_signal_noise_masks(n=192, s0=95, ds=5, n0=60, dn=15):
    inds = np.arange(n)

    signal = abs(inds - s0) <= ds
    noise = (abs(inds - (s0-n0)) <= dn) + (abs(inds - (s0+n0)) <= dn)
    return signal, noise


def denoise_SSPC(X, mask=None, alpha=0.5e-3,
                 signalmask=None, noisemask=None,
                 print_rank=False):
    '''
    This function implements the selection of signal-related principal
    components (SSPC) denoising method for MRSI data.

    Abdoli, Abas, Stoyanova, Radka, Maudsley Andrew, "Denoising of MR
    spectroscopic imaging data using statistical selection of principal
    components," Magnetic Resonance Materials in Physics, Biology and
    Medicine, Vol 29, Issue 6, 2016.


    Parameters
    _________
    data : array
        N x w spectroscopic MR data, where N is the number of
        voxels and w is the number of frequency points. If the array
        is higher dimensional, use mask to select voxels.
    mask : array
        Used to exclude background voxels if data is not already masked

    '''

    # Check if data has empty dimension
    X = np.squeeze(X)
    shape = X.shape

    # Mask if applicable
    if (X.ndim > 2) & (mask is not None):
        mask = np.squeeze(mask)
        if mask.dtype != bool:
            mask = mask.astype(bool)

        X = X[mask]  # form data matrix

    _, Nw = X.shape

    # Fits data to PCs and returns coefficients
    pca = PCA(n_components=Nw)
    Z = pca.fit_transform(X)

    if (signalmask is None) & (noisemask is None):
        signalmask, noisemask = default_signal_noise_masks()

    # Levene test for equal variance between custom-defined signal
    # and noise frequency regions. Only keep potentially non-consecutive
    # PCs that meach significance alpha.
    is_sig = np.array([levene(pc[signalmask], pc[noisemask])
                       [1] <= alpha for pc in pca.components_])

    # Reconstruct spectra as weighted sum of significant PCs
    if mask is not None:
        denoised_data = np.zeros(shape)
        denoised_data[mask] = np.dot(Z[:, is_sig], pca.components_[
            is_sig, :]) + pca.mean_
    else:
        denoised_data = np.dot(Z[:, is_sig], pca.components_[
            is_sig, :]) + pca.mean_

    if print_rank:
        print(f'Rank: {is_sig.sum()}')

    return denoised_data
