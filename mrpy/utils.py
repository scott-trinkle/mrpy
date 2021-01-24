import pickle
import nibabel as nib
import numpy as np
from skimage.io import imsave


def loadnii(fn):
    '''
    Loads NifTi files and returns the affine matrix and data array.

    Parameters
    __________
    fn : str
        Image filename

    Returns
    _______
    aff : ndarray
        Affine matrix
    data : ndarray
        Image data
    '''

    if fn.split('.')[-1] not in ['.nii', '.nii.gz']:
        raise ValueError('Filename is not .nii or .nii.gz')

    img = nib.load(fn)
    return img.affine, img.get_data()


def savenii(data, aff, fn):
    '''
    Saves NifTi files.

    Parameters
    __________
    data : ndarray
        Image data
    aff : ndarray
        Affine matrix
    fn : str
        Image filename

    '''

    if fn.split('.')[-1] not in ['.nii', '.nii.gz']:
        raise ValueError('Filename is not .nii or .nii.gz')

    nib.save(nib.Nifti1Image(data, aff), fn)


def dec2tif(fn):
    '''
    Converts MRTrix3 DEC image to .tif for easier visualization in ImageJ

    Parameters
    __________
    fn : str
        Image filename
    '''
    img, data = loadnii(fn)
    data = np.moveaxis(data, [0, 1, 2, 3], [2, 0, 1, 3])
    data = np.flip(data, axis=1)
    data = (255 * (data - data.min()) /
            (data - data.min()).max()).astype(np.uint8)
    imsave(fn.split('.nii')[0] + '.tif', data)
