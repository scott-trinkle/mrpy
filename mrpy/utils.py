import nibabel as nib
import numpy as np
from skimage.io import imsave


def loadnii(fn):
    img = nib.load(fn)
    return img.affine, img.get_data()


def savenii(data, aff, fn):
    nib.save(nib.Nifti1Image(data, aff), fn)


def dec2tif(fn):
    img, data = loadnii(fn)
    data = np.moveaxis(data, [0, 1, 2, 3], [2, 0, 1, 3])
    data = np.flip(data, axis=1)
    data = (255 * (data - data.min()) /
            (data - data.min()).max()).astype(np.uint8)
    imsave(fn.split('.nii')[0] + '.tif', data)
