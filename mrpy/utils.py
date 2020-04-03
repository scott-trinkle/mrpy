import nibabel as nib


def loadnii(fn):
    img = nib.load(fn)
    return img.affine, img.get_data()


def savenii(data, aff, fn):
    nib.save(nib.Nifti1Image(data, aff), fn)
