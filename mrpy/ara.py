'''
Functions for interacting with the Allen Mouse Brain Connectivity Atlas
'''

import nrrd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.grid_data_api import GridDataApi
from os.path import expanduser

# Default manifest file location is user's home directory
home = expanduser('~')
manifest_file = home + '/mouse_connectivity/manifest.json'


def make_aff(val):
    '''
    Makes an affine matrix assuming isotropic voxels of length "val". For
    use in saving .nii files

    Parameters
    __________
    val : int or float
        Voxel size

    Returns
    _______
    aff : ndarray
        Affine matrix
    '''
    aff = np.eye(4)
    aff[:3, :3] = val * np.eye(3)
    return aff


def set_manifest(fn):
    '''
    Changes manifest file location from default (home directory) to custom

    Parameters
    __________
    fn : str
        Location of manifest file
    '''
    global manifest_file
    manifest_file = fn


def get_mcc(res=50, manifest_file=manifest_file):
    '''
    Fetches MouseConnectivityCache object to interact with Allen data. Defaults
    to setting manifest file to user's home directory.

    Parameters
    __________
    res : int
        Sets voxel size for Allen data. Must be 100, 50, 25, or 10.
        Default is 50.
    manifest_file : str
        Sets location of manifest file

    Returns
    _______
    MouseConnectivityCache object.
    '''
    if res not in [100, 50, 25, 10]:
        raise ValueError('Res must be 100, 50, 25, or 10')

    return MouseConnectivityCache(resolution=res, manifest_file=manifest_file)


def reorient_ara_data(data):
    '''
    Utility function for reorienting allen data to neurological display
    convention.

    Parameters
    __________
    data : ndarray
        3D image data

    Returns
    _______
    data_reoriented : ndarray
        Data reoriented to neurological display convention

    '''
    data_reoriented = np.flip(
        np.flip(np.moveaxis(data, [0, 1, 2], [1, 2, 0]), axis=2), axis=1)
    return data_reoriented


def get_structure_mask(acronym=None, res=50):
    '''
    Wraps allensdk to fetch structure mask given structure acronym.

    Parameters
    __________
    acronym : str
        Allen structure acronym
    res : int
        Sets voxel size for Allen data. Must be 100, 50, 25, or 10. Default is
        50.

    Returns
    _______
    aff : ndarray
        Affine matrix
    mask : ndarray
        Binary mask for given structure, in Allen CCF
    '''

    if res not in [100, 50, 25, 10]:
        raise ValueError('Res must be 100, 50, 25, or 10')

    mcc = get_mcc(res)
    tree = mcc.get_structure_tree()
    ID = tree.get_structures_by_acronym([acronym])[0]['id']
    mask, _ = mcc.get_structure_mask(ID)
    mask = reorient_ara_data(mask)
    aff = make_aff(res / 1000)
    return aff, mask.astype(float)


def get_injection_density(exp_id, res=50):
    '''
    Wraps allensdk to fetch injection density given experiment ID

    Parameters
    __________
    exp_id : int
        Allen experiment ID
    res : int
        Sets voxel size for Allen data. Must be 100, 50, 25, or 10. Default is
        50.

    Returns
    _______
    aff : ndarray
        Affine matrix
    inj : ndarray
        Injection density image for given experiment, in Allen CCF
    '''

    if res not in [100, 50, 25, 10]:
        raise ValueError('Res must be 100, 50, 25, or 10')

    mcc = get_mcc(res)
    inj, _ = mcc.get_injection_density(exp_id)
    inj = reorient_ara_data(inj)
    aff = make_aff(res / 1000)
    return aff, inj


def get_projection_density(exp_id, res=50):
    '''
    Wraps allensdk to fetch projection density given experiment ID

    Parameters
    __________
    exp_id : int
        Allen experiment ID
    res : int
        Sets voxel size for Allen data. Must be 100, 50, 25, or 10. Default is
        50.

    Returns
    _______
    aff : ndarray
        Affine matrix
    proj : ndarray
        Projection density image for given experiment, in Allen CCF
    '''

    if res not in [100, 50, 25, 10]:
        raise ValueError('Res must be 100, 50, 25, or 10')

    mcc = get_mcc(res)
    proj, _ = mcc.get_projection_density(exp_id)
    proj = reorient_ara_data(proj)
    aff = make_aff(res / 1000)
    return aff, proj


def get_projection_energy(exp_id, res=50):
    '''
    Wraps allensdk to fetch projection energy given experiment ID

    Parameters
    __________
    exp_id : int
        Allen experiment ID
    res : int
        Sets voxel size for Allen data. Must be 100, 50, 25, or 10. Default is
        50.

    Returns
    _______
    aff : ndarray
        Affine matrix
    energy : ndarray
        Projection energy image for given experiment, in Allen CCF
    '''

    if res not in [100, 50, 25, 10]:
        raise ValueError('Res must be 100, 50, 25, or 10')

    fn = manifest_file.split('manifest.json')[
        0] + f'/experiment_{exp_id}/projection_energy_{res}.nrrd'
    gda = GridDataApi(res)
    gda.download_projection_grid_data(exp_id, image=['projection_energy'],
                                      resolution=res,
                                      save_file_path=fn)
    energy, _ = nrrd.read(fn)
    energy = reorient_ara_data(energy)
    aff = make_aff(res / 1000)
    return aff, energy
