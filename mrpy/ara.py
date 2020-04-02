import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pkg_resources
data_path = pkg_resources.resource_filename('mrpy', 'data/')


def make_aff(val):
    aff = np.eye(4)
    aff[:3, :3] = val * np.eye(3)
    return aff


def get_mcc(res=50):
    return MouseConnectivityCache(resolution=res,
                                  manifest_file=data_path + 'mouse_connectivity/manifest.json')


def reorient_ara_data(data):
    data = np.flip(
        np.flip(np.moveaxis(data, [0, 1, 2], [1, 2, 0]), axis=2), axis=1)
    return data


def get_structure_mask(acronym=None, res=50):
    mcc = get_mcc(res)
    tree = mcc.get_structure_tree()
    ID = tree.get_structures_by_acronym([acronym])[0]['id']
    mask, _ = mcc.get_structure_mask(ID)
    mask = reorient_ara_data(mask)
    aff = make_aff(res / 1000)
    return aff, mask.astype(float)


def get_injection_density(exp_id, res=50):
    mcc = get_mcc(res)
    inj, _ = mcc.get_injection_density(exp_id)
    inj = reorient_ara_data(inj)
    aff = make_aff(res / 1000)
    return aff, inj


def get_projection_density(exp_id, res=50):
    mcc = get_mcc(res)
    proj, _ = mcc.get_projection_density(exp_id)
    proj = reorient_ara_data(proj)
    aff = make_aff(res / 1000)
    return aff, proj
