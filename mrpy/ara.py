import nrrd
import numpy as np
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from allensdk.api.queries.grid_data_api import GridDataApi
from os.path import expanduser
home = expanduser('~')


def make_aff(val):
    aff = np.eye(4)
    aff[:3, :3] = val * np.eye(3)
    return aff


def get_mcc(res=50):
    return MouseConnectivityCache(resolution=res,
                                  manifest_file=home + '/mouse_connectivity/manifest.json')


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


def get_projection_energy(exp_id, res=50):
    fn = home + \
        f'/mouse_connectivity/experiment_{exp_id}/projection_energy_{res}.nrrd'
    gda = GridDataApi(res)
    gda.download_projection_grid_data(exp_id, image=['projection_energy'],
                                      resolution=res,
                                      save_file_path=fn)
    energy, _ = nrrd.read(fn)
    energy = reorient_ara_data(energy)
    aff = make_aff(res / 1000)
    return aff, energy
