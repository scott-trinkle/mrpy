'''
Wrapper functions for the following MRTrix3 command line tools:

tckgen
tckmap
tcksift2
tck2connectome
'''

import numpy as np
import subprocess
import multiprocessing

n_cpus = multiprocessing.cpu_count()


def tckgen(source, tracks, mel=False, **options):
    '''
    Wrapper for the tckgen command line function in MRTri3. For details, see
    https://mrtrix.readthedocs.io/en/latest/reference/commands/tckgen.html

    Parameters
    __________
    source : str
        The image containing the source data. The type of image data
        required depends on the algorithm used.
    tracks : str
        The output file containing the tracks generated.
    mel : bool
        Appends 'assign-nodes' for multiprocessing on MEL workstation
    options : dict
        Additional arguments. Set flag options to bool values.

    '''

    # Set default options
    default = {'algorithm': 'iFOD2',
               'select': 500000,
               'step': 0.0375,
               'curvature': 22.5,
               'cutoff': 0.08,
               'minlength': 0.4,
               'nthreads': n_cpus}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    # List of bool parameters
    bool_params = ['noprecomputed', 'rk4', 'stop', 'seed_unidirectional',
                   'backtrack', 'crop_at_gmwmi', 'info', 'quiet', 'debug',
                   'force', 'help', 'version']

    mrtrix_call = 'assign-nodes tckgen' if mel else 'tckgen'

    # Constructs command line call string
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        elif (param == 'curvature'):
            # Allows for curvature input instead of angle
            angle_from_curve = 2 * np.arcsin(float(options['step'])
                                             * 1000 / (2 * val)) * 180 / np.pi
            if 'angle' not in options.keys():
                mrtrix_call += f' -angle {angle_from_curve}'
        else:
            mrtrix_call += f' -{param} {val}'
    mrtrix_call += f' {source} {tracks}'

    # Calls function
    subprocess.run(mrtrix_call.split(' '))


def tckmap(tracks, output, mel=False, **options):
    '''
    Wrapper for the tckmap command line function in MRTri3. For details, see
    https://mrtrix.readthedocs.io/en/latest/reference/commands/tckmap.html

    Parameters
    __________
    tracks : str
        The input track file.
    output : str
        The output track-weighted image.
    mel : bool
        Appends 'assign-nodes' for multiprocessing on MEL workstation
    options : dict
        Additional arguments. Set flag options to bool values.

    '''

    # Set default options
    default = {'nthreads': n_cpus,
               'precise': True}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    # List of bool parameters
    bool_params = ['dec', 'map_zero', 'backtrack', 'precise',
                   'ends_only', 'info', 'quiet', 'debug', 'force', 'help',
                   'version']

    mrtrix_call = 'assign-nodes tckmap' if mel else 'tckmap'

    # Constructs command line call string
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'
    mrtrix_call += f' {tracks} {output}'

    # Calls function
    subprocess.run(mrtrix_call.split(' '))


def tcksift2(in_tracks, in_fod, out_weights, mel=False, **options):
    '''
    Wrapper for the tcksift2 command line function in MRTri3. For details, see
    https://mrtrix.readthedocs.io/en/latest/reference/commands/tcksift2.html

    Parameters
    __________
    in_tracks : str
        The input track file.
    in_fod : str
        Input image containing the spherical harmonics of the fibre
        orientation distributions
    out_weights : str
        Output text file containing the weighting factor for each streamline
    mel : bool
        Appends 'assign-nodes' for multiprocessing on MEL workstation
    options : dict
        Additional arguments. Set flag options to bool values.

    '''

    # Set default options
    default = {'nthreads': n_cpus}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    # List of bool parameters
    bool_params = ['fd_scale_gm', 'no_dilate_lut', 'make_null_lobes',
                   'remove_untracked', 'output_debug', 'linear',
                   'info', 'quiet', 'debug', 'force', 'help', 'version']

    mrtrix_call = 'assign-nodes tcksift2' if mel else 'tcksift2'

    # Constructs command line call string
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'
    mrtrix_call += f' {in_tracks} {in_fod} {out_weights}'

    # Calls function
    subprocess.run(mrtrix_call.split(' '))


def tck2connectome(tracks_in, nodes_in, connectome_out, mel=False, **options):
    '''
    Wrapper for the tck2connectome command line function in MRTri3. For details,
    see https://mrtrix.readthedocs.io/en/latest/reference/commands/tck2connectome.html

    Parameters
    __________
    tracks_in : str
        The input track file.
    nodes_in : str
        The input node parcellation image
    connectome_out : str
        The output .csv containing edge weights
    mel : bool
        Appends 'assign-nodes' for multiprocessing on MEL workstation
    options : dict
        Additional arguments. Set flag options to bool values.

    '''

    # Set default options
    default = {'nthreads': n_cpus,
               'vector': True}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    # List of bool parameters
    bool_params = ['assignment_end_voxels', 'assignment_all_voxels',
                   'scale_length', 'scale_invlength', 'scale_invnodevol',
                   'symmetric', 'zero_diagonal', 'keep_unassigned', 'vector',
                   'info', 'quiet', 'debug', 'force', 'help', 'version']

    mrtrix_call = 'assign-nodes tck2connectome' if mel else 'tck2connectome'

    # Constructs command line call string
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'
    mrtrix_call += f' {tracks_in} {nodes_in} {connectome_out}'

    # Calls function
    subprocess.run(mrtrix_call.split(' '))
