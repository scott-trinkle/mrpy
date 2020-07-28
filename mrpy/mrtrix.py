import numpy as np
import subprocess


def tckgen(source, tracks, mel=False, **options):
    default = {'algorithm': 'iFOD2',
               'select': 500000,
               'step': 0.0375,
               'curvature': 22.5,
               'cutoff': 0.08,
               'minlength': 0.4,
               'nthreads': 32}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    bool_params = ['noprecomputed', 'rk4', 'stop', 'seed_unidirectional',
                   'backtrack', 'crop_at_gmwmi', 'info', 'quiet', 'debug',
                   'force', 'help', 'version']

    mrtrix_call = 'assign-nodes tckgen' if mel else 'tckgen'
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        elif (param == 'curvature'):
            angle_from_curve = 2 * np.arcsin(float(options['step'])
                                             * 1000 / (2 * val)) * 180 / np.pi
            if 'angle' not in options.keys():
                mrtrix_call += f' -angle {angle_from_curve}'
        else:
            mrtrix_call += f' -{param} {val}'

    mrtrix_call += f' {source} {tracks}'
    subprocess.run(mrtrix_call.split(' '))


def tckmap(tracks, output, mel=False, **options):
    default = {'nthreads': 32,
               'precise': True}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    bool_params = ['dec', 'map_zero', 'backtrack', 'precise',
                   'ends_only', 'info', 'quiet', 'debug', 'force', 'help',
                   'version']

    mrtrix_call = 'assign-nodes tckmap' if mel else 'tckmap'
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'

    mrtrix_call += f' {tracks} {output}'
    subprocess.run(mrtrix_call.split(' '))


def tcksift2(in_tracks, in_fod, out_weights, mel=False, **options):
    default = {'nthreads': 16}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    bool_params = ['fd_scale_gm', 'no_dilate_lut', 'make_null_lobes',
                   'remove_untracked', 'output_debug', 'linear',
                   'info', 'quiet', 'debug', 'force', 'help', 'version']

    mrtrix_call = 'assign-nodes tcksift2' if mel else 'tcksift2'
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'

    mrtrix_call += f' {in_tracks} {in_fod} {out_weights}'
    subprocess.run(mrtrix_call.split(' '))


def tck2connectome(tracks_in, nodes_in, connectome_ouit, mel=False, **options):
    default = {'nthreads': 16,
               'assignment_radial_search': 0.15,
               'vector': True}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    bool_params = ['assignment_end_voxels', 'assignment_all_voxels',
                   'scale_length', 'scale_invlength', 'scale_invnodevol',
                   'symmetric', 'zero_diagonal', 'keep_unassigned', 'vector',
                   'info', 'quiet', 'debug', 'force', 'help', 'version']

    mrtrix_call = 'assign-nodes tck2connectome' if mel else 'tck2connectome'
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'

    mrtrix_call += f' {tracks_in} {nodes_in} {connectome_out}'
    subprocess.run(mrtrix_call.split(' '))
