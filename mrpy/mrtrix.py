import numpy as np
import subprocess


def tckgen(source, tracks, **options):
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

    mrtrix_call = 'tckgen'
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        elif param == 'curvature':
            angle = 2 * np.arcsin(float(options['step'])
                                  * 1000 / (2 * val)) * 180 / np.pi
            mrtrix_call += f' -angle {angle}'
        else:
            mrtrix_call += f' -{param} {val}'

    mrtrix_call += f' {source} {tracks}'
    subprocess.run(mrtrix_call.split(' '))


def tckmap(tracks, output, **options):
    default = {'nthreads': 32,
               'precise': True}
    for default_param, default_val in default.items():
        if default_param not in options.keys():
            options[default_param] = default_val

    bool_params = ['dec', 'map_zero', 'backtrack', 'precise',
                   'ends_only', 'info', 'quiet', 'debug', 'force', 'help',
                   'version']

    mrtrix_call = 'tckmap'
    for param, val in options.items():
        if param in bool_params:
            mrtrix_call += f' -{param}' if val else ''
        else:
            mrtrix_call += f' -{param} {val}'

    mrtrix_call += f' {tracks} {output}'
    subprocess.run(mrtrix_call.split(' '))
