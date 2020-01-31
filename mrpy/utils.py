

import numpy as np

for i, num in enumerate(np.arange(.08, .28, .03 )):
            print(i, num)
            cutoff = num
def mrtrix(injection_rois, nstreamlines, cutoff, step, curvature ,odf_fn, output_fn, force=FALSE, algorithm = 'iFOD2'):
    angle = (2*np.arcsin(step/(2*curvature))*(180/np.pi))
    output = f'tckgen -seed_rejection {injection_rois} -select {nstreamlines} -cutoff {cutoff} -angle {angle} -step {step} -algorithm {algorithm} {odf_fn} {output_fn}'
    return output          

roi = 'C:/Users/laugh/Downloads/dataforchineze/mri/seed_images/I1.nii.gz'
streamlines = 1000
mystep = .15/4  # mm
mycurvature = 50 / 1000
myodf = 'C:/Users/laugh/Downloads/dataforchineze/mri/mri_odfs.nii.gz'
fn_list = ['various_cutoffs_{}.tck'.format(i) for i in range(6)]
my_output_fn = 'various_cutoffs_{}.tck'.format(i)
num = 0.1
forcestring = if '-force' is false else ''            


tckgen_track = mrtrix( injection_rois = roi, nstreamlines = streamlines, cutoff = num, step = mystep, curvature = mycurvature, odf_fn = myodf, output_fn = my_output_fn)
print(tckgen_track)
