


import numpy as np


for i, num in enumerate(np.arange(.08, .28, .03 )):
            print(i, num)
            cutoff = num
def mrtrix(injection_rois, nstreamlines, cutoff, step, curvature ,odf_fn, output_fn, force=False, algorithm = 'iFOD2', nthreads = 2):
    angle = (2*np.arcsin(step/(2*curvature))*(180/np.pi))
    output = f'tckgen -seed_rejection {injection_rois} -select {nstreamlines} -cutoff {cutoff} -angle {angle} -step {step} -force {} -algorithm {algorithm} -nthreads {} {odf_fn} {output_fn}'
    return output       



roi = 'C:/Users/laugh/Downloads/dataforchineze/mri/seed_images/I1.nii.gz'
streamlines = 1000
mystep = .15/4  # mm
mycurvature = 50 / 1000
myodf = 'C:/Users/laugh/Downloads/dataforchineze/mri/mri_odfs.nii.gz'
fn_list = ['various_cutoffs_{}.tck'.format(i) for i in range(6)]
my_output_fn = 'various_cutoffs_{}.tck'.format(i)
num = 0.1
nthreads = 2
forcestring = '-force' if force else ''  
minlength = 0.3
maxlength = 35


tckgen_track = mrtrix(injection_rois = roi, nstreamlines = streamlines, cutoff = num, step = mystep, curvature = mycurvature, odf_fn = myodf, output_fn = my_output_fn)
print(tckgen_track)


def tckmap(template_image):
    output = f'tckmap -template {template_image} '
    return output    

temp_img = 'test_track.tck'

tckmap_img = mrtrix(template_image = temp_img)
print(tckmap_img)
