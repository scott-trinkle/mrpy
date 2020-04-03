'''
This tutorial will show you how to download and work with data from the Allen
Institute.

It will be helpful to reference the web-based Allen Brain Explorer for
reference:
https://connectivity.brain-map.org/3d-viewer?v=1
'''

# First we import the package
import mrpy as mr

# There are three main types of data we will need from the Allen Atlas. The
# first is a "structure mask". This gives a binary mask of a specific brain
# structure, for use in both tracking and in evaluating connectivity between
# different structures. To get a structure mask, we just need a structure
# acronym. These can be found in the Brain Explorer link above. Let's
# make a mask of the whole cerebellum. In the lower left panel of the website,
# I searched around and found that the abbreviation for the cerebellum is
# "CB", so here is how I call the function:
aff, CB_mask = mr.get_structure_mask(acronym='CB', res=50)

# A few comments. "res" is a keyword for all of the Allen functions that sets
# the spatial resolution (in microns) of the output. The default is 50, and
# that should be fine for almost everything that we do, so you can almost
# always ignore that argument and call the function without it. If you do want
# to change the resolution, the options are 10, 25, 50, and 100. The
# function is called like so:
aff, CB_mask = mr.get_structure_mask('CB')

# This function returns two things. CB_mask is the mask itself, it is a
# NumPy array with 1s where the cerebellum is, and 0s everywhere else. When
# we save these masks and load them into MRTrix for tracking, we will need to
# tell the software the spatial dimensions and resolution of the data. That
# information is contained in a separate NumPy array called "aff", which stands
# for "affine". This needs to be saved along with the data in a .nii file,
# which you can do with the following function:
mr.savenii(CB_mask, aff, 'CB_mask.nii.gz')

# Likewise, when we load a .nii file, it also returns the affine:
aff, CB_mask = mr.loadnii('CB_mask.nii.gz')


# The next type of data we will need is the injection density. This is a
# spatial map of where the tracer was injected into the mouse for a given
# Allen experiment. For this, we just need the Allen experiment ID, the
# 9-digit number that the Allen Institute uses to organize experiments. We will
# be getting this from the supplementary document of the Wu paper. You can
# also find the experiment ID for injections to a given structure online
# here: http://connectivity.brain-map.org. Here I am loading the injection
# density from the first experiment in the Wu document:
aff, injection = mr.get_injection_density(exp_id=100147861)

# "injection" returns a NumPy array where the value in each voxel is the percent
# of that voxel that contained tracer in the injection region. This can be used
# for creating tractography seeds using the "seed_rejection" option in
# mr.tckgen. Note that just like mr.get_structure_mask, we can also change
# the resolution if we want with the "res" keyword.

# The final function returns the projection density. This is the full tracer map
# from each experiment, which we will be comparing tractography to. It is
# called exactly the same as mr.get_injection_density:
aff, projection = mr.get_projection_density(exp_id=100147861)
