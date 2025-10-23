# 22/05/2024, B. Giblin, Hawking Fellow, Edinburgh
# This script creates a shuffle matrix to use when you
# have lots of noise realisations. It uses the outputs from:
# Shuffle_SLICS-K1000-Cov.py 

# How does it work?
# If you have 10 noise realisations,
# it reads in 10 of the [217,18] shuffle  matrices,
# where each one has a different shuffle config,
# and stacks them making something 3D.

# SO... if the first (of 10 layers) looks like this:

# [ 74, 75, 76, ...., 91 ]           (18 columns)
# [ 75, 76, 77, ...., 92 ]
# [ ................ ]
# [ ................ ]
# [ 291,292, 74, ......, 89 ]
# [ 292, 74, 75, ......, 90 ] 

# it means survey realn 0 (first row) consists of:
# LOS74R1n0, LOS75R2n0, LOS76R3n0, ...., LOS91R18n0
# And the next layer down will all have n=1, but different LOS for each R.
# Dont worry: n=0 maps do not have same shape noise! LOS, R, and zbin all used
# to seed the shape noise realn.

import numpy as np

Nnoise = 10 # Num. of noise realns (how many matrices to stack)
config = 0  # What shuffle config is used for first layer (n=0 surveys)

# SM stands for shuffle matrix (here we start stacking).
m0 = np.load('Shuffled_SLICS-K1000-Mosaic_Configs/Shuffled_Config%s.npy' %config )
SM = np.zeros([ Nnoise, m0.shape[0], m0.shape[1] ])
for nn in range(Nnoise):
    SM[nn] = np.load('Shuffled_SLICS-K1000-Mosaic_Configs/Shuffled_Config%s.npy' %(config+nn) )

np.save('Shuffled_SLICS-K1000-Mosaic_Configs/Shuffled_Config%s_Nnoise%s.npy' %(config,Nnoise), SM)

