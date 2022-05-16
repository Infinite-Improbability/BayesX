import matplotlib
import numpy as np
from matplotlib import pyplot
import os
import sys

# Needs tidying up and generalizing.  Currently works on Coma data.

matplotlib.use('TkAgg')

try:
    mode = sys.argv[1]
except IndexError:
    mode = ''

while mode not in ['evts', 'bg', 'mask']:
    mode = input('Enter mode (evts/bg/mask): ')

if mode == 'evts':
    data1 = np.loadtxt('13996_evts_07_7keV.txt', usecols=(10, 11, 40))
    outfilename = '13996data_256by256by433.txt'
elif mode == 'mask':
    data1 = np.load('mask.npz')['arr_0']
    outfilename = 'mask_binned.txt'
    PI_chan_max = int(input('Enter number of energy channels: '))
    # This should case the energy mask file to match the structure of the binned data
    # The file is accordingly PI_chan_max times longer, but we avoid having to expand the
    # mask across all channels in BayesX, which could confusion in the future.
else:
    data1 = np.loadtxt('13996_bg_07_7keV.txt', usecols=(
        6, 7, 10))  # For the background file
    outfilename = '13996BG_256by256by433.txt'

sky_x = data1[:, 0]  # spatial coordinate
sky_y = data1[:, 1]  # spatial coordinate
# For data this is energy channel
# For mask this is 1 (masked) or 0 (unmasked)
PI_chan = data1[:, 2]

del data1

x_point = 4096.5
y_point = 4096.5

dx = sky_x-x_point
dy = sky_y-y_point
nbins = 256  # 1024

cellsize = 10.
xyrange = cellsize*(nbins - 1.)/2.

if mode == 'mask':
    PI_chan_min = 1
else:
    PI_chan_min = int(np.min(PI_chan))
    PI_chan_max = int(np.max(PI_chan))

PI_counts = np.zeros((nbins, nbins, PI_chan_max)) # [x, y, channel]

# Could also use scipy.stats.binned_statistic here but it's not too slow as is
icentre = nbins/2
jcentre = nbins/2
m = len(dx)
for ii in range(m):
    # Get cell coordinates, relative to centre
    u = dx[ii]
    v = dy[ii]
    # Get energy channel and convert to index by subtracting 1
    # Not used for masking.
    vr = int(PI_chan[ii])-1
    # Get bin indices
    signx = 1.
    if(u < 0.0):
        signx = -1
    signy = 1.
    if(v < 0.0):
        signy = -1
    i = int(np.floor(u/cellsize + signx*0.5) + icentre)
    j = int(np.floor(v/cellsize + signy*0.5) + jcentre)
    # If bin indices out of range then abort
    if((i < 0) | (i >= nbins)):
        continue
    if((j < 0) | (j >= nbins)):
        continue

    if mode == 'mask':
        # Masks all channels at spatial coordinates, 1 if masked, 0 otherwise
        # For masking, PI_chan[ii] = 1 implies that region (u,v) is masked.
        # We take the maximum compared with existing value so we can't accidentally overwrite the mask
        # The entire bin is masked if any pixel in it is masked.
        PI_counts[i, j, :] = np.maximum(PI_counts[i, j, :], PI_chan[ii])
    else:
        # Increments count for point
        PI_counts[i, j, vr] += 1

# Select inner 256x256 pixels for output
npix_out = 256
blc = int(icentre-npix_out/2)
trc = int(icentre+npix_out/2)
# PI_counts_out=PI_counts[PI_chan_min-1:,blc:trc,blc:trc]
PI_counts_out = PI_counts[blc:trc, blc:trc, PI_chan_min-1:]

pyplot.imshow(np.sum(PI_counts_out, axis=-1))
pyplot.show()

np.savetxt(outfilename, PI_counts_out.ravel(), '%5d')
#fprintf(fid,' %5d\n',count_obs);
