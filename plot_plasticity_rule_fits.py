__author__ = 'milsteina'
from plot_results import *
import sys


if len(sys.argv) > 1:
    cell_id = str(sys.argv[1])
else:
    cell_id = '1'

output_filename = '121116 plasticity rule optimization summary'
with h5py.File(data_dir+output_filename+'.hdf5', 'r') as f:
    ramp = f['long'][cell_id]['ramp'][:]
    long_ramp = f['long'][cell_id]['model_ramp'][:]
    short_ramp = f['short'][cell_id]['model_ramp'][:]
    for this_ramp in ramp, long_ramp, short_ramp:
        baseline = np.mean(this_ramp[np.where(this_ramp <= np.percentile(this_ramp, 10.))[0]])
        this_ramp -= baseline
    track_length = f['long'][cell_id].attrs['track_length']
    induction_loc = np.mean(f['long'][cell_id].attrs['induction_loc'])
    x = np.arange(track_length / len(this_ramp) / 2., track_length, track_length / len(this_ramp))
    ylim = max(np.max(ramp), np.max(long_ramp), np.max(short_ramp))
    fig, axes = plt.subplots(2, 2)
    fig.set_size_inches(5.2, 3.9)
    axes[0][0].plot(x, ramp, color='k', label='Experiment')
    axes[0][0].plot(x, long_ramp, color='b', label='Long model')
    axes[0][0].plot(x, short_ramp, color='r', label='Short model')
    axes[0][0].axhline(y=ylim + 0.25, xmin=x_start, xmax=x_end, linewidth=2, c='k')
    axes[0][0].set_ylabel('Depolarization (mV)')
    axes[0][0].set_xlabel('Location (cm)')
    axes[0][0].set_xlim(0., track_length)
    axes[0][0].set_ylim(-0.5, ylim + 0.5)
    axes[0][0].set_title('Induced Vm ramp', fontsize=12.)
    clean_axes(axes[0])
    plt.tight_layout()
    plt.show()
    plt.close()