__author__ = 'milsteina'
from plot_results import *
import sys
import matplotlib.gridspec as gridspec


output_filename = '121316 plasticity rule optimization summary'

dt = 1.
back = 10000.
forward = 30000.
full_t = np.arange(-back, forward, dt)

fig = plt.figure()
gs = gridspec.GridSpec(3, 3)
ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[0, 1])

mean_vals = {}
mean_vals['long'] = [2.066E+02, 1.033E+03, 2.861E+01, 4.317E+02]
mean_vals['short'] = [10., 100., 2.301E+01, 1.240E+02]

for rule in ['short', 'long']:
    kernel_t = np.arange(0., 30000., dt)
    local_kernel = np.exp(-kernel_t/mean_vals[rule][1]) - np.exp(-kernel_t/mean_vals[rule][0])
    peak_index = np.where(local_kernel == np.max(local_kernel))[0][0]
    decay_indexes = np.where(local_kernel[peak_index:] < 0.005*np.max(local_kernel))[0]
    if np.any(decay_indexes):
        local_kernel = local_kernel[:peak_index+decay_indexes[0]]
    local_kernel /= np.max(local_kernel)
    local_peak_index = np.where(local_kernel == 1.)[0][0]
    mean_local_kernel = np.zeros_like(full_t)
    mean_local_kernel[int(back/dt)-local_peak_index:int(back/dt)-local_peak_index+len(local_kernel)] += local_kernel

    global_kernel = np.exp(-kernel_t/mean_vals[rule][3]) - np.exp(-kernel_t/mean_vals[rule][2])
    peak_index = np.where(global_kernel == np.max(global_kernel))[0][0]
    decay_indexes = np.where(global_kernel[peak_index:] < 0.005*np.max(global_kernel))[0]
    if np.any(decay_indexes):
        global_kernel = global_kernel[:peak_index+decay_indexes[0]]
    global_kernel /= np.max(global_kernel)
    global_peak_index = np.where(global_kernel == 1.)[0][0]
    mean_global_kernel = np.zeros_like(full_t)
    mean_global_kernel[int(back/dt)-global_peak_index:int(back/dt)-global_peak_index+len(global_kernel)] += global_kernel

    this_color = 'b' if rule == 'long' else 'r'

    with h5py.File(data_dir+output_filename+'.hdf5', 'r') as f:
        for cell in f[rule].values():
            local_kernel = cell['local_kernel'][:]
            local_kernel /= np.max(local_kernel)
            local_peak_index = np.where(local_kernel == 1.)[0][0]
            this_local_kernel = np.zeros_like(full_t)
            this_local_kernel[int(back/dt)-local_peak_index:int(back/dt)-local_peak_index+len(local_kernel)] += \
                local_kernel
            global_kernel = cell['global_kernel'][:]
            global_kernel /= np.max(global_kernel)
            global_peak_index = np.where(global_kernel == 1.)[0][0]
            this_global_kernel = np.zeros_like(full_t)
            this_global_kernel[int(back / dt) - global_peak_index:int(back / dt) - global_peak_index +
                                                                  len(global_kernel)] += global_kernel
            ax0.plot(full_t/1000., this_local_kernel, color=this_color, alpha=0.1)
            ax1.plot(full_t/1000., this_global_kernel, color=this_color, alpha=0.1)
    ax0.plot(full_t/1000., mean_local_kernel, color=this_color)
    ax1.plot(full_t/1000., mean_global_kernel, color=this_color)
ax0.set_xlim(-1, 5.)
ax1.set_xlim(-0.5, 2.)
ax0.set_xlabel('Time (s)')
ax0.set_ylabel('Normalized kernel amplitude')
ax0.set_title('Optimized local kernels', fontsize=mpl.rcParams['font.size'])
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Normalized kernel amplitude')
ax1.set_title('Optimized global kernels', fontsize=mpl.rcParams['font.size'])
clean_axes(ax0)
clean_axes(ax1)
gs.tight_layout(fig)
plt.show()
plt.close()
