__author__ = 'milsteina'
from plot_results import *
import matplotlib.gridspec as gridspec


output_filename = '032017 discrete plasticity summary'

dt = 1.
back = 10000.
forward = 30000.
full_t = np.arange(-back, forward, dt)

track_length = 187.  # cm

binned_dx = track_length / 100.  # cm
binned_x = np.arange(0., track_length+binned_dx/2., binned_dx)[:100] + binned_dx/2.
generic_dx = binned_dx / 100.  # cm
generic_x = np.arange(0., track_length, generic_dx)

default_run_vel = 30.  # cm/s
generic_position_dt = generic_dx / default_run_vel * 1000.
generic_t = np.arange(0., len(generic_x)*generic_position_dt, generic_position_dt)[:len(generic_x)]

default_interp_t = np.arange(0., generic_t[-1], dt)
default_interp_x = np.interp(default_interp_t, generic_t, generic_x)


orig_font_size = mpl.rcParams['font.size']
orig_fig_size = mpl.rcParams['figure.figsize']
mpl.rcParams['font.size'] = 8.
mpl.rcParams['figure.figsize'] = 5.0236, 4.0236
fig = plt.figure()
gs = gridspec.GridSpec(2, 2)
ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[0, 1])

mean_vals = {}
mean_vals['long'] = [2.129E+02, 1.196E+03, 1.015E+02, 4.548E+02]
mean_vals['short'] = [10., 100., 10., 100.]

ramp = {'induced': {}, 'spont': {}}
induction_loc = {}
induction_dur = {}

model_ramp_type = 'model_ramp_simulation'

def subtract_baseline(ramp):
    """

    :param ramp: array
    :return: array
    """
    baseline_indexes = np.where(ramp <= np.percentile(ramp, 10.))[0]
    baseline = np.mean(ramp[baseline_indexes])
    return ramp - baseline


with h5py.File(data_dir+output_filename+'.hdf5', 'r') as f:
    rule = 'long'
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

    local_kernel = f['short']['1']['local_kernel'][:]
    local_kernel /= np.max(local_kernel)
    local_peak_index = np.where(local_kernel == 1.)[0][0]
    short_local_kernel = np.zeros_like(full_t)
    short_local_kernel[int(back / dt) - local_peak_index:int(back / dt) - local_peak_index + len(local_kernel)] += \
        local_kernel
    global_kernel = f['short']['1']['global_kernel'][:]
    global_kernel /= np.max(global_kernel)
    global_peak_index = np.where(global_kernel == 1.)[0][0]
    short_global_kernel = np.zeros_like(full_t)
    short_global_kernel[int(back / dt) - global_peak_index:int(back / dt) - global_peak_index +
                                                          len(global_kernel)] += global_kernel
    for rule in ['long', 'long_spont']:
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
            ax0.plot(full_t/1000., this_local_kernel, color='k', alpha=0.1)
            ax1.plot(full_t/1000., this_global_kernel, color='k', alpha=0.1)
    for rule, cell_id in [('induced', '5'), ('spont', '6')]:
        if rule == 'induced':
            long_group = 'long'
            short_group = 'short'
        else:
            long_group = 'long_spont'
            short_group = 'short_spont'
        induction_loc[rule] = f[long_group][cell_id].attrs['induction_loc']
        induction_dur[rule] = f[long_group][cell_id].attrs['induction_dur']
        ramp[rule]['exp'] = f[long_group][cell_id]['exp_ramp'][:]
        ramp[rule]['long'] = subtract_baseline(f[long_group][cell_id][model_ramp_type][:])
        ramp[rule]['short'] = subtract_baseline(f[short_group][cell_id][model_ramp_type][:])


ax0.plot(full_t/1000., short_local_kernel, color='r', label='standard rule')
ax1.plot(full_t/1000., short_global_kernel, color='r')
ax0.plot(full_t/1000., mean_local_kernel, color='k', label='fit to data')
ax1.plot(full_t/1000., mean_global_kernel, color='k')
ax0.set_xlim(-1, 5.)
ax1.set_xlim(-0.5, 2.)
ax0.set_xlabel('time (s)')
ax0.set_ylabel('kernel ampl. (norm.)')
ax0.set_title('local signal kernels', fontsize=mpl.rcParams['font.size'])
ax0.legend(loc='best', frameon=False, framealpha=0.5)
ax1.set_xlabel('time (s)')
ax1.set_ylabel('kernel ampl. (norm.)')
ax1.set_title('global signal kernels', fontsize=mpl.rcParams['font.size'])
clean_axes(ax0)
clean_axes(ax1)

rule = 'induced'
ymin = np.min(ramp[rule].values())
ymax = np.max(ramp[rule].values())
start_index = np.where(default_interp_x >= induction_loc[rule])[0][0]
end_index = start_index + int(induction_dur[rule] / dt)
x_start = induction_loc[rule] / track_length
x_end = default_interp_x[end_index] / track_length
axes2 = plt.subplot(gs[1, 0])
axes2.plot(binned_x, ramp[rule]['exp'], c='b', label='exp')
axes2.plot(binned_x, ramp[rule]['long'], c='k', label='fit to data')
axes2.plot(binned_x, ramp[rule]['short'], c='r', label='standard rule')
axes2.axhline(y=ymax*1.05, xmin=x_start, xmax=x_end, linewidth=1, c='k')
axes2.set_xlabel('position (cm)')
axes2.set_ylabel('Vm ramp ampl. (mV)')
axes2.set_ylim(ymin, ymax*1.1)
axes2.set_title('experimentally induced field', fontsize=mpl.rcParams['font.size'])
axes2.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes2)

rule = 'spont'
ymin = np.min(ramp[rule].values())
ymax = np.max(ramp[rule].values())
start_index = np.where(default_interp_x >= induction_loc[rule])[0][0]
end_index = start_index + int(induction_dur[rule] / dt)
x_start = induction_loc[rule] / track_length
x_end = default_interp_x[end_index] / track_length
axes3 = plt.subplot(gs[1, 1])
axes3.plot(binned_x, ramp[rule]['exp'], c='b')
axes3.plot(binned_x, ramp[rule]['long'], c='k')
axes3.plot(binned_x, ramp[rule]['short'], c='r')
axes3.axhline(y=ymax*1.05, xmin=x_start, xmax=x_end, linewidth=1, c='k')
axes3.set_xlabel('position (cm)')
axes3.set_ylabel('Vm ramp ampl. (mV)')
axes3.set_ylim(ymin, ymax*1.1)
axes3.set_title('naturally occurring field', fontsize=mpl.rcParams['font.size'])
clean_axes(axes3)



gs.tight_layout(fig)
plt.show()
plt.close()
