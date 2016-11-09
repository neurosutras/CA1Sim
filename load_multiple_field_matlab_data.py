__author__ = 'milsteina'
from plot_results import *
import scipy.io as io
import scipy.signal as signal
import sys

"""
Loads .mat file from MATLAB containing:
1) 'pos_cell': waves of timepoint (20 KHz) vs. position for each place field induction trial - need to be cleaned of
    nans, and used to create a wave of position vs. mean elapsed time.
2) 'vm_cell': waves of position vs. vm (V) - need to be scaled to mV, and plotted vs. either position or cumulative
    elapsed time

Induction locations:
newcell1: 125 - 7
newcell2: 128 - 32
newcell3: 90 - 0
newcell4: 150 - 65
newcell5: 88 - 45

"""

if len(sys.argv) > 1:
    filename1 = sys.argv[1]
if len(sys.argv) > 2:
    filename2 = sys.argv[2]
else:
    filename2 = None

mat_dir = 'katie/'

track_length = 190.
dx = track_length/100.
dt = 1./20000.

x = np.arange(0., track_length, dx)
t_x = np.arange(0., track_length+dx/2., dx)

field = []
vm = []
filtered_vm = []
t = []
delta_t = []
mean_delta_t = []
mean_vel = []
peak_loc = []
field.append(io.loadmat(data_dir+mat_dir+filename1+'.mat'))
if filename2 is None:
    field.append(io.loadmat(data_dir+mat_dir+filename1+'b.mat'))
else:
    field.append(io.loadmat(data_dir+mat_dir+filename2+'.mat'))

if filename1 == 'newcell1':
    induction_locs = [125., 7.]
elif filename1 == 'newcell2':
    induction_locs = [128., 32.]
elif filename1 == 'newcell3':
    induction_locs = [90., 0.]
elif filename1 == 'newcell4':
    induction_locs = [150., 65.]
elif filename1 == 'newcell5':
    induction_locs = [88., 45.]


def get_delta_t(pos):
    """
    Convert timepoint vs. position to position vs. elapsed time in bin. Renormalize the track length to 190 cm.
    :param pos: array
    :return: array
    """
    this_delta_t = []
    pos = np.array(pos[~np.isnan(pos)])
    pos /= max(pos)
    pos *= track_length
    for i, this_x in enumerate(x):
        start_x = this_x
        stop_x = this_x + dx
        start_index = np.where(np.array(pos) >= start_x)[0][0]
        stop_index = np.where(np.array(pos) < stop_x)[0][-1]
        this_dur = (stop_index - start_index) * dt
        this_delta_t.append(this_dur)
    return this_delta_t


for this_field in field:
    this_vm = this_field['vm_cell'][:,0]
    this_vm = this_vm[~np.isnan(this_vm)]
    this_vm = this_vm[:len(x)]
    this_vm *= 1000.
    vm.append(this_vm)
    # some peaks around periods of slow running velocity need less filtering
    if filename1 in ['newcell2']:
        filtered_vm.append(signal.savgol_filter(this_vm, 19, 4, mode='wrap'))
    elif filename1 in ['newcell1']:
        filtered_vm.append(signal.savgol_filter(this_vm, 21, 3, mode='wrap'))
    else:
        filtered_vm.append(signal.savgol_filter(this_vm, 51, 3, mode='wrap'))

if filename2 is None:
    difference = np.subtract(filtered_vm[1], filtered_vm[0])
else:
    difference = None

for this_field in field:
    delta_t.append([])
    for i in range(len(this_field['pos_cell'][0,:])):
        this_delta_t = get_delta_t(this_field['pos_cell'][:,i])
        delta_t[-1].append(this_delta_t)
    mean_delta_t.append(np.mean(delta_t[-1], axis=0))
    t.append(np.insert(np.cumsum(mean_delta_t[-1]), 0, 0.))
    mean_vel.append(np.mean(np.divide(dx, mean_delta_t[-1])))
print 'Mean vel: [%s]' % ('%.1f '*len(mean_vel) % tuple(mean_vel))

vel = []
smoothed_vel = []
for i in range(len(field)):
    this_vel = np.divide(dx, mean_delta_t[i])
    this_smoothed_vel = signal.savgol_filter(this_vel, 19, 4 , mode='wrap')
    vel.append(this_vel)
    smoothed_vel.append(this_smoothed_vel)

colors = ['c', 'g', 'r', 'k', 'purple']
x_start = [induction_loc/track_length for induction_loc in induction_locs]

for i in range(2):
    plt.plot(x, vm[i], color=colors[i])
    plt.plot(x, filtered_vm[i], color=colors[i+2])
    plt.show()
    plt.close()

fig, axes = plt.subplots(1)
ylim = max(np.max(vm), np.max(filtered_vm))
for i in range(2):
    axes.plot(x, vm[i], label='Mean vel: %.1f cm/s' % mean_vel[i], color=colors[i])
    axes.plot(x, filtered_vm[i], color=colors[i+2])
    axes.axhline(y=ylim + 0.3, xmin=x_start[i], xmax=x_start[i] + 0.02, c=colors[i], linewidth=3., zorder=0)
axes.set_ylabel('Vm (mV)')
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_title('Depolarization')
if difference is not None:
    axes.plot(x, np.add(difference, np.min(filtered_vm[0])), label='Difference', color=colors[-1])
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
for i in range(2):
    axes.plot(t[i], t_x, color=colors[i], label='Field induction '+str(i+1))
    axes.axvline(ymin=x_start[i], ymax=x_start[i]+0.02, x=0., c=colors[i], linewidth=3., zorder=0)
axes.set_ylim(-3., track_length)
axes.set_ylabel('Location (cm)')
axes.set_xlabel('Time (s)')
axes.set_xlim(-1., np.max(t))
axes.set_title('Position vs. Time')
clean_axes(axes)
axes.legend(loc='best', frameon=False, framealpha=0.5)
plt.show()
plt.close()

ylim = np.max(smoothed_vel)
max_vel = math.ceil(ylim/10.) * 10.
fig, axes = plt.subplots(1)
for i in range(len(vel)):
    # axes.plot(x, vel[i])
    axes.plot(x, smoothed_vel[i], label='Mean vel: %.1f cm/s' % mean_vel[i], color=colors[i])
    axes.axhline(y=ylim + 0.3, xmin=x_start[i], xmax=x_start[i] + 0.02, c=colors[i], linewidth=3., zorder=0)
axes.set_ylabel('Velocity (cm/s)')
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_ylim(0., max(40., max_vel))
axes.set_title('Running Speed')
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()


saved_filename = '110716 katie '+filename1+' saved output'
saved = {'t': t, 'ramp': filtered_vm, 'difference': difference, 'induction_locs': induction_locs}
write_to_pkl(data_dir+saved_filename+'.pkl', saved)
