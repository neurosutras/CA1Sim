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
    filtered_vm.append(signal.savgol_filter(this_vm, 51, 3))
    this_peak_loc = x[np.where(this_vm == max(this_vm))[0][0]]
    peak_loc.append(this_peak_loc)

if filename2 is None:
    difference = np.add(np.subtract(filtered_vm[1], filtered_vm[0]), np.min(filtered_vm[0]))
else:
    difference = None

for j, this_field in enumerate(field):
    delta_t.append([])
    for i in range(len(this_field['pos_cell'][0,:])):
        this_delta_t = get_delta_t(this_field['pos_cell'][:,i])
        delta_t[-1].append(this_delta_t)
    mean_delta_t.append(np.mean(delta_t[-1], axis=0))
    t.append(np.insert(np.cumsum(mean_delta_t[-1]), 0, 0.))
    mean_vel.append(np.mean(np.divide(dx, mean_delta_t[-1])))
print 'Mean vel:', mean_vel

fig, axes = plt.subplots(1)
for i in range(2):
    axes.plot(x, vm[i], label='Vel: %.1f; Loc: %.1f' % (mean_vel[i], peak_loc[i]))
    axes.plot(x, filtered_vm[i])
axes.set_ylabel('Vm (mV)')
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_title('Depolarization')
if difference is not None:
    axes.plot(x, difference, label='Difference')
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()

for i in range(2):
    plt.plot(t[i], t_x, label=str(i+1))
plt.legend(loc='best', frameon=False, framealpha=0.5)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
for i in range(2):
    axes.plot(x, np.divide(dx, mean_delta_t[i]), label='Vel: %.1f; Loc: %.1f' % (mean_vel[i], peak_loc[i]))
axes.set_ylabel('Velocity (cm/s)')
axes.set_xlabel('Location (cm)')
axes.set_xlim(0., track_length)
axes.set_title('Running Speed')
axes.legend(loc='best', frameon=False, framealpha=0.5)
clean_axes(axes)
plt.show()
plt.close()

