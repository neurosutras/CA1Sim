from plot_results import *

svg_title = '060716'
#svg_title = None
track_equilibrate = 300.
track_duration = 7500.
dt = 0.1

if svg_title is not None:
    remember_font_size = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 8

t = np.arange(0., track_duration, dt)
offset = int(track_equilibrate/dt)

peak_locs, ampa_forces, ampa_forces_sum, gaba_forces_sum = {}, {}, {}, {}
# length of ampa_forces is 2.*track_equilibrate + track_duration (extra padding is for filtering)
# length of ampa_forces_sum is track_duration (pre-filtered)
"""
for filename, condition in zip(['052416 - E_I distributions - biased density',
                                '052416 - E_I distributions - uniform density'], ['biased', 'uniform']):
    peak_locs[condition], ampa_forces[condition], ampa_forces_sum[condition] = read_from_pkl(data_dir+filename+'.pkl')

filtered = low_pass_filter(ampa_forces['uniform'][np.where((np.array(peak_locs['uniform']) >= 1000.) &
                                                           (np.array(peak_locs['uniform']) <= 1500.))[0][0]], 2.,
                           2.*track_equilibrate+track_duration, dt)
unit_baseline = np.max(filtered)
population_baseline = np.mean(ampa_forces_sum['uniform'][int(600./dt):int(1200./dt)])
"""

# length of gaba_forces is track_duration
for filename, condition in zip(['060716 - E_I distributions - subtr0_shape_inh1.0_modinh0',
                                '060716 - E_I distributions - subtr0_shape_inh1.0_modinh3',
                                '060716 - E_I distributions - subtr50_shape_inh0.5_modinh0',
                                '060716 - E_I distributions - subtr50_shape_inh0.5_modinh3',
                                '060716 - E_I distributions - subtr40_shape_inh2.0_modinh0',
                                '060716 - E_I distributions - subtr40_shape_inh2.0_modinh3'],
                               ['uniform0', 'uniform3', 'decreased0', 'decreased3', 'increased0', 'increased3']):
    gaba_forces_sum[condition] = read_from_pkl(data_dir+filename+'.pkl')

gaba_baseline = np.mean(gaba_forces_sum['uniform0'][int(600./dt):int(1200./dt)])

"""
fig, axes = plt.subplots(1)
for ampa_force in ampa_forces['uniform'][::100]:
    filtered = low_pass_filter(ampa_force, 2., 2.*track_equilibrate+track_duration, dt)
    axes.plot(t, filtered[offset:-offset]/unit_baseline)
axes.set_ylim(0., 3.)
axes.set_ylabel('Normalized Conductance')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
axes.set_title('Uniform input density', fontsize=mpl.rcParams['font.size'])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    fig.set_size_inches(2.05, 1.40)
    fig.savefig(data_dir+svg_title+' - uniform E - individual inputs.svg', format='svg', transparent=True)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
for ampa_force in ampa_forces['biased'][::100]:
    filtered = low_pass_filter(ampa_force, 2., 2.*track_equilibrate + track_duration, dt)
    axes.plot(t, filtered[offset:-offset]/unit_baseline)
for i, start in zip([1, 3, 5, 3, 1], np.arange(3000., 6000., 600.)):
    end = start + 600.
    interval = 600./float(i+1)
    for j in range(1, i+1):
        target = start + interval * float(j)
        index = np.where((np.array(peak_locs['biased'])>=target-10.)&(np.array(peak_locs['biased'])<=target+10.))[0][0]
        ampa_force = ampa_forces['biased'][index]
        filtered = low_pass_filter(ampa_force, 2., 2.*track_equilibrate+track_duration, dt)
        axes.plot(t, filtered[offset:-offset]/unit_baseline)
axes.set_ylim(0., 3.)
axes.set_ylabel('Normalized Conductance')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
axes.set_title('Biased input density', fontsize=mpl.rcParams['font.size'])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    fig.set_size_inches(2.05, 1.40)
    fig.savefig(data_dir+svg_title+' - biased E - individual inputs.svg', format='svg', transparent=True)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
axes.plot(t, ampa_forces_sum['uniform']/population_baseline, c='k')
axes.set_ylim(0., 2.5)
axes.set_ylabel('Normalized Conductance')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
#axes.set_title('Biased input density', fontsize=mpl.rcParams['font.size'])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    fig.set_size_inches(2.05, 1.40)
    fig.savefig(data_dir+svg_title+' - uniform E - total input.svg', format='svg', transparent=True)
plt.show()
plt.close()

fig, axes = plt.subplots(1)
axes.plot(t, ampa_forces_sum['biased']/population_baseline, c='r')
axes.set_ylim(0., 2.5)
axes.set_ylabel('Normalized Conductance')
axes.set_xlim(0., 7500.)
axes.set_xlabel('Time (s)')
axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
#axes.set_title('Biased input density', fontsize=mpl.rcParams['font.size'])
#axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
clean_axes(axes)
axes.tick_params(direction='out')
if svg_title is not None:
    fig.set_size_inches(2.05, 1.40)
    fig.savefig(data_dir+svg_title+' - biased E - total input.svg', format='svg', transparent=True)
plt.show()
plt.close()
"""
for (control, reduced), title in zip([('uniform0', 'uniform3'), ('increased0', 'increased3'),
                                      ('decreased0', 'decreased3')], ['Uniform inhibition',
                                                                      'Increased inhibition in field',
                                                                      'Decreased inhibition in field']):
    fig, axes = plt.subplots(1)
    axes.plot(t, gaba_forces_sum[control]/gaba_baseline, c='k', label='Control')
    axes.plot(t, gaba_forces_sum[reduced]/gaba_baseline, c='orange', label='Reduced inhibition')
    axes.set_ylim(0., 2.5)
    axes.set_ylabel('Normalized conductance')
    axes.set_xlim(0., 7500.)
    axes.set_xlabel('Time (s)')
    axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
    axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
    axes.set_title(title, fontsize=mpl.rcParams['font.size'])
    axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if svg_title is not None:
        fig.set_size_inches(2.05, 1.40)
        fig.savefig(data_dir+svg_title+' - I distribution - '+title+'.svg', format='svg', transparent=True)
    plt.show()
    plt.close()


"""
dt = 0.02
rec_filenames, residuals, mean_theta_envelope, scatter_vm_mean, scatter_vm_var, binned_t, mean_binned_vm, \
mean_binned_var, mean_ramp, mean_output = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}

for rec_filename, group in zip(['060716 - cell140 - subtr_inh0_shape_inh1.0 - saved parameters',
                                '060716 - cell128 - subtr_inh50_shape_inh0.5 - saved parameters',
                                '060716 - cell129 - subtr_inh40_shape_inh2.0 - saved parameters',
                                '060716 - cell130 - subtr_inh60_density_shape_inh1.0 - saved parameters',
                                '060716 - cell139 - subtr_inh80_density_shape_inh0.5 - saved parameters',
                                '060716 - cell138 - subtr_inh70_density_shape_inh2.0 - saved parameters'],
                               ['shape_inh1.0', 'shape_inh0.5', 'shape_inh2.0', 'density_shape_inh1.0',
                                'density_shape_inh0.5', 'density_shape_inh2.0']):
    rec_filenames[group], rec_t, residuals[group], mean_theta_envelope[group], scatter_vm_mean[group], \
    scatter_vm_var[group], binned_t[group], mean_binned_vm[group], mean_binned_var[group], mean_ramp[group], \
    mean_output[group] = read_from_pkl(data_dir+rec_filename+'.pkl')

for group in ['shape_inh1.0', 'shape_inh0.5', 'shape_inh2.0', 'density_shape_inh1.0', 'density_shape_inh0.5',
                  'density_shape_inh2.0']:
    baseline = np.mean(mean_ramp[group]['modinh0'][:int(600. / dt)])
    fig, axes = plt.subplots(1)
    axes.plot(rec_t, np.subtract(mean_ramp[group]['modinh0'], baseline), color='k')
    axes.plot(rec_t, np.subtract(mean_ramp[group]['modinh3'], baseline), color='orange')
    axes.set_ylim(-0.8, 16.5)
    axes.set_ylabel('DVm (mV)')
    axes.set_xlim(0., 7500.)
    axes.set_xlabel('Time (s)')
    axes.set_xticks([0., 1500., 3000., 4500., 6000., 7500.])
    axes.set_xticklabels([0, 1.5, 3, 4.5, 6, 7.5])
    axes.set_title(group, fontsize=mpl.rcParams['font.size'])
    #axes.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if svg_title is not None:
        fig.set_size_inches(2.05, 1.40)
        fig.savefig(data_dir + svg_title + ' - no_na ramp - ' + group + '.svg', format='svg', transparent=True)
    plt.show()
    plt.close()
"""

if svg_title is not None:
    mpl.rcParams['font.size'] = remember_font_size
