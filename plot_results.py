__author__ = 'milsteina'
from function_lib import *
import numpy as np


def plot_AR(rec_filename):
    """
    Plots a graph of Amplitude Ratio (spine to branch) vs. distance from primary branch for each dendritic sec_type.
    :param rec_filename: str
    """
    f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
    sec_types = []
    distances = {}
    AR = {}
    dendR = {}
    neckR = {}
    amp = f['0'].attrs['amp']
    equilibrate = f['0'].attrs['equilibrate']
    duration = f['0'].attrs['duration']
    simiter = 0
    while simiter < len(f):
        if f[str(simiter)].attrs['stim_loc'] == 'spine':
            spine_stim = f[str(simiter)]['rec']
            spine_tvec = f[str(simiter)]['time']
            branch_stim = f[str(simiter+1)]['rec']
            branch_tvec = f[str(simiter+1)]['time']
        elif f[str(simiter)].attrs['stim_loc'] == 'branch':
            spine_stim = f[str(simiter+1)]['rec']
            spine_tvec = f[str(simiter+1)]['time']
            branch_stim = f[str(simiter)]['rec']
            branch_tvec = f[str(simiter)]['time']
        for rec in spine_stim.itervalues():
            if rec.attrs['description'] == 'branch':
                branch_rec = rec
                sec_type = rec.attrs['type']
            elif rec.attrs['description'] == 'spine':
                spine_rec = rec
        if not sec_type in sec_types:
            sec_types.append(sec_type)
            distances[sec_type] = []
            AR[sec_type] = []
            dendR[sec_type] = []
            neckR[sec_type] = []
        distances[sec_type].append(branch_rec.attrs['branch_distance'])
        left, right = time2index(branch_tvec[:], equilibrate-3.0, equilibrate-1.0)
        baseline_branch = np.average(branch_rec[left:right])
        left, right = time2index(branch_tvec[:], equilibrate, duration)
        peak_branch = np.max(branch_rec[left:right]) - baseline_branch
        left, right = time2index(spine_tvec[:], equilibrate-3.0, equilibrate-1.0)
        baseline_spine = np.average(spine_rec[left:right])
        left, right = time2index(spine_tvec[:], equilibrate, duration)
        peak_spine = np.max(spine_rec[left:right]) - baseline_spine
        this_AR = peak_spine / peak_branch
        AR[sec_type].append(this_AR)
        for rec in branch_stim.itervalues():
            if rec.attrs['description'] == 'branch':
                branch_rec = rec
                break
        left, right = time2index(branch_tvec[:], equilibrate-3.0, equilibrate-1.0)
        baseline_branch = np.average(branch_rec[left:right])
        left, right = time2index(branch_tvec[:], equilibrate, duration)
        peak_branch = np.max(branch_rec[left:right]) - baseline_branch
        this_dendR = peak_branch / amp
        dendR[sec_type].append(this_dendR)
        this_neckR = (this_AR - 1) * this_dendR
        neckR[sec_type].append(this_neckR)
        simiter += 2
    fig, axes = plt.subplots(3, len(sec_types))
    colors = ['b', 'g', 'r', 'c']
    sorted_distances = {}
    for i, sec_type in enumerate(sec_types):
        indexes = range(len(distances[sec_type]))
        indexes.sort(key=distances[sec_type].__getitem__)
        sorted_distances[sec_type] = map(distances[sec_type].__getitem__, indexes)
        color = colors[i]
        axes[0][i].scatter(sorted_distances[sec_type], map(AR[sec_type].__getitem__, indexes), marker='o', color=color)
        axes[0][i].set_xlabel('Location (um)')
        axes[0][i].set_ylabel('Amplitude Ratio')
        axes[0][i].set_title(sec_type)
        axes[1][i].scatter(sorted_distances[sec_type], map(dendR[sec_type].__getitem__, indexes), marker='s',
                                                                                            color=color)
        axes[1][i].set_xlabel('Location (um)')
        axes[1][i].set_ylabel('R_Dend (MOhm)')
        axes[1][i].set_title(sec_type)
        axes[2][i].scatter(sorted_distances[sec_type], map(neckR[sec_type].__getitem__, indexes), marker='v',
                                                                                            color=color)
        axes[2][i].set_xlabel('Location (um)')
        axes[2][i].set_ylabel('R_Neck (MOhm)')
        axes[2][i].set_title(sec_type)
    plt.subplots_adjust(hspace=0.5)
    plt.show()
    plt.close()
    f.close()


def plot_spine_amp(rec_filename):
    """
    Plots a graph of spine EPSP amp and branch EPSP amp vs. distance from primary branch for each dendritic sec_type.
    :param rec_filename: str
    """
    f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
    sec_types = []
    distances = {}
    spine_amp = {}
    branch_amp = {}
    amp = f['0'].attrs['amp']
    equilibrate = f['0'].attrs['equilibrate']
    duration = f['0'].attrs['duration']
    simiter = 0
    while simiter < len(f):
        if f[str(simiter)].attrs['stim_loc'] == 'spine':
            spine_stim = f[str(simiter)]['rec']
            spine_tvec = f[str(simiter)]['time']
        elif f[str(simiter)].attrs['stim_loc'] == 'branch':
            spine_stim = f[str(simiter+1)]['rec']
            spine_tvec = f[str(simiter+1)]['time']
        for rec in spine_stim.itervalues():
            if rec.attrs['description'] == 'branch':
                branch_rec = rec
                sec_type = rec.attrs['type']
            elif rec.attrs['description'] == 'spine':
                spine_rec = rec
        if not sec_type in sec_types:
            sec_types.append(sec_type)
            distances[sec_type] = []
            spine_amp[sec_type] = []
            branch_amp[sec_type] = []
        distances[sec_type].append(branch_rec.attrs['branch_distance'])
        left, right = time2index(spine_tvec[:], equilibrate-3.0, equilibrate-1.0)
        baseline_branch = np.average(branch_rec[left:right])
        left, right = time2index(spine_tvec[:], equilibrate, duration)
        peak_branch = np.max(branch_rec[left:right]) - baseline_branch
        left, right = time2index(spine_tvec[:], equilibrate-3.0, equilibrate-1.0)
        baseline_spine = np.average(spine_rec[left:right])
        left, right = time2index(spine_tvec[:], equilibrate, duration)
        peak_spine = np.max(spine_rec[left:right]) - baseline_spine
        spine_amp[sec_type].append(peak_spine)
        branch_amp[sec_type].append(peak_branch)
        simiter += 2
    fig, axes = plt.subplots(2, len(sec_types))
    colors = ['b', 'g', 'r', 'c']
    sorted_distances = {}
    for i, sec_type in enumerate(sec_types):
        indexes = range(len(distances[sec_type]))
        indexes.sort(key=distances[sec_type].__getitem__)
        sorted_distances[sec_type] = map(distances[sec_type].__getitem__, indexes)
        color = colors[i]
        axes[0][i].scatter(sorted_distances[sec_type], map(spine_amp[sec_type].__getitem__, indexes), marker='o',
                           color=color)
        axes[0][i].set_xlabel('Location (um)')
        axes[0][i].set_ylabel('Spine Amp (mV)')
        axes[0][i].set_title(sec_type)
        axes[1][i].scatter(sorted_distances[sec_type], map(branch_amp[sec_type].__getitem__, indexes), marker='s',
                                                                                            color=color)
        axes[1][i].set_xlabel('Location (um)')
        axes[1][i].set_ylabel('Branch Amp (mV)')
        axes[1][i].set_title(sec_type)
    plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.05, right=0.98, top=0.95, bottom=0.05)
    plt.show()
    plt.close()
    f.close()


def plot_Rinp(rec_filename):
    """
    Produces a plot of input resistance vs. distance from primary branch for each dendritic subtype.
    :param rec_filename: str
    """
    f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
    sec_types = []
    distances = {}
    peak = {}
    steady = {}
    sag = {}
    amp = f['0']['stim']['0'].attrs['amp']
    start = f['0']['stim']['0'].attrs['delay']
    stop = start + f['0']['stim']['0'].attrs['dur']
    simiter = 0
    for sim in f.itervalues():
        rec = sim['rec']['0']
        sec_type = rec.attrs['type']
        if not sec_type in sec_types:
            sec_types.append(sec_type)
            distances[sec_type] = []
            peak[sec_type] = []
            steady[sec_type] = []
            sag[sec_type] = []
        distances[sec_type].append(rec.attrs['branch_distance'])
        tvec = sim['time']
        this_peak, this_steady = get_Rinp(tvec[:], rec[:], start, stop, amp)
        peak[sec_type].append(this_peak)
        steady[sec_type].append(this_steady)
        sag[sec_type].append(100*(1-this_steady/this_peak))
    rowlen = len(sec_types)
    fig, axes = plt.subplots(3, rowlen)
    colors = ['b', 'g', 'r', 'c']
    sorted_distances = {}
    for i, sec_type in enumerate(sec_types):
        indexes = range(len(distances[sec_type]))
        indexes.sort(key=distances[sec_type].__getitem__)
        sorted_distances[sec_type] = map(distances[sec_type].__getitem__, indexes)
        color = colors[i]
        axes[0][i].scatter(sorted_distances[sec_type], map(peak[sec_type].__getitem__, indexes), marker='o',
                           color=color)
        axes[0][i].set_xlabel('Location (um)')
        axes[0][i].set_ylabel('Input Resistance - Peak (MOhm)')
        axes[0][i].set_title(sec_type)
        axes[1][i].scatter(sorted_distances[sec_type], map(steady[sec_type].__getitem__, indexes), marker='s',
                                                                                            color=color)
        axes[1][i].set_xlabel('Location (um)')
        axes[1][i].set_ylabel('Input Resistance - Steady-state (MOhm)')
        axes[1][i].set_title(sec_type)
        axes[2][i].scatter(sorted_distances[sec_type], map(sag[sec_type].__getitem__, indexes), marker='v',
                                                                                            color=color)
        axes[2][i].set_xlabel('Location (um)')
        axes[2][i].set_ylabel('% Sag')
        axes[2][i].set_title(sec_type)
    plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.05, right=0.98, top=0.95, bottom=0.05)
    plt.show()
    plt.close()
    f.close()


def plot_superimpose_conditions(rec_filename):
    """
    File contains simulation results from iterating through some changes in parameters or stimulation conditions.
    This function produces one plot per recorded vector. Each plot superimposes the recordings from each of the
    simulation iterations.
    :param rec_filename: str
    """
    f = h5py.File(data_dir+rec_filename+'.hdf5', 'r')
    rec_ids = []
    sim_ids = []
    for sim in f.itervalues():
        if ('description' in sim.attrs) and not (sim.attrs['description'] in sim_ids):
            sim_ids.append(sim.attrs['description'])
    for rec in f['0']['rec'].itervalues():
        if ('description' in rec.attrs):
            rec_id = rec.attrs['description']
        else:
            rec_id = rec.attrs['type']+str(rec.attrs['index'])
        if not rec_id in (id['id'] for id in rec_ids):
            rec_ids.append({'id': rec_id, 'ylabel': rec.attrs['ylabel']+' ('+rec.attrs['units']+')'})
    if len(rec_ids) > 1:
        fig, axes = plt.subplots(1, len(rec_ids))
    else:
        fig, ax = plt.subplots(1, 1)
        axes = [ax]
    for i in range(len(rec_ids)):
        axes[i].set_xlabel('Time (ms)')
        axes[i].set_ylabel(rec_ids[i]['ylabel'])
        axes[i].set_title(rec_ids[i]['id'])
    for simiter in f:
        if 'description' in f[simiter].attrs:
            sim_id = f[simiter].attrs['description']
        else:
            sim_id = ''
        tvec = f[simiter]['time']
        for rec in f[simiter]['rec'].itervalues():
            if ('description' in rec.attrs):
                rec_id = rec.attrs['description']
            else:
                rec_id = rec.attrs['type']+str(rec.attrs['index'])
            i = [index for index, id in enumerate(rec_ids) if id['id'] == rec_id][0]
            axes[i].plot(tvec[:], rec[:], label=sim_id)
    for i in range(len(rec_ids)):
        axes[i].legend(loc='best')
    plt.subplots_adjust(hspace=0.3, wspace=0.3, left=0.05, right=0.98, top=0.95, bottom=0.05)
    plt.show()
    plt.close()
    f.close()