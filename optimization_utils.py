from scipy.optimize import minimize, basinhopping
from scipy.signal import find_peaks
import numpy as np
import csv

def calculate_cycle_features(trace, dt):
    """
    Returns areas (list), spikes (list), and troughs (list) for all 5 theta cycles.
    - Areas:   AUC (mV*s) within each 150ms burst window
    - Spikes:  AP count within each 150ms burst window
    - Troughs: minimum mV in the 50ms window just before the next burst onset
 
    FIX: Now returns three plain lists instead of (list, dict-of-lists, dict-of-lists).
    The dict-of-single-element-lists was unnecessarily awkward and forced [0] unwrapping
    at every call site.
 
    FIX: Converts trace to ndarray once at the top, then slices into views — previously
    called np.array(..., copy=True) on every cycle, creating multiple redundant copies.
    """
    # Single conversion; zero-copy if trace is already an ndarray
    t = np.asarray(trace) if trace is not None else None
 
    areas   = []
    spikes  = []
    troughs = []
 
    for cycle_num, start in enumerate(BURST_STARTS_MS, start=1):
        s_idx = int(start / dt)
        e_idx = int((start + BURST_DURATION_MS) / dt)
 
        # --- AREA ---
        if t is None or len(t) <= s_idx:
            areas.append(0.0)
        else:
            segment      = t[s_idx:min(len(t), e_idx)]          # view, no copy
            segment_rect = np.clip(segment, 0, None)
            areas.append(np.trapz(segment_rect, dx=(dt / 1000.0)))
 
        # --- SPIKES ---
        if t is None or len(t) <= s_idx:
            spikes.append(0)
        else:
            segment  = t[s_idx:min(len(t), e_idx)]              # view, no copy
            AP_peaks = find_peaks(segment, height=-50, width=(10, 1000), distance=50, prominence=20)[0]
            spikes.append(len(AP_peaks))
 
        # --- TROUGH ---
        if cycle_num < 5:
            next_burst = BURST_STARTS_MS[cycle_num]
            t_end      = int(next_burst / dt)
            t_start    = int((next_burst - TROUGH_WINDOW_MS) / dt)
        else:
            t_start = int(1100.0 / dt)
            t_end   = int(1250.0 / dt)
 
        if t is None or len(t) <= t_start:
            troughs.append(0.0)
        else:
            trough_seg = t[t_start:min(len(t), t_end)]          # view, no copy
            troughs.append(float(np.min(trough_seg)))
 
    return areas, spikes, troughs

def clear_all_synapses():
    for node in cell.tree:
        node.content['synapses'] = []
    for pathway in stim_exc_syns: stim_exc_syns[pathway] = []
    for pathway in stim_inh_syns: stim_inh_syns[pathway] = []



def update_calcium_and_potassium(m_cal, m_calH, m_car, m_cat, m_mykca, m_kca):
    cell.modify_mech_param('soma', 'cal', 'gcalbar', value=0.007 * m_cal) #optimize soma cal
    for sec in ['trunk', 'basal']:
        cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='soma')
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='trunk')

    cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.000031635 * m_calH, max_loc=50.0, origin='soma') #optimize soma calH
    cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.001455 * m_calH, min_loc=50.0, origin='soma', replace=False) #optimize trunk calH
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'calH', 'gcalbar', origin='trunk')

    cell.modify_mech_param('soma', 'car', 'gcabar', value=0.003 * m_car) #optimize soma car
    cell.modify_mech_param('trunk', 'car', 'gcabar', value=0.00003 * m_car) #optimize trunk car
    cell.modify_mech_param('basal', 'car', 'gcabar', origin='soma')
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'car', 'gcabar', origin='trunk')

    #The max value of the trunk is the free param - then compute the slope from max value
    cell.modify_mech_param('soma', 'cat', 'gcatbar', value=0.00005 * m_cat)
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0 * m_cat, max_loc=100.0, origin='soma') 
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0, slope=1.143e-06 * m_cat, min_loc=100.0, max_loc=350.0, origin='soma', replace=False) 
    cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0004 * m_cat, min_loc=350.0, origin='soma', replace=False) 
    cell.modify_mech_param('basal', 'cat', 'gcatbar', origin='soma')
    for sec in ['apical', 'tuft']:
        cell.modify_mech_param(sec, 'cat', 'gcatbar', origin='trunk')

    #The max value of the trunk is the free param - then compute the slope from max value
    cell.modify_mech_param('soma', 'mykca', 'gkbar', value=0.09075 * m_mykca) #optimize soma mykca
    cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, max_loc=50.0, origin='soma')
    cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, slope=-0.0005543 * m_mykca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.007585 * m_mykca, min_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('tuft', 'mykca', 'gkbar', value=0.007585 * m_mykca)
    cell.modify_mech_param('basal', 'mykca', 'gkbar', origin='soma')
    cell.modify_mech_param('apical', 'mykca', 'gkbar', origin='trunk')

    #The max value of the trunk is the free param - then compute the slope from max value
    cell.modify_mech_param('soma', 'kca', 'gbar', value=0.0005 * m_kca)
    cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005 * m_kca, max_loc=50.0, origin='soma')
    cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005 * m_kca, slope=-3.056e-6 * m_kca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('trunk', 'kca', 'gbar', value=4.167e-05 * m_kca, min_loc=200.0, origin='soma', replace=False)
    cell.modify_mech_param('tuft', 'kca', 'gbar', value=4.167e-05 * m_kca)
    cell.modify_mech_param('basal', 'kca', 'gbar', origin='soma')
    cell.modify_mech_param('apical', 'kca', 'gbar', origin='trunk')

    all_secs  = ['soma', 'trunk', 'basal', 'apical', 'tuft']
    all_mechs = ['cal', 'calH', 'car', 'cat', 'mykca', 'kca', 'pas']
    for sec in all_secs:
        for mech in all_mechs:
            cell.reinitialize_subset_mechanisms(sec, mech)