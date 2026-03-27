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
