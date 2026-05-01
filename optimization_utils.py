from scipy.optimize import minimize, basinhopping
from scipy.signal import find_peaks
import numpy as np
import csv

def calculate_cycle_features(trace, dt, burst_starts=None, burst_duration=None):
    """
    Returns areas (list), spikes (list), and troughs (list) for theta cycles.
    - Areas:   AUC (mV*s) within each burst window
    - Spikes:  AP count within each burst window
    - Troughs: minimum mV in the 50ms window just before the next burst onset
    """
    if burst_starts is None:
        burst_starts = [350.0, 550.0, 750.0, 950.0, 1150.0]
    if burst_duration is None:
        burst_duration = 150.0

    t = np.asarray(trace) if trace is not None else None
 
    areas   = []
    spikes  = []
    troughs = []
 
    for cycle_num, start in enumerate(burst_starts, start=1):
        s_idx = int(start / dt)
        e_idx = int((start + burst_duration) / dt)
 
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
        trough_window = 50.0  # ms before next burst to measure trough
        if cycle_num < len(burst_starts):
            next_burst = burst_starts[cycle_num]
            t_end      = int(next_burst / dt)
            t_start    = int((next_burst - trough_window) / dt)
        else:
            # Last cycle: measure trough in 50ms window after burst ends
            t_start = int((start + burst_duration) / dt)
            t_end   = int((start + burst_duration + trough_window) / dt)
 
        if t is None or len(t) <= t_start:
            troughs.append(0.0)
        else:
            trough_seg = t[t_start:min(len(t), t_end)]          # view, no copy
            troughs.append(float(np.min(trough_seg)) if len(trough_seg) > 0 else 0.0)
 
    return areas, spikes, troughs

def clear_all_synapses():
    for node in cell.tree:
        node.content['synapses'] = []
    for pathway in stim_exc_syns: stim_exc_syns[pathway] = []
    for pathway in stim_inh_syns: stim_inh_syns[pathway] = []





    