"""
generate_from_eval.py — Recover a full biophysics YAML from a specific optimization log entry.

Purpose:
  After identifying a promising evaluation (e.g. via plot_error_log.py), this script:
    1. Reads the CSV log to extract the 6 channel multipliers for that eval
    2. Builds a NEURON cell model using the base biophysics from the config YAML
    3. Applies the channel multipliers (cal, calH, car, cat, mykca, kca) to scale
       the corresponding ion channel conductances across all compartments
    4. Exports the resulting biophysics as a new YAML file (data/recovered_biophysics_eval_*.yaml)
    5. Saves a metadata YAML with the original log file, error, and feature values

  The exported YAML can then be loaded directly as mech_filename in the notebook or optimizer
  to reproduce that exact model configuration.

Usage:
  python generate_from_eval.py \\
    --log-file optimization_log_20260324_2200.csv \\
    --eval-num 232 \\
    --optimizer BasinHopping \\
    --config-path optimization_config_BasinHopping_20260324_5th_Run_Overnight.yaml
"""

import pandas as pd
import click
import os
import yaml
import datetime
from specify_cells import CA1_Pyr


# ========================================================================================
# CLI INTERFACE
# ========================================================================================
@click.command()
@click.option('--log-file', default='optimization_log_20260321_2046.csv',
              help='Path to the CSV log file.')
@click.option('--eval-num',
              help='The Evaluation Number (Eval column) you want to recover.')
@click.option('--optimizer', default='L-BFGS-B',
              help='The optimizer used for this run (e.g. L-BFGS-B, BasinHopping). Metadata only.')
@click.option('--config-path', default='optimization_config_L-BFGS-B_20260321.yaml',
              help='Path to the optimization config file (provides morph/mech filenames).')

# python generate_from_eval.py --log-file optimization_log_20260324_2200.csv --eval-num 232 --optimizer BasinHopping --config-path optimization_config_BasinHopping_20260324_5th_Run_Overnight.yaml
def recover_from_eval(log_file, eval_num, optimizer, config_path):
    """
    Reconstruct a full biophysics YAML from a single row of the optimization CSV.

    The 6 channel multipliers (cal, calH, car, cat, mykca, kca) are applied to the
    base conductance values exactly as update_calcium_and_potassium() does in the
    optimizer — ensuring the exported YAML matches what that evaluation actually ran.
    """

    # ── 1. SAFETY CHECKS ─────────────────────────────────────────────────────
    if not os.path.exists(log_file):
        print(f"Error: File '{log_file}' not found.")
        return
    if eval_num is None:
        print("Error: Please provide --eval-num [integer].")
        return

    # ── 2. LOAD CSV ───────────────────────────────────────────────────────────
    df = pd.read_csv(log_file)

    # ── 3. LOCATE THE SPECIFIC EVALUATION ─────────────────────────────────────
    row_data = df[df['Eval'] == int(eval_num)]
    if row_data.empty:
        print(f"Error: Evaluation {eval_num} not found in {log_file}.")
        return

    row = row_data.iloc[0]   # first matching row (should be exactly one)

    # ── 4. EXTRACT PARAMETERS & FEATURES ──────────────────────────────────────
    # Column groups — must match the CSV headers from optimize_soma_plateaus.py
    channel_cols = ['cal', 'calH', 'car', 'cat', 'mykca', 'kca']
    synapse_cols = ['n_exc_CA3', 'n_exc_ECIII', 'n_inh_CA3', 'n_inh_ECIII']
    area_cols    = ['Area_C1', 'Area_C2', 'Area_C3', 'Area_C4', 'Area_C5']
    spike_cols   = ['Spks_C1', 'Spks_C2', 'Spks_C3', 'Spks_C4', 'Spks_C5']
    trough_cols  = ['Trough_C1', 'Trough_C2', 'Trough_C3', 'Trough_C4', 'Trough_C5']

    # Build dicts — skip columns that don't exist in this CSV or are NaN
    m       = {k: float(row[k])          for k in channel_cols if k in row.index and not pd.isna(row[k])}
    syns    = {k: int(round(row[k]))     for k in synapse_cols if k in row.index and not pd.isna(row[k])}
    areas   = {k: round(float(row[k]), 4) for k in area_cols   if k in row.index and not pd.isna(row[k])}
    spikes  = {k: int(row[k])             for k in spike_cols  if k in row.index and not pd.isna(row[k])}
    troughs = {k: round(float(row[k]), 4) for k in trough_cols if k in row.index and not pd.isna(row[k])}

    # ── 5. PRINT RECOVERED FEATURES FOR CONFIRMATION ──────────────────────────
    print(f"\n--- Recovering Eval #{int(row['Eval'])} | Total Error: {row['Total_Error']:.4f} ---")
    print(f"Channels : {m}")
    print(f"Areas    : {[areas.get(k, '?') for k in area_cols]}")
    print(f"Spikes   : {[spikes.get(k, '?') for k in spike_cols]}")
    print(f"Troughs  : {[troughs.get(k, '?') for k in trough_cols]}")
    print(f"\nTargets:")
    print(f"  Areas  : [3.337, 3.5297, 4.356, 5.2235, 5.6729] mV*s")
    print(f"  Spikes : [6, 6, 7, 7, 7]")
    print(f"  Troughs: [1.2946, 3.1996, 8.0606, 17.4871, 22.5451] mV")

    # ── 6. BUILD NEURON MODEL & APPLY CHANNEL MULTIPLIERS ─────────────────────
    # This replicates exactly what update_calcium_and_potassium() does in the optimizer.
    # Each multiplier scales the BASE conductance for that channel family.
    print(f"\n--- Constructing Biophysical Model ---")
    try:
        with open(config_path, 'r') as f:
            sim_config = yaml.safe_load(f)

        morph_filename = sim_config['simulation']['morph_filename']
        mech_filename  = sim_config['simulation']['mech_filename']

        # Build cell from the BASE biophysics (before optimization)
        cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)

        # Extract multipliers (default to 1.0 = no change if missing)
        m_cal   = m.get('cal',   1.0)
        m_calH  = m.get('calH',  1.0)
        m_car   = m.get('car',   1.0)
        m_cat   = m.get('cat',   1.0)
        m_mykca = m.get('mykca', 1.0)
        m_kca   = m.get('kca',   1.0)

        # ─── cal (L-type calcium, Cav1.2) ──────────────────────────────
        # Soma sets the base; trunk/basal inherit via 'origin'; apical/tuft inherit from trunk
        cell.modify_mech_param('soma', 'cal', 'gcalbar', value=0.007 * m_cal)
        for sec in ['trunk', 'basal']:
            cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='soma')
        for sec in ['apical', 'tuft']:
            cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='trunk')

        # ─── calH (high-threshold L-type calcium) ─────────────────────
        # Piecewise: proximal trunk (< 50μm) vs. distal trunk (> 50μm)
        cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.000031635 * m_calH, max_loc=50.0, origin='soma')
        cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.001455 * m_calH, min_loc=50.0, origin='soma', replace=False)
        for sec in ['apical', 'tuft']:
            cell.modify_mech_param(sec, 'calH', 'gcalbar', origin='trunk')

        # ─── car (R-type calcium, Cav2.3) ─────────────────────────────
        cell.modify_mech_param('soma', 'car', 'gcabar', value=0.003 * m_car)
        cell.modify_mech_param('trunk', 'car', 'gcabar', value=0.00003 * m_car)
        cell.modify_mech_param('basal', 'car', 'gcabar', origin='soma')
        for sec in ['apical', 'tuft']:
            cell.modify_mech_param(sec, 'car', 'gcabar', origin='trunk')

        # ─── cat (T-type calcium, Cav3.x) ─────────────────────────────
        # Piecewise gradient: zero near soma, linear ramp 100-350μm, plateau beyond 350μm
        cell.modify_mech_param('soma', 'cat', 'gcatbar', value=0.00005 * m_cat)
        cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0 * m_cat, max_loc=100.0, origin='soma')
        cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0, slope=1.143e-06 * m_cat,
                               min_loc=100.0, max_loc=350.0, origin='soma', replace=False)
        cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0004 * m_cat,
                               min_loc=350.0, origin='soma', replace=False)
        cell.modify_mech_param('basal', 'cat', 'gcatbar', origin='soma')
        for sec in ['apical', 'tuft']:
            cell.modify_mech_param(sec, 'cat', 'gcatbar', origin='trunk')

        # ─── mykca (BK-type calcium-activated K+) ─────────────────────
        # High at soma, decreasing along trunk, low plateau distally
        cell.modify_mech_param('soma', 'mykca', 'gkbar', value=0.09075 * m_mykca)
        cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, max_loc=50.0, origin='soma')
        cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca,
                               slope=-0.0005543 * m_mykca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
        cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.007585 * m_mykca,
                               min_loc=200.0, origin='soma', replace=False)
        cell.modify_mech_param('tuft', 'mykca', 'gkbar', value=0.007585 * m_mykca)
        cell.modify_mech_param('basal', 'mykca', 'gkbar', origin='soma')
        cell.modify_mech_param('apical', 'mykca', 'gkbar', origin='trunk')

        # ─── kca (SK-type calcium-activated K+) ───────────────────────
        # Similar spatial gradient to mykca
        cell.modify_mech_param('soma', 'kca', 'gbar', value=0.0005 * m_kca)
        cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005 * m_kca, max_loc=50.0, origin='soma')
        cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005 * m_kca,
                               slope=-3.056e-6 * m_kca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
        cell.modify_mech_param('trunk', 'kca', 'gbar', value=4.167e-05 * m_kca,
                               min_loc=200.0, origin='soma', replace=False)
        cell.modify_mech_param('tuft', 'kca', 'gbar', value=4.167e-05 * m_kca)
        cell.modify_mech_param('basal', 'kca', 'gbar', origin='soma')
        cell.modify_mech_param('apical', 'kca', 'gbar', origin='trunk')

        # ─── REINITIALIZE mechanisms to apply all modify_mech_param changes ──
        # This walks each section×mechanism pair and calls NEURON's internal
        # re-initialization to commit the new conductance values.
        all_secs  = ['soma', 'trunk', 'basal', 'apical', 'tuft']
        all_mechs = ['cal', 'calH', 'car', 'cat', 'mykca', 'kca', 'pas']
        for sec in all_secs:
            for mech in all_mechs:
                cell.reinitialize_subset_mechanisms(sec, mech)

        # ── 7. EXPORT BIOPHYSICS YAML ─────────────────────────────────────────
        # This captures the FULL mech_dict (all channels, all compartments) into a YAML
        # that can be loaded directly by CA1_Pyr(mech_filename=...) in the notebook.
        timestamp    = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        biophys_name = f"recovered_biophysics_eval_{eval_num}_{timestamp}.yaml"
        cell.export_mech_dict(mech_filename=biophys_name)
        print(f"\nSUCCESS: {biophys_name} saved.")

        # ── 8. SAVE METADATA YAML ────────────────────────────────────────────
        # Separate provenance file linking the YAML back to the original log + eval
        metadata = {
            'recovered_on':      datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'optimizer_used':    optimizer,
            'original_log':      log_file,
            'evaluation_number': int(eval_num),
            'total_error':       float(row['Total_Error']),
            'channels':          m,
            'features': {
                'areas':   [areas.get(k, None)   for k in area_cols],
                'spikes':  [spikes.get(k, None)  for k in spike_cols],
                'troughs': [troughs.get(k, None) for k in trough_cols]
            },
            'targets': {
                'areas':   [3.337, 3.5297, 4.356, 5.2235, 5.6729],
                'spikes':  [6, 6, 7, 7, 7],
                'troughs': [1.2946, 3.1996, 8.0606, 17.4871, 22.5451]
            }
        }
        meta_name = f"recovered_metadata_eval_{eval_num}_{timestamp}.yaml"
        with open(meta_name, 'w') as f:
            yaml.dump(metadata, f, default_flow_style=False)
        print(f"Metadata saved: {meta_name}")

    except Exception as e:
        print(f"Error during recovery: {e}")
        raise  # re-raise so the full traceback is visible


if __name__ == "__main__":
    recover_from_eval()