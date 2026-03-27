"""
plot_error_log.py — Marder-group (parameter degeneracy) analysis for optimization logs.

Purpose:
  After running optimize_soma_plateaus.py, this script reads the CSV log and identifies
  parameter sets whose total error is within a threshold of the global minimum.
  These "degenerate" solutions (Marder groups) show that multiple conductance combinations
  can produce equally good theta-burst plateau behaviour — a hallmark of biological robustness.

Outputs:
  1. A YAML summary (<log_file>_marder_summary.yaml) listing all degenerate solutions
  2. A terminal table comparing channels, areas, spikes, and troughs
  3. A PNG plot (<log_file>_degeneracy_plot.png) showing error vs. evaluation number

Usage:
  python plot_error_log.py --log-file optimization_log_20260324_2200.csv --threshold 0.5 --optimizer BasinHopping
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import click
import os
import yaml

# ========================================================================================
# CLI INTERFACE
# ========================================================================================
@click.command()
@click.option('--log-file', default='optimization_log_20260324_2200.csv',
              help='Path to the CSV log file produced by optimize_soma_plateaus.py.')
@click.option('--threshold', default=0.5,
              help='Degeneracy threshold: fraction above global minimum. '
                   '0.5 → all solutions with error ≤ 1.5× best error.')
@click.option('--optimizer', default='BasinHopping',
              help='Name of the optimizer (metadata only — saved in the YAML summary).')
@click.option('--num-evals', default=84,
              help='Total number of evaluations (metadata only).')
@click.option('--num-iters', default=5,
              help='Number of major iterations (metadata only).')

#python plot_error_log.py --log-file optimization_log_20260324_2200.csv --threshold 0.5 --optimizer BasinHopping
def analyze_marder_groups(log_file, threshold, optimizer, num_evals, num_iters):
    """
    Main analysis pipeline:
      1. Load CSV  →  2. Find global minimum  →  3. Identify degenerate solutions
      →  4. Build YAML summary  →  5. Print terminal report  →  6. Save plot
    """

    # ── 1. SAFETY CHECK ───────────────────────────────────────────────────────
    if not os.path.exists(log_file):
        print(f"Error: File '{log_file}' not found.")
        return

    # ── 2. LOAD DATA ──────────────────────────────────────────────────────────
    df = pd.read_csv(log_file)

    # ── 3. FIND GLOBAL MINIMUM ────────────────────────────────────────────────
    best_error = df['Total_Error'].min()           # lowest error across all evals
    best_idx   = df['Total_Error'].idxmin()        # row index of that eval
    best_eval  = df.loc[best_idx]                  # full row as a Series

    print("-" * 70)
    print(f"MARDER GROUP ANALYSIS: {log_file}")
    print(f"Global Minimum Found at Eval {int(best_eval['Eval'])}: {best_error:.6f}")
    print("-" * 70)

    # ── 4. IDENTIFY MARDER GROUP (degenerate solutions) ───────────────────────
    # Any evaluation with error ≤ (1 + threshold) × best_error is "degenerate"
    threshold_val = best_error * (1 + threshold)
    marder_df     = df[df['Total_Error'] <= threshold_val].sort_values(by='Total_Error')

    print(f"Found {len(marder_df)} Degenerate Solutions (Error <= {threshold_val:.4f})")

    # ── 5. BUILD DETAILED SUMMARY (channel params + features per eval) ────────
    # Column groups — must match the CSV headers from optimize_soma_plateaus.py
    channel_cols = ['cal', 'calH', 'car', 'cat', 'mykca', 'kca']
    synapse_cols = ['n_exc_CA3', 'n_exc_ECIII', 'n_inh_CA3', 'n_inh_ECIII']
    area_cols    = ['Area_C1', 'Area_C2', 'Area_C3', 'Area_C4', 'Area_C5']
    spike_cols   = ['Spks_C1', 'Spks_C2', 'Spks_C3', 'Spks_C4', 'Spks_C5']
    trough_cols  = ['Trough_C1', 'Trough_C2', 'Trough_C3', 'Trough_C4', 'Trough_C5']

    marder_group_details = []

    for _, row in marder_df.iterrows():
        # Safely extract each column group — skip columns that don't exist or are NaN
        found_channels = {k: float(row[k])          for k in channel_cols if k in df.columns and not pd.isna(row[k])}
        found_synapses = {k: int(round(row[k]))     for k in synapse_cols if k in df.columns and not pd.isna(row[k])}
        found_areas    = {k: round(float(row[k]), 4) for k in area_cols   if k in df.columns and not pd.isna(row[k])}
        found_spikes   = {k: int(row[k])             for k in spike_cols  if k in df.columns and not pd.isna(row[k])}
        found_troughs  = {k: round(float(row[k]), 4) for k in trough_cols if k in df.columns and not pd.isna(row[k])}

        entry = {
            'Eval':        int(row['Eval']),
            'Total_Error': float(row['Total_Error']),
            'params': {
                'channels': found_channels,
                'synapses': found_synapses
            },
            'features': {
                'areas':   found_areas,
                'spikes':  found_spikes,
                'troughs': found_troughs
            }
        }
        marder_group_details.append(entry)

    # ── 6. SAVE YAML SUMMARY ─────────────────────────────────────────────────
    summary_data = {
        'metadata': {
            'log_file':             log_file,
            'optimizer':            optimizer,
            'num_evaluations':      num_evals,
            'num_iterations':       num_iters,
            'threshold_percentage': f"{threshold * 100}%",
            'best_error_baseline':  float(best_error)
        },
        'marder_groups': marder_group_details
    }

    summary_name = log_file.replace('.csv', '_marder_summary.yaml')
    with open(summary_name, 'w') as f:
        yaml.dump(summary_data, f, default_flow_style=False)
    print(f"Enhanced Marder Group Summary saved as: {summary_name}")

    # ── 7. TERMINAL REPORT ────────────────────────────────────────────────────
    sep = "-" * 120
    print(sep)
    print(f"{'Eval':<8} | {'Error':<14} | {'cal':>6} {'calH':>6} {'car':>6} {'cat':>6} {'mykca':>6} {'kca':>6} | {'Areas (C1-5)':<38} | {'Spikes':<18} | Troughs (C1-5)")
    print(sep)

    for entry in marder_group_details:
        c  = entry['params']['channels']
        a  = entry['features']['areas']
        s  = entry['features']['spikes']
        tr = entry['features']['troughs']

        chan_str   = (f"{c.get('cal',0):>6.3f} {c.get('calH',0):>6.3f} {c.get('car',0):>6.3f} "
                      f"{c.get('cat',0):>6.3f} {c.get('mykca',0):>6.3f} {c.get('kca',0):>6.3f}")
        area_str   = str([a.get(k, 0) for k in area_cols])
        spike_str  = str([s.get(k, 0) for k in spike_cols])
        trough_str = str([tr.get(k, 0) for k in trough_cols])

        print(f"{entry['Eval']:<8} | {entry['Total_Error']:<14.2f} | {chan_str} | {area_str:<38} | {spike_str:<18} | {trough_str}")
    print(sep)

    # Reference targets for quick comparison
    print(f"\nTARGETS:")
    print(f"  Areas:   [3.337, 3.5297, 4.356, 5.2235, 5.6729] mV*s")
    print(f"  Spikes:  [6, 6, 7, 7, 7]")
    print(f"  Troughs: [1.2946, 3.1996, 8.0606, 17.4871, 22.5451] mV")

    # ── 8. PLOT ───────────────────────────────────────────────────────────────
    plt.figure(figsize=(12, 6))

    # Full search history (light grey line)
    plt.plot(df['Eval'].values, df['Total_Error'].values,
             color='#bdc3c7', alpha=0.5, label='Search History')

    # Marder group members (orange dots)
    plt.scatter(marder_df['Eval'], marder_df['Total_Error'],
                color='#e67e22', s=50, edgecolors='black',
                label=f'Marder Group (Top {len(marder_df)})')

    # Global minimum (gold star)
    plt.scatter(best_eval['Eval'], best_eval['Total_Error'],
                color='#f1c40f', s=120, edgecolors='black', zorder=5,
                label='Global Minimum')

    plt.yscale('log')
    plt.title(f'Degeneracy Analysis (Marder Groups): {log_file}', fontsize=14)
    plt.xlabel('Evaluation Number')
    plt.ylabel('Total Error (Log Scale)')
    plt.grid(True, which="both", ls="-", alpha=0.1)
    plt.legend()
    plt.tight_layout()

    plot_name = log_file.replace('.csv', '_degeneracy_plot.png')
    plt.savefig(plot_name, dpi=300)
    print(f"\nPlot saved as: {plot_name}")
    plt.show()


if __name__ == "__main__":
    analyze_marder_groups()