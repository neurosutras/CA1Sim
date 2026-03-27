"""
Export the YAML biophysics for Eval #1 of optimization_log_20260323_0042.csv.

Source:     optimization_log_20260323_0042.csv  (overnight Basin-Hopping, started 2026-03-23 00:42)
Eval:       1  (also evals 2-7 are identical — all share error = 117277.80)
Base mech:  recovered_biophysics_Best_Model_so_Far_20260322_2311.yaml

Multipliers applied on top of Best_Model:
  cal   = 0.599473
  calH  = 1.033540
  car   = 1.119427
  cat   = 1.136134
  mykca = 1.300000
  kca   = 1.300000

Features (from log):
  Areas:    [1.528, 1.734, 2.004, 2.132, 4.756]  — increasing
  Spikes:   [5, 5, 5, 5, 6]                       — all 5 cycles active
  Troughs:  [-1.057, -1.846, 3.192, 3.936, 3.928] — positive in C3-C5
  Error:    117277.80
"""
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from specify_cells import CA1_Pyr

# ── Load base biophysics ──────────────────────────────────────────────────────
cell = CA1_Pyr(
    morph_filename='EB2-late-bifurcation.swc',
    mech_filename='recovered_biophysics_Best_Model_so_Far_20260322_2311.yaml',
    full_spines=False
)

# ── Apply Eval #1 multipliers ─────────────────────────────────────────────────
m_cal   = 0.599473
m_calH  = 1.033540
m_car   = 1.119427
m_cat   = 1.136134
m_mykca = 1.300000
m_kca   = 1.300000

# --- cal ---
cell.modify_mech_param('soma',  'cal', 'gcalbar', value=0.007 * m_cal)
for sec in ['trunk', 'basal']:
    cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='soma')
for sec in ['apical', 'tuft']:
    cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='trunk')

# --- calH ---
cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.000031635 * m_calH, max_loc=50.0, origin='soma')
cell.modify_mech_param('trunk', 'calH', 'gcalbar', value=0.001455 * m_calH, min_loc=50.0, origin='soma', replace=False)
for sec in ['apical', 'tuft']:
    cell.modify_mech_param(sec, 'calH', 'gcalbar', origin='trunk')

# --- car ---
cell.modify_mech_param('soma',  'car', 'gcabar', value=0.003   * m_car)
cell.modify_mech_param('trunk', 'car', 'gcabar', value=0.00003 * m_car)
cell.modify_mech_param('basal', 'car', 'gcabar', origin='soma')
for sec in ['apical', 'tuft']:
    cell.modify_mech_param(sec, 'car', 'gcabar', origin='trunk')

# --- cat ---
cell.modify_mech_param('soma',  'cat', 'gcatbar', value=0.00005 * m_cat)
cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0 * m_cat, max_loc=100.0, origin='soma')
cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0, slope=1.143e-06 * m_cat, min_loc=100.0, max_loc=350.0, origin='soma', replace=False)
cell.modify_mech_param('trunk', 'cat', 'gcatbar', value=0.0004 * m_cat, min_loc=350.0, origin='soma', replace=False)
cell.modify_mech_param('basal', 'cat', 'gcatbar', origin='soma')
for sec in ['apical', 'tuft']:
    cell.modify_mech_param(sec, 'cat', 'gcatbar', origin='trunk')

# --- mykca ---
cell.modify_mech_param('soma',  'mykca', 'gkbar', value=0.09075 * m_mykca)
cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, max_loc=50.0, origin='soma')
cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * m_mykca, slope=-0.0005543 * m_mykca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.007585 * m_mykca, min_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('tuft',  'mykca', 'gkbar', value=0.007585 * m_mykca)
cell.modify_mech_param('basal', 'mykca', 'gkbar', origin='soma')
cell.modify_mech_param('apical','mykca', 'gkbar', origin='trunk')

# --- kca ---
cell.modify_mech_param('soma',  'kca', 'gbar', value=0.0005    * m_kca)
cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005    * m_kca, max_loc=50.0, origin='soma')
cell.modify_mech_param('trunk', 'kca', 'gbar', value=0.0005    * m_kca, slope=-3.056e-6 * m_kca, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('trunk', 'kca', 'gbar', value=4.167e-05 * m_kca, min_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('tuft',  'kca', 'gbar', value=4.167e-05 * m_kca)
cell.modify_mech_param('basal', 'kca', 'gbar', origin='soma')
cell.modify_mech_param('apical','kca', 'gbar', origin='trunk')

# --- reinitialize ---
for sec in ['soma', 'trunk', 'basal', 'apical', 'tuft']:
    for mech in ['cal', 'calH', 'car', 'cat', 'mykca', 'kca', 'pas']:
        cell.reinitialize_subset_mechanisms(sec, mech)

# ── Export ────────────────────────────────────────────────────────────────────
out_filename = 'recovered_biophysics_eval_1_20260323_0042.yaml'
cell.export_mech_dict(out_filename)
print(f"\n✅ Saved: data/{out_filename}")
print()
print("=" * 60)
print("TO USE IN NOTEBOOK (Cell 2):")
print("=" * 60)
print(f"mech = '{out_filename}'")
print("cell = CA1_Pyr(morph_filename=morph, mech_filename=mech, full_spines=False)")
print()
print("No Cell 3 multiplier block needed — conductances are baked in.")
print()
print("Source: optimization_log_20260323_0042.csv")
print("Run:    Overnight Basin-Hopping started 2026-03-23 00:42")
print("Eval:   1  (error = 117277.80)")
print("  cal=0.599473  calH=1.03354  car=1.119427")
print("  cat=1.136134  mykca=1.3     kca=1.3")
