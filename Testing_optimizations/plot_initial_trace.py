"""
Run one simulation with Eval #1 biophysics and plot the soma Vm trace.
Saves initial_eval_soma_trace.svg and .png
"""
import sys, os
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from neuron import h
from specify_cells import CA1_Pyr, QuickSim

# ── Config (must match optimize_soma_plateaus.py exactly) ──────────────────────
MORPH        = 'EB2-late-bifurcation.swc'
MECH         = 'recovered_biophysics_eval_1_20260323_0042.yaml'   # eval_1 overnight 0042 run
MECH_LABEL   = 'eval1_0042_GIRK_x0p9'

# GIRK_Yim insertion — matches notebook Cell 2
# Default gbar (girk_yim.mod) = 1.44e-05 S/cm2
# CHANGE GIRK_MULTIPLIER here and in notebook Cell 2 together
GIRK_MULTIPLIER   = 0.9                       # ← CHANGE THIS (0.9 / 0.8 / 0.7 / 0.5)
GIRK_GBAR         = 1.44e-05 * GIRK_MULTIPLIER

# Channel multipliers — all 1.0 (conductances baked into YAML)
M_CAL   = 1.0; M_CALH = 1.0; M_CAR  = 1.0
M_CAT   = 1.0; M_MYKCA = 1.0; M_KCA = 1.0

SYNAPSE_SEED = 0
EQUILIBRATE  = 250.0
SIM_DURATION = 2500.0
DT           = 0.025
V_INIT       = -70.0
DURATION     = EQUILIBRATE + SIM_DURATION

TBS_TIMES    = [500., 510., 520., 530., 540.,
                650., 660., 670., 680., 690.,
                800., 810., 820., 830., 840.,
                950., 960., 970., 980., 990.,
                1100., 1110., 1120., 1130., 1140.]

N_EXC = {'CA3': 100, 'ECIII': 160}
N_INH = {'CA3': 120, 'ECIII': 120}

EXC_SYN_TYPES = ['AMPA_KIN', 'NMDA_KIN5']
INH_SYN_TYPES = ['GABAb']

# ── Build cell ─────────────────────────────────────────────────────────────────
print("--- Initializing cell ---")
# specify_cells looks for mech in data_dir; we pass the path without data_dir prefix
# Testing_optimizations/ is not data_dir, so we change CWD trick:
os.chdir(os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd())

cell = CA1_Pyr(morph_filename=MORPH, mech_filename=MECH, full_spines=False)

# ── Apply channel multipliers ──────────────────────────────────────────────────
cell.modify_mech_param('soma',  'cal',   'gcalbar', value=0.007 * M_CAL)
for sec in ['trunk', 'basal']: cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='soma')
for sec in ['apical', 'tuft']: cell.modify_mech_param(sec, 'cal', 'gcalbar', origin='trunk')

cell.modify_mech_param('trunk', 'calH',  'gcalbar', value=0.000031635 * M_CALH, max_loc=50.0, origin='soma')
cell.modify_mech_param('trunk', 'calH',  'gcalbar', value=0.001455 * M_CALH, min_loc=50.0, origin='soma', replace=False)
for sec in ['apical', 'tuft']: cell.modify_mech_param(sec, 'calH', 'gcalbar', origin='trunk')

cell.modify_mech_param('soma',  'car',   'gcabar', value=0.003   * M_CAR)
cell.modify_mech_param('trunk', 'car',   'gcabar', value=0.00003 * M_CAR)
cell.modify_mech_param('basal', 'car',   'gcabar', origin='soma')
for sec in ['apical', 'tuft']: cell.modify_mech_param(sec, 'car', 'gcabar', origin='trunk')

cell.modify_mech_param('soma',  'cat',   'gcatbar', value=0.00005 * M_CAT)
cell.modify_mech_param('trunk', 'cat',   'gcatbar', value=0.0 * M_CAT, max_loc=100.0, origin='soma')
cell.modify_mech_param('trunk', 'cat',   'gcatbar', value=0.0, slope=1.143e-06 * M_CAT, min_loc=100.0, max_loc=350.0, origin='soma', replace=False)
cell.modify_mech_param('trunk', 'cat',   'gcatbar', value=0.0004 * M_CAT, min_loc=350.0, origin='soma', replace=False)
cell.modify_mech_param('basal', 'cat',   'gcatbar', origin='soma')
for sec in ['apical', 'tuft']: cell.modify_mech_param(sec, 'cat', 'gcatbar', origin='trunk')

cell.modify_mech_param('soma',  'mykca', 'gkbar', value=0.09075 * M_MYKCA)
cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * M_MYKCA, max_loc=50.0, origin='soma')
cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.09075 * M_MYKCA, slope=-0.0005543 * M_MYKCA, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('trunk', 'mykca', 'gkbar', value=0.007585 * M_MYKCA, min_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('tuft',  'mykca', 'gkbar', value=0.007585 * M_MYKCA)
cell.modify_mech_param('basal', 'mykca', 'gkbar', origin='soma')
cell.modify_mech_param('apical','mykca', 'gkbar', origin='trunk')

cell.modify_mech_param('soma',  'kca',   'gbar', value=0.0005    * M_KCA)
cell.modify_mech_param('trunk', 'kca',   'gbar', value=0.0005    * M_KCA, max_loc=50.0, origin='soma')
cell.modify_mech_param('trunk', 'kca',   'gbar', value=0.0005    * M_KCA, slope=-3.056e-6 * M_KCA, min_loc=50.0, max_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('trunk', 'kca',   'gbar', value=4.167e-05 * M_KCA, min_loc=200.0, origin='soma', replace=False)
cell.modify_mech_param('tuft',  'kca',   'gbar', value=4.167e-05 * M_KCA)
cell.modify_mech_param('basal', 'kca',   'gbar', origin='soma')
cell.modify_mech_param('apical','kca',   'gbar', origin='trunk')

for sec in ['soma','trunk','basal','apical','tuft']:
    for mech in ['cal','calH','car','cat','mykca','kca','pas']:
        cell.reinitialize_subset_mechanisms(sec, mech)

# ── Insert GIRK_Yim × 0.5 ─────────────────────────────────────────────────
print(f'Inserting GIRK_Yim into all compartments: gbar = {GIRK_GBAR:.3e} S/cm2 ({GIRK_MULTIPLIER}x default)')
for node in cell.tree:
    if node.type in ['soma', 'trunk', 'basal', 'apical', 'tuft']:
        node.sec.insert('GIRK_Yim')
        node.sec.gbar_GIRK_Yim = GIRK_GBAR

# ── Synapse locations ──────────────────────────────────────────────────────────
exc_locs = cell.get_excitatory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
inh_locs = cell.get_inhibitory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
local_random = random.Random(SYNAPSE_SEED)

stim_exc_syns = {'CA3': [], 'ECIII': []}
stim_inh_syns = {'CA3': [], 'ECIII': []}

def assign_synapses():
    global stim_exc_syns, stim_inh_syns
    stim_exc_syns = {'CA3': [], 'ECIII': []}
    stim_inh_syns = {'CA3': [], 'ECIII': []}

    for pathway, n in N_EXC.items():
        locs = exc_locs['tuft'] if pathway == 'ECIII' else exc_locs['trunk'] + exc_locs['apical']
        sel  = local_random.sample(locs, min(n, len(locs)))
        stim_exc_syns[pathway].extend(cell.insert_synapses_at_syn_locs(sel, EXC_SYN_TYPES, stochastic=True))

    for pathway, n in N_INH.items():
        locs = inh_locs['tuft'] if pathway == 'ECIII' else inh_locs['trunk'] + inh_locs['apical']
        sel  = local_random.sample(locs, min(n, len(locs)))
        stim_inh_syns[pathway].extend(cell.insert_synapses_at_syn_locs(sel, INH_SYN_TYPES, stochastic=False))

local_random.seed(SYNAPSE_SEED)
assign_synapses()
cell.init_synaptic_mechanisms()

# ── Sim setup ──────────────────────────────────────────────────────────────────
print("--- Setting up simulation ---")
sim = QuickSim(DURATION, cvode=False, dt=DT, verbose=0)
sim.append_rec(cell, cell.tree.root, description='soma_v')
spike_vec = h.Vector()
cell.spike_detector.record(spike_vec)

tbs_vec = h.Vector(TBS_TIMES)
for stim_dict in (stim_exc_syns, stim_inh_syns):
    for pathway in stim_dict:
        for syn in stim_dict[pathway]:
            syn.source.play(tbs_vec)

# ── Run ────────────────────────────────────────────────────────────────────────
print("--- Running simulation ---")
sim.run(v_init=V_INIT)
print("--- Done ---")

soma_v   = np.array(sim.get_rec('soma_v')['vec'])
time_ms  = np.arange(0, DURATION, DT)
n        = min(len(soma_v), len(time_ms))
soma_v   = soma_v[:n]
time_ms  = time_ms[:n]

# ── Publication-style plot with scale bar ─────────────────────────────────────
def add_scale_bar(ax, x_size=300, y_size=50, x_label='300 ms', y_label='50 mV',
                  x_pos=None, y_pos=None, lw=2.5):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    if x_pos is None:
        x_pos = xlim[1] - 0.05 * (xlim[1] - xlim[0]) - x_size
    if y_pos is None:
        y_pos = ylim[0] + 0.05 * (ylim[1] - ylim[0])
    ax.plot([x_pos, x_pos + x_size], [y_pos, y_pos], color='black', lw=lw, clip_on=False)
    ax.plot([x_pos + x_size, x_pos + x_size], [y_pos, y_pos + y_size], color='black', lw=lw, clip_on=False)
    ax.text(x_pos + x_size/2, y_pos - 0.04*(ylim[1]-ylim[0]), x_label,
            ha='center', va='top', fontsize=12, fontweight='bold')
    ax.text(x_pos + x_size + 0.01*(xlim[1]-xlim[0]), y_pos + y_size/2, y_label,
            ha='left', va='center', fontsize=12, fontweight='bold', rotation=90)

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 11})

fig, ax = plt.subplots(figsize=(10, 5), facecolor='white')
ax.plot(time_ms, soma_v, color='black', lw=1.0)
ax.set_axis_off()
ax.set_ylim(-90, 40)
ax.set_xlim(300, 1350)

add_scale_bar(ax)

ax.text(-0.04, 0.5, 'Soma', transform=ax.transAxes,
        rotation=90, verticalalignment='center', fontsize=22, fontweight='bold')

plt.tight_layout()

for ext in ('svg', 'png'):
    out = f'{MECH_LABEL}_soma_trace.{ext}'
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f'Saved: {out}')

plt.close()
