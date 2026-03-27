"""
GIRK LOF (I80T mutation) comparison.
WT:  recovered_biophysics_eval_1_20260323_0042.yaml     → GABAb gmax unmodified
I80T: recovered_biophysics_I80T_mutation.yaml           → GABAb gmax × 0.5

The I80T YAML contains:
  gabab_mutation:
    gmax_mult: 0.5
Applied directly to syn.target.gmax after init_synaptic_mechanisms().

Output: Soma_Dendrite_LOF.svg / .png
"""
import os, sys, random
import yaml
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from neuron import h
from specify_cells import CA1_Pyr, QuickSim

os.chdir(os.path.dirname(os.path.abspath(__file__)))

SYNAPSE_SEED = 0
EQUILIBRATE  = 250.0
SIM_DURATION = 2500.0
DT           = 0.025
V_INIT       = -70.0
DURATION     = EQUILIBRATE + SIM_DURATION
XLIM         = (400, 1500)
MORPH        = 'EB2-late-bifurcation.swc'

TBS_TIMES = [500., 510., 520., 530., 540.,
             650., 660., 670., 680., 690.,
             800., 810., 820., 830., 840.,
             950., 960., 970., 980., 990.,
             1100., 1110., 1120., 1130., 1140.]

N_EXC = {'CA3': 100,  'ECIII': 160}
N_INH = {'CA3': 120,  'ECIII': 120}
EXC_SYN_TYPES = ['AMPA_KIN', 'NMDA_KIN5']
INH_SYN_TYPES = ['GABAb']

CONDITIONS = [
    ('WT',    'recovered_biophysics_eval_1_20260323_0042.yaml'),
    ('I80T',  'recovered_biophysics_I80T_mutation.yaml'),
]

def run_sim(mech_yaml):
    # Load YAML and check for gabab_mutation field
    mech_dict = yaml.safe_load(open(f'data/{mech_yaml}'))
    gabab_mult = 1.0
    if 'gabab_mutation' in mech_dict:
        gabab_mult = mech_dict['gabab_mutation']['gmax_mult']
        print(f"  gabab_mutation: {mech_dict['gabab_mutation']['label']}  gmax × {gabab_mult}")

    cell = CA1_Pyr(morph_filename=MORPH, mech_filename=mech_yaml, full_spines=False)

    # Trunk recording site
    tbif = [t for t in cell.trunk if cell.is_bifurcation(t, 'trunk')]
    if tbif:
        tbranches = [b for b in tbif[0].children if b.type == 'trunk']
        trunk = max(tbranches, key=lambda n: n.sec(0.).diam)
        trunk = next((n for n in cell.trunk
                      if cell.node_in_subtree(trunk, n)
                      and 'tuft' in (c.type for c in n.children)))
    else:
        trunk = [n for n in cell.trunk if 'tuft' in (c.type for c in n.children)][0]

    # Synapses
    exc_locs = cell.get_excitatory_syn_locs(sec_type_list=['trunk','apical','tuft'])
    inh_locs = cell.get_inhibitory_syn_locs(sec_type_list=['trunk','apical','tuft'])
    lr = random.Random(SYNAPSE_SEED); lr.seed(SYNAPSE_SEED)
    exc_syns = {'CA3':[], 'ECIII':[]}
    inh_syns = {'CA3':[], 'ECIII':[]}
    for pw, n in N_EXC.items():
        locs = exc_locs['tuft'] if pw=='ECIII' else exc_locs['trunk']+exc_locs['apical']
        sel  = lr.sample(locs, min(n, len(locs)))
        exc_syns[pw].extend(cell.insert_synapses_at_syn_locs(sel, EXC_SYN_TYPES, stochastic=True))
    for pw, n in N_INH.items():
        locs = inh_locs['tuft'] if pw=='ECIII' else inh_locs['trunk']+inh_locs['apical']
        sel  = lr.sample(locs, min(n, len(locs)))
        inh_syns[pw].extend(cell.insert_synapses_at_syn_locs(sel, INH_SYN_TYPES, stochastic=False))
    cell.init_synaptic_mechanisms()

    # Apply GABAb gmax mutation from YAML — s._syn['GABAb']['target'].gmax
    if gabab_mult != 1.0:
        count = 0
        for pw in inh_syns:
            for s in inh_syns[pw]:
                try:
                    pp = s._syn['GABAb']['target']   # NEURON GABAb point process
                    pp.gmax *= gabab_mult
                    count += 1
                except (KeyError, AttributeError):
                    pass
        print(f"  Applied gmax × {gabab_mult} to {count} GABAb point processes")

    # Sim
    sim = QuickSim(DURATION, cvode=False, dt=DT, verbose=0)
    sim.parameters['equilibrate'] = EQUILIBRATE
    sim.append_rec(cell, cell.tree.root, description='soma',        loc=0.)
    sim.append_rec(cell, tbif[0],        description='proximal_trunk', loc=1.)
    sim.append_rec(cell, trunk,          description='distal_trunk', loc=1.)
    sv = h.Vector(); cell.spike_detector.record(sv)

    tvec = h.Vector(TBS_TIMES)
    for sd in (exc_syns, inh_syns):
        for pw in sd:
            for syn in sd[pw]:
                syn.source.play(tvec)

    sim.run(v_init=V_INIT)

    t     = np.array(sim.tvec.to_python())
    soma_v = distal_v = None
    for rec in sim.rec_list:
        desc = rec['description'].lower()
        v    = np.array(rec['vec'].to_python())
        if 'soma'   in desc:                        soma_v   = v
        if 'distal' in desc and 'trunk' in desc:    distal_v = v
    return t, soma_v, distal_v

# ── Scale bar ──────────────────────────────────────────────────────────────────
def add_scalebar(ax, scale_x=100, scale_y=50, lw=2.0):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    x_pos = xlim[1] - scale_x - 40
    y_pos = ylim[0] + (ylim[1]-ylim[0])*0.06
    ax.plot([x_pos, x_pos+scale_x], [y_pos, y_pos],         'k-', lw=lw)
    ax.plot([x_pos, x_pos],         [y_pos, y_pos+scale_y], 'k-', lw=lw)
    ax.text(x_pos+scale_x/2, y_pos-(ylim[1]-ylim[0])*0.04, f'{scale_x} ms',
            ha='center', va='top', fontsize=10, fontweight='bold')
    ax.text(x_pos-(xlim[1]-xlim[0])*0.01, y_pos+scale_y/2, f'{scale_y} mV',
            ha='right', va='center', fontsize=10, fontweight='bold', rotation=90)

# ── Run ────────────────────────────────────────────────────────────────────────
results = {}
for label, mech in CONDITIONS:
    print(f"\n{'='*50}\n{label}  ({mech})\n{'='*50}")
    results[label] = run_sim(mech)

# ── 2×2 plot ──────────────────────────────────────────────────────────────────
plt.rcParams.update({'font.family':'sans-serif','font.size':11})
fig, axes = plt.subplots(2, 2, figsize=(16, 8), facecolor='white')
fig.patch.set_facecolor('white')

for col_i, (label, mech) in enumerate(CONDITIONS):
    t, sv, dv = results[label]
    n = min(len(t), *(len(v) for v in [sv, dv] if v is not None))
    for row_i, (vec, row_label) in enumerate([(sv,'Soma'), (dv,'Dendrite')]):
        ax = axes[row_i][col_i]
        if vec is not None:
            ax.plot(t[:n], vec[:n], color='black', lw=1.2)
        ax.set_xlim(*XLIM)
        ax.set_ylim(-90, 40)
        ax.set_axis_off()
        if col_i == 0:
            ax.text(-0.06, 0.5, row_label, transform=ax.transAxes,
                    rotation=90, va='center', fontsize=16, fontweight='bold')
        if row_i == 0:
            ax.set_title(label, fontsize=14, fontweight='bold', pad=6)
        if row_i == 1:
            add_scalebar(ax)

plt.tight_layout()
for ext in ('svg', 'png'):
    fig.savefig(f'Soma_Dendrite_LOF.{ext}', dpi=200, bbox_inches='tight')
    print(f"Saved: Soma_Dendrite_LOF.{ext}")
plt.close()
