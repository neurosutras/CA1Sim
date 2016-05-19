__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import sys
import random
"""
Synaptic mechanisms were implemented with a kinetic scheme that allows temporal summation to be tuned to exhibit
facilitation. This script generates plots to demonstrate the normalized summation for single trunk synapses stimulated
5 x 100 Hz.
"""

if len(sys.argv) > 1:
    syn_type = str(sys.argv[1])
else:
    syn_type = 'NMDA_KIN5'
if len(sys.argv) > 2:
    title = str(sys.argv[2])
else:
    title = None
if len(sys.argv) > 3:
    svg_title = str(sys.argv[3])
else:
    svg_title = None

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'


def plot_facilitation(syn_type, title=None, svg_title=None):
    """

    :param syn_type: str
    :param svg_title: str
    """
    syn_list = []
    syn_list.extend([syn for syn in trunk.synapses if syn_type in syn._syn])
    for spine in trunk.spines:
        syn_list.extend([syn for syn in spine.synapses if syn_type in syn._syn])
    syn = local_random.choice(syn_list)
    syn.source.play(spike_times)
    sim.append_rec(cell, syn.node, param='_ref_g', object=syn.target(syn_type), description='g_' + syn_type,
                       ylabel='Conductance', units='uS')
    sim.append_rec(cell, syn.node, param='_ref_Ro', object=syn.target(syn_type), description='O_' + syn_type,
                   ylabel='Occupancy', units=' ')
    sim.run(v_init)
    syn.source.play(h.Vector([]))

    if svg_title is not None:
        remember_font_size = mpl.rcParams['font.size']
        mpl.rcParams['font.size'] = 20

    dt = 0.01
    t = np.arange(0., duration, dt)
    left = int(equilibrate / dt)
    right = int((equilibrate + ISI) / dt)
    g = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
    O = np.interp(t, sim.tvec, sim.rec_list[1]['vec'])
    unit_gmax = np.max(g[left:right])
    g /= unit_gmax
    O_max = np.max(O[left:right])
    print '%s 1st pulse occupancy: %.4f' % (syn_type, O_max)
    t -= equilibrate
    start = int((equilibrate - 5.)/dt)
    fig, axes = plt.subplots(1)
    axes.plot(t[start:], g[start:], color='k')
    axes.set_xlabel('Time (ms)')
    if title is not None:
        axes.set_title(title, fontsize=mpl.rcParams['font.size'])
    axes.set_ylabel('Normalized conductance')
    axes.set_ylim(-0.1, 1.5)
    axes.set_xlim(-5., 140.)
    axes.set_xticks([0., 25., 50., 75., 100., 125.])
    clean_axes(axes)
    axes.tick_params(direction='out')
    if not svg_title is None:
        fig.set_size_inches(5.27, 4.37)
        fig.savefig(data_dir + svg_title + ' - ' + syn_type + ' facilitation.svg', format='svg', transparent=True)
    plt.show()
    plt.close()
    if svg_title is not None:
        mpl.rcParams['font.size'] = remember_font_size


equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -67.
num_stims = 5
ISI = 10.
spike_times = [equilibrate+i*ISI for i in range(num_stims)]
spike_times = h.Vector(spike_times)

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.zero_na()

sim = QuickSim(duration)

# look for a trunk bifurcation
trunk_bifurcation = [trunk for trunk in cell.trunk if len(trunk.children) > 1 and trunk.children[0].type == 'trunk' and
                     trunk.children[1].type == 'trunk']

# get where the thickest trunk branch gives rise to the tuft
if trunk_bifurcation:  # follow the thicker trunk
    trunk = max(trunk_bifurcation[0].children[:2], key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type for child in
                                                                                             node.children)).next()
else:
    trunk = (node for node in cell.trunk if 'tuft' in (child.type for child in node.children)).next()
tuft = (child for child in trunk.children if child.type == 'tuft').next()
trunk = trunk_bifurcation[0]

spine_list = []
if syn_type in ['AMPA_KIN', 'NMDA_KIN5']:
    for branch in cell.trunk:
        for spine in branch.spines:
            syn = Synapse(cell, spine, [syn_type], stochastic=0)
    spine_list.extend(trunk.spines)
cell.init_synaptic_mechanisms()
cell.insert_inhibitory_synapses_in_subset(['trunk'])

local_random = random.Random()
local_random.seed(0)

plot_facilitation(syn_type, title, svg_title)

# AMPAR unit occupancy: 0.3258
# NMDAR unit occupancy: 0.5211
# GABAR unit occupancy: 0.4478