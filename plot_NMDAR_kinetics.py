__author__ = 'Aaron D. Milstein'
from specify_cells import *
from plot_results import *
import random
import sys
"""
This simulation checks the voltage dependence of NMDAR kinetics while clamping the trunk voltage.
"""
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '043016 Type A - km2_NMDA_KIN5_Pr'

if len(sys.argv) > 1:
    svg_title = str(sys.argv[1])
else:
    svg_title = None

equilibrate = 250.  # time to steady-state
duration = 450.
v_init = -67.
NMDA_type = 'NMDA_KIN5'
syn_types = ['AMPA_KIN', NMDA_type]

spike_times = h.Vector([equilibrate])

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

syn_list = []
for spine in trunk.spines:
    syn = Synapse(cell, spine, syn_types, stochastic=0)
    syn_list.append(syn)
cell.init_synaptic_mechanisms()

local_random = random.Random()
local_random.seed(0)

syn = local_random.choice(syn_list)
syn.source.play(spike_times)
sim.append_rec(cell, syn.node, object=syn.target(NMDA_type), param='_ref_g', description='g')

clamp = h.SEClamp(trunk.sec(1.))
clamp.dur1=duration
clamp.rs = 0.015
vc_range = np.arange(-70., 0., 10.)

results = {}
dt = 0.02
t = np.arange(0., duration, dt)

for vc in vc_range:
    print 'Holding Voltage:', vc
    clamp.amp1 = vc
    sim.run(v_init)
    result = np.interp(t, sim.tvec, sim.rec_list[0]['vec'])
    result /= np.max(result)
    results[vc] = result
syn.source.play(h.Vector([]))

if svg_title is not None:
    remember_font_size = mpl.rcParams['font.size']
    mpl.rcParams['font.size'] = 20

fig, axes = plt.subplots(1)
colors = ['k', 'r', 'c', 'y', 'm', 'g', 'b'][::-1]
left = int(equilibrate / dt)
t -= equilibrate
start = int((equilibrate - 5.) / dt)
for i, vc in enumerate(vc_range[::-1]):
    axes.plot(t[start:], results[vc][start:], color=colors[i], label=vc)
axes.set_xlabel('Time (ms)')
axes.set_ylabel('Normalized conductance')
axes.set_ylim(-0.1, 1.1)
axes.set_xlim(-5., 200.)
axes.set_xticks([0., 50., 100., 150., 200.])
clean_axes(axes)
axes.tick_params(direction='out')
plt.legend(loc='best', frameon=False, framealpha=0.5, fontsize=mpl.rcParams['font.size'])
if not svg_title is None:
    fig.set_size_inches(5.27, 4.37)
    fig.savefig(data_dir + svg_title + ' - NMDAR kinetics.svg', format='svg', transparent=True)
plt.show()
plt.close()

gamma = syn.target(NMDA_type).gamma
Kd = syn.target(NMDA_type).Kd
mg = syn.target(NMDA_type).mg
v = np.arange(-80., 40., 1.)
g = 1. / (1. + np.exp(gamma * (-v)) * (mg / Kd))

fig, axes = plt.subplots(1)
axes.plot(v, g, color='k')
axes.set_xlabel('Voltage (mV)')
axes.set_ylabel('Normalized conductance')
axes.set_ylim(0., 1.)
axes.set_xlim(-80., 40.)
axes.set_xticks([-80., -60., -40., -20., 0., 20., 40.])
clean_axes(axes)
axes.tick_params(direction='out')
if not svg_title is None:
    fig.set_size_inches(5.27, 4.37)
    fig.savefig(data_dir + svg_title + ' - NMDAR g-V.svg', format='svg', transparent=True)
plt.show()
plt.close()

if svg_title is not None:
    mpl.rcParams['font.size'] = remember_font_size