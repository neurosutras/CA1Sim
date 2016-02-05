from plot_results import *
from specify_cells import *

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '020516 altered km2 rinp - ampa'

equilibrate = 250.  # time to steady-state
duration = 400.
v_init = -67.
NMDA_type = 'NMDA_KIN3'
syn_types = ['AMPA_KIN', NMDA_type]
excitatory_stochastic = False

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
cell.insert_inhibitory_synapses_in_subset()

trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = (node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                                    for child in node.children)).next()
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]
#trunk = trunk_bifurcation[0]

for node in cell.trunk:
    for spine in node.spines:
        syn = Synapse(cell, spine, syn_types, stochastic=excitatory_stochastic)
cell.init_synaptic_mechanisms()

sim = QuickSim(duration)
sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
sim.append_rec(cell, trunk, description='proximal_trunk', loc=1.)
#AMPA_syn = trunk.spines[0].synapses[0]
#GABA_syn = trunk.synapses[0]
GABA_syn = cell.tree.root.synapses[0]
GABA_syn.target('GABA_A_KIN').Erev = 0.
GABA_syn.target('GABA_A_KIN').gmax = 0.000492
#sim.append_rec(cell, trunk.spines[0], object=AMPA_syn.target('AMPA_KIN'), param='_ref_i', description='AMPA_i')
sim.append_rec(cell, cell.tree.root, object=GABA_syn.target('GABA_A_KIN'), param='_ref_i', description='GABA_i')

GABA_syn.source.play(h.Vector([equilibrate]))
sim.run(v_init)
sim.plot()
