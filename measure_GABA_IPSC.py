from plot_results import *
from specify_cells import *

morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '20220808_default_biophysics.yaml'

equilibrate = 250.  # time to steady-state
duration = 400.
v_init = -67.
NMDA_type = 'NMDA_KIN5'
syn_types = ['AMPA_KIN', NMDA_type]
GABA_syn_type = 'GABA_A_KIN'
excitatory_stochastic = False

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
inh_syn_locs_by_sec_type = cell.get_inhibitory_syn_locs(sec_type_list=['soma'])
inh_syn_locs = [inh_syn_locs_by_sec_type['soma'][0]]
syn_list = cell.insert_synapses_at_syn_locs(inh_syn_locs, [GABA_syn_type])

trunk_bifurcation = [trunk for trunk in cell.trunk if cell.is_bifurcation(trunk, 'trunk')]
if trunk_bifurcation:
    trunk_branches = [branch for branch in trunk_bifurcation[0].children if branch.type == 'trunk']
    # get where the thickest trunk branch gives rise to the tuft
    trunk = max(trunk_branches, key=lambda node: node.sec(0.).diam)
    trunk = next(node for node in cell.trunk if cell.node_in_subtree(trunk, node) and 'tuft' in (child.type
                                                                                    for child in node.children))
else:
    trunk_bifurcation = [node for node in cell.trunk if 'tuft' in (child.type for child in node.children)]
    trunk = trunk_bifurcation[0]

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
GABA_syn.target('GABA_A_KIN').gmax = 0.000492 * 1.2
sim.append_rec(cell, cell.tree.root, object=GABA_syn.target('GABA_A_KIN'), param='_ref_i', description='GABA_i')

GABA_syn.source.play(h.Vector([equilibrate]))
sim.run(v_init)
sim.plot()
