__author__ = 'milsteina'
from specify_cells import *
from plot_results import *

#morph_filename = 'EB1-early-bifurcation.swc'
morph_filename = 'EB2-late-bifurcation.swc'
mech_filename = '041515 soma_pas spines - EB2.pkl'

cell = CA1_Pyr(morph_filename, mech_filename, full_spines=True)
for sec_type in ['spine_head', 'spine_neck']:
    cell.modify_mech_param(sec_type, 'pas', 'g', origin='parent')
cell.modify_mech_param('soma', 'kdr', 'gkdrbar', 0.005)
cell.modify_mech_param('soma', 'kap', 'gkabar', 0.022)
cell.modify_mech_param('soma', 'ions', 'ek', -90.)
cell.modify_mech_param('basal', 'kap', 'gkabar', origin='soma')
for sec_type in ['axon', 'ais']:
    cell.modify_mech_param(sec_type, 'kap', 'gkabar', cell.mech_dict['soma']['kap']['gkabar']['value']*0.2)
cell.modify_mech_param('trunk', 'kap', 'gkabar', origin='soma', slope=0.00022, max_loc=75.)
cell.modify_mech_param('trunk', 'kap', 'gkabar', value=0.0385, origin='soma', slope=-0.00077, min_loc=75., max_loc=125.,
                       replace=False)
cell.modify_mech_param('trunk', 'kap', 'gkabar', value=0., origin='soma', min_loc=125., replace=False)
cell.modify_mech_param('trunk', 'kad', 'gkabar', value=0., origin='soma', max_loc=75.)
cell.modify_mech_param('trunk', 'kad', 'gkabar', value=0., origin='soma', slope=0.00099, min_loc=75., max_loc=125.,
                       replace=False)
cell.modify_mech_param('trunk', 'kad', 'gkabar', value=0.0495, origin='soma', slope=0.00022, min_loc=125.,
                       replace=False)
cell.modify_mech_param('apical', 'kap', 'gkabar', origin='trunk')
cell.modify_mech_param('apical', 'kad', 'gkabar', origin='trunk')
cell.modify_mech_param('tuft', 'kad', 'gkabar', origin='trunk')
for sec_type in ['basal', 'axon', 'ais', 'trunk', 'apical', 'tuft']:
    cell.modify_mech_param(sec_type, 'kdr', 'gkdrbar', origin='soma')
    cell.modify_mech_param(sec_type, 'ions', 'ek', origin='soma')
#cell.reinit_mechanisms()
plot_mech_param_distribution(cell, 'kad', 'gkabar')