import os
import sys
import yaml
import numpy as np
import random
import time
import datetime
import copy
import matplotlib.pyplot as plt
from neuron import h, gui
from specify_cells import CA1_Pyr, QuickSim, data_dir
from plot_utils import *
import click


def assign_exc_and_inh_synapse_stims(cell, num_exc_syns, num_inh_syns, exc_syn_locs_by_sec_type, 
                                     inh_syn_locs_by_sec_type, exc_syn_types, 
                                     inh_syn_types, local_random, excitatory_stochastic):
    """
    Standard assignment logic:
    CA3 -> trunk, apical
    ECIII -> tuft
    """
    stim_exc_syns = {'CA3': [], 'ECIII': []}
    stim_inh_syns = {'CA3': [], 'ECIII': []}

    # 1. Place Excitatory Synapses
    for pathway, num_to_insert in num_exc_syns.items():
        list_of_valid_syn_locs = []
        if pathway == 'ECIII':
            list_of_valid_syn_locs.extend(exc_syn_locs_by_sec_type['tuft'])
        else:
            for sec_type in ['trunk', 'apical']:
                list_of_valid_syn_locs.extend(exc_syn_locs_by_sec_type[sec_type])

        actual_num = min(int(num_to_insert), len(list_of_valid_syn_locs))
        selected_locs = local_random.sample(list_of_valid_syn_locs, actual_num)
        syn_list = cell.insert_synapses_at_syn_locs(selected_locs, exc_syn_types, stochastic=excitatory_stochastic)
        stim_exc_syns[pathway].extend(syn_list)
    
    # 2. Place Inhibitory Synapses
    for pathway, num_to_insert in num_inh_syns.items():
        list_of_valid_syn_locs = []
        if pathway == 'ECIII':
            list_of_valid_syn_locs.extend(inh_syn_locs_by_sec_type['tuft'])
        else:
            for sec_type in ['trunk', 'apical']:
                list_of_valid_syn_locs.extend(inh_syn_locs_by_sec_type[sec_type])

        actual_num = min(int(num_to_insert), len(list_of_valid_syn_locs))
        selected_locs = local_random.sample(list_of_valid_syn_locs, actual_num)
        syn_list = cell.insert_synapses_at_syn_locs(selected_locs, inh_syn_types, stochastic=False)
        stim_inh_syns[pathway].extend(syn_list)

    return stim_exc_syns, stim_inh_syns


@click.command()
@click.option("--config-filename", type=str, default='default_sim_config.yaml')
@click.option("--mech-filename", type=str, default=None)
@click.option("--synapses-seed", type=int, default=0)
@click.option("--trial-seed", type=int, default=0)
@click.option("--plot", is_flag=True)
@click.option("--interactive", is_flag=True)
def main(config_filename, mech_filename, synapses_seed, trial_seed, plot, interactive):
    # 1. LOAD CONFIG
    with open(config_filename, 'r') as f:
        config = yaml.safe_load(f)
    
    sim_params = config['simulation']
    if mech_filename is None:
        mech_filename = sim_params['mech_filename']
    
    duration = sim_params['equilibrate'] + sim_params['sim_duration']
    dt = sim_params['dt']
    v_init = sim_params['v_init']
    tbs_times = sim_params['tbs_times']

    # 2. INITIALIZE CELL
    print(f"--- Initializing Cell with {mech_filename} ---")
    cell = CA1_Pyr(morph_filename=sim_params['morph_filename'],
                   mech_filename=mech_filename, full_spines=False)
    
    # 3. ASSIGN SYNAPSES
    local_random = random.Random(synapses_seed)
    exc_syn_types = ['AMPA_KIN', 'NMDA_KIN5']
    inh_syn_types = ['GABAb']
    
    exc_syn_locs = cell.get_excitatory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
    inh_syn_locs = cell.get_inhibitory_syn_locs(sec_type_list=['trunk', 'apical', 'tuft'])
    
    # Pull synapse counts and settings from config
    num_exc = sim_params['num_exc_syns']
    num_inh = sim_params['num_inh_syns']
    excitatory_stochastic = sim_params.get('excitatory_stochastic', True)
    
    stim_exc_syns, stim_inh_syns = assign_exc_and_inh_synapse_stims(
        cell, num_exc, num_inh, exc_syn_locs, inh_syn_locs, exc_syn_types, 
        inh_syn_types, local_random, excitatory_stochastic=excitatory_stochastic)
    
    cell.init_synaptic_mechanisms()

    # 4. APPLY MUTATION LOGIC
    with open(data_dir + mech_filename, 'r') as f:
        biophys = yaml.safe_load(f)
    gabab_mut = biophys.get('gabab_mutation', None)
    if gabab_mut:
        gmax_mult = gabab_mut['gmax_mult']
        for group in stim_inh_syns.values():
            for syn in group:
                if 'GABAb' in syn._syn:
                    syn._syn['GABAb']['target'].gmax *= gmax_mult
        print(f"I80T mutation applied: GABAb gmax x {gmax_mult}")

    # 5. SETUP RECORDINGS & STIMULATION
    trunk_bif = [n for n in cell.trunk if cell.is_bifurcation(n, 'trunk')]
    if trunk_bif:
        branches = [b for b in trunk_bif[0].children if b.type == 'trunk']
        target_branch = max(branches, key=lambda n: n.sec(0.).diam)
        distal_trunk = next((n for n in cell.trunk if cell.node_in_subtree(target_branch, n) 
                             and 'tuft' in (c.type for c in n.children)))
    else:
        distal_trunk = [n for n in cell.trunk if 'tuft' in (c.type for c in n.children)][0]

    sim = QuickSim(duration, cvode=False, dt=dt, verbose=0)
    sim.append_rec(cell, cell.tree.root, description='soma', loc=0.)
    sim.append_rec(cell, distal_trunk, description='distal_trunk', loc=1.)
    
    tbs_vec = h.Vector(tbs_times)
    for group in list(stim_exc_syns.values()) + list(stim_inh_syns.values()):
        for syn in group:
            syn.source.play(tbs_vec)

    # 6. RUN
    print(f"Starting simulation for {duration} ms (Trial Seed: {trial_seed})...")
    local_random.seed(trial_seed) 
    sim.run(v_init=v_init)
    print("Simulation complete.")

    # 7. PLOT
    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        t = sim.tvec.to_python()
        for rec in sim.rec_list:
            if 'soma' in rec['description']:
                ax1.plot(t, rec['vec'].to_python(), color='black', label='Soma')
            elif 'distal' in rec['description']:
                ax2.plot(t, rec['vec'].to_python(), color='red', label='Distal Trunk')
        
        for ax in [ax1, ax2]:
            ax.set_xlim(400, 1500)
            add_scalebar(ax)
        plt.tight_layout()
        plt.show()

    if interactive:
        h.gui()


if __name__ == "__main__":
    main()