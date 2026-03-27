from neuron import h, gui
from specify_cells import CA1_Pyr, QuickSim
from plot_utils import *

import click


@click.command()
@click.option("--mech-filename", type=str, default='20260319_default_biophysics_with_Calcium_Channels.yaml') #Calcium, KCa, and GABAb from Poirazi et al. 2003
@click.option("--synapses-seed", type=int, default=0)
@click.option("--trial-seed", type=int, default=0)
@click.option("--label", type=str, default=None)
@click.option("--sim_duration", type=float, default=None)
@click.option("--spines", type=bool, default=False)
@click.option("--export", is_flag=True)
@click.option("--plot", is_flag=True)
@click.option("--interactive", is_flag=True)
@click.option("--debug", is_flag=True)

def main(mech_filename, synapses_seed, trial_seed, label, sim_duration, spines, export, plot, interactive, debug):

    """
    :param mech_filename:
        .yaml file must be located in the data subdirectory
    :param synapses_seed:
        a unique random seed can be used to shuffle the number and locations of synapses, and the place field locs of
        the presynaptic CA3 inputs (like simulating a different cell with the same morphology)
    :param trial_seed:
        a unique random seed shuffles the input spike times and synaptic release probabilities to allows simulation of
        multiple independent trials
    :param label: append a label when exporting data to .hdf5
    :param sim_duration: float (ms) - sim duration can be truncated during testing
    :param spines: bool, whether to include explicit spine neck and head compartments for every excitatory synapse
    :param export: bool, whether to export to .hdf5
    :param plot: bool, whether to plot
    :param interactive: bool, whether to enable live object inspection after simulation
    :param debug: bool, does not run simulation in debug mode
    """

    morph_filename = 'EB2-late-bifurcation.swc'

    if label is None:
        rec_filename = 'output'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'-pid'+str(os.getpid())+'-seed'+\
                       str(synapses_seed)+'-trial'+str(trial_seed) + '_theta_burst_stim.hdf5'
    else:
        rec_filename = 'output' + datetime.datetime.today().strftime('%m%d%Y%H%M') + '-pid' + str(os.getpid()) + \
                       '-' + label + '-seed' + str(synapses_seed) + '-trial' + str(trial_seed) + '_theta_burst_stim.hdf5'

    def run_trial(simiter):
        """

        :param simiter: int
        """
        



    

    