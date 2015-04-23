__author__ = 'milsteina'
from IPython.parallel import Client
from IPython.display import clear_output
from plot_results import *
import parallel_spine_attenuation_ratio_engine
import sys
"""
Parallel version: Iterates through spines injecting an EPSC-shaped current and recording from the spine and branch to
calculate the amplitude attenuation ratio, spine neck, and branch impedance. Mechanism dictionary is specified in
parallel_spine_attenuation_ratio_engine.py

Assumes a controller is already running in another process with:
ipcluster start -n num_cores
"""
#new_rec_filename = '042015 soma_pas spine_adjusted - EB2 - spine AR'
#new_rec_filename = '042015 soma_pas kdr ka_scale - adjusted - EB2 - spine AR'
new_rec_filename = '042015 pas_ka_scale kdr - EB2 - spine AR'

num_syns = len(parallel_spine_attenuation_ratio_engine.spine_syn_list)
c = Client()
dv = c[:]
dv.clear()
dv.block = True
start_time = time.time()
dv.execute('from parallel_spine_attenuation_ratio_engine import *')
v = c.load_balanced_view()
result = v.map_async(parallel_spine_attenuation_ratio_engine.calculate_single_attenuation_ratio, range(num_syns))
while not result.ready():
    clear_output()
    for stdout in [stdout for stdout in result.stdout if stdout][-len(c):]:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
    sys.stdout.flush()
    time.sleep(60)
for stdout in result.stdout:
    if stdout:
        lines = stdout.split('\n')
        if lines[-2]:
            print lines[-2]
print 'Parallel execution took:', time.time()-start_time, 's'
rec_file_list = dv['rec_filename']
combine_output_files(rec_file_list, new_rec_filename)
plot_AR(new_rec_filename)