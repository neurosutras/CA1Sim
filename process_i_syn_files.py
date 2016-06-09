__author__ = 'milsteina'
from function_lib import *
import sys

if len(sys.argv) > 1:
    syn_type = str(sys.argv[1])
if len(sys.argv) > 2:
    seed = int(sys.argv[2])
if len(sys.argv) > 3:
    condition = int(sys.argv[3])

file_list, rec_filenames = {}, {}
for this_syn_type in ['i_AMPA', 'i_NMDA', 'i_GABA']:
    for parameter in file_list, rec_filenames:
        parameter[this_syn_type] = {}
        for this_seed in range(5):
            parameter[this_syn_type][this_seed] = {}

file_list['i_AMPA'][0]['modinh0'] = []
file_list['i_AMPA'][0]['modinh3'] = []
file_list['i_AMPA'][1]['modinh0'] = []
file_list['i_AMPA'][1]['modinh3'] = []
file_list['i_AMPA'][2]['modinh0'] = []
file_list['i_AMPA'][2]['modinh3'] = []
file_list['i_AMPA'][3]['modinh0'] = []
file_list['i_AMPA'][3]['modinh3'] = []
file_list['i_AMPA'][4]['modinh0'] = []
file_list['i_AMPA'][4]['modinh3'] = []

file_list['i_NMDA'][0]['modinh0'] = []
file_list['i_NMDA'][0]['modinh3'] = []
file_list['i_NMDA'][1]['modinh0'] = []
file_list['i_NMDA'][1]['modinh3'] = []
file_list['i_NMDA'][2]['modinh0'] = []
file_list['i_NMDA'][2]['modinh3'] = []
file_list['i_NMDA'][3]['modinh0'] = []
file_list['i_NMDA'][3]['modinh3'] = []
file_list['i_NMDA'][4]['modinh0'] = []
file_list['i_NMDA'][4]['modinh3'] = []

file_list['i_GABA'][0]['modinh0'] = []
file_list['i_GABA'][0]['modinh3'] = []
file_list['i_GABA'][1]['modinh0'] = []
file_list['i_GABA'][1]['modinh3'] = []
file_list['i_GABA'][2]['modinh0'] = []
file_list['i_GABA'][2]['modinh3'] = []
file_list['i_GABA'][3]['modinh0'] = []
file_list['i_GABA'][3]['modinh3'] = []
file_list['i_GABA'][4]['modinh0'] = []
file_list['i_GABA'][4]['modinh3'] = []

rec_filenames['i_AMPA'][0]['modinh0'] = 'output060916 - cell141 - e3200-i600-subtr_0-i_AMPA-inh0 - 10 trials'
rec_filenames['i_AMPA'][0]['modinh3'] = 'output060916 - cell141 - e3200-i600-subtr_0-i_AMPA-inh3 - 10 trials'
rec_filenames['i_AMPA'][1]['modinh0'] = 'output060916 - cell143 - e3200-i600-subtr_1-i_AMPA-inh0 - 10 trials'
rec_filenames['i_AMPA'][1]['modinh3'] = 'output060916 - cell143 - e3200-i600-subtr_1-i_AMPA-inh3 - 10 trials'
rec_filenames['i_AMPA'][2]['modinh0'] = 'output060916 - cell144 - e3200-i600-subtr_2-i_AMPA-inh0 - 10 trials'
rec_filenames['i_AMPA'][2]['modinh3'] = 'output060916 - cell144 - e3200-i600-subtr_2-i_AMPA-inh3 - 10 trials'
rec_filenames['i_AMPA'][3]['modinh0'] = 'output060916 - cell145 - e3200-i600-subtr_3-i_AMPA-inh0 - 10 trials'
rec_filenames['i_AMPA'][3]['modinh3'] = 'output060916 - cell145 - e3200-i600-subtr_3-i_AMPA-inh3 - 10 trials'
rec_filenames['i_AMPA'][4]['modinh0'] = 'output060916 - cell146 - e3200-i600-subtr_4-i_AMPA-inh0 - 10 trials'
rec_filenames['i_AMPA'][4]['modinh3'] = 'output060916 - cell146 - e3200-i600-subtr_4-i_AMPA-inh3 - 10 trials'

rec_filenames['i_NMDA'][0]['modinh0'] = 'output060916 - cell141 - e3200-i600-subtr_0-i_NMDA-inh0 - 10 trials'
rec_filenames['i_NMDA'][0]['modinh3'] = 'output060916 - cell141 - e3200-i600-subtr_0-i_NMDA-inh3 - 10 trials'
rec_filenames['i_NMDA'][1]['modinh0'] = 'output060916 - cell143 - e3200-i600-subtr_1-i_NMDA-inh0 - 10 trials'
rec_filenames['i_NMDA'][1]['modinh3'] = 'output060916 - cell143 - e3200-i600-subtr_1-i_NMDA-inh3 - 10 trials'
rec_filenames['i_NMDA'][2]['modinh0'] = 'output060916 - cell144 - e3200-i600-subtr_2-i_NMDA-inh0 - 10 trials'
rec_filenames['i_NMDA'][2]['modinh3'] = 'output060916 - cell144 - e3200-i600-subtr_2-i_NMDA-inh3 - 10 trials'
rec_filenames['i_NMDA'][3]['modinh0'] = 'output060916 - cell145 - e3200-i600-subtr_3-i_NMDA-inh0 - 10 trials'
rec_filenames['i_NMDA'][3]['modinh3'] = 'output060916 - cell145 - e3200-i600-subtr_3-i_NMDA-inh3 - 10 trials'
rec_filenames['i_NMDA'][4]['modinh0'] = 'output060916 - cell146 - e3200-i600-subtr_4-i_NMDA-inh0 - 10 trials'
rec_filenames['i_NMDA'][4]['modinh3'] = 'output060916 - cell146 - e3200-i600-subtr_4-i_NMDA-inh3 - 10 trials'

rec_filenames['i_GABA'][0]['modinh0'] = 'output060916 - cell141 - e3200-i600-subtr_0-i_GABA-inh0 - 10 trials'
rec_filenames['i_GABA'][0]['modinh3'] = 'output060916 - cell141 - e3200-i600-subtr_0-i_GABA-inh3 - 10 trials'
rec_filenames['i_GABA'][1]['modinh0'] = 'output060916 - cell143 - e3200-i600-subtr_1-i_GABA-inh0 - 10 trials'
rec_filenames['i_GABA'][1]['modinh3'] = 'output060916 - cell143 - e3200-i600-subtr_1-i_GABA-inh3 - 10 trials'
rec_filenames['i_GABA'][2]['modinh0'] = 'output060916 - cell144 - e3200-i600-subtr_2-i_GABA-inh0 - 10 trials'
rec_filenames['i_GABA'][2]['modinh3'] = 'output060916 - cell144 - e3200-i600-subtr_2-i_GABA-inh3 - 10 trials'
rec_filenames['i_GABA'][3]['modinh0'] = 'output060916 - cell145 - e3200-i600-subtr_3-i_GABA-inh0 - 10 trials'
rec_filenames['i_GABA'][3]['modinh3'] = 'output060916 - cell145 - e3200-i600-subtr_3-i_GABA-inh3 - 10 trials'
rec_filenames['i_GABA'][4]['modinh0'] = 'output060916 - cell146 - e3200-i600-subtr_4-i_GABA-inh0 - 10 trials'
rec_filenames['i_GABA'][4]['modinh3'] = 'output060916 - cell146 - e3200-i600-subtr_4-i_GABA-inh3 - 10 trials'

compress_i_syn_rec_files(file_list[syn_type][seed]['modinh'+str(condition)])
combine_output_files(file_list[syn_type][seed]['modinh'+str(condition)],
                     rec_filenames[syn_type][seed]['modinh'+str(condition)])