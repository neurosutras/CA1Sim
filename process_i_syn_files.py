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

file_list['i_GABA'][0]['modinh3'] = ['output051120161658-pid15614-seed0-e3200-i500-mod_inh2-i_GABA_2.5_20',
 'output051120161658-pid15612-seed0-e3200-i500-mod_inh2-i_GABA_2.5_21',
 'output051120161658-pid15613-seed0-e3200-i500-mod_inh2-i_GABA_2.5_22',
 'output051120161658-pid10369-seed0-e3200-i500-mod_inh2-i_GABA_2.5_23',
 'output051120161658-pid10374-seed0-e3200-i500-mod_inh2-i_GABA_2.5_24',
 'output051120161658-pid57559-seed0-e3200-i500-mod_inh2-i_GABA_2.5_25',
 'output051120161658-pid57558-seed0-e3200-i500-mod_inh2-i_GABA_2.5_26',
 'output051120161658-pid57557-seed0-e3200-i500-mod_inh2-i_GABA_2.5_27',
 'output051120161658-pid18383-seed0-e3200-i500-mod_inh2-i_GABA_2.5_28',
 'output051120161658-pid18385-seed0-e3200-i500-mod_inh2-i_GABA_2.5_29']
file_list['i_GABA'][1]['modinh3'] = ['output051120161727-pid13169-seed1-e3200-i500-mod_inh2-i_GABA_2.5_50',
 'output051120161727-pid28442-seed1-e3200-i500-mod_inh2-i_GABA_2.5_51',
 'output051120161727-pid24919-seed1-e3200-i500-mod_inh2-i_GABA_2.5_52',
 'output051120161727-pid24929-seed1-e3200-i500-mod_inh2-i_GABA_2.5_53',
 'output051120161727-pid4338-seed1-e3200-i500-mod_inh2-i_GABA_2.5_54',
 'output051120161727-pid22673-seed1-e3200-i500-mod_inh2-i_GABA_2.5_55',
 'output051120161727-pid22672-seed1-e3200-i500-mod_inh2-i_GABA_2.5_56',
 'output051120161727-pid13504-seed1-e3200-i500-mod_inh2-i_GABA_2.5_57',
 'output051120161727-pid9878-seed1-e3200-i500-mod_inh2-i_GABA_2.5_58',
 'output051120161727-pid9880-seed1-e3200-i500-mod_inh2-i_GABA_2.5_59']
file_list['i_GABA'][2]['modinh3'] = ['output051120161728-pid25163-seed2-e3200-i500-mod_inh2-i_GABA_2.5_80',
 'output051120161728-pid25162-seed2-e3200-i500-mod_inh2-i_GABA_2.5_81',
 'output051120161728-pid5001-seed2-e3200-i500-mod_inh2-i_GABA_2.5_82',
 'output051120161728-pid32399-seed2-e3200-i500-mod_inh2-i_GABA_2.5_83',
 'output051120161728-pid47428-seed2-e3200-i500-mod_inh2-i_GABA_2.5_84',
 'output051120161728-pid13990-seed2-e3200-i500-mod_inh2-i_GABA_2.5_85',
 'output051120161728-pid42950-seed2-e3200-i500-mod_inh2-i_GABA_2.5_86',
 'output051120161728-pid42951-seed2-e3200-i500-mod_inh2-i_GABA_2.5_87',
 'output051120161728-pid25445-seed2-e3200-i500-mod_inh2-i_GABA_2.5_88',
 'output051120161728-pid45406-seed2-e3200-i500-mod_inh2-i_GABA_2.5_89']
file_list['i_GABA'][3]['modinh3'] = ['output051120161733-pid84204-seed3-e3200-i500-mod_inh2-i_GABA_2.5_110',
 'output051120161733-pid84220-seed3-e3200-i500-mod_inh2-i_GABA_2.5_111',
 'output051120161734-pid29341-seed3-e3200-i500-mod_inh2-i_GABA_2.5_112',
 'output051120161734-pid29342-seed3-e3200-i500-mod_inh2-i_GABA_2.5_113',
 'output051120161734-pid29343-seed3-e3200-i500-mod_inh2-i_GABA_2.5_114',
 'output051120161733-pid41382-seed3-e3200-i500-mod_inh2-i_GABA_2.5_115',
 'output051120161733-pid1477-seed3-e3200-i500-mod_inh2-i_GABA_2.5_116',
 'output051120161733-pid10828-seed3-e3200-i500-mod_inh2-i_GABA_2.5_117',
 'output051120161734-pid1553-seed3-e3200-i500-mod_inh2-i_GABA_2.5_118',
 'output051120161734-pid1552-seed3-e3200-i500-mod_inh2-i_GABA_2.5_119']
file_list['i_GABA'][4]['modinh3'] = ['output051120161734-pid81169-seed4-e3200-i500-mod_inh2-i_GABA_2.5_140',
 'output051120161734-pid101521-seed4-e3200-i500-mod_inh2-i_GABA_2.5_141',
 'output051120161734-pid101522-seed4-e3200-i500-mod_inh2-i_GABA_2.5_142',
 'output051120161734-pid32248-seed4-e3200-i500-mod_inh2-i_GABA_2.5_143',
 'output051120161734-pid32249-seed4-e3200-i500-mod_inh2-i_GABA_2.5_144',
 'output051120161734-pid4600-seed4-e3200-i500-mod_inh2-i_GABA_2.5_145',
 'output051120161734-pid9613-seed4-e3200-i500-mod_inh2-i_GABA_2.5_146',
 'output051120161734-pid46598-seed4-e3200-i500-mod_inh2-i_GABA_2.5_147',
 'output051120161734-pid10884-seed4-e3200-i500-mod_inh2-i_GABA_2.5_148',
 'output051120161734-pid3583-seed4-e3200-i500-mod_inh2-i_GABA_2.5_149']

file_list['i_AMPA'][0]['modinh3'] = ['output051120161658-pid17200-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_20',
 'output051120161658-pid17201-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_21',
 'output051120161658-pid15629-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_22',
 'output051120161658-pid10386-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_23',
 'output051120161658-pid10398-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_24',
 'output051120161658-pid57583-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_25',
 'output051120161658-pid18397-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_26',
 'output051120161658-pid109405-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_27',
 'output051120161658-pid31762-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_28',
 'output051120161658-pid31763-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_29']
file_list['i_AMPA'][1]['modinh3'] = ['output051120161729-pid39452-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_50',
 'output051120161729-pid15539-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_51',
 'output051120161729-pid15540-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_52',
 'output051120161729-pid15541-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_53',
 'output051120161729-pid82782-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_54',
 'output051120161729-pid82790-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_55',
 'output051120161729-pid82802-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_56',
 'output051120161729-pid82781-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_57',
 'output051120161732-pid40938-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_58',
 'output051120161732-pid10821-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_59']
file_list['i_AMPA'][2]['modinh3'] = ['output051120161744-pid37799-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_80',
 'output051120161745-pid38315-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_81',
 'output051120161746-pid41234-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_82',
 'output051120161747-pid5351-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_83',
 'output051120161747-pid5340-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_84',
 'output051120161747-pid5326-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_85',
 'output051120161749-pid30823-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_86',
 'output051120161749-pid36892-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_87',
 'output051120161749-pid36899-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_88',
 'output051120161749-pid2304-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_89']
file_list['i_AMPA'][3]['modinh3'] = ['output051120161940-pid44407-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_110',
 'output051120161940-pid31191-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_111',
 'output051120161940-pid34277-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_112',
 'output051120161940-pid34277-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_113',
 'output051120161940-pid23581-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_114',
 'output051120161940-pid10663-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_115',
 'output051120161940-pid14437-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_116',
 'output051120161940-pid44598-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_117',
 'output051120161941-pid45082-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_118',
 'output051120161941-pid34716-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_119']
file_list['i_AMPA'][4]['modinh3'] = ['output051120162113-pid39101-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_140',
 'output051120162113-pid33030-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_141',
 'output051120162113-pid24264-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_142',
 'output051120162114-pid5128-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_143',
 'output051120162114-pid44422-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_144',
 'output051120162114-pid53626-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_145',
 'output051120162114-pid23843-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_146',
 'output051120162115-pid16704-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_147',
 'output051120162115-pid10675-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_148',
 'output051120162115-pid10453-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_149']

file_list['i_NMDA'][0]['modinh3'] = ['output051120161658-pid11797-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_20',
 'output051120161658-pid11798-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_21',
 'output051120161658-pid859-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_22',
 'output051120161658-pid868-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_23',
 'output051120161658-pid844-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_24',
 'output051120161658-pid38164-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_25',
 'output051120161658-pid38181-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_26',
 'output051120161658-pid38172-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_27',
 'output051120161658-pid5978-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_28',
 'output051120161658-pid5976-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_29']
file_list['i_NMDA'][1]['modinh3'] = ['output051120161734-pid33025-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_50',
 'output051120161734-pid29377-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_51',
 'output051120161734-pid29382-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_52',
 'output051120161734-pid14127-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_53',
 'output051120161734-pid1587-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_54',
 'output051120161734-pid1566-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_55',
 'output051120161734-pid1591-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_56',
 'output051120161734-pid48006-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_57',
 'output051120161735-pid35865-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_58',
 'output051120161735-pid35874-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_59']
file_list['i_NMDA'][2]['modinh3'] = ['output051120161749-pid31006-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_80',
 'output051120161749-pid2370-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_81',
 'output051120161750-pid25261-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_82',
 'output051120161750-pid31141-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_83',
 'output051120161750-pid5247-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_84',
 'output051120161750-pid5272-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_85',
 'output051120161750-pid5269-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_86',
 'output051120161751-pid42557-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_87',
 'output051120161751-pid42554-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_88',
 'output051120161751-pid46889-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_89']
file_list['i_NMDA'][3]['modinh3'] = ['output051120161941-pid28603-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_110',
 'output051120161941-pid44449-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_111',
 'output051120161941-pid27039-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_112',
 'output051120161942-pid37966-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_113',
 'output051120161942-pid44890-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_114',
 'output051120161942-pid24988-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_115',
 'output051120161943-pid23250-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_116',
 'output051120161944-pid44622-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_117',
 'output051120162107-pid10070-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_118',
 'output051120162107-pid24992-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_119']
file_list['i_NMDA'][4]['modinh3'] = ['output051120162115-pid47320-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_140',
 'output051120162115-pid4170-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_141',
 'output051120162115-pid24870-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_142',
 'output051120162116-pid12146-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_143',
 'output051120162218-pid135020-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_144',
 'output051120162219-pid36464-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_145',
 'output051120162219-pid110737-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_146',
 'output051120162219-pid54583-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_147',
 'output051120162219-pid91321-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_148',
 'output051120162220-pid110772-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_149']

rec_filenames['i_GABA'][0]['modinh3'] = 'output052616 - cell107 - e3200-i500-i_GABA_0-inh3 - 10 trials'
rec_filenames['i_AMPA'][0]['modinh3'] = 'output052616 - cell107 - e3200-i500-i_AMPA_0-inh3 - 10 trials'
rec_filenames['i_NMDA'][0]['modinh3'] = 'output052616 - cell107 - e3200-i500-i_NMDA_0-inh3 - 10 trials'
rec_filenames['i_GABA'][1]['modinh3'] = 'output052616 - cell108 - e3200-i500-i_GABA_1-inh3 - 10 trials'
rec_filenames['i_AMPA'][1]['modinh3'] = 'output052616 - cell108 - e3200-i500-i_AMPA_1-inh3 - 10 trials'
rec_filenames['i_NMDA'][1]['modinh3'] = 'output052616 - cell108 - e3200-i500-i_NMDA_1-inh3 - 10 trials'
rec_filenames['i_GABA'][2]['modinh3'] = 'output052616 - cell109 - e3200-i500-i_GABA_2-inh3 - 10 trials'
rec_filenames['i_AMPA'][2]['modinh3'] = 'output052616 - cell109 - e3200-i500-i_AMPA_2-inh3 - 10 trials'
rec_filenames['i_NMDA'][2]['modinh3'] = 'output052616 - cell109 - e3200-i500-i_NMDA_2-inh3 - 10 trials'
rec_filenames['i_GABA'][3]['modinh3'] = 'output052616 - cell113 - e3200-i500-i_GABA_3-inh3 - 10 trials'
rec_filenames['i_AMPA'][3]['modinh3'] = 'output052616 - cell113 - e3200-i500-i_AMPA_3-inh3 - 10 trials'
rec_filenames['i_NMDA'][3]['modinh3'] = 'output052616 - cell113 - e3200-i500-i_NMDA_3-inh3 - 10 trials'
rec_filenames['i_GABA'][4]['modinh3'] = 'output052616 - cell114 - e3200-i500-i_GABA_4-inh3 - 10 trials'
rec_filenames['i_AMPA'][4]['modinh3'] = 'output052616 - cell114 - e3200-i500-i_AMPA_4-inh3 - 10 trials'
rec_filenames['i_NMDA'][4]['modinh3'] = 'output052616 - cell114 - e3200-i500-i_NMDA_4-inh3 - 10 trials'

compress_i_syn_rec_files(file_list[syn_type][seed]['modinh'+str(condition)])
combine_output_files(file_list[syn_type][seed]['modinh'+str(condition)],
                     rec_filenames[syn_type][seed]['modinh'+str(condition)])