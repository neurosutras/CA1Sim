__author__ = 'milsteina'
from function_lib import *

file_list = {}
rec_filenames = {}
for parameter in file_list, rec_filenames:
    for this_syn_type in ['i_AMPA', 'i_NMDA', 'i_GABA']:
        parameter[this_syn_type] = {}
        for this_seed in range(5):
            parameter[this_syn_type][this_seed] = {}
file_list['i_GABA'][0]['modinh0'] = ['output051120161658-pid10295-seed0-e3200-i500-mod_inh0-i_GABA_2.5_0',
 'output051120161658-pid10296-seed0-e3200-i500-mod_inh0-i_GABA_2.5_1',
 'output051120161658-pid100288-seed0-e3200-i500-mod_inh0-i_GABA_2.5_2',
 'output051120161658-pid100289-seed0-e3200-i500-mod_inh0-i_GABA_2.5_3',
 'output051120161658-pid100290-seed0-e3200-i500-mod_inh0-i_GABA_2.5_4',
 'output051120161658-pid39314-seed0-e3200-i500-mod_inh0-i_GABA_2.5_5',
 'output051120161658-pid39312-seed0-e3200-i500-mod_inh0-i_GABA_2.5_6',
 'output051120161658-pid39313-seed0-e3200-i500-mod_inh0-i_GABA_2.5_7',
 'output051120161658-pid34649-seed0-e3200-i500-mod_inh0-i_GABA_2.5_8',
 'output051120161658-pid34650-seed0-e3200-i500-mod_inh0-i_GABA_2.5_9']
file_list['i_GABA'][0]['modinh1'] = ['output051120161658-pid34675-seed0-e3200-i500-mod_inh1-i_GABA_2.5_10',
 'output051120161658-pid420-seed0-e3200-i500-mod_inh1-i_GABA_2.5_11',
 'output051120161658-pid463-seed0-e3200-i500-mod_inh1-i_GABA_2.5_12',
 'output051120161658-pid25485-seed0-e3200-i500-mod_inh1-i_GABA_2.5_13',
 'output051120161658-pid41025-seed0-e3200-i500-mod_inh1-i_GABA_2.5_14',
 'output051120161658-pid5332-seed0-e3200-i500-mod_inh1-i_GABA_2.5_15',
 'output051120161658-pid17177-seed0-e3200-i500-mod_inh1-i_GABA_2.5_16',
 'output051120161658-pid17172-seed0-e3200-i500-mod_inh1-i_GABA_2.5_17',
 'output051120161658-pid17176-seed0-e3200-i500-mod_inh1-i_GABA_2.5_18',
 'output051120161658-pid17171-seed0-e3200-i500-mod_inh1-i_GABA_2.5_19']
file_list['i_GABA'][0]['modinh2'] = ['output051120161658-pid15614-seed0-e3200-i500-mod_inh2-i_GABA_2.5_20',
 'output051120161658-pid15612-seed0-e3200-i500-mod_inh2-i_GABA_2.5_21',
 'output051120161658-pid15613-seed0-e3200-i500-mod_inh2-i_GABA_2.5_22',
 'output051120161658-pid10369-seed0-e3200-i500-mod_inh2-i_GABA_2.5_23',
 'output051120161658-pid10374-seed0-e3200-i500-mod_inh2-i_GABA_2.5_24',
 'output051120161658-pid57559-seed0-e3200-i500-mod_inh2-i_GABA_2.5_25',
 'output051120161658-pid57558-seed0-e3200-i500-mod_inh2-i_GABA_2.5_26',
 'output051120161658-pid57557-seed0-e3200-i500-mod_inh2-i_GABA_2.5_27',
 'output051120161658-pid18383-seed0-e3200-i500-mod_inh2-i_GABA_2.5_28',
 'output051120161658-pid18385-seed0-e3200-i500-mod_inh2-i_GABA_2.5_29']
file_list['i_GABA'][1]['modinh0'] = ['output051120161658-pid109393-seed1-e3200-i500-mod_inh0-i_GABA_2.5_30',
 'output051120161658-pid109392-seed1-e3200-i500-mod_inh0-i_GABA_2.5_31',
 'output051120161658-pid38142-seed1-e3200-i500-mod_inh0-i_GABA_2.5_32',
 'output051120161658-pid5954-seed1-e3200-i500-mod_inh0-i_GABA_2.5_33',
 'output051120161658-pid5953-seed1-e3200-i500-mod_inh0-i_GABA_2.5_34',
 'output051120161700-pid19072-seed1-e3200-i500-mod_inh0-i_GABA_2.5_35',
 'output051120161701-pid2264-seed1-e3200-i500-mod_inh0-i_GABA_2.5_36',
 'output051120161706-pid110669-seed1-e3200-i500-mod_inh0-i_GABA_2.5_37',
 'output051120161706-pid10299-seed1-e3200-i500-mod_inh0-i_GABA_2.5_38',
 'output051120161707-pid13775-seed1-e3200-i500-mod_inh0-i_GABA_2.5_39']
file_list['i_GABA'][1]['modinh1'] = ['output051120161707-pid10378-seed1-e3200-i500-mod_inh1-i_GABA_2.5_40',
 'output051120161707-pid10376-seed1-e3200-i500-mod_inh1-i_GABA_2.5_41',
 'output051120161727-pid46464-seed1-e3200-i500-mod_inh1-i_GABA_2.5_42',
 'output051120161727-pid46463-seed1-e3200-i500-mod_inh1-i_GABA_2.5_43',
 'output051120161727-pid64344-seed1-e3200-i500-mod_inh1-i_GABA_2.5_44',
 'output051120161727-pid64345-seed1-e3200-i500-mod_inh1-i_GABA_2.5_45',
 'output051120161727-pid2522-seed1-e3200-i500-mod_inh1-i_GABA_2.5_46',
 'output051120161727-pid2521-seed1-e3200-i500-mod_inh1-i_GABA_2.5_47',
 'output051120161727-pid82380-seed1-e3200-i500-mod_inh1-i_GABA_2.5_48',
 'output051120161727-pid33773-seed1-e3200-i500-mod_inh1-i_GABA_2.5_49']
file_list['i_GABA'][1]['modinh2'] = ['output051120161727-pid13169-seed1-e3200-i500-mod_inh2-i_GABA_2.5_50',
 'output051120161727-pid28442-seed1-e3200-i500-mod_inh2-i_GABA_2.5_51',
 'output051120161727-pid24919-seed1-e3200-i500-mod_inh2-i_GABA_2.5_52',
 'output051120161727-pid24929-seed1-e3200-i500-mod_inh2-i_GABA_2.5_53',
 'output051120161727-pid4338-seed1-e3200-i500-mod_inh2-i_GABA_2.5_54',
 'output051120161727-pid22673-seed1-e3200-i500-mod_inh2-i_GABA_2.5_55',
 'output051120161727-pid22672-seed1-e3200-i500-mod_inh2-i_GABA_2.5_56',
 'output051120161727-pid13504-seed1-e3200-i500-mod_inh2-i_GABA_2.5_57',
 'output051120161727-pid9878-seed1-e3200-i500-mod_inh2-i_GABA_2.5_58',
 'output051120161727-pid9880-seed1-e3200-i500-mod_inh2-i_GABA_2.5_59']
file_list['i_GABA'][2]['modinh0'] = ['output051120161727-pid2854-seed2-e3200-i500-mod_inh0-i_GABA_2.5_60',
 'output051120161727-pid48725-seed2-e3200-i500-mod_inh0-i_GABA_2.5_61',
 'output051120161727-pid64380-seed2-e3200-i500-mod_inh0-i_GABA_2.5_62',
 'output051120161727-pid13294-seed2-e3200-i500-mod_inh0-i_GABA_2.5_63',
 'output051120161727-pid13297-seed2-e3200-i500-mod_inh0-i_GABA_2.5_64',
 'output051120161727-pid29029-seed2-e3200-i500-mod_inh0-i_GABA_2.5_65',
 'output051120161727-pid4537-seed2-e3200-i500-mod_inh0-i_GABA_2.5_66',
 'output051120161727-pid4536-seed2-e3200-i500-mod_inh0-i_GABA_2.5_67',
 'output051120161727-pid4540-seed2-e3200-i500-mod_inh0-i_GABA_2.5_68',
 'output051120161727-pid41316-seed2-e3200-i500-mod_inh0-i_GABA_2.5_69']
file_list['i_GABA'][2]['modinh1'] = ['output051120161727-pid31545-seed2-e3200-i500-mod_inh1-i_GABA_2.5_70',
 'output051120161727-pid41885-seed2-e3200-i500-mod_inh1-i_GABA_2.5_71',
 'output051120161727-pid4443-seed2-e3200-i500-mod_inh1-i_GABA_2.5_72',
 'output051120161727-pid45006-seed2-e3200-i500-mod_inh1-i_GABA_2.5_73',
 'output051120161728-pid42118-seed2-e3200-i500-mod_inh1-i_GABA_2.5_74',
 'output051120161728-pid42117-seed2-e3200-i500-mod_inh1-i_GABA_2.5_75',
 'output051120161728-pid32309-seed2-e3200-i500-mod_inh1-i_GABA_2.5_76',
 'output051120161728-pid13620-seed2-e3200-i500-mod_inh1-i_GABA_2.5_77',
 'output051120161728-pid98349-seed2-e3200-i500-mod_inh1-i_GABA_2.5_78',
 'output051120161728-pid41570-seed2-e3200-i500-mod_inh1-i_GABA_2.5_79']
file_list['i_GABA'][2]['modinh2'] = ['output051120161728-pid25163-seed2-e3200-i500-mod_inh2-i_GABA_2.5_80',
 'output051120161728-pid25162-seed2-e3200-i500-mod_inh2-i_GABA_2.5_81',
 'output051120161728-pid5001-seed2-e3200-i500-mod_inh2-i_GABA_2.5_82',
 'output051120161728-pid32399-seed2-e3200-i500-mod_inh2-i_GABA_2.5_83',
 'output051120161728-pid47428-seed2-e3200-i500-mod_inh2-i_GABA_2.5_84',
 'output051120161728-pid13990-seed2-e3200-i500-mod_inh2-i_GABA_2.5_85',
 'output051120161728-pid42950-seed2-e3200-i500-mod_inh2-i_GABA_2.5_86',
 'output051120161728-pid42951-seed2-e3200-i500-mod_inh2-i_GABA_2.5_87',
 'output051120161728-pid25445-seed2-e3200-i500-mod_inh2-i_GABA_2.5_88',
 'output051120161728-pid45406-seed2-e3200-i500-mod_inh2-i_GABA_2.5_89']
file_list['i_GABA'][3]['modinh0'] = ['output051120161729-pid5313-seed3-e3200-i500-mod_inh0-i_GABA_2.5_90',
 'output051120161729-pid45591-seed3-e3200-i500-mod_inh0-i_GABA_2.5_91',
 'output051120161729-pid43508-seed3-e3200-i500-mod_inh0-i_GABA_2.5_92',
 'output051120161729-pid25673-seed3-e3200-i500-mod_inh0-i_GABA_2.5_93',
 'output051120161729-pid25675-seed3-e3200-i500-mod_inh0-i_GABA_2.5_94',
 'output051120161729-pid9158-seed3-e3200-i500-mod_inh0-i_GABA_2.5_95',
 'output051120161733-pid31948-seed3-e3200-i500-mod_inh0-i_GABA_2.5_96',
 'output051120161733-pid1369-seed3-e3200-i500-mod_inh0-i_GABA_2.5_97',
 'output051120161733-pid18803-seed3-e3200-i500-mod_inh0-i_GABA_2.5_98',
 'output051120161733-pid38460-seed3-e3200-i500-mod_inh0-i_GABA_2.5_99']
file_list['i_GABA'][3]['modinh1'] = ['output051120161733-pid84178-seed3-e3200-i500-mod_inh1-i_GABA_2.5_100',
 'output051120161733-pid17507-seed3-e3200-i500-mod_inh1-i_GABA_2.5_101',
 'output051120161733-pid41333-seed3-e3200-i500-mod_inh1-i_GABA_2.5_102',
 'output051120161733-pid41335-seed3-e3200-i500-mod_inh1-i_GABA_2.5_103',
 'output051120161733-pid52398-seed3-e3200-i500-mod_inh1-i_GABA_2.5_104',
 'output051120161733-pid52397-seed3-e3200-i500-mod_inh1-i_GABA_2.5_105',
 'output051120161733-pid52399-seed3-e3200-i500-mod_inh1-i_GABA_2.5_106',
 'output051120161733-pid128977-seed3-e3200-i500-mod_inh1-i_GABA_2.5_107',
 'output051120161733-pid1436-seed3-e3200-i500-mod_inh1-i_GABA_2.5_108',
 'output051120161733-pid18840-seed3-e3200-i500-mod_inh1-i_GABA_2.5_109']
file_list['i_GABA'][3]['modinh2'] = ['output051120161733-pid84204-seed3-e3200-i500-mod_inh2-i_GABA_2.5_110',
 'output051120161733-pid84220-seed3-e3200-i500-mod_inh2-i_GABA_2.5_111',
 'output051120161734-pid29341-seed3-e3200-i500-mod_inh2-i_GABA_2.5_112',
 'output051120161734-pid29342-seed3-e3200-i500-mod_inh2-i_GABA_2.5_113',
 'output051120161734-pid29343-seed3-e3200-i500-mod_inh2-i_GABA_2.5_114',
 'output051120161733-pid41382-seed3-e3200-i500-mod_inh2-i_GABA_2.5_115',
 'output051120161733-pid1477-seed3-e3200-i500-mod_inh2-i_GABA_2.5_116',
 'output051120161733-pid10828-seed3-e3200-i500-mod_inh2-i_GABA_2.5_117',
 'output051120161734-pid1553-seed3-e3200-i500-mod_inh2-i_GABA_2.5_118',
 'output051120161734-pid1552-seed3-e3200-i500-mod_inh2-i_GABA_2.5_119']
file_list['i_GABA'][4]['modinh0'] = ['output051120161734-pid1554-seed4-e3200-i500-mod_inh0-i_GABA_2.5_120',
 'output051120161734-pid18873-seed4-e3200-i500-mod_inh0-i_GABA_2.5_121',
 'output051120161734-pid32198-seed4-e3200-i500-mod_inh0-i_GABA_2.5_122',
 'output051120161734-pid5217-seed4-e3200-i500-mod_inh0-i_GABA_2.5_123',
 'output051120161734-pid33490-seed4-e3200-i500-mod_inh0-i_GABA_2.5_124',
 'output051120161734-pid38828-seed4-e3200-i500-mod_inh0-i_GABA_2.5_125',
 'output051120161734-pid45938-seed4-e3200-i500-mod_inh0-i_GABA_2.5_126',
 'output051120161734-pid2275-seed4-e3200-i500-mod_inh0-i_GABA_2.5_127',
 'output051120161734-pid2272-seed4-e3200-i500-mod_inh0-i_GABA_2.5_128',
 'output051120161734-pid2273-seed4-e3200-i500-mod_inh0-i_GABA_2.5_129']
file_list['i_GABA'][4]['modinh1'] = ['output051120161734-pid43208-seed4-e3200-i500-mod_inh1-i_GABA_2.5_130',
 'output051120161734-pid1514-seed4-e3200-i500-mod_inh1-i_GABA_2.5_131',
 'output051120161734-pid46430-seed4-e3200-i500-mod_inh1-i_GABA_2.5_132',
 'output051120161734-pid46436-seed4-e3200-i500-mod_inh1-i_GABA_2.5_133',
 'output051120161734-pid46439-seed4-e3200-i500-mod_inh1-i_GABA_2.5_134',
 'output051120161734-pid10860-seed4-e3200-i500-mod_inh1-i_GABA_2.5_135',
 'output051120161734-pid26593-seed4-e3200-i500-mod_inh1-i_GABA_2.5_136',
 'output051120161734-pid38869-seed4-e3200-i500-mod_inh1-i_GABA_2.5_137',
 'output051120161734-pid47223-seed4-e3200-i500-mod_inh1-i_GABA_2.5_138',
 'output051120161734-pid81157-seed4-e3200-i500-mod_inh1-i_GABA_2.5_139']
file_list['i_GABA'][4]['modinh2'] = ['output051120161734-pid81169-seed4-e3200-i500-mod_inh2-i_GABA_2.5_140',
 'output051120161734-pid101521-seed4-e3200-i500-mod_inh2-i_GABA_2.5_141',
 'output051120161734-pid101522-seed4-e3200-i500-mod_inh2-i_GABA_2.5_142',
 'output051120161734-pid32248-seed4-e3200-i500-mod_inh2-i_GABA_2.5_143',
 'output051120161734-pid32249-seed4-e3200-i500-mod_inh2-i_GABA_2.5_144',
 'output051120161734-pid4600-seed4-e3200-i500-mod_inh2-i_GABA_2.5_145',
 'output051120161734-pid9613-seed4-e3200-i500-mod_inh2-i_GABA_2.5_146',
 'output051120161734-pid46598-seed4-e3200-i500-mod_inh2-i_GABA_2.5_147',
 'output051120161734-pid10884-seed4-e3200-i500-mod_inh2-i_GABA_2.5_148',
 'output051120161734-pid3583-seed4-e3200-i500-mod_inh2-i_GABA_2.5_149']

file_list['i_AMPA'][0]['modinh0'] = ['output051120161658-pid100314-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_0',
 'output051120161658-pid100305-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_1',
 'output051120161658-pid44512-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_2',
 'output051120161658-pid44513-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_3',
 'output051120161658-pid2461-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_4',
 'output051120161658-pid2460-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_5',
 'output051120161658-pid25526-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_6',
 'output051120161658-pid25538-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_7',
 'output051120161658-pid39335-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_8',
 'output051120161658-pid39338-seed0-e3200-i500-mod_inh0-i_AMPA_2.5_9']
file_list['i_AMPA'][0]['modinh1'] = ['output051120161658-pid397-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_10',
 'output051120161658-pid22021-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_11',
 'output051120161658-pid33763-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_12',
 'output051120161658-pid33760-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_13',
 'output051120161658-pid4372-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_14',
 'output051120161658-pid25961-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_15',
 'output051120161658-pid46532-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_16',
 'output051120161658-pid33317-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_17',
 'output051120161658-pid46758-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_18',
 'output051120161658-pid30027-seed0-e3200-i500-mod_inh1-i_AMPA_2.5_19']
file_list['i_AMPA'][0]['modinh2'] = ['output051120161658-pid17200-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_20',
 'output051120161658-pid17201-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_21',
 'output051120161658-pid15629-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_22',
 'output051120161658-pid10386-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_23',
 'output051120161658-pid10398-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_24',
 'output051120161658-pid57583-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_25',
 'output051120161658-pid18397-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_26',
 'output051120161658-pid109405-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_27',
 'output051120161658-pid31762-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_28',
 'output051120161658-pid31763-seed0-e3200-i500-mod_inh2-i_AMPA_2.5_29']
file_list['i_AMPA'][1]['modinh0'] = ['output051120161714-pid115531-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_30',
 'output051120161714-pid45337-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_31',
 'output051120161714-pid62495-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_32',
 'output051120161714-pid62496-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_33',
 'output051120161714-pid115572-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_34',
 'output051120161714-pid115596-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_35',
 'output051120161714-pid115592-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_36',
 'output051120161714-pid38650-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_37',
 'output051120161714-pid86519-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_38',
 'output051120161714-pid86520-seed1-e3200-i500-mod_inh0-i_AMPA_2.5_39']
file_list['i_AMPA'][1]['modinh1'] = ['output051120161726-pid35258-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_40',
 'output051120161726-pid35257-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_41',
 'output051120161726-pid35259-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_42',
 'output051120161726-pid24127-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_43',
 'output051120161727-pid44882-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_44',
 'output051120161727-pid9705-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_45',
 'output051120161727-pid46566-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_46',
 'output051120161727-pid2578-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_47',
 'output051120161727-pid48676-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_48',
 'output051120161727-pid48675-seed1-e3200-i500-mod_inh1-i_AMPA_2.5_49']
file_list['i_AMPA'][1]['modinh2'] = ['output051120161729-pid39452-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_50',
 'output051120161729-pid15539-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_51',
 'output051120161729-pid15540-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_52',
 'output051120161729-pid15541-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_53',
 'output051120161729-pid82782-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_54',
 'output051120161729-pid82790-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_55',
 'output051120161729-pid82802-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_56',
 'output051120161729-pid82781-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_57',
 'output051120161732-pid40938-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_58',
 'output051120161732-pid10821-seed1-e3200-i500-mod_inh2-i_AMPA_2.5_59']
file_list['i_AMPA'][2]['modinh0'] = ['output051120161735-pid35871-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_60',
 'output051120161735-pid81771-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_61',
 'output051120161735-pid37671-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_62',
 'output051120161735-pid37674-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_63',
 'output051120161735-pid37652-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_64',
 'output051120161739-pid34487-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_65',
 'output051120161739-pid73718-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_66',
 'output051120161739-pid73747-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_67',
 'output051120161739-pid73746-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_68',
 'output051120161739-pid73730-seed2-e3200-i500-mod_inh0-i_AMPA_2.5_69']
file_list['i_AMPA'][2]['modinh1'] = ['output051120161739-pid33436-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_70',
 'output051120161739-pid33433-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_71',
 'output051120161739-pid34739-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_72',
 'output051120161739-pid34742-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_73',
 'output051120161739-pid42040-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_74',
 'output051120161739-pid38117-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_75',
 'output051120161739-pid38129-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_76',
 'output051120161739-pid33490-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_77',
 'output051120161739-pid33478-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_78',
 'output051120161739-pid34783-seed2-e3200-i500-mod_inh1-i_AMPA_2.5_79']
file_list['i_AMPA'][2]['modinh2'] = ['output051120161744-pid37799-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_80',
 'output051120161745-pid38315-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_81',
 'output051120161746-pid41234-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_82',
 'output051120161747-pid5351-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_83',
 'output051120161747-pid5340-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_84',
 'output051120161747-pid5326-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_85',
 'output051120161749-pid30823-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_86',
 'output051120161749-pid36892-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_87',
 'output051120161749-pid36899-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_88',
 'output051120161749-pid2304-seed2-e3200-i500-mod_inh2-i_AMPA_2.5_89']
file_list['i_AMPA'][3]['modinh0'] = ['output051120161751-pid46910-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_90',
 'output051120161751-pid46907-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_91',
 'output051120161819-pid115964-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_92',
 'output051120161819-pid115957-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_93',
 'output051120161819-pid115961-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_94',
 'output051120161819-pid12237-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_95',
 'output051120161931-pid40712-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_96',
 'output051120161932-pid23562-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_97',
 'output051120161933-pid14193-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_98',
 'output051120161933-pid44423-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_99']
file_list['i_AMPA'][3]['modinh1'] = ['output051120161934-pid23372-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_100',
 'output051120161934-pid47736-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_101',
 'output051120161934-pid28605-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_102',
 'output051120161934-pid10044-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_103',
 'output051120161934-pid44262-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_104',
 'output051120161934-pid33775-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_105',
 'output051120161934-pid44467-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_106',
 'output051120161935-pid21165-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_107',
 'output051120161936-pid18911-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_108',
 'output051120161936-pid16108-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_109']
file_list['i_AMPA'][3]['modinh2'] = ['output051120161940-pid44407-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_110',
 'output051120161940-pid31191-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_111',
 'output051120161940-pid34277-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_112',
 'output051120161940-pid34277-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_113',
 'output051120161940-pid23581-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_114',
 'output051120161940-pid10663-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_115',
 'output051120161940-pid14437-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_116',
 'output051120161940-pid44598-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_117',
 'output051120161941-pid45082-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_118',
 'output051120161941-pid34716-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_119']
file_list['i_AMPA'][4]['modinh0'] = ['output051120162107-pid22198-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_120',
 'output051120162108-pid9740-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_121',
 'output051120162108-pid38985-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_122',
 'output051120162108-pid32881-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_123',
 'output051120162108-pid46726-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_124',
 'output051120162108-pid3425-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_125',
 'output051120162108-pid11537-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_126',
 'output051120162109-pid23555-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_127',
 'output051120162109-pid33307-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_128',
 'output051120162109-pid47913-seed4-e3200-i500-mod_inh0-i_AMPA_2.5_129']
file_list['i_AMPA'][4]['modinh1'] = ['output051120162110-pid40568-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_130',
 'output051120162110-pid26087-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_131',
 'output051120162112-pid10155-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_132',
 'output051120162110-pid32963-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_133',
 'output051120162110-pid46785-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_134',
 'output051120162111-pid3496-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_135',
 'output051120162111-pid11600-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_136',
 'output051120162111-pid47973-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_137',
 'output051120162111-pid26127-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_138',
 'output051120162112-pid1987-seed4-e3200-i500-mod_inh1-i_AMPA_2.5_139']
file_list['i_AMPA'][4]['modinh2'] = ['output051120162113-pid39101-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_140',
 'output051120162113-pid33030-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_141',
 'output051120162113-pid24264-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_142',
 'output051120162114-pid5128-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_143',
 'output051120162114-pid44422-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_144',
 'output051120162114-pid53626-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_145',
 'output051120162114-pid23843-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_146',
 'output051120162115-pid16704-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_147',
 'output051120162115-pid10675-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_148',
 'output051120162115-pid10453-seed4-e3200-i500-mod_inh2-i_AMPA_2.5_149']

file_list['i_NMDA'][0]['modinh0'] = ['output051120161658-pid401-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_0',
 'output051120161658-pid402-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_1',
 'output051120161658-pid2436-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_2',
 'output051120161658-pid2418-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_3',
 'output051120161658-pid2442-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_4',
 'output051120161658-pid34660-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_5',
 'output051120161658-pid16103-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_6',
 'output051120161658-pid16102-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_7',
 'output051120161658-pid22269-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_8',
 'output051120161658-pid22270-seed0-e3200-i500-mod_inh0-i_NMDA_2.5_9']
file_list['i_NMDA'][0]['modinh1'] = ['output051120161658-pid14927-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_10',
 'output051120161658-pid22162-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_11',
 'output051120161658-pid48372-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_12',
 'output051120161658-pid37982-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_13',
 'output051120161658-pid2867-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_14',
 'output051120161658-pid33600-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_15',
 'output051120161658-pid33599-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_16',
 'output051120161658-pid33601-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_17',
 'output051120161658-pid39790-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_18',
 'output051120161658-pid5344-seed0-e3200-i500-mod_inh1-i_NMDA_2.5_19']
file_list['i_NMDA'][0]['modinh2'] = ['output051120161658-pid11797-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_20',
 'output051120161658-pid11798-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_21',
 'output051120161658-pid859-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_22',
 'output051120161658-pid868-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_23',
 'output051120161658-pid844-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_24',
 'output051120161658-pid38164-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_25',
 'output051120161658-pid38181-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_26',
 'output051120161658-pid38172-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_27',
 'output051120161658-pid5978-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_28',
 'output051120161658-pid5976-seed0-e3200-i500-mod_inh2-i_NMDA_2.5_29']
file_list['i_NMDA'][1]['modinh0'] = ['output051120161714-pid86518-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_30',
 'output051120161714-pid30475-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_31',
 'output051120161719-pid32922-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_32',
 'output051120161719-pid2604-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_33',
 'output051120161725-pid42889-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_34',
 'output051120161725-pid40925-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_35',
 'output051120161725-pid48590-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_36',
 'output051120161725-pid43041-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_37',
 'output051120161725-pid2963-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_38',
 'output051120161726-pid12323-seed1-e3200-i500-mod_inh0-i_NMDA_2.5_39']
file_list['i_NMDA'][1]['modinh1'] = ['output051120161727-pid48679-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_40',
 'output051120161727-pid41278-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_41',
 'output051120161727-pid44940-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_42',
 'output051120161727-pid9914-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_43',
 'output051120161727-pid46953-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_44',
 'output051120161727-pid38873-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_45',
 'output051120161727-pid45000-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_46',
 'output051120161728-pid98278-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_47',
 'output051120161728-pid3413-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_48',
 'output051120161728-pid35910-seed1-e3200-i500-mod_inh1-i_NMDA_2.5_49']
file_list['i_NMDA'][1]['modinh2'] = ['output051120161734-pid33025-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_50',
 'output051120161734-pid29377-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_51',
 'output051120161734-pid29382-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_52',
 'output051120161734-pid14127-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_53',
 'output051120161734-pid1587-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_54',
 'output051120161734-pid1566-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_55',
 'output051120161734-pid1591-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_56',
 'output051120161734-pid48006-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_57',
 'output051120161735-pid35865-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_58',
 'output051120161735-pid35874-seed1-e3200-i500-mod_inh2-i_NMDA_2.5_59']
file_list['i_NMDA'][2]['modinh0'] = ['output051120161739-pid30842-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_60',
 'output051120161739-pid30841-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_61',
 'output051120161739-pid39088-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_62',
 'output051120161739-pid41862-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_63',
 'output051120161739-pid41861-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_64',
 'output051120161739-pid38078-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_65',
 'output051120161739-pid38077-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_66',
 'output051120161739-pid37926-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_67',
 'output051120161739-pid37924-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_68',
 'output051120161739-pid37938-seed2-e3200-i500-mod_inh0-i_NMDA_2.5_69']
file_list['i_NMDA'][2]['modinh1'] = ['output051120161739-pid42185-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_70',
 'output051120161740-pid30910-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_71',
 'output051120161740-pid30920-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_72',
 'output051120161743-pid73939-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_73',
 'output051120161743-pid73958-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_74',
 'output051120161743-pid73956-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_75',
 'output051120161744-pid25942-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_76',
 'output051120161744-pid25941-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_77',
 'output051120161744-pid25991-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_78',
 'output051120161744-pid26003-seed2-e3200-i500-mod_inh1-i_NMDA_2.5_79']
file_list['i_NMDA'][2]['modinh2'] = ['output051120161749-pid31006-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_80',
 'output051120161749-pid2370-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_81',
 'output051120161750-pid25261-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_82',
 'output051120161750-pid31141-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_83',
 'output051120161750-pid5247-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_84',
 'output051120161750-pid5272-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_85',
 'output051120161750-pid5269-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_86',
 'output051120161751-pid42557-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_87',
 'output051120161751-pid42554-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_88',
 'output051120161751-pid46889-seed2-e3200-i500-mod_inh2-i_NMDA_2.5_89']
file_list['i_NMDA'][3]['modinh0'] = ['output051120161933-pid14228-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_90',
 'output051120161933-pid22862-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_91',
 'output051120161933-pid41477-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_92',
 'output051120161933-pid28484-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_93',
 'output051120161933-pid11674-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_94',
 'output051120161933-pid40641-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_95',
 'output051120161934-pid18119-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_96',
 'output051120161934-pid31062-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_97',
 'output051120161934-pid42956-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_98',
 'output051120161934-pid29421-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_99']
file_list['i_NMDA'][3]['modinh1'] = ['output051120161936-pid33368-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_100',
 'output051120161936-pid36297-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_101',
 'output051120161936-pid15330-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_102',
 'output051120161937-pid11508-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_103',
 'output051120161937-pid22929-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_104',
 'output051120161938-pid34221-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_105',
 'output051120161939-pid23426-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_106',
 'output051120161939-pid23474-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_107',
 'output051120161939-pid67043-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_108',
 'output051120161940-pid28825-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_109']
file_list['i_NMDA'][3]['modinh2'] = ['output051120161941-pid28603-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_110',
 'output051120161941-pid44449-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_111',
 'output051120161941-pid27039-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_112',
 'output051120161942-pid37966-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_113',
 'output051120161942-pid44890-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_114',
 'output051120161942-pid24988-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_115',
 'output051120161943-pid23250-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_116',
 'output051120161944-pid44622-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_117',
 'output051120162107-pid10070-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_118',
 'output051120162107-pid24992-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_119']
file_list['i_NMDA'][4]['modinh0'] = ['output051120162109-pid16042-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_120',
 'output051120162109-pid16043-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_121',
 'output051120162109-pid18934-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_122',
 'output051120162109-pid23365-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_123',
 'output051120162109-pid52954-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_124',
 'output051120162109-pid23591-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_125',
 'output051120162109-pid23724-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_126',
 'output051120162110-pid53020-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_127',
 'output051120162110-pid9917-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_128',
 'output051120162110-pid1916-seed4-e3200-i500-mod_inh0-i_NMDA_2.5_129']
file_list['i_NMDA'][4]['modinh1'] = ['output051120162112-pid53274-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_130',
 'output051120162112-pid32956-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_131',
 'output051120162112-pid9852-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_132',
 'output051120162112-pid38856-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_133',
 'output051120162113-pid18999-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_134',
 'output051120162113-pid22289-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_135',
 'output051120162113-pid40903-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_136',
 'output051120162113-pid24206-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_137',
 'output051120162113-pid40639-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_138',
 'output051120162113-pid23660-seed4-e3200-i500-mod_inh1-i_NMDA_2.5_139']
file_list['i_NMDA'][4]['modinh2'] = ['output051120162115-pid47320-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_140',
 'output051120162115-pid4170-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_141',
 'output051120162115-pid24870-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_142',
 'output051120162116-pid12146-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_143',
 'output051120162218-pid135020-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_144',
 'output051120162219-pid36464-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_145',
 'output051120162219-pid110737-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_146',
 'output051120162219-pid54583-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_147',
 'output051120162219-pid91321-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_148',
 'output051120162220-pid110772-seed4-e3200-i500-mod_inh2-i_NMDA_2.5_149']

rec_filenames['i_AMPA'][0]['modinh0'] = 'output050916 - cell107 - e3200-i500-i_AMPA_0-inh0 - 10 trials'
rec_filenames['i_AMPA'][0]['modinh1'] = 'output050916 - cell107 - e3200-i500-i_AMPA_0-inh1 - 10 trials'
rec_filenames['i_AMPA'][0]['modinh2'] = 'output050916 - cell107 - e3200-i500-i_AMPA_0-inh2 - 10 trials'

rec_filenames['i_NMDA'][0]['modinh0'] = 'output050916 - cell107 - e3200-i500-i_NMDA_0-inh0 - 10 trials'
rec_filenames['i_NMDA'][0]['modinh1'] = 'output050916 - cell107 - e3200-i500-i_NMDA_0-inh1 - 10 trials'
rec_filenames['i_NMDA'][0]['modinh2'] = 'output050916 - cell107 - e3200-i500-i_NMDA_0-inh2 - 10 trials'

rec_filenames['i_GABA'][0]['modinh0'] = 'output050916 - cell107 - e3200-i500-i_GABA_0-inh0 - 10 trials'
rec_filenames['i_GABA'][0]['modinh1'] = 'output050916 - cell107 - e3200-i500-i_GABA_0-inh1 - 10 trials'
rec_filenames['i_GABA'][0]['modinh2'] = 'output050916 - cell107 - e3200-i500-i_GABA_0-inh2 - 10 trials'

rec_filenames['i_AMPA'][1]['modinh0'] = 'output050916 - cell108 - e3200-i500-i_AMPA_1-inh0 - 10 trials'
rec_filenames['i_AMPA'][1]['modinh1'] = 'output050916 - cell108 - e3200-i500-i_AMPA_1-inh1 - 10 trials'
rec_filenames['i_AMPA'][1]['modinh2'] = 'output050916 - cell108 - e3200-i500-i_AMPA_1-inh2 - 10 trials'

rec_filenames['i_NMDA'][1]['modinh0'] = 'output050916 - cell108 - e3200-i500-i_NMDA_1-inh0 - 10 trials'
rec_filenames['i_NMDA'][1]['modinh1'] = 'output050916 - cell108 - e3200-i500-i_NMDA_1-inh1 - 10 trials'
rec_filenames['i_NMDA'][1]['modinh2'] = 'output050916 - cell108 - e3200-i500-i_NMDA_1-inh2 - 10 trials'

rec_filenames['i_GABA'][1]['modinh0'] = 'output050916 - cell108 - e3200-i500-i_GABA_1-inh0 - 10 trials'
rec_filenames['i_GABA'][1]['modinh1'] = 'output050916 - cell108 - e3200-i500-i_GABA_1-inh1 - 10 trials'
rec_filenames['i_GABA'][1]['modinh2'] = 'output050916 - cell108 - e3200-i500-i_GABA_1-inh2 - 10 trials'

rec_filenames['i_AMPA'][2]['modinh0'] = 'output050916 - cell109 - e3200-i500-i_AMPA_2-inh0 - 10 trials'
rec_filenames['i_AMPA'][2]['modinh1'] = 'output050916 - cell109 - e3200-i500-i_AMPA_2-inh1 - 10 trials'
rec_filenames['i_AMPA'][2]['modinh2'] = 'output050916 - cell109 - e3200-i500-i_AMPA_2-inh2 - 10 trials'

rec_filenames['i_NMDA'][2]['modinh0'] = 'output050916 - cell109 - e3200-i500-i_NMDA_2-inh0 - 10 trials'
rec_filenames['i_NMDA'][2]['modinh1'] = 'output050916 - cell109 - e3200-i500-i_NMDA_2-inh1 - 10 trials'
rec_filenames['i_NMDA'][2]['modinh2'] = 'output050916 - cell109 - e3200-i500-i_NMDA_2-inh2 - 10 trials'

rec_filenames['i_GABA'][2]['modinh0'] = 'output050916 - cell109 - e3200-i500-i_GABA_2-inh0 - 10 trials'
rec_filenames['i_GABA'][2]['modinh1'] = 'output050916 - cell109 - e3200-i500-i_GABA_2-inh1 - 10 trials'
rec_filenames['i_GABA'][2]['modinh2'] = 'output050916 - cell109 - e3200-i500-i_GABA_2-inh2 - 10 trials'

rec_filenames['i_AMPA'][3]['modinh0'] = 'output050916 - cell113 - e3200-i500-i_AMPA_3-inh0 - 10 trials'
rec_filenames['i_AMPA'][3]['modinh1'] = 'output050916 - cell113 - e3200-i500-i_AMPA_3-inh1 - 10 trials'
rec_filenames['i_AMPA'][3]['modinh2'] = 'output050916 - cell113 - e3200-i500-i_AMPA_3-inh2 - 10 trials'

rec_filenames['i_NMDA'][3]['modinh0'] = 'output050916 - cell113 - e3200-i500-i_NMDA_3-inh0 - 10 trials'
rec_filenames['i_NMDA'][3]['modinh1'] = 'output050916 - cell113 - e3200-i500-i_NMDA_3-inh1 - 10 trials'
rec_filenames['i_NMDA'][3]['modinh2'] = 'output050916 - cell113 - e3200-i500-i_NMDA_3-inh2 - 10 trials'

rec_filenames['i_GABA'][3]['modinh0'] = 'output050916 - cell113 - e3200-i500-i_GABA_3-inh0 - 10 trials'
rec_filenames['i_GABA'][3]['modinh1'] = 'output050916 - cell113 - e3200-i500-i_GABA_3-inh1 - 10 trials'
rec_filenames['i_GABA'][3]['modinh2'] = 'output050916 - cell113 - e3200-i500-i_GABA_3-inh2 - 10 trials'

rec_filenames['i_AMPA'][4]['modinh0'] = 'output050916 - cell114 - e3200-i500-i_AMPA_4-inh0 - 10 trials'
rec_filenames['i_AMPA'][4]['modinh1'] = 'output050916 - cell114 - e3200-i500-i_AMPA_4-inh1 - 10 trials'
rec_filenames['i_AMPA'][4]['modinh2'] = 'output050916 - cell114 - e3200-i500-i_AMPA_4-inh2 - 10 trials'

rec_filenames['i_NMDA'][4]['modinh0'] = 'output050916 - cell114 - e3200-i500-i_NMDA_4-inh0 - 10 trials'
rec_filenames['i_NMDA'][4]['modinh1'] = 'output050916 - cell114 - e3200-i500-i_NMDA_4-inh1 - 10 trials'
rec_filenames['i_NMDA'][4]['modinh2'] = 'output050916 - cell114 - e3200-i500-i_NMDA_4-inh2 - 10 trials'

rec_filenames['i_GABA'][4]['modinh0'] = 'output050916 - cell114 - e3200-i500-i_GABA_4-inh0 - 10 trials'
rec_filenames['i_GABA'][4]['modinh1'] = 'output050916 - cell114 - e3200-i500-i_GABA_4-inh1 - 10 trials'
rec_filenames['i_GABA'][4]['modinh2'] = 'output050916 - cell114 - e3200-i500-i_GABA_4-inh2 - 10 trials'

for syn_type in ['i_AMPA', 'i_NMDA', 'i_GABA']:
    for seed in range(5):
        for condition in ['modinh0', 'modinh1', 'modinh2']:
            combine_output_files(file_list[syn_type][seed][condition], rec_filenames[syn_type][seed][condition])