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

file_list['i_AMPA'][0]['modinh0'] = ['output060820161304-pid47126-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_0',
 'output060820161303-pid47756-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_1',
 'output060920161653-pid100108-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_2',
 'output060920161653-pid15101-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_3',
 'output060820161303-pid47758-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_4',
 'output061020161207-pid20403-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_5',
 'output060820161303-pid26172-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_6',
 'output060820161303-pid26171-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_7',
 'output060820161303-pid26175-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_8',
 'output060920161654-pid58898-seed0-e3200-i600-subtr_mod_inh0-i_AMPA_9']
file_list['i_AMPA'][0]['modinh3'] = ['output061020161207-pid839-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_20',
 'output060820161303-pid19190-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_21',
 'output060820161303-pid19194-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_22',
 'output060920161649-pid59644-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_23',
 'output061020161207-pid840-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_24',
 'output060820161303-pid16527-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_25',
 'output060920161649-pid112428-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_26',
 'output061020161207-pid12887-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_27',
 'output060820161303-pid23447-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_28',
 'output060920161649-pid109360-seed0-e3200-i600-subtr_mod_inh3-i_AMPA_29']
file_list['i_AMPA'][1]['modinh0'] = ['output060820161304-pid39439-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_30',
 'output060820161304-pid39440-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_31',
 'output060920161655-pid100607-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_32',
 'output060820161303-pid22791-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_33',
 'output060920161658-pid13819-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_34',
 'output061020161207-pid28718-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_35',
 'output060920161659-pid131694-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_36',
 'output060820161303-pid34490-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_37',
 'output060820161303-pid34491-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_38',
 'output060920161659-pid131780-seed1-e3200-i600-subtr_mod_inh0-i_AMPA_39']
file_list['i_AMPA'][1]['modinh3'] = ['output060820161303-pid32514-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_50',
 'output060820161303-pid32498-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_51',
 'output060820161303-pid32491-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_52',
 'output061020161207-pid29377-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_53',
 'output060920161650-pid12212-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_54',
 'output060920161650-pid54812-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_55',
 'output061020161207-pid36560-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_56',
 'output060820161303-pid48472-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_57',
 'output060820161303-pid48479-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_58',
 'output060820161303-pid48475-seed1-e3200-i600-subtr_mod_inh3-i_AMPA_59']
file_list['i_AMPA'][2]['modinh0'] = ['output060920161708-pid125840-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_60',
 'output060920161708-pid13922-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_61',
 'output060920161708-pid37136-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_62',
 'output060920161708-pid29659-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_63',
 'output060920161709-pid26035-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_64',
 'output060920161709-pid29684-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_65',
 'output060920161709-pid37254-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_66',
 'output060920161709-pid29752-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_67',
 'output060920161709-pid25103-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_68',
 'output060920161709-pid26091-seed2-e3200-i600-subtr_mod_inh0-i_AMPA_69']
file_list['i_AMPA'][2]['modinh3'] = ['output060820161303-pid17744-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_80',
 'output060920161651-pid55546-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_81',
 'output060820161303-pid3527-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_82',
 'output060920161651-pid8053-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_83',
 'output060920161651-pid97232-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_84',
 'output061520161622-pid34558-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_85',
 'output060920161651-pid39565-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_86',
 'output060820161303-pid12408-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_87',
 'output060820161303-pid12410-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_88',
 'output060820161303-pid29425-seed2-e3200-i600-subtr_mod_inh3-i_AMPA_89']
file_list['i_AMPA'][3]['modinh0'] = ['output060820161304-pid29822-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_90',
 'output060820161304-pid29832-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_91',
 'output060820161304-pid29831-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_92',
 'output060920161718-pid14449-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_93',
 'output061020161207-pid30755-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_94',
 'output060920161718-pid48060-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_95',
 'output060920161719-pid42833-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_96',
 'output061020161207-pid32317-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_97',
 'output060920161719-pid26821-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_98',
 'output060820161303-pid36166-seed3-e3200-i600-subtr_mod_inh0-i_AMPA_99']
file_list['i_AMPA'][3]['modinh3'] = ['output060920161651-pid4611-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_110',
 'output060920161651-pid4607-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_111',
 'output061020161207-pid32360-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_112',
 'output060820161303-pid43295-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_113',
 'output060920161652-pid17798-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_114',
 'output060820161303-pid26749-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_115',
 'output060920161652-pid14230-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_116',
 'output060820161303-pid26750-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_117',
 'output060920161652-pid14231-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_118',
 'output060820161303-pid30936-seed3-e3200-i600-subtr_mod_inh3-i_AMPA_119']
file_list['i_AMPA'][4]['modinh0'] = ['output060820161304-pid20342-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_120',
 'output061020161207-pid689-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_121',
 'output061020161207-pid82185-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_122',
 'output060920161724-pid27076-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_123',
 'output061020161207-pid82182-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_124',
 'output060820161303-pid48006-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_125',
 'output060820161303-pid48000-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_126',
 'output060920161724-pid27382-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_127',
 'output060820161303-pid48004-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_128',
 'output060920161724-pid27139-seed4-e3200-i600-subtr_mod_inh0-i_AMPA_129']
file_list['i_AMPA'][4]['modinh3'] = ['output060820161304-pid28407-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_140',
 'output060820161304-pid28401-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_141',
 'output060820161303-pid13482-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_142',
 'output060820161303-pid13478-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_143',
 'output060820161303-pid13479-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_144',
 'output060920161653-pid5051-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_145',
 'output060820161303-pid3126-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_146',
 'output060820161303-pid3130-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_147',
 'output060920161653-pid5054-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_148',
 'output060920161653-pid58516-seed4-e3200-i600-subtr_mod_inh3-i_AMPA_149']

file_list['i_NMDA'][0]['modinh0'] = ['output060920161654-pid30570-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_0',
 'output060920161654-pid113085-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_1',
 'output060820161303-pid13522-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_2',
 'output061020161207-pid20404-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_3',
 'output060920161654-pid13327-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_4',
 'output061020161207-pid30259-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_5',
 'output060920161655-pid13548-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_6',
 'output060820161304-pid36265-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_7',
 'output060820161304-pid36264-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_8',
 'output060820161304-pid36253-seed0-e3200-i600-subtr_mod_inh0-i_NMDA_9']
file_list['i_NMDA'][0]['modinh3'] = ['output060920161649-pid42227-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_20',
 'output060920161650-pid39134-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_21',
 'output061020161207-pid12886-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_22',
 'output061020161207-pid12885-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_23',
 'output060820161303-pid739-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_24',
 'output060820161303-pid757-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_25',
 'output060820161303-pid738-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_26',
 'output060820161303-pid40755-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_27',
 'output060820161303-pid40757-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_28',
 'output061020161207-pid39969-seed0-e3200-i600-subtr_mod_inh3-i_NMDA_29']
file_list['i_NMDA'][1]['modinh0'] = ['output060820161303-pid34488-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_30',
 'output060920161701-pid61439-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_31',
 'output060920161703-pid46586-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_32',
 'output060920161704-pid25479-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_33',
 'output061020161207-pid28717-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_34',
 'output060920161706-pid54471-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_35',
 'output060920161708-pid111366-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_36',
 'output060920161708-pid26023-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_37',
 'output060820161304-pid47580-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_38',
 'output060820161304-pid47575-seed1-e3200-i600-subtr_mod_inh0-i_NMDA_39']
file_list['i_NMDA'][1]['modinh3'] = ['output060920161650-pid109542-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_50',
 'output061020161207-pid36559-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_51',
 'output060920161650-pid96830-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_52',
 'output060920161650-pid55035-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_53',
 'output060920161650-pid59785-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_54',
 'output060920161650-pid59784-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_55',
 'output061020161207-pid48108-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_56',
 'output060820161303-pid17875-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_57',
 'output060820161303-pid17876-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_58',
 'output060820161303-pid17760-seed1-e3200-i600-subtr_mod_inh3-i_NMDA_59']
file_list['i_NMDA'][2]['modinh0'] = ['output060920161711-pid5516-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_60',
 'output060820161303-pid24243-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_61',
 'output060820161303-pid15174-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_62',
 'output060920161717-pid47971-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_63',
 'output060920161717-pid1988-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_64',
 'output061020161207-pid30754-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_65',
 'output060920161718-pid26628-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_66',
 'output060820161304-pid17655-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_67',
 'output060820161304-pid17660-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_68',
 'output060820161304-pid17653-seed2-e3200-i600-subtr_mod_inh0-i_NMDA_69']
file_list['i_NMDA'][2]['modinh3'] = ['output061020161207-pid32361-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_80',
 'output061020161207-pid32362-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_81',
 'output060920161651-pid109605-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_82',
 'output060920161651-pid42407-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_83',
 'output060920161651-pid17436-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_84',
 'output060920161651-pid55218-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_85',
 'output060920161651-pid99140-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_86',
 'output060820161304-pid37928-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_87',
 'output060920161651-pid4608-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_88',
 'output060820161304-pid37932-seed2-e3200-i600-subtr_mod_inh3-i_NMDA_89']
file_list['i_NMDA'][3]['modinh0'] = ['output060820161303-pid36168-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_90',
 'output060920161719-pid5913-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_91',
 'output061020161207-pid32316-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_92',
 'output061020161207-pid41557-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_93',
 'output060920161719-pid33617-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_94',
 'output060820161303-pid21427-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_95',
 'output060920161720-pid21842-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_96',
 'output061020161207-pid41556-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_97',
 'output060920161721-pid32459-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_98',
 'output060820161304-pid20330-seed3-e3200-i600-subtr_mod_inh0-i_NMDA_99']
file_list['i_NMDA'][3]['modinh3'] = ['output060820161303-pid30931-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_110',
 'output060920161652-pid112803-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_111',
 'output060920161652-pid56462-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_112',
 'output060920161652-pid97777-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_113',
 'output061020161207-pid3191-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_114',
 'output060920161652-pid58427-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_115',
 'output060920161652-pid112984-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_116',
 'output061020161207-pid3192-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_117',
 'output060820161303-pid16150-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_118',
 'output060920161653-pid113250-seed3-e3200-i600-subtr_mod_inh3-i_NMDA_119']
file_list['i_NMDA'][4]['modinh0'] = ['output060920161734-pid28704-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_120',
 'output060920161734-pid28736-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_121',
 'output060920161734-pid28748-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_122',
 'output060920161734-pid28660-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_123',
 'output060920161734-pid28703-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_124',
 'output060920161734-pid28700-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_125',
 'output060920161734-pid25589-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_126',
 'output061020161207-pid82183-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_127',
 'output060820161304-pid23832-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_128',
 'output061020161207-pid82184-seed4-e3200-i600-subtr_mod_inh0-i_NMDA_129']
file_list['i_NMDA'][4]['modinh3'] = ['output060920161653-pid43177-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_140',
 'output060820161303-pid38641-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_141',
 'output060920161653-pid14663-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_142',
 'output060920161653-pid30228-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_143',
 'output060820161303-pid29747-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_144',
 'output060820161303-pid29749-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_145',
 'output060920161653-pid12722-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_146',
 'output060820161303-pid29743-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_147',
 'output060820161303-pid37585-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_148',
 'output060820161303-pid37587-seed4-e3200-i600-subtr_mod_inh3-i_NMDA_149']

file_list['i_GABA'][0]['modinh0'] = ['output060820161304-pid47117-seed0-e3200-i600-subtr_mod_inh0-i_GABA_0',
 'output060820161304-pid47119-seed0-e3200-i600-subtr_mod_inh0-i_GABA_1',
 'output060820161304-pid47123-seed0-e3200-i600-subtr_mod_inh0-i_GABA_2',
 'output060820161304-pid36242-seed0-e3200-i600-subtr_mod_inh0-i_GABA_3',
 'output060920161649-pid8906-seed0-e3200-i600-subtr_mod_inh0-i_GABA_4',
 'output060820161304-pid36254-seed0-e3200-i600-subtr_mod_inh0-i_GABA_5',
 'output060820161304-pid36252-seed0-e3200-i600-subtr_mod_inh0-i_GABA_6',
 'output060820161304-pid39434-seed0-e3200-i600-subtr_mod_inh0-i_GABA_7',
 'output060820161304-pid39438-seed0-e3200-i600-subtr_mod_inh0-i_GABA_8',
 'output060920161649-pid8913-seed0-e3200-i600-subtr_mod_inh0-i_GABA_9']
file_list['i_GABA'][0]['modinh3'] = ['output060920161649-pid44728-seed0-e3200-i600-subtr_mod_inh3-i_GABA_20',
 'output060920161649-pid34792-seed0-e3200-i600-subtr_mod_inh3-i_GABA_21',
 'output060820161303-pid19173-seed0-e3200-i600-subtr_mod_inh3-i_GABA_22',
 'output060920161649-pid130248-seed0-e3200-i600-subtr_mod_inh3-i_GABA_23',
 'output060920161649-pid130249-seed0-e3200-i600-subtr_mod_inh3-i_GABA_24',
 'output060820161303-pid40753-seed0-e3200-i600-subtr_mod_inh3-i_GABA_25',
 'output060820161303-pid40747-seed0-e3200-i600-subtr_mod_inh3-i_GABA_26',
 'output060920161649-pid39347-seed0-e3200-i600-subtr_mod_inh3-i_GABA_27',
 'output060820161303-pid32494-seed0-e3200-i600-subtr_mod_inh3-i_GABA_28',
 'output061020161207-pid39966-seed0-e3200-i600-subtr_mod_inh3-i_GABA_29']
file_list['i_GABA'][1]['modinh0'] = ['output060820161304-pid39425-seed1-e3200-i600-subtr_mod_inh0-i_GABA_30',
 'output060920161649-pid45937-seed1-e3200-i600-subtr_mod_inh0-i_GABA_31',
 'output060920161649-pid4164-seed1-e3200-i600-subtr_mod_inh0-i_GABA_32',
 'output060820161304-pid47576-seed1-e3200-i600-subtr_mod_inh0-i_GABA_33',
 'output060820161304-pid47572-seed1-e3200-i600-subtr_mod_inh0-i_GABA_34',
 'output060820161304-pid47578-seed1-e3200-i600-subtr_mod_inh0-i_GABA_35',
 'output060820161304-pid47577-seed1-e3200-i600-subtr_mod_inh0-i_GABA_36',
 'output060820161304-pid47574-seed1-e3200-i600-subtr_mod_inh0-i_GABA_37',
 'output060820161304-pid47579-seed1-e3200-i600-subtr_mod_inh0-i_GABA_38',
 'output060820161304-pid47571-seed1-e3200-i600-subtr_mod_inh0-i_GABA_39']
file_list['i_GABA'][1]['modinh3'] = ['output060920161649-pid36121-seed1-e3200-i600-subtr_mod_inh3-i_GABA_50',
 'output060920161649-pid20550-seed1-e3200-i600-subtr_mod_inh3-i_GABA_51',
 'output060820161303-pid17732-seed1-e3200-i600-subtr_mod_inh3-i_GABA_52',
 'output060820161303-pid17741-seed1-e3200-i600-subtr_mod_inh3-i_GABA_53',
 'output060820161303-pid17750-seed1-e3200-i600-subtr_mod_inh3-i_GABA_54',
 'output060820161303-pid17757-seed1-e3200-i600-subtr_mod_inh3-i_GABA_55',
 'output060820161303-pid17756-seed1-e3200-i600-subtr_mod_inh3-i_GABA_56',
 'output060820161303-pid17755-seed1-e3200-i600-subtr_mod_inh3-i_GABA_57',
 'output060820161303-pid17747-seed1-e3200-i600-subtr_mod_inh3-i_GABA_58',
 'output060820161303-pid17742-seed1-e3200-i600-subtr_mod_inh3-i_GABA_59']
file_list['i_GABA'][2]['modinh0'] = ['output060820161304-pid47573-seed2-e3200-i600-subtr_mod_inh0-i_GABA_60',
 'output060820161303-pid15175-seed2-e3200-i600-subtr_mod_inh0-i_GABA_61',
 'output060820161303-pid15171-seed2-e3200-i600-subtr_mod_inh0-i_GABA_62',
 'output060920161649-pid109297-seed2-e3200-i600-subtr_mod_inh0-i_GABA_63',
 'output060820161304-pid17647-seed2-e3200-i600-subtr_mod_inh0-i_GABA_64',
 'output060820161304-pid17661-seed2-e3200-i600-subtr_mod_inh0-i_GABA_65',
 'output060820161304-pid17657-seed2-e3200-i600-subtr_mod_inh0-i_GABA_66',
 'output060820161304-pid29829-seed2-e3200-i600-subtr_mod_inh0-i_GABA_67',
 'output060920161649-pid109296-seed2-e3200-i600-subtr_mod_inh0-i_GABA_68',
 'output060820161304-pid29819-seed2-e3200-i600-subtr_mod_inh0-i_GABA_69']
file_list['i_GABA'][2]['modinh3'] = ['output060920161649-pid12049-seed2-e3200-i600-subtr_mod_inh3-i_GABA_80',
 'output060820161304-pid37929-seed2-e3200-i600-subtr_mod_inh3-i_GABA_81',
 'output060820161304-pid37926-seed2-e3200-i600-subtr_mod_inh3-i_GABA_82',
 'output060920161649-pid39996-seed2-e3200-i600-subtr_mod_inh3-i_GABA_83',
 'output060820161303-pid35857-seed2-e3200-i600-subtr_mod_inh3-i_GABA_84',
 'output060820161303-pid35866-seed2-e3200-i600-subtr_mod_inh3-i_GABA_85',
 'output060820161303-pid35873-seed2-e3200-i600-subtr_mod_inh3-i_GABA_86',
 'output060920161649-pid34463-seed2-e3200-i600-subtr_mod_inh3-i_GABA_87',
 'output060820161303-pid35877-seed2-e3200-i600-subtr_mod_inh3-i_GABA_88',
 'output060920161649-pid7705-seed2-e3200-i600-subtr_mod_inh3-i_GABA_89']
file_list['i_GABA'][3]['modinh0'] = ['output060920161649-pid146376-seed3-e3200-i600-subtr_mod_inh0-i_GABA_90',
 'output060920161649-pid146374-seed3-e3200-i600-subtr_mod_inh0-i_GABA_91',
 'output060920161649-pid29702-seed3-e3200-i600-subtr_mod_inh0-i_GABA_92',
 'output060820161304-pid20311-seed3-e3200-i600-subtr_mod_inh0-i_GABA_93',
 'output060820161304-pid20318-seed3-e3200-i600-subtr_mod_inh0-i_GABA_94',
 'output060820161304-pid20312-seed3-e3200-i600-subtr_mod_inh0-i_GABA_95',
 'output060820161304-pid20313-seed3-e3200-i600-subtr_mod_inh0-i_GABA_96',
 'output060820161304-pid20317-seed3-e3200-i600-subtr_mod_inh0-i_GABA_97',
 'output060820161304-pid20314-seed3-e3200-i600-subtr_mod_inh0-i_GABA_98',
 'output060820161304-pid20315-seed3-e3200-i600-subtr_mod_inh0-i_GABA_99']
file_list['i_GABA'][3]['modinh3'] = ['output060920161649-pid7706-seed3-e3200-i600-subtr_mod_inh3-i_GABA_110',
 'output060820161303-pid35861-seed3-e3200-i600-subtr_mod_inh3-i_GABA_111',
 'output060820161303-pid35879-seed3-e3200-i600-subtr_mod_inh3-i_GABA_112',
 'output060820161304-pid28406-seed3-e3200-i600-subtr_mod_inh3-i_GABA_113',
 'output060820161304-pid28403-seed3-e3200-i600-subtr_mod_inh3-i_GABA_114',
 'output060820161304-pid28405-seed3-e3200-i600-subtr_mod_inh3-i_GABA_115',
 'output060820161304-pid28404-seed3-e3200-i600-subtr_mod_inh3-i_GABA_116',
 'output060820161304-pid28400-seed3-e3200-i600-subtr_mod_inh3-i_GABA_117',
 'output060820161304-pid28408-seed3-e3200-i600-subtr_mod_inh3-i_GABA_118',
 'output060820161304-pid28409-seed3-e3200-i600-subtr_mod_inh3-i_GABA_119']
file_list['i_GABA'][4]['modinh0'] = ['output060920161649-pid112082-seed4-e3200-i600-subtr_mod_inh0-i_GABA_120',
 'output060820161304-pid23831-seed4-e3200-i600-subtr_mod_inh0-i_GABA_121',
 'output060820161304-pid23834-seed4-e3200-i600-subtr_mod_inh0-i_GABA_122',
 'output060820161304-pid23833-seed4-e3200-i600-subtr_mod_inh0-i_GABA_123',
 'output060920161649-pid112085-seed4-e3200-i600-subtr_mod_inh0-i_GABA_124',
 'output060820161303-pid40626-seed4-e3200-i600-subtr_mod_inh0-i_GABA_125',
 'output060820161303-pid40638-seed4-e3200-i600-subtr_mod_inh0-i_GABA_126',
 'output060820161303-pid40634-seed4-e3200-i600-subtr_mod_inh0-i_GABA_127',
 'output060820161303-pid40633-seed4-e3200-i600-subtr_mod_inh0-i_GABA_128',
 'output060820161303-pid40636-seed4-e3200-i600-subtr_mod_inh0-i_GABA_129']
file_list['i_GABA'][4]['modinh3'] = ['output060820161304-pid28402-seed4-e3200-i600-subtr_mod_inh3-i_GABA_140',
 'output060820161304-pid47116-seed4-e3200-i600-subtr_mod_inh3-i_GABA_141',
 'output060820161304-pid47118-seed4-e3200-i600-subtr_mod_inh3-i_GABA_142',
 'output060820161304-pid47121-seed4-e3200-i600-subtr_mod_inh3-i_GABA_143',
 'output060820161304-pid47125-seed4-e3200-i600-subtr_mod_inh3-i_GABA_144',
 'output060820161304-pid47120-seed4-e3200-i600-subtr_mod_inh3-i_GABA_145',
 'output060820161304-pid47122-seed4-e3200-i600-subtr_mod_inh3-i_GABA_146',
 'output060820161304-pid47115-seed4-e3200-i600-subtr_mod_inh3-i_GABA_147',
 'output060820161304-pid47124-seed4-e3200-i600-subtr_mod_inh3-i_GABA_148',
 'output060820161304-pid47127-seed4-e3200-i600-subtr_mod_inh3-i_GABA_149']

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