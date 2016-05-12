__author__ = 'milsteina'
from function_lib import *

target_file_list = ['output051120161734-pid45981-seed3-e3200-i500-mod_inh0-track_spines_2.5_90',
 'output051120161734-pid45996-seed3-e3200-i500-mod_inh0-track_spines_2.5_91',
 'output051120161734-pid45994-seed3-e3200-i500-mod_inh0-track_spines_2.5_92',
 'output051120161734-pid101565-seed3-e3200-i500-mod_inh0-track_spines_2.5_93',
 'output051120161734-pid17557-seed3-e3200-i500-mod_inh0-track_spines_2.5_94',
 'output051120161734-pid17569-seed3-e3200-i500-mod_inh0-track_spines_2.5_95',
 'output051120161734-pid11739-seed3-e3200-i500-mod_inh0-track_spines_2.5_96',
 'output051120161734-pid11750-seed3-e3200-i500-mod_inh0-track_spines_2.5_97',
 'output051120161734-pid32298-seed3-e3200-i500-mod_inh0-track_spines_2.5_98',
 'output051120161734-pid43234-seed3-e3200-i500-mod_inh0-track_spines_2.5_99',
 'output051120161734-pid129036-seed3-e3200-i500-mod_inh1-track_spines_2.5_100',
 'output051120161734-pid1542-seed3-e3200-i500-mod_inh1-track_spines_2.5_101',
 'output051120161734-pid14536-seed3-e3200-i500-mod_inh1-track_spines_2.5_102',
 'output051120161735-pid46556-seed3-e3200-i500-mod_inh1-track_spines_2.5_103',
 'output051120161734-pid46495-seed3-e3200-i500-mod_inh1-track_spines_2.5_104',
 'output051120161734-pid1330-seed3-e3200-i500-mod_inh1-track_spines_2.5_105',
 'output051120161734-pid1326-seed3-e3200-i500-mod_inh1-track_spines_2.5_106',
 'output051120161734-pid1329-seed3-e3200-i500-mod_inh1-track_spines_2.5_107',
 'output051120161734-pid17608-seed3-e3200-i500-mod_inh1-track_spines_2.5_108',
 'output051120161734-pid41419-seed3-e3200-i500-mod_inh1-track_spines_2.5_109',
 'output051120161734-pid66514-seed3-e3200-i500-mod_inh2-track_spines_2.5_110',
 'output051120161734-pid32344-seed3-e3200-i500-mod_inh2-track_spines_2.5_111',
 'output051120161734-pid19638-seed3-e3200-i500-mod_inh2-track_spines_2.5_112',
 'output051120161734-pid46521-seed3-e3200-i500-mod_inh2-track_spines_2.5_113',
 'output051120161735-pid35835-seed3-e3200-i500-mod_inh2-track_spines_2.5_114',
 'output051120161735-pid35834-seed3-e3200-i500-mod_inh2-track_spines_2.5_115',
 'output051120161735-pid29594-seed3-e3200-i500-mod_inh2-track_spines_2.5_116',
 'output051120161735-pid29595-seed3-e3200-i500-mod_inh2-track_spines_2.5_117',
 'output051120161735-pid29598-seed3-e3200-i500-mod_inh2-track_spines_2.5_118',
 'output051120161735-pid37414-seed3-e3200-i500-mod_inh2-track_spines_2.5_119']

i_AMPA_files = ['output051120161751-pid46910-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_90',
 'output051120161751-pid46907-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_91',
 'output051120161819-pid115964-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_92',
 'output051120161819-pid115957-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_93',
 'output051120161819-pid115961-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_94',
 'output051120161819-pid12237-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_95',
 'output051120161931-pid40712-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_96',
 'output051120161932-pid23562-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_97',
 'output051120161933-pid14193-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_98',
 'output051120161933-pid44423-seed3-e3200-i500-mod_inh0-i_AMPA_2.5_99',
 'output051120161934-pid23372-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_100',
 'output051120161934-pid47736-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_101',
 'output051120161934-pid28605-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_102',
 'output051120161934-pid10044-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_103',
 'output051120161934-pid44262-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_104',
 'output051120161934-pid33775-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_105',
 'output051120161934-pid44467-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_106',
 'output051120161935-pid21165-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_107',
 'output051120161936-pid18911-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_108',
 'output051120161936-pid16108-seed3-e3200-i500-mod_inh1-i_AMPA_2.5_109',
 'output051120161940-pid44407-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_110',
 'output051120161940-pid31191-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_111',
 'output051120161940-pid34277-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_112',
 'output051120161940-pid34277-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_113',
 'output051120161940-pid23581-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_114',
 'output051120161940-pid10663-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_115',
 'output051120161940-pid14437-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_116',
 'output051120161940-pid44598-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_117',
 'output051120161941-pid45082-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_118',
 'output051120161941-pid34716-seed3-e3200-i500-mod_inh2-i_AMPA_2.5_119']

i_NMDA_files = ['output051120161933-pid14228-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_90',
 'output051120161933-pid22862-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_91',
 'output051120161933-pid41477-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_92',
 'output051120161933-pid28484-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_93',
 'output051120161933-pid11674-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_94',
 'output051120161933-pid40641-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_95',
 'output051120161934-pid18119-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_96',
 'output051120161934-pid31062-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_97',
 'output051120161934-pid42956-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_98',
 'output051120161934-pid29421-seed3-e3200-i500-mod_inh0-i_NMDA_2.5_99',
 'output051120161936-pid33368-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_100',
 'output051120161936-pid36297-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_101',
 'output051120161936-pid15330-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_102',
 'output051120161937-pid11508-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_103',
 'output051120161937-pid22929-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_104',
 'output051120161938-pid34221-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_105',
 'output051120161939-pid23426-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_106',
 'output051120161939-pid23474-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_107',
 'output051120161939-pid67043-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_108',
 'output051120161940-pid28825-seed3-e3200-i500-mod_inh1-i_NMDA_2.5_109',
 'output051120161941-pid28603-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_110',
 'output051120161941-pid44449-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_111',
 'output051120161941-pid27039-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_112',
 'output051120161942-pid37966-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_113',
 'output051120161942-pid44890-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_114',
 'output051120161942-pid24988-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_115',
 'output051120161943-pid23250-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_116',
 'output051120161944-pid44622-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_117',
 'output051120162107-pid10070-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_118',
 'output051120162107-pid24992-seed3-e3200-i500-mod_inh2-i_NMDA_2.5_119']

index_list = []
key_list = []
target_filename = target_file_list[0]
with h5py.File(data_dir+target_filename+'.hdf5', 'r') as target_file:
    target_trial = target_file.itervalues().next()
    for target_rec in target_trial['rec'].itervalues():
        if 'description' in target_rec.attrs and target_rec.attrs['description'] == 'spine_vm':
            index_list.append(target_rec.attrs['index'])
with h5py.File(data_dir+i_AMPA_files[0]+'.hdf5', 'r') as source_file:
    source_trial = source_file.itervalues().next()
    for source_key in source_trial['rec']:
        if ('index' in source_trial['rec'][source_key].attrs) and \
                (source_trial['rec'][source_key].attrs['index'] in index_list):
            key_list.append(source_key)
for i, target_filename in enumerate(target_file_list):
    with h5py.File(data_dir + target_filename + '.hdf5', 'a') as target_file:
        for trial_key in target_file:
            for source_file_list in i_AMPA_files, i_NMDA_files:
                with h5py.File(data_dir+source_file_list[i]+'.hdf5', 'r') as source_file:
                    target_key_int = len(target_file[trial_key]['rec'])
                    for source_key in key_list:
                        target_file.copy(source_file[trial_key]['rec'][source_key], target_file[trial_key]['rec'],
                                         str(target_key_int))
                        target_key_int += 1
