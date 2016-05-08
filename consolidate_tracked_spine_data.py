__author__ = 'milsteina'
from function_lib import *

target_file_list = ['output050720162252-pid59164-seed3-e3200-i500-mod_inh0-track_spines_90.0_0',
 'output050720162252-pid30124-seed3-e3200-i500-mod_inh0-track_spines_91.0_0',
 'output050720162252-pid42149-seed3-e3200-i500-mod_inh0-track_spines_92.0_0',
 'output050720162252-pid30068-seed3-e3200-i500-mod_inh0-track_spines_93.0_0',
 'output050720162252-pid30067-seed3-e3200-i500-mod_inh0-track_spines_94.0_0',
 'output050720162252-pid30066-seed3-e3200-i500-mod_inh0-track_spines_95.0_0',
 'output050720162252-pid25271-seed3-e3200-i500-mod_inh0-track_spines_96.0_0',
 'output050720162252-pid3834-seed3-e3200-i500-mod_inh0-track_spines_97.0_0',
 'output050720162252-pid2616-seed3-e3200-i500-mod_inh0-track_spines_98.0_0',
 'output050720162252-pid30138-seed3-e3200-i500-mod_inh0-track_spines_99.0_0',
 'output050720162252-pid9794-seed3-e3200-i500-mod_inh1-track_spines_100.0_0',
 'output050720162252-pid18707-seed3-e3200-i500-mod_inh1-track_spines_101.0_0',
 'output050720162252-pid18708-seed3-e3200-i500-mod_inh1-track_spines_102.0_0',
 'output050720162252-pid81846-seed3-e3200-i500-mod_inh1-track_spines_103.0_0',
 'output050720162252-pid42914-seed3-e3200-i500-mod_inh1-track_spines_104.0_0',
 'output050720162252-pid42912-seed3-e3200-i500-mod_inh1-track_spines_105.0_0',
 'output050720162252-pid42913-seed3-e3200-i500-mod_inh1-track_spines_106.0_0',
 'output050720162252-pid13559-seed3-e3200-i500-mod_inh1-track_spines_107.0_0',
 'output050720162252-pid13558-seed3-e3200-i500-mod_inh1-track_spines_108.0_0',
 'output050720162252-pid36651-seed3-e3200-i500-mod_inh1-track_spines_109.0_0',
 'output050720162252-pid36649-seed3-e3200-i500-mod_inh2-track_spines_110.0_0',
 'output050720162252-pid36650-seed3-e3200-i500-mod_inh2-track_spines_111.0_0',
 'output050720162252-pid30082-seed3-e3200-i500-mod_inh2-track_spines_112.0_0',
 'output050720162252-pid25285-seed3-e3200-i500-mod_inh2-track_spines_113.0_0',
 'output050720162252-pid13407-seed3-e3200-i500-mod_inh2-track_spines_114.0_0',
 'output050720162252-pid23159-seed3-e3200-i500-mod_inh2-track_spines_115.0_0',
 'output050720162252-pid48178-seed3-e3200-i500-mod_inh2-track_spines_116.0_0',
 'output050720162252-pid118249-seed3-e3200-i500-mod_inh2-track_spines_117.0_0',
 'output050720162252-pid118248-seed3-e3200-i500-mod_inh2-track_spines_118.0_0',
 'output050720162252-pid39547-seed3-e3200-i500-mod_inh2-track_spines_119.0_0']

i_AMPA_files = ['output050720162233-pid13641-seed3-e3200-i500-mod_inh0-i_AMPA_90.0_0',
 'output050720162234-pid80106-seed3-e3200-i500-mod_inh0-i_AMPA_91.0_0',
 'output050720162234-pid871-seed3-e3200-i500-mod_inh0-i_AMPA_92.0_0',
 'output050720162234-pid27563-seed3-e3200-i500-mod_inh0-i_AMPA_93.0_0',
 'output050720162234-pid14335-seed3-e3200-i500-mod_inh0-i_AMPA_94.0_0',
 'output050720162234-pid23807-seed3-e3200-i500-mod_inh0-i_AMPA_95.0_0',
 'output050720162234-pid143937-seed3-e3200-i500-mod_inh0-i_AMPA_96.0_0',
 'output050720162235-pid5830-seed3-e3200-i500-mod_inh0-i_AMPA_97.0_0',
 'output050720162235-pid5851-seed3-e3200-i500-mod_inh0-i_AMPA_98.0_0',
 'output050720162235-pid11104-seed3-e3200-i500-mod_inh0-i_AMPA_99.0_0',
 'output050720162235-pid14375-seed3-e3200-i500-mod_inh1-i_AMPA_100.0_0',
 'output050720162235-pid18300-seed3-e3200-i500-mod_inh1-i_AMPA_101.0_0',
 'output050720162235-pid144006-seed3-e3200-i500-mod_inh1-i_AMPA_102.0_0',
 'output050720162236-pid37527-seed3-e3200-i500-mod_inh1-i_AMPA_103.0_0',
 'output050720162236-pid34105-seed3-e3200-i500-mod_inh1-i_AMPA_104.0_0',
 'output050720162236-pid16409-seed3-e3200-i500-mod_inh1-i_AMPA_105.0_0',
 'output050720162236-pid1832-seed3-e3200-i500-mod_inh1-i_AMPA_106.0_0',
 'output050720162236-pid20065-seed3-e3200-i500-mod_inh1-i_AMPA_107.0_0',
 'output050720162236-pid30902-seed3-e3200-i500-mod_inh1-i_AMPA_108.0_0',
 'output050720162236-pid27609-seed3-e3200-i500-mod_inh1-i_AMPA_109.0_0',
 'output050720162237-pid131009-seed3-e3200-i500-mod_inh2-i_AMPA_110.0_0',
 'output050720162237-pid80257-seed3-e3200-i500-mod_inh2-i_AMPA_111.0_0',
 'output050720162238-pid80288-seed3-e3200-i500-mod_inh2-i_AMPA_112.0_0',
 'output050720162238-pid12016-seed3-e3200-i500-mod_inh2-i_AMPA_113.0_0',
 'output050720162238-pid71134-seed3-e3200-i500-mod_inh2-i_AMPA_114.0_0',
 'output050720162238-pid42795-seed3-e3200-i500-mod_inh2-i_AMPA_115.0_0',
 'output050720162238-pid90822-seed3-e3200-i500-mod_inh2-i_AMPA_116.0_0',
 'output050720162238-pid131111-seed3-e3200-i500-mod_inh2-i_AMPA_117.0_0',
 'output050720162238-pid13835-seed3-e3200-i500-mod_inh2-i_AMPA_118.0_0',
 'output050720162238-pid145099-seed3-e3200-i500-mod_inh2-i_AMPA_119.0_0']

i_NMDA_files = ['output050720162235-pid98226-seed3-e3200-i500-mod_inh0-i_NMDA_90.0_0',
 'output050720162235-pid28269-seed3-e3200-i500-mod_inh0-i_NMDA_91.0_0',
 'output050720162235-pid22311-seed3-e3200-i500-mod_inh0-i_NMDA_92.0_0',
 'output050720162235-pid11591-seed3-e3200-i500-mod_inh0-i_NMDA_93.0_0',
 'output050720162235-pid90644-seed3-e3200-i500-mod_inh0-i_NMDA_94.0_0',
 'output050720162235-pid5051-seed3-e3200-i500-mod_inh0-i_NMDA_95.0_0',
 'output050720162235-pid27250-seed3-e3200-i500-mod_inh0-i_NMDA_96.0_0',
 'output050720162235-pid30521-seed3-e3200-i500-mod_inh0-i_NMDA_97.0_0',
 'output050720162235-pid42612-seed3-e3200-i500-mod_inh0-i_NMDA_98.0_0',
 'output050720162235-pid19283-seed3-e3200-i500-mod_inh0-i_NMDA_99.0_0',
 'output050720162236-pid22350-seed3-e3200-i500-mod_inh1-i_NMDA_100.0_0',
 'output050720162236-pid47776-seed3-e3200-i500-mod_inh1-i_NMDA_101.0_0',
 'output050720162236-pid70764-seed3-e3200-i500-mod_inh1-i_NMDA_102.0_0',
 'output050720162236-pid19326-seed3-e3200-i500-mod_inh1-i_NMDA_103.0_0',
 'output050720162236-pid129615-seed3-e3200-i500-mod_inh1-i_NMDA_104.0_0',
 'output050720162236-pid114782-seed3-e3200-i500-mod_inh1-i_NMDA_105.0_0',
 'output050720162236-pid23276-seed3-e3200-i500-mod_inh1-i_NMDA_106.0_0',
 'output050720162236-pid27228-seed3-e3200-i500-mod_inh1-i_NMDA_107.0_0',
 'output050720162237-pid20207-seed3-e3200-i500-mod_inh1-i_NMDA_108.0_0',
 'output050720162237-pid9405-seed3-e3200-i500-mod_inh1-i_NMDA_109.0_0',
 'output050720162238-pid35210-seed3-e3200-i500-mod_inh2-i_NMDA_110.0_0',
 'output050720162238-pid80324-seed3-e3200-i500-mod_inh2-i_NMDA_111.0_0',
 'output050720162238-pid10870-seed3-e3200-i500-mod_inh2-i_NMDA_112.0_0',
 'output050720162238-pid112716-seed3-e3200-i500-mod_inh2-i_NMDA_113.0_0',
 'output050720162238-pid20453-seed3-e3200-i500-mod_inh2-i_NMDA_114.0_0',
 'output050720162239-pid23333-seed3-e3200-i500-mod_inh2-i_NMDA_115.0_0',
 'output050720162239-pid10964-seed3-e3200-i500-mod_inh2-i_NMDA_116.0_0',
 'output050720162239-pid36811-seed3-e3200-i500-mod_inh2-i_NMDA_117.0_0',
 'output050720162239-pid96605-seed3-e3200-i500-mod_inh2-i_NMDA_118.0_0',
 'output050720162239-pid1900-seed3-e3200-i500-mod_inh2-i_NMDA_119.0_0']

index_list = []
key_list = []
target_filename = target_file_list[0]
with h5py.File(data_dir+target_filename+'.hdf5', 'r') as target_file:
    target_trial = target_file.itervalues().next()
    for target_rec in target_trial['rec'].itervalues():
        if 'description' in target_rec.attrs and target_rec.attrs['description'] == 'spine_vm':
            index_list.append(target_rec.attrs['index'])
with h5py.File(data_dir+i_AMPA_files[i]+'.hdf5', 'r') as source_file:
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
