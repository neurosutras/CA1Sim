#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_30 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 30'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_31 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 31'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_32 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 32'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_33 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 33'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_34 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 34'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_35 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 35'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_36 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 36'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_37 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 37'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_38 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 38'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh0_39 -m e -M milsteina@janelia.hhmi.org -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 0 0.8 39'

qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_60 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 60'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_61 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 61'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_62 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 62'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_63 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 63'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_64 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 64'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_65 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 65'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_66 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 66'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_67 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 67'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_68 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 68'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh0_69 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 0 0.8 69'

qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_90 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 90'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_91 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 91'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_92 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 92'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_93 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 93'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_94 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 94'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_95 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 95'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_96 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 96'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_97 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 97'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_98 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 98'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh0_99 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 0 0.8 99'

qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_120 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 120'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_121 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 121'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_122 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 122'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_123 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 123'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_124 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 124'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_125 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 125'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_126 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 126'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_127 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 127'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_128 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 128'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh0_129 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 0 0.8 129'

qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_20 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 20'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_21 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 21'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_22 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 22'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_23 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 23'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_24 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 24'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_25 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 25'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_26 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 26'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_27 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 27'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_28 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 28'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_inh3_29 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 0 3200 600 3 0.8 29'

qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_50 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 50'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_51 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 51'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_52 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 52'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_53 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 53'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_54 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 54'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_55 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 55'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_56 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 56'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_57 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 57'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_58 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 58'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_inh3_59 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 1 3200 600 3 0.8 59'

qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_80 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 80'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_81 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 81'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_82 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 82'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_83 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 83'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_84 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 84'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_85 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 85'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_86 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 86'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_87 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 87'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_88 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 88'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_inh3_89 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 2 3200 600 3 0.8 89'

qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_110 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 110'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_111 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 111'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_112 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 112'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_113 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 113'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_114 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 114'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_115 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 115'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_116 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 116'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_117 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 117'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_118 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 118'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_inh3_119 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 3 3200 600 3 0.8 119'

qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_140 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 140'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_141 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 141'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_142 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 142'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_143 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 143'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_144 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 144'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_145 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 145'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_146 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 146'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_147 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 147'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_148 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 148'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_inh3_149 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise.py 4 3200 600 3 0.8 149'

qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_0 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 0'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_1 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 1'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_2 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 2'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_3 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 3'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_4 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 4'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_5 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 5'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_6 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 6'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_7 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 7'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_8 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 8'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh0_9 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 0 0.8 9'

qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_30 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 30'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_31 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 31'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_32 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 32'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_33 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 33'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_34 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 34'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_35 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 35'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_36 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 36'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_37 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 37'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_38 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 38'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh0_39 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 0 0.8 39'

qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_60 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 60'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_61 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 61'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_62 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 62'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_63 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 63'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_64 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 64'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_65 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 65'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_66 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 66'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_67 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 67'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_68 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 68'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh0_69 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 0 0.8 69'

qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_90 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 90'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_91 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 91'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_92 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 92'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_93 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 93'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_94 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 94'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_95 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 95'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_96 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 96'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_97 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 97'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_98 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 98'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh0_99 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 0 0.8 99'

qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_120 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 120'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_121 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 121'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_122 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 122'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_123 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 123'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_124 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 124'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_125 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 125'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_126 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 126'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_127 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 127'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_128 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 128'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh0_129 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 0 0.8 129'

qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_20 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 20'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_21 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 21'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_22 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 22'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_23 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 23'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_24 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 24'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_25 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 25'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_26 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 26'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_27 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 27'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_28 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 28'
qsub -N job_092816_e3200_i600_subtr0_noise0.8_no_na_inh3_29 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 0 3200 600 3 0.8 29'

qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_50 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 50'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_51 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 51'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_52 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 52'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_53 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 53'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_54 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 54'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_55 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 55'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_56 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 56'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_57 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 57'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_58 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 58'
qsub -N job_092816_e3200_i600_subtr1_noise0.8_no_na_inh3_59 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 1 3200 600 3 0.8 59'

qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_80 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 80'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_81 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 81'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_82 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 82'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_83 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 83'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_84 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 84'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_85 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 85'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_86 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 86'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_87 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 87'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_88 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 88'
qsub -N job_092816_e3200_i600_subtr2_noise0.8_no_na_inh3_89 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 2 3200 600 3 0.8 89'

qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_110 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 110'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_111 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 111'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_112 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 112'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_113 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 113'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_114 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 114'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_115 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 115'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_116 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 116'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_117 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 117'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_118 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 118'
qsub -N job_092816_e3200_i600_subtr3_noise0.8_no_na_inh3_119 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 3 3200 600 3 0.8 119'

qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_140 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 140'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_141 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 141'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_142 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 142'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_143 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 143'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_144 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 144'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_145 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 145'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_146 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 146'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_147 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 147'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_148 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 148'
qsub -N job_092816_e3200_i600_subtr4_noise0.8_no_na_inh3_149 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh_add_noise_no_na.py 4 3200 600 3 0.8 149'