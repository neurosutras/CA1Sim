import btmorph
import numpy
import matplotlib.pyplot as plt

#--------analyzing and visualizing morphology--------

swc_tree = btmorph.STree2()

morph_dir = 'morphologies/'
cell_path = 'EB022715-stitched-proofread.swc'
#cell_path = 'Erik_Bloss_CA1_0215_Stitched_Proofread.swc'

swc_tree.read_SWC_tree_from_file(morph_dir+cell_path)
stats = btmorph.BTStats(swc_tree)

# get the total length
total_length = stats.total_length()
print 'total_length = %f' % total_length

# get the max degree
max_degree = stats.degree_of_node(swc_tree.get_root())
print 'max_degree = %f' % max_degree

#generate dendrogram
btmorph.plot_dendrogram(morph_dir+cell_path)

#generate 2D_projection
btmorph.plot_2D_SWC(morph_dir+cell_path, show_axis=False, color_scheme='default', depth='Y')

#generate population_2D_density_projection
#btmorph.population_density_projection(destination="morphologies",filter='*.swc', precision=[10, 10, 10],depth='Y')

plt.show()