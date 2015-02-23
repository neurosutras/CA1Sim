import btmorph
import numpy
import matplotlib.pyplot as plt

#--------analyzing and visualizing morphology--------

swc_tree= btmorph.STree2()

#cell_path="morphologies/020910_B09Lpy_JK.swc" #from Jinny but has no tuft
cell_path="morphologies/MC120914100xC3-tweaked.swc"
#cell_path="/Users/milsteina/neuron/CA3project/SWCfiles/DH070213C3-.Edit.scaled.swc.swc"

swc_tree.read_SWC_tree_from_file(cell_path)
stats = btmorph.BTStats(swc_tree)

# get the total length
total_length = stats.total_length()
print 'total_length = %f' % total_length

# get the max degree
max_degree = stats.degree_of_node(swc_tree.get_root())
print 'max_degree = %f' % max_degree

#generate dendrogram
btmorph.plot_dendrogram(cell_path)

#generate 2D_projection
btmorph.plot_2D_SWC(cell_path,show_axis=False,color_scheme='default',depth='Y')

#generate population_2D_density_projection
#btmorph.population_density_projection(destination="morphologies",filter='*.swc', precision=[10, 10, 10],depth='Y')

plt.show()