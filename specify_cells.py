__author__ = 'Aaron D. Milstein'
# Includes modification of an early version of SWC_neuron.py by Daniele Linaro
import copy
import pprint
import scipy as sp
import btmorph
from neuron import h
from function_lib import *


# SWC files must use this nonstandard convention to exploit trunk and tuft categorization
swc_types = [soma_type, axon_type, basal_type, apical_type, trunk_type, tuft_type] = [1, 2, 3, 4, 5, 6]
sec_types = ['soma', 'axon', 'ais', 'basal', 'apical', 'trunk', 'tuft', 'spine_neck', 'spine_head']

verbose = 0  # Turn on for text reporting during model initialization and simulation
gid_max = 0  # Every new HocCell will receive a global identifier for network simulation, and increment this counter

#-------Wrapper for converting SWC --> BtMorph --> Skeleton morphology in NEURON hoc------------------------


class HocCell(object):
    def __init__(self, morph_filename=None, mech_filename=None):
        """
        :param morph_filename: str : path to .swc file containing morphology
        :param mech_filename: str : path to .pkl file specifying cable parameters and membrane mechanisms
        """
        global gid_max
        self._gid = gid_max
        gid_max += 1
        self.tree = btmorph.STree2()  # Builds a simple tree to store nodes of type 'SHocNode'
        self.index = 0  # Keep track of number of nodes
        self._node_dict = {'soma': [], 'axon': [], 'basal': [], 'apical': [], 'trunk': [], 'tuft': [], 'spine': []}
        self.mech_dict = self.load_mech_dict(mech_filename)  # Refer to function_lib for description of structure of
                                                             # mechanism dictionary. loads from .pkl or
                                                             # default_mech_dict in function_lib
        if not morph_filename is None:
            self.load_morphology_from_swc(morph_filename)
            self.reinit_mechanisms()  # Membrane mechanisms must be reinitialized whenever cable properties (Ra, cm) or
                                      # spatial resolution (nseg) changes.

    def load_morphology_from_swc(self, morph_filename):
        """
        This method builds an STree2 comprised of SHocNode nodes associated with hoc sections, connects the hoc
        sections, and initializes various parameters: Ra, cm, L, diam, nseg
        This method implements a standardized soma and axon:
        The soma consists of two cylindrical hoc sections of equal length and diameter, connected (0) to (0).
        The basal dendritic tree is connected to soma[1](1), and the apical tree is connected to soma[0](1).
        The axon is attached to soma[0](0), and consists of three sequential hoc sections:
            1) axon[0] : a tapered cylindrical 'axon hillock' section connected to soma[0](0)
            2) axon[1] : a tapered cylindrical 'axon initial segment' section connected to axon[0](1)
            3) axon[2] : a cylindrical 'axon' section connected to axon[1](1)
        """
        raw_tree = btmorph.STree2() # import the full tree from an SWC file
        raw_tree.read_SWC_tree_from_file(morph_dir+morph_filename, types=range(10))
        soma_length = 14
        soma_diam = 9
        for index in range(2):
            node = self.make_section('soma')
            node.sec.L = soma_length/2
            node.sec.diam = soma_diam
            self._init_cable(node)  # consults the mech_dict to initialize Ra, cm, and nseg
        self.tree.root = self.soma[0]
        self.soma[1].connect(self.soma[0], 0, 0)
        for index in range(3):
            self.make_section('axon')
        self.axon[0].sec.L = 20
        self.axon[0].set_diam_bounds(3, 2)  # stores the diameter boundaries for a tapered cylindrical section
        self.axon[1].type = 'ais'
        self.axon[1].sec.L = 25
        self.axon[1].set_diam_bounds(2, 1)
        self.axon[2].sec.L = 500
        self.axon[2].sec.diam = 1
        self.axon[0].connect(self.soma[0], 0, 0)
        self.axon[1].connect(self.axon[0], 1, 0)
        self.axon[2].connect(self.axon[1], 1, 0)
        for node in self.axon:
            self._init_cable(node)
        for child in raw_tree.root.children:
            self.make_skeleton(child, self.tree.root)

    def make_section(self, sec_type):
        """
        Create a new hoc section to associate with this node, and this cell, and store information about it in the
        node's content dictionary.
        :param sec_type: str
        :return node: :class:'SHocNode'
        """
        node = SHocNode(self.index)
        if self.index == 0:
            self.tree.root = node
        self.index += 1
        node.type = sec_type
        if sec_type in ['spine_head', 'spine_neck']:
            self._node_dict['spine'].append(node)
        else:
            self._node_dict[sec_type].append(node)
        node.sec = h.Section(name=node.name, cell=self)
        return node

    def test_sec_properties(self, node=None):
        """
        Used for debugging and validating model specification.
        :param node:
        """
        if node is None:
            node = self.tree.root
#        node.sec.push()
#        h.psection()
#        h.pop_section()
        print '{} [nseg: {}, Ra: {}]'.format(node.name, node.sec.nseg, node.sec.Ra)
#        h('for (x) print (x), diam(x)', sec=node.sec)
        for child in node.children:
            self.test_sec_properties(child)

    def make_skeleton(self, raw_node, parent, length=0, diams=None):
        """
        Following construction of soma and axon nodes of type 'SHocNode' in the tree of type 'STree2', this method
        recursively converts dendritic 'SNode2' nodes into 'SHocNode' nodes, and connects them to the appropriate
        somatic nodes. Skeletonized dendritic nodes have only one hoc section for each stretch of unbranched dendrite,
        with length equal to the sum of the lengths of the converted SNode2 nodes.
        Nodes that taper more than 1 um remain tapered, otherwise they are converted into untapered cylinders with
        diameter equal to the mean diameter of the the converted SNode2 nodes.
        Dendrite types that are pre-categorized as basal, apical, trunk, or tuft in the input .swc file are preserved.
        :param raw_node: :class:'SNode2'
        :param parent: :class:'SHocNode'
        :param length: int or float
        :param diams: None or (list: float)
        """
        dend_types=([basal_type, apical_type, trunk_type, tuft_type], ['basal', 'apical', 'trunk', 'tuft'])
        swc = raw_node.get_content()['p3d']
        swc_type = swc.type
        if swc_type in dend_types[0]:
            diam = swc.radius*2
            length += self.get_node_length_swc(raw_node)
            leaves = len(raw_node.children)
            """
            print 'raw node {}: [type:{}, length:{}, diam:{}] has {} children, parent is {}'.format(raw_node.index,
                                                        swc_type, self.get_node_length_swc(raw_node), diam, leaves,
                                                        raw_node.parent.index)
            """
            if (leaves == 0) | (leaves > 1):
                sec_type = dend_types[1][dend_types[0].index(swc_type)]
                new_node = self.make_section(sec_type)
                new_node.sec.L = length
                if (self.tree.is_root(parent)) & (sec_type == 'basal'):
                    parent = self.soma[1]
                new_node.connect(parent)
                if diams is None:
                    new_node.sec.diam = diam
                    self._init_cable(new_node)
                    if verbose:
                        print '{} [nseg: {}, diam: {}, length: {}, parent: {}]'.format(new_node.name, new_node.sec.nseg,
                                                                                    diam, length, new_node.parent.name)
                else:
                    diams.append(diam)
                    if len(diams) > 2:
                        mean = sp.mean(diams)
                        stdev = sp.std(diams)
                        if stdev*2 > 1:  # If 95% of the values are within 1 um, don't taper
                            new_node.set_diam_bounds(mean+stdev, mean-stdev)
                            self._init_cable(new_node)
                            if verbose:
                                print '{} [nseg: {}, diam: ({}:{}), length: {}, parent: {}]'.format(new_node.name,
                                                new_node.sec.nseg, mean+stdev, mean-stdev, length, new_node.parent.name)
                        else:
                            new_node.sec.diam = mean
                            self._init_cable(new_node)
                            if verbose:
                                print '{} [nseg: {}, diam: {}, length: {}, parent: {}]'.format(new_node.name,
                                                                new_node.sec.nseg, mean, length, new_node.parent.name)
                    elif abs(diams[0]-diams[1]) > 1:
                        new_node.set_diam_bounds(diams[0], diams[1])
                        self._init_cable(new_node)
                        if verbose:
                            print '{} [diam: ({}:{}), length: {}, parent: {}]'.format(new_node.name, new_node.sec.nseg,
                                                                    diams[0], diams[1], length, new_node.parent.name)
                    else:
                        mean = sp.mean(diams)
                        new_node.sec.diam = mean
                        self._init_cable(new_node)
                        if verbose:
                            print '{} [nseg: {}, diam: {}, length: {}, parent: {}]'.format(new_node.name,
                                                                new_node.sec.nseg, mean, length, new_node.parent.name)
                if leaves > 1:  # Follow all branches from this fork
                    for child in raw_node.children:
                        self.make_skeleton(child, new_node)
            else:  # Follow unbranched dendrite
                if diams is None:
                    diams = [diam]
                else:
                    diams.append(diam)
                self.make_skeleton(raw_node.children[0], parent, length, diams)

    def get_nodes_of_subtype(self, sec_type):
        """
        This method searches the node dictionary for nodes of a given sec_type and returns them in a list. Used during
        specification of membrane mechanisms.
        :param sec_type: str
        :return: list of :class:'SHocNode'
        """
        if sec_type == 'ais':
            return [node for node in self.axon if node.type == sec_type]
        elif sec_type == 'axon':
            return [node for node in self.axon if not node.type == 'ais']
        elif sec_type in ['spine_head', 'spine_neck']:
            return [node for node in self.spine if node.type == sec_type]
        else:
            return self._node_dict[sec_type]

    def load_mech_dict(self, mech_filename=None):
        """
        This method loads the dictionary specifying membrane mechanism parameters. If a .pkl file is not provided, a
        global variable default_mech_dict from function_lib is used.
        :param mech_filename: str
        """
        if not mech_filename is None:
            return read_from_pkl(data_dir+mech_filename)
        else:
            local_mech_dict = copy.deepcopy(default_mech_dict)
            return local_mech_dict

    def _init_cable(self, node):
        """
        If the mechanism dictionary specifies the cable properties 'Ra' or 'cm', then _modify_mechanism() properly sets
        those parameters, and reinitializes the number of segments per section. To avoid redundancy, this
        method passes _modify_mechanism() a copy of the dictionary with the spatial_res parameter removed, since this is
        consulted in setting nseg. However, if spatial_res is the only parameter being specified, it is passed to
        _modify_mechanism()
        :param node: :class:'SHocNode'
        """
        sec_type = node.type
        if sec_type in self.mech_dict and 'cable' in self.mech_dict[sec_type]:
            mech_content = copy.deepcopy(self.mech_dict[sec_type]['cable'])
            if ('Ra' in mech_content) | ('cm' in mech_content):
                if 'spatial_res' in mech_content:
                    del mech_content['spatial_res']
                self._modify_mechanism(node, 'cable', mech_content)
            elif 'spatial_res' in mech_content:
                self._modify_mechanism(node, 'cable', mech_content)
        else:
            node.init_nseg()
            node.reinit_diam()

    def reinit_mechanisms(self, reset_cable=0):
        """
        Once a mechanism dictionary has been loaded, and a morphology has been specified, this method traverses through
        the tree of SHocNode nodes breadth-first and properly sets membrane mechanism parameters, including gradients
        and inheritance of parameters from nodes along the path from root. Since cable parameters are set during
        specification of morphology, it is not necessary to immediately reinitialize these parameters again. However,
        they can be manually reinitialized with the reset_cable flag.
        :param reset_cable: boolean
        """
        nodes = [self.tree.root]
        for node in nodes:
            nodes.extend(node.children)
        self._reinit_mech(nodes, reset_cable)

    def _reinit_mech(self, nodes, reset_cable=0):
        """
        Given a list of nodes, this method loops through all the mechanisms specified in the mechanism dictionary for
        the hoc section type of each node and updates their associated parameters. If the reset_cable flag is set to 1,
        cable parameters are modified first, then the parameters for all other mechanisms are reinitialized.
        :param nodes: list of :class:'SHocNode'
        :param reset_cable: boolean
        """
        for node in nodes:
            sec_type = node.type
            if sec_type in self.mech_dict:
                if ('cable' in self.mech_dict[sec_type]) & reset_cable:  # cable properties must be set first, as they
                                                                         # can change nseg, which will affect insertion
                                                                         # of membrane mechanism gradients
                    self._init_cable(node)
                for mech_name in (mech_name for mech_name in self.mech_dict[sec_type]
                                  if not mech_name in ['cable', 'ions']):
                    node.sec.insert(mech_name)
                    self._modify_mechanism(node, mech_name, self.mech_dict[sec_type][mech_name])
                # ion-related parameters do not exist until after membrane mechanisms have been inserted
                if 'ions' in self.mech_dict[sec_type]:
                    self._modify_mechanism(node, 'ions', self.mech_dict[sec_type]['ions'])

    def _modify_mechanism(self, node, mech_name, mech_content):
        """
        This method loops through all the parameters for a single mechanism specified in the mechanism dictionary and
        calls self._parse_mech_content to interpret the rules and set the values for the given node.
        :param node: :class:'SHocNode'
        :param mech_name: str
        :param mech_content: dict
        """
        if not mech_content is None:
            for param_name in mech_content:
                self._parse_mech_content(node, mech_name, param_name, mech_content[param_name])

    def _parse_mech_content(self, node, mech_name, param_name, rules):
        """
        This method loops through all the segments in a node and sets the value(s) for a single mechanism parameter by
        interpreting the rules specified in the mechanism dictionary. Properly handles ion channel gradients and
        inheritance of values from the closest segment of a specified type of section along the path from root.
        :param node: :class:'SHocNode'
        :param mech_name: str
        :param param_name: str
        :param rules: dict
        """
        if rules['origin'] == 'self':
            baseline = rules['value']
            donor = None
        elif rules['origin'] in sec_types:
            sec_type = rules['origin']
            donor = self._get_node_along_path_to_root(node, sec_type)
            if (mech_name == 'cable') & (param_name == 'spatial_res'):
                baseline = self._get_spatial_res(donor)
            else:
                baseline = self._inherit_mech_param(node, donor, mech_name, param_name)
        else:
            raise Exception('Cannot inherit mechanism: {} parameter: {} from unknown sec_type: {}'.format(mech_name,
                                                                                                param_name, node.type))
        if mech_name == 'cable':  # cable properties can be inherited, but cannot be specified as gradients
            if param_name == 'spatial_res':
                node.init_nseg(baseline)
            else:
                setattr(node.sec, param_name, baseline)
                node.init_nseg(self._get_spatial_res(node))
            node.reinit_diam()
        elif 'slope' in rules:
            if not donor is None:
                distance = self.get_distance_to_node(donor, node)
                baseline += rules['slope'] * distance
            for seg in node.sec:
                value = baseline + rules['slope'] * seg.x * node.sec.L
                if (not rules['min'] is None) & (value < rules['min']):
                    value = rules['min']
                elif (rules['min'] is None) & (value < 0):
                    value = 0
                elif (not rules['max'] is None) & (value > rules['max']):
                    value = rules['max']
                if mech_name == 'ions':
                    setattr(seg, param_name, value)
                else:
                    setattr(getattr(seg, mech_name), param_name, value)
        else:
            if mech_name == 'ions':
                setattr(node.sec, param_name, baseline)
            else:
                setattr(node.sec, param_name+"_"+mech_name, baseline)

    def get_dendrite_origin(self, node):
        """
        This method determines the section type of the given node, and returns the node representing the primary branch
        point for the given section type. Basal and trunk sections originate at the soma, and apical and tuft dendrites
        originate at the trunk. For spines, recursively calls with parent node to identify the parent branch first.
        :param node: :class:'SHocNode'
        :return: :class:'SHocNode'
        """
        sec_type = node.type
        if sec_type in ['spine_head', 'spine_neck']:
            return self.get_dendrite_origin(node.parent)
        elif sec_type in ['basal', 'trunk', 'axon', 'ais']:
            return self._get_node_along_path_to_root(node, 'soma')
        elif sec_type in ['apical', 'tuft']:
            return self._get_node_along_path_to_root(node, 'trunk')
        elif sec_type == 'soma':
            return node

    def _get_node_along_path_to_root(self, node, sec_type):
        """
        This method follows the path from the given node to the root node, and returns the first node with section type
        sec_type.
        :param node: :class:'SHocNode'
        :param sec_type: str
        :return: :class:'SHocNode'
        """
        parent = node
        while not parent is None:
            if parent in self.soma and not sec_type == 'soma':
                parent = None
            elif parent.type == sec_type:
                return parent
            else:
                parent = parent.parent
        raise Exception('The path from node: {} to root does not contain sections of type: {}'.format(node.name,
                                                                                    sec_type))

    def _inherit_mech_param(self, node, donor, mech_name, param_name):
        """
        When the mechanism dictionary specifies that a node inherit a parameter value from a donor node, this method
        returns the value of that parameter found in the section or final segment of the donor node.
        :param node: :class:'SHocNode'
        :param donor: :class:'SHocNode'
        :param mech_name: str
        :param param_name: str
        :return: list: float
        """
        try:
            if mech_name in ['cable', 'ions']:
                return getattr(donor.sec, param_name)
            else:
                loc = donor.sec.nseg/(donor.sec.nseg + 1.)  # accesses the last segment of the section
                return getattr(getattr(donor.sec(loc), mech_name), param_name)
        except (AttributeError, NameError):
            print 'Exception: Mechanism: {} parameter: {} cannot be inherited from sec_type: {}'.format(mech_name,
                                                                                                param_name, donor.type)
            raise KeyError

    def _get_spatial_res(self, node):
        """
        Checks the mechanism dictionary if the section type of this node has a specified spatial resolution factor.
        Used to scale the number of segments per section in the hoc model by a factor of an exponent of 3.
        :param node: :class:'SHocNode
        :return: int
        """
        try:  # if spatial_res has not been specified for the origin type of section, it defaults to 0
            rules = self.mech_dict[node.type]['cable']['spatial_res']
        except KeyError:
            return 0
        if rules['origin'] == 'self':
            return rules['value']
        elif rules['origin'] in sec_types:  # if this sec_type also inherits the value, continue following the path
            return self._get_spatial_res(self._get_node_along_path_to_root(node, rules['origin']))

    def modify_mech_param(self, sec_type, mech_name, param_name=None, value=None, origin=None, slope=None, max=None,
                                                                                                    min=None):
        """
        Modifies or inserts new membrane mechanisms into hoc sections of type sec_type. First updates the mechanism
        dictionary, the sets the corresponding hoc parameters. This method is meant to be called manually during initial
        model specification, or during parameter optimization. For modifications to persist across simulations, the
        mechanism dictionary must be saved to a file using self.export_mech_dict() and re-imported during HocCell
        initialization.
        :param sec_type: str
        :param mech_name: str
        :param param_name: str
        :param value: float
        :param origin: str
        :param slope: float
        :param max: float
        :param min: float
        """
        insert = False
        if param_name is None:  # In this case, the mechanism will be inserted, but no parameters will be specified.
            if mech_name in ['cable', 'ions']:
                raise Exception('No parameter specified for mechanism: {}'.format(mech_name))
            else:
                mech_content = None
                insert = True
        else:
            if origin is None:
                origin = 'self'
            elif not origin in sec_types:
                raise Exception('Cannot inherit mechanism: {} parameter: {} from unknown sec_type: {}'.format(mech_name,
                                                                                                    param_name, origin))
            rules = {'origin': origin}
            if origin == 'self':
                if not value is None:
                    rules['value'] = value
                else:
                    raise Exception('No value provided for mechanism: {} parameter: {}'.format(mech_name, param_name))
            if not slope is None:
                rules['slope'] = slope
                rules['max'] = max
                rules['min'] = min
            mech_content = {param_name: rules}
        if not sec_type in self.mech_dict:  # No mechanisms have been inserted into this type of section yet
            """
            if (self.get_nodes_of_subtype(sec_type)):
                backup_content = None
                self.mech_dict[sec_type] = {mech_name: mech_content}
                insert = True
            else:
                # raise Exception('This cell does not yet contain sections of sec_type: {}'.format(sec_type))
                print 'This cell does not yet contain sections of sec_type: {}'.format(sec_type)
            """
            backup_content = None
            self.mech_dict[sec_type] = {mech_name: mech_content}
            insert = True
        elif not mech_name in self.mech_dict[sec_type]:
            backup_content = copy.deepcopy(self.mech_dict[sec_type])
            self.mech_dict[sec_type][mech_name] = mech_content  # This mechanism has not yet been inserted into this
            insert = True                                       # type of section
        elif self.mech_dict[sec_type][mech_name] is None:
            backup_content = copy.deepcopy(self.mech_dict[sec_type])
            self.mech_dict[sec_type][mech_name] = mech_content  # This mechanism has been inserted, but no parameters
                                                                # have been specified
        elif not param_name is None:
            backup_content = copy.deepcopy(self.mech_dict[sec_type])
            self.mech_dict[sec_type][mech_name][param_name] = rules  # This mechanism has been inserted, but this
                                                                     # parameter has not yet been specified
        try:
            if mech_name == 'cable':  # all membrane mechanisms in sections of type sec_type must be reinitialized after
                                      # changing cable properties
                if param_name in ['Ra', 'cm', 'spatial_res']:
                    self._reinit_mech(self.get_nodes_of_subtype(sec_type), reset_cable=1)
                else:
                    print 'Exception: Unknown cable property: {}'.format(param_name)
                    raise KeyError
            else:
                for node in self.get_nodes_of_subtype(sec_type):
                    if not mech_name == 'ions':
                        if insert:
                            try:
                                node.sec.insert(mech_name)
                            except ValueError:
                                print 'Exception: Unknown mechanism: {}'.format(mech_name)
                                raise KeyError
                        try:
                            self._modify_mechanism(node, mech_name, mech_content)
                        except AttributeError:
                            if insert:
                                self.mech_dict[sec_type][mech_name] = None
                                raise Exception('Mechanism: {} inserted but parameter: {} is unknown'.format(mech_name,
                                                                                                            param_name))
                            else:
                                self.mech_dict[sec_type] = copy.deepcopy(backup_content)
                                raise Exception('Mechanism: {} parameter: {} is unknown'.format(mech_name, param_name))
                    else:
                        try:
                            self._modify_mechanism(node, mech_name, mech_content)
                        except AttributeError:
                            print 'Exception: Unknown ion parameter: {}'.format(param_name)
                            raise KeyError
                        except NameError:
                            print 'Exception: Ion parameter: {} has not yet been defined for section type: {}'.format(
                                                                                                param_name, node.type)
                            raise KeyError
        except KeyError:
            if backup_content is None:
                del self.mech_dict[sec_type]
            else:
                self.mech_dict[sec_type] = copy.deepcopy(backup_content)
        finally:
            if verbose:
                pprint.pprint(self.mech_dict)

    def export_mech_dict(self, mech_filename=None):
        """
        Following modifications to the mechanism dictionary either during model specification or parameter optimization,
        this method stores the current mech_dict to a pickle file stamped with the date and time. This allows the
        current set of mechanism parameters to be recalled later.
        """
        if mech_filename is None:
            mech_filename = 'mech_dict_'+datetime.datetime.today().strftime('%m%d%Y%H%M')+'.pkl'
        write_to_pkl(data_dir+mech_filename, self.mech_dict)
        print "Exported mechanism dictionary to "+mech_filename

    def get_node_by_distance_to_soma(self, distance, sec_type):
        """
        Gets the first node of the given section type at least the given distance from a soma node.
        Not particularly useful, since it will always return the same node.
        :param distance: int or float
        :param sec_type: str
        :return: :class:'SHocNode'
        """
        nodes = self._node_dict[sec_type]
        for node in nodes:
            if self.get_distance_to_node(self.tree.root, node) >= distance:
                return node
        raise Exception('No node is {} um from a soma node.'.format(distance))

    def get_distance_to_node(self, root, node, loc=None):
        """
        Returns the distance from the given location on the given node to its connection with a root node.
        :param root: :class:'SHocNode'
        :param node: :class:'SHocNode'
        :param loc: float
        :return: int or float
        """
        length = 0.
        if node in self.soma:
            return length
        if not loc is None:
            length += loc*node.sec.L
        if root in self.soma:
            while not node.parent in self.soma:
                node.sec.push()
                loc = h.parent_connection()
                h.pop_section()
                node = node.parent
                length += loc*node.sec.L
        elif self.node_in_subtree(root, node):
            while not node.parent is root:
                node.sec.push()
                loc = h.parent_connection()
                h.pop_section()
                node = node.parent
                length += loc*node.sec.L
        else:
            return None  # node is not connected to root
        return length

    def node_in_subtree(self, root, node):
        """
        Checks if a node is contained within a subtree of root.
        :param root: 'class':SNode2 or SHocNode
        :param node: 'class':SNode2 or SHocNode
        :return: boolean
        """
        nodelist = []
        self.tree._gather_nodes(root, nodelist)
        if node in nodelist:
            return True
        else:
            return False

    def get_path_length_swc(self, path):
        """
        Calculates the distance between nodes given a list of SNode2 nodes connected in a path.
        :param path: list : :class:'SNode2'
        :return: int or float
        """
        distance = 0
        for i in range(len(path)-1):
            distance += np.sqrt(np.sum((path[i].get_content()['p3d'].xyz - path[i+1].get_content()['p3d'].xyz)**2))
        return distance

    def get_node_length_swc(self, node):
        """
        Calculates the distance between the center points of an SNode2 node and its parent.
        :param node: :class:'SNode2'
        :return: float
        """
        if not self.tree.is_root(node):
            return np.sqrt(np.sum((node.get_content()['p3d'].xyz - node.parent.get_content()['p3d'].xyz)**2))
        else:
            return np.sqrt(np.sum(node.get_content()['p3d'].xyz**2))

    def get_branch_order(self, node):
        """
        Calculates the branch order of a SHocNode node. The order is defined as 0 for all soma, axon, and apical trunk
        dendrite nodes, but defined as 1 for basal dendrites that branch from the soma, and apical and tuft dendrites
        that branch from the trunk. Increases by 1 after each additional branch point. Makes sure not to count spines.
        :param node: :class:'SHocNode'
        :return: int
        """
        order = 0
        while not node in self.soma+self.axon+self.trunk:
            if len([child for child in node.parent.children if not child.type == 'spine_neck']) > 1:
                order += 1
                node = node.parent
        return order

    def is_terminal(self, node):
        """
        Calculates if a node is a terminal branch.
        :param node: :class:'SHocNode'
        :return: Boolean
        """
        while node.children:
            if len([child for child in node.children if not child.type == 'spine_neck']) > 1:
                return False
            else:
                node = node.children[0]
        return True

    def set_stochastic_synapses(self, value):
        """
        This method turns stochastic filtering of release probability on or off for all synapses contained in this cell.
        :param value: int in [0, 1]
        """
        for nodelist in self._node_dict.itervalues():
            for node in nodelist:
                for syn in node.synapses:
                    syn.stochastic = value

    @property
    def gid(self):
        return self._gid

    @property
    def soma(self):
        return self._node_dict['soma']
            
    @property
    def axon(self):
        return self._node_dict['axon']
            
    @property
    def basal(self):
        return self._node_dict['basal']

    @property
    def apical(self):
        return self._node_dict['apical']

    @property
    def trunk(self):
        return self._node_dict['trunk']

    @property
    def tuft(self):
        return self._node_dict['tuft']

    @property
    def spine(self):
        return self._node_dict['spine']

#------------------------------Extend SNode2 to interact with NEURON hoc sections------------------------


class SHocNode(btmorph.btstructs2.SNode2):
    """
    Extends SNode2 with some methods for storing and retrieving additional information in the node's content
    dictionary related to running NEURON models specified in the hoc language.
    """

    def __init__(self, index=0):
        """
        :param index: int : unique node identifier
        """
        btmorph.btstructs2.SNode2.__init__(self, index)
        self.content['spines'] = []
        self.content['synapses'] = []

    def get_sec(self):
        """
        Returns the hoc section associated with this node, stored in the node's content dictionary.
        :return: :class:'neuron.h.Section'
        """
        if 'sec' in self.content:
            return self.content['sec']
        else:
            raise Exception('This node does not yet have an associated hoc section.')

    def set_sec(self, sec):
        """
        Stores the hoc section associated with this node in the node's content dictionary.
        :param sec: :class:'neuron.h.Section'
        """
        self.content['sec'] = sec

    sec = property(get_sec, set_sec)

    def init_nseg(self, spatial_res=0):
        """
        Initializes the number of hoc segments in this node's hoc section (nseg) based on the AC length constant.
        Must be re-initialized whenever basic cable properties Ra or cm are changed. If the node is a tapered cylinder,
        it should contain at least 3 segments. The spatial resolution parameter increases the number of segments per
        section by a factor of an exponent of 3.
        If a section's nseg has been manually increased beyond the suggestion of the mechanism dictionary, this method
        does not decrease it.
        :param spatial_res: int
        """
        sugg_nseg = d_lambda_nseg(self.sec)
        if not self.get_diam_bounds() is None:
            sugg_nseg = max(sugg_nseg, 3)
        sugg_nseg *= 3**spatial_res
        if self.sec.nseg < sugg_nseg:
            self.sec.nseg = sugg_nseg

    def reinit_diam(self):
        """
        For a node associated with a hoc section that is a tapered cylinder, every time the spatial resolution
        of the section (nseg) is changed, the section diameters must be reinitialized. This method checks the
        node's content dictionary for diameter boundaries and recalibrates the hoc section associated with this node.
        """
        if not self.get_diam_bounds() is None:
            [diam1, diam2] = self.get_diam_bounds()
            h('diam(0:1)={}:{}'.format(diam1, diam2), sec=self.sec)

    def get_diam_bounds(self):
        """
        If the hoc section associated with this node is a tapered cylinder, this method returns a list containing
        the values of the diameters at the 0 and 1 ends of the section, stored in the node's content dictionary.
        Otherwise, it returns None (for non-conical cylinders).
        :return: (list: int) or None
        """
        if 'diam' in self.content:
            return self.content['diam']
        else:
            return None

    def set_diam_bounds(self, diam1, diam2):
        """
        For a node associated with a hoc section that is a tapered cylinder, this stores a list containing the values
        of the diameters at the 0 and 1 ends of the section in the node's content dictionary.
        :param diam1: int
        :param diam2: int
        """
        self.content['diam'] = [diam1, diam2]
        self.reinit_diam()

    def get_type(self):
        """
        NEURON sections are assigned a node type for convenience in order to later specify membrane mechanisms and
        properties for each type of compartment.
        :return: str
        """
        if 'type' in self.content:
            return self.content['type']
        else:
            raise Exception('This node does not yet have a defined type.')

    def set_type(self, type):
        """
        Checks that type is a string in the list of defined section types, and stores the value in the node's content
        dictionary.
        :param type: str
        """
        if type in sec_types:
            self.content['type'] = type
        else:
            raise Exception('That is not a defined type of section.')

    type = property(get_type, set_type)

    def connect(self, parent, pindex=1, cindex=0):
        """
        Connects this SHocNode node to a parent node, and establishes a connection between their associated
        hoc sections.
        :param parent: :class:'SHocNode'
        :param pindex: int in [0,1] Connect to this end of the parent hoc section.
        :param cindex: int in [0,1] Connect this end of the child hoc section
        """
        self.parent=parent
        parent.add_child(self)
        self.sec.connect(parent.sec, pindex, cindex)

    @property
    def name(self):
        """
        Returns a str containing the name of the hoc section associated with this node. Consists of a type descriptor
        and an index identifier.
        :return: str
        """
        if 'type' in self.content:
            return '{0.type}{0.index}'.format(self)
        else:
            raise Exception('This node does not yet have a defined type.')

    @property
    def spines(self):
        """
        Returns a list of the spine head sections attached to the hoc section associated with this node.
        :return: list of :class:'SHocNode' of sec_type == 'spine_head'
        """
        return self.content['spines']

    @property
    def synapses(self):
        """
        Returns a list of the objects of :class:'Synapse' associated with this node.
        :return: list of hoc objects, type depends on .mod file(s) used to implement synapses
        """
        return self.content['synapses']


class CA1_Pyr(HocCell):
    def __init__(self, morph_filename=None, mech_filename=None, full_spines=True):
        HocCell.__init__(self, morph_filename, mech_filename)
        if full_spines:
            self.insert_spines()

    def insert_spines(self):
        """
        This method populates the cell tree with spines following spine density information from Erk Bloss &
        Nelson Spruston. Basal dendrites have no spines until the first branch point, and a higher density beyond the
        second branch point. Trunk dendrites have no spines until the first branch point, and an increasing density
        until the tuft branch point(s). Apical dendrites have a density that varies with the distance from the soma of
        their original branch point from the trunk. Terminal tuft branches have a higher density than their parents.
        Should standardize the implementation of the rules for each type of dendrite and import the
        density dictionary from a file, similar to the implementation of the membrane mechanism dictionary.
        """
        np.random.seed(self.gid)  # This cell will always have the same spine locations
        densities = {'trunk': {'min': 0.2418, 'max': 3.8,
                               'start': self.get_distance_to_node(self.tree.root, self.apical[0]),
                               'end': self.get_distance_to_node(self.tree.root, self.tuft[0])},
                     'basal': {'1': 0., '2': 0.4428, '3': 1.891},
                     'apical': {'min': 2.273, 'max': 2.688,
                                'start': self.get_distance_to_node(self.tree.root, self.apical[0]),
                                'end': max([self.get_distance_to_node(self.tree.root, branch)
                                            for branch in self.apical if self.get_branch_order(branch) == 1])},
                     'tuft': {'1': 1.354, 'terminal': 0.7157}
                    }
        for node in self.basal:
            order = self.get_branch_order(node)
            if order == 2:
                self.insert_spines_every(node, densities['basal']['2'])
            elif order > 2:
                self.insert_spines_every(node, densities['basal']['3'])
        for node in self.trunk:
            distance = self.get_distance_to_node(self.tree.root, node)
            if distance >= densities['trunk']['start']:
                slope = (densities['trunk']['max'] - densities['trunk']['min']) / \
                        (densities['trunk']['end'] - densities['trunk']['start'])
                density = densities['trunk']['min'] + slope * (distance - densities['trunk']['start'])
                self.insert_spines_every(node, density)
        for node in self.apical:
            distance = self.get_distance_to_node(self.tree.root, self.get_dendrite_origin(node), loc=1)
            slope = (densities['apical']['max'] - densities['apical']['min']) / \
                    (densities['apical']['end'] - densities['apical']['start'])
            density = densities['apical']['min'] + slope * (distance - densities['apical']['start'])
            self.insert_spines_every(node, density)
        for node in self.tuft:
            if self.is_terminal(node):
                self.insert_spines_every(node, densities['tuft']['terminal'])
            else:
                self.insert_spines_every(node, densities['tuft']['1'])
        self._reinit_mech(self.spine)

    def insert_spines_every(self, node, lam):
        """
        Given a mean inter-spine interval in um, insert spines in the node at the specified density.
        :param node: :class:'SHocNode'
        :param lam: float: mean interval
        """
        L = node.sec.L
        interval = np.random.poisson(10000.*lam)/10000.  # random intervals with correct significant digits
        while interval < L:
            self.insert_spine(node, interval/L)
            interval += np.random.poisson(10000.*lam)/10000.

    def insert_spine(self, node, parent_loc, child_loc=0):
        """
        Spines consist of two hoc sections: a cylindrical spine head and a cylindrical spine neck.
        :param node: :class:'SHocNode'
        :param parent_loc: float
        :param child_loc: int
        """
        neck = self.make_section('spine_neck')
        neck.connect(node, parent_loc, child_loc)
        neck.sec.L = 1.58
        neck.sec.diam = 0.077
        self._init_cable(neck)
        head = self.make_section('spine_head')
        head.connect(neck)
        node.spines.append(head)
        head.sec.L = 0.408  # matches surface area of sphere with diam = 0.5
        head.sec.diam = 0.408
        self._init_cable(head)


class Synapse(object):
    """
    The implementation in hoc of synaptic mechanisms that can be triggered is complicated. This container is an attempt
    to wrap all the objects required to deliver synaptic events to a section, and have separable synaptic mechanisms
    (e.g. GluA-Rs and GluN-Rs) respond with individually specifiable weights and kinetics.
    To make model specification and simulation implementation straightforward, synapses are not meant to be moved once
    they are initialized.
    """
    def __init__(self, cell, node, type_list=None, stochastic=1, loc=0.5, delay=0, source=None):
        """
        Design goals: A source (like a spike detector in a presynaptic neuron) can be specified. If not, a VecStim
        object is used a source, which can be played events at specified times using its .play method. If stochastic,
        all spikes are intercepted by a point process with release probability dynamics and its own unique and
        independent random variable from a uniform distribution. If not, the specified synaptic mechanisms are connected
        directly to the source of spikes.
        :param cell: :class:'HocCell'
        :param node: :class:'SHoCNode'
        :param type_list: list of str
        :param stochastic: int in [0, 1]
        :param loc: int or float
        :param delay: int or float
        :param source: :class:'h.VecStim' or other source of spike events
        """
        self._cell = cell
        self._node = node
        self._stochastic = stochastic
        self._loc = loc
        self._delay = delay
        self._syn = {}
        self.randObj = None
        if not source is None:
            self.source = source
        else:
            self.source = h.VecStim()
        if type_list is None:
            type_list = ['AMPA_S']
        if self.stochastic:
            self._init_stochastic()
        for type in type_list:
            syn = getattr(h, type)(self.node.sec(self.loc))
            self._syn[type] = {'target': syn}
            if self.stochastic:
                self._syn[type]['netcon'] = h.NetCon(self.target('Pr'), syn)
            else:
                self._syn[type]['netcon'] = h.NetCon(self.source, syn)
            self.netcon(type).delay = self.delay
            self.netcon(type).weight[0] = 1
        self._node.synapses.append(self)

    def _init_stochastic(self):
        """
        This method constructs and initializes a stochastic filtering mechanism that intercepts spikes delivered to this
        synapse and calculates whether or not to pass a spike to the rest of the specified synaptic mechanisms.
        """
        if self.randObj is None:  # if this synapse has never been stochastic, it needs a new random number generator
            self.randObj = h.Random()
            self.randObj.MCellRan4(1, self.cell.gid*1e10+self.node.index*1e4+len(self.node.synapses)+1)
            # a unique seed for up to 10,000 synapses per node and 1,000,000 sections per cell
            self.randObj.uniform(0,1)
        else:  # if this synapse has already been stochastic before, this restarts its random number generator
            self.randObj.seq(1)
        syn = getattr(h, 'Pr')(self.node.sec(self.loc))
        self._syn['Pr'] = {'target': syn}
        self._syn['Pr']['netcon'] = h.NetCon(self.source, syn)
        self.netcon('Pr').delay = 0
        self.netcon('Pr').weight[0] = 1
        self.target('Pr').setRandObjRef(self.randObj)

    def target(self, type):
        """
        Returns the hoc object for the synaptic mechanism of the specified type
        :param type: str
        :return: :class:'h.HocObject'
        """
        return self._syn[type]['target']

    def netcon(self, type):
        """
        Returns the hoc network connection linking the synaptic mechanism of the specified type to a source of spikes.
        :param type: str
        :return: :class:'h.NetCon'
        """
        return self._syn[type]['netcon']

    def get_stochastic(self):
        """
        Returns the value of an internal variable indicating if this synapse has a stochastic filter for spikes.
        :return: int in [0, 1]
        """
        return self._stochastic

    def set_stochastic(self, value):
        """
        Turns on or off stochastic filtering of spikes.
        :param value: int in [0, 1]
        """
        if not (value == self._stochastic):
            self._stochastic = value
            if value:
                self._init_stochastic()
                for type in (type for type in self._syn if not type == 'Pr'):
                    del self._syn[type]['netcon']
                    self._syn[type]['netcon'] = h.NetCon(self.target('Pr'), self.target(type))
                    self.netcon(type).delay = self._delay
                    self.netcon(type).weight[0] = 1
            else:
                for type in (type for type in self._syn if not type == 'Pr'):
                    del self._syn[type]['netcon']
                    self._syn[type]['netcon'] = h.NetCon(self.source, self.target(type))
                    self.netcon(type).delay = self._delay
                    self.netcon(type).weight[0] = 1
                del self._syn['Pr']

    stochastic = property(get_stochastic, set_stochastic)

    def get_delay(self):
        """
        Returns the value of the time delay (ms) between spike and activation for the specified synaptic mechanisms.
        :return: int or float
        """
        return self._delay

    def set_delay(self, value):
        """
        Changes the value of the time delay (ms) between spike and activation for the specified synaptic mechanisms.
        :param value: int or float
        """
        self._delay = value
        for type in (type for type in self._syn if not type == 'Pr'):
            self.netcon(type).delay = value

    delay = property(get_delay, set_delay)

    @property
    def cell(self):
        """
        Returns the cell containing this synapse.
        :return: :class:'HocCell'
        """
        return self._cell

    @property
    def node(self):
        """
        Returns the node containing this synapse.
        :return: :class:'SHocNode'
        """
        return self._node

    @property
    def loc(self):
        """
        Returns the location along the hoc section containing this synapse.
        :return: int or float
        """
        return self._loc
