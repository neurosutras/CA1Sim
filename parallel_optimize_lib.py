__author__ = 'milsteina'

class EPSP_Amp_Custom_Callback(object):
    """
    Used by parallel_optimize_EPSP_amp to prints a status report during iterations of optimize.minimize
    :param xk: list
    """
    def __init__(self, syn, pid):
        self.spine = syn.node.index
        self.node = syn.node.parent.parent.name
        self.pid = pid
        self.simiter = 0

    def __call__(self, xk):
        self.simiter += 1
        print 'Process:', self.pid, ', Node:', self.node, 'Spine:', self.spine, ', Current x:', xk, 'after', \
            self.simiter, 'iterations'