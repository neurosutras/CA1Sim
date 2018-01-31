"""

"""
from nested.utils import *


class StateMachine(object):
    """

    """

    def __init__(self, ti=0., dt=1., states=None, rates=None):
        """

        :param ti: float
        :param dt: float (milliseconds)
        :param states: dict
        :param rates: dict (proportion of total/second)
        """
        self.dt = dt
        self.init_states = dict()
        self.states = dict()
        self.states_history = dict()
        self.rates = dict()
        if states is not None:
            self.update_states(states)
        if rates is not None:
            self.update_rates(rates)  # {'A': {'B': constant or iterable}}
        self.ti = ti
        self.t = self.ti
        self.t_history = np.array([self.t])
        self.i = 0

    def reset(self):
        """

        :param ti:
        :return:
        """
        self.t = self.ti
        self.i = 0
        self.t_history = np.array([self.t])
        self.states = dict(self.init_states)
        for s0 in self.states:
            self.states_history[s0] = np.array([self.states[s0]])

    def get_current_rates(self):
        """

        :return: dict
        """
        current = {}
        for s0 in self.rates:
            if s0 not in current:
                current[s0] = {}
            for s1 in self.rates[s0]:
                r = self.rates[s0][s1]
                if hasattr(r, '__iter__'):
                    if len(r) - 1 < self.i:
                        raise Exception('StateMachine: Insufficient array length for non-stationary rate: %s to %s ' %
                                        (s0, s1))
                    this_r = r[self.i]
                else:
                    this_r = r
                current[s0][s1] = this_r
        return current

    def update_transition(self, s0, s1, r):
        """

        :param s0: str
        :param s1: str
        :param r: float or array
        """
        if s0 not in self.states:
            raise Exception('StateMachine: Cannot update transition from invalid state: %s' % s0)
        if s1 not in self.states:
            raise Exception('StateMachine: Cannot update transition to invalid state: %s' % s1)
        if s0 not in self.rates:
            self.rates[s0] = {}
        self.rates[s0][s1] = r

    def update_rates(self, rates):
        """

        :param rates: dict
        """
        for s0 in rates:
            for s1, r in rates[s0].iteritems():
                self.update_transition(s0, s1, r)

    def update_states(self, states):
        """

        :param states: dict
        """
        for s, v in states.iteritems():
            self.init_states[s] = v
            self.states[s] = v
            self.states_history[s] = np.array([v])

    def get_out_rate(self, state):
        """

        :param state: str
        :return: float
        """
        if state not in self.states:
            raise Exception('StateMachine: Invalid state: %s' % state)
        if state not in self.rates:
            return 0.
        out_rate = 0.
        for s1 in self.rates[state]:
            r = self.rates[state][s1]
            if hasattr(r, '__iter__'):
                if len(r) - 1 < self.i:
                    raise Exception('StateMachine: Insufficient array length for non-stationary rate: %s to %s ' %
                                    (state, s1))
                this_r = r[self.i]
            else:
                this_r = r
            out_rate += this_r
        return out_rate

    def step(self, n=1):
        """
        Need to think about how to deal with weights that vary from 0.8 to unknown amount....
        :param n: int
        """
        for i in range(n):
            next_states = dict(self.states)
            for s0 in self.rates:
                factor = 1.
                if self.states[s0] > 0.:
                    total_out = self.get_out_rate(s0) * self.dt / 1000. * self.states[s0]
                    if total_out > self.states[s0]:
                        factor = self.states[s0] / total_out
                for s1 in self.rates[s0]:
                    r = self.rates[s0][s1]
                    if hasattr(r, '__iter__'):
                        if len(r) - 1 < self.i:
                            raise Exception('StateMachine: Insufficient array length for non-stationary rate: %s to '
                                            '%s ' % (s0, s1))
                        this_r = r[self.i]
                    else:
                        this_r = r
                    # print 'this_r: %.4E, factor: %.4E, %s: %.4E' % (this_r, factor, s0, self.states[s0])
                    this_delta = this_r * self.dt / 1000. * factor * self.states[s0]
                    next_states[s0] -= this_delta
                    next_states[s1] += this_delta
            self.states = dict(next_states)
            for s0 in self.states:
                self.states_history[s0] = np.append(self.states_history[s0], self.states[s0])
            self.i += 1
            self.t += self.dt
            self.t_history = np.append(self.t_history, self.t)

    def run(self):
        """

        """
        self.reset()
        min_steps = None
        for s0 in self.rates:
            for s1 in self.rates[s0]:
                r = self.rates[s0][s1]
                if hasattr(r, '__iter__'):
                    if min_steps is None:
                        min_steps = len(r)
                    else:
                        min_steps = min(min_steps, len(r))
        if min_steps is None:
            raise Exception('StateMachine: Use step method to specify number of steps for stationary process.')
        self.step(min_steps)

    def plot(self, states=None):
        """

        :param states:
        """
        if states is None:
            states = self.states.keys()
        elif not hasattr(states, '__iter__'):
            states = [states]
        fig, axes = plt.subplots(1)
        for state in states:
            if state in self.states:
                axes.plot(self.t_history, self.states_history[state], label=state)
            else:
                print 'StateMachine: Not including invalid state: %s' % state
        axes.set_xlabel('Time (ms)')
        axes.set_ylabel('Occupancy')
        axes.legend(loc='best', frameon=False, framealpha=0.5)
        clean_axes(axes)
        plt.show()
        plt.close()


def sigmoid_segment(slope, th, xlim=None, ylim=None):
    """
    Transform a sigmoid to intersect x and y range limits.
    :param slope: float
    :param th: float
    :param xlim: pair of float
    :param ylim: pair of float
    :return: callable
    """
    if xlim is None:
        xlim = (0., 1.)
    if ylim is None:
        ylim = (0., 1.)
    if th < xlim[0] or th > xlim[1]:
        raise ValueError('sigmoid_segment: th: %.2E is out of range for xlim: %s' % (th, str(xlim)))
    y = lambda x: 1. / (1. + np.exp(-slope * (x - th)))
    target_amp = ylim[1] - ylim[0]
    y0 = y(xlim[0])
    y1 = y(xlim[1])
    current_amp = y1 - y0
    return lambda x: (target_amp / current_amp) * (1. / (1. + np.exp(-slope * (x - th))) - y0) + ylim[0]


def subtract_baseline(waveform, baseline=None):
    """

    :param waveform: array
    :param baseline: float
    :return: array
    """
    new_waveform = np.array(waveform)
    if baseline is None:
        baseline = np.mean(new_waveform[np.where(new_waveform <= np.percentile(new_waveform, 10.))[0]])
    new_waveform -= baseline
    return new_waveform, baseline


def wrap_around_and_compress(waveform, interp_x):
    """

    :param waveform: array of len(3 * interp_x)
    :param interp_x: array
    :return: array of len(interp_x)
    """
    before = np.array(waveform[:len(interp_x)])
    after = np.array(waveform[2 * len(interp_x):])
    within = np.array(waveform[len(interp_x):2 * len(interp_x)])
    waveform = within[:len(interp_x)] + before[:len(interp_x)] + after[:len(interp_x)]
    return waveform
