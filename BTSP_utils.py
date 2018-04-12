"""

"""
from nested.utils import *
from scipy.optimize import minimize


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


def get_indexes_from_ramp_bounds_with_wrap(x, start, peak, end, min):
    """

    :param x: array
    :param start: float
    :param peak: float
    :param end: float
    :param min: float
    :return: tuple of float: (start_index, peak_index, end_index, min_index)
    """
    peak_index = np.where(x >= peak)[0]
    if np.any(peak_index):
        peak_index = peak_index[0]
    else:
        peak_index = len(x) - 1
    min_index = np.where(x >= min)[0]
    if np.any(min_index):
        min_index = min_index[0]
    else:
        min_index = len(x) - 1
    if start < peak:
        start_index = np.where(x[:peak_index] <= start)[0][-1]
    else:
        start_index = peak_index + np.where(x[peak_index:] <= start)[0]
        if np.any(start_index):
            start_index = start_index[-1]
        else:
            start_index = len(x) - 1
    if end < peak:
        end_index = np.where(x > end)[0][0]
    else:
        end_index = peak_index + np.where(x[peak_index:] > end)[0]
        if np.any(end_index):
            end_index = end_index[0]
        else:
            end_index = len(x) - 1
    return start_index, peak_index, end_index, min_index


def calculate_ramp_features(local_context, ramp, induction_loc, offset=False, smooth=False):
    """

    :param local_context: :class:'Context'
    :param ramp: array
    :param induction_loc: float
    :param offset: bool
    :param smooth: bool
    :return tuple of float
    """
    binned_x = local_context.binned_x
    track_length = local_context.track_length
    default_interp_x = local_context.default_interp_x
    extended_binned_x = np.concatenate([binned_x - track_length, binned_x, binned_x + track_length])
    if smooth:
        local_ramp = signal.savgol_filter(ramp, 21, 3, mode='wrap')
    else:
        local_ramp = np.array(ramp)
    extended_binned_ramp = np.concatenate([local_ramp] * 3)
    extended_interp_x = np.concatenate([default_interp_x - track_length, default_interp_x,
                                        default_interp_x + track_length])
    extended_ramp = np.interp(extended_interp_x, extended_binned_x, extended_binned_ramp)
    interp_ramp = extended_ramp[len(default_interp_x):2 * len(default_interp_x)]
    baseline_indexes = np.where(interp_ramp <= np.percentile(interp_ramp, 10.))[0]
    baseline = np.mean(interp_ramp[baseline_indexes])
    if offset:
        interp_ramp -= baseline
        extended_ramp -= baseline
    peak_index = np.where(interp_ramp == np.max(interp_ramp))[0][0] + len(interp_ramp)
    peak_val = extended_ramp[peak_index]
    peak_x = extended_interp_x[peak_index]
    start_index = np.where(extended_ramp[:peak_index] <=
                           0.15 * (peak_val - baseline) + baseline)[0][-1]
    end_index = peak_index + np.where(extended_ramp[peak_index:] <= 0.15 *
                                      (peak_val - baseline) + baseline)[0][0]
    start_loc = float(start_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    end_loc = float(end_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    peak_loc = float(peak_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    min_index = np.where(interp_ramp == np.min(interp_ramp))[0][0] + len(interp_ramp)
    min_val = extended_ramp[min_index]
    min_loc = float(min_index % len(default_interp_x)) / float(len(default_interp_x)) * track_length
    peak_shift = peak_x - induction_loc
    if peak_shift > track_length / 2.:
        peak_shift = -(track_length - peak_shift)
    elif peak_shift < -track_length / 2.:
        peak_shift += track_length
    # ramp_width = extended_interp_x[end_index] - extended_interp_x[start_index]
    scaled_ramp = np.divide(local_ramp, peak_val)
    ramp_width = np.trapz(scaled_ramp, binned_x)
    before_width = induction_loc - start_loc
    if induction_loc < start_loc:
        before_width += track_length
    after_width = end_loc - induction_loc
    if induction_loc > end_loc:
        after_width += track_length
    ratio = before_width / after_width
    return peak_val, ramp_width, peak_shift, ratio, start_loc, peak_loc, end_loc, min_val, min_loc


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