__author__ = 'Aaron D. Milstein'
import numpy as np
import matplotlib.pyplot as plt
import os
import time

dt = 1.  # ms
const_vel = 20.  # cm/s
space_res = const_vel * dt / 1000.  # cm
space_res_low = 5.  # cm


def get_stepwise_array(a, axis=0):
    """
    Duplicates values within an array to produce a stepwise version of the array. axis=0 simply duplicates values to
    create the appropriate time or x-axis values. axis=1 repeats the previous value to create the appropriate output or
    y-axis values.
    :param a: array
    :param axis: int in [0, 1]
    :return: array
    """
    a_steps = np.zeros(2 * (len(a) - 2) + 2)
    a_steps[0] = a[0]
    if axis == 0:
        a_steps[-1] = a[-1]
    elif axis == 1:
        a_steps[-1] = a[-2]
    i = 1
    for j in range(1, len(a) - 1):
        if axis == 0:
            a_steps[i] = a[j]
            a_steps[i + 1] = a[j]
        elif axis == 1:
            a_steps[i] = a[j - 1]
            a_steps[i + 1] = a[j]
        i += 2
    return a_steps


X = np.arange(-100., 100., space_res)
Y = np.zeros_like(X)
# X, Y = np.meshgrid(np.arange(-120., 120., res), np.arange(-120., 120., res))
# r = G(X, Y)
V = const_vel * np.ones(len(X) - 1)
X_diff = np.diff(X)
Y_diff = np.diff(Y)
D_diff = np.sqrt(X_diff**2.+Y_diff**2.)
t_diff = np.divide(D_diff, V) * 1000.
D = np.cumsum(np.insert(D_diff, 0, 0.))
t = np.cumsum(np.insert(t_diff, 0, 0.))
t_interp = np.arange(0., t[-1], dt)

X_low = np.arange(-100., 100., space_res_low)
Y_low = np.zeros_like(X_low)
V_low = const_vel * np.ones(len(X_low) - 1)
X_low_diff = np.diff(X_low)
Y_low_diff = np.diff(Y_low)
D_low_diff = np.sqrt(X_low_diff**2.+Y_low_diff**2.)
t_low_diff = np.divide(D_low_diff, V_low) * 1000.
D_low = np.cumsum(np.insert(D_low_diff, 0, 0.))
t_low = np.cumsum(np.insert(t_low_diff, 0, 0.))

D_steps = get_stepwise_array(D_low)
t_steps = get_stepwise_array(t_low)


a = 0.3
b = -3. / 2.
u = lambda theta: (np.cos(theta), np.sin(theta))
ori = 2. * np.pi * np.array([-30., 30., 90.]) / 360.  # rads
g = lambda x: np.exp(a * (x - b)) - 1.
grid_peak_rate = 10.  # Hz


def get_1D_rate_map(xoff, yoff, lam, theta_off, X, Y, plot=True):
    """

    :param xoff: float: x offset for this grid cell (cm)
    :param yoff: float: y offset for this grid cell (cm)
    :param lam: float: grid spacing for this grid cell (cm)
    :param theta_off: float: orientation offset for this grid cell (rads)
    :param X: array: x coordinates of linear trajectory (n elements; cm)
    :param Y: array: y coordinates of linear trajectory (n elements; cm)
    :param plot: bool
    :return: array, array, array: linear distance, time, rate
    """
    G = lambda x, y: g(np.sum([np.cos(4. * np.pi / np.sqrt(3.) / lam * np.dot(u(theta - theta_off),
                                                                              (x - xoff, y - yoff))) for theta in ori]))
    G = np.vectorize(G)
    peak = G(xoff, yoff)
    r = G(X, Y) / peak * grid_peak_rate
    r_interp = np.interp(t_interp, t, r)
    if plot:
        fig, axes = plt.subplots(2)
        axes[0].plot(D, r)
        axes[0].set_ylabel('Firing Rate (Hz)')
        axes[0].set_xlabel('Distance (cm)')
        axes[1].plot(t_interp, r_interp)
        axes[1].set_ylabel('Firing Rate (Hz)')
        axes[1].set_xlabel('Time (ms)')
        plt.show()
        plt.close()
    return r_interp


def get_1D_rate_map_discrete(xoff, yoff, lam, theta_off, X, Y, stepwise=True, plot=True):
    """

    :param xoff: float: x offset for this grid cell (cm)
    :param yoff: float: y offset for this grid cell (cm)
    :param lam: float: grid spacing for this grid cell (cm)
    :param theta_off: float: orientation offset for this grid cell (rads)
    :param X: array: x coordinates of linear trajectory (n elements; cm)
    :param Y: array: y coordinates of linear trajectory (n elements; cm)
    :param plot: bool
    :return: array, array, array: linear distance, time, rate
    """
    G = lambda x, y: g(np.sum([np.cos(4. * np.pi / np.sqrt(3.) / lam * np.dot(u(theta - theta_off),
                                                                              (x - xoff, y - yoff))) for theta in ori]))
    G = np.vectorize(G)
    peak = G(xoff, yoff)
    r = G(X, Y) / peak * grid_peak_rate
    if stepwise:
        r = get_stepwise_array(r, axis=1)
        this_D = D_steps
        this_t = t_steps
    else:
        this_D = D_low
        this_t = t_low
    if plot:
        fig, axes = plt.subplots(2)
        axes[0].plot(this_D, r)
        axes[1].plot(this_t, r)
        axes[0].set_ylabel('Firing Rate (Hz)')
        axes[0].set_xlabel('Distance (cm)')
        axes[1].set_ylabel('Firing Rate (Hz)')
        axes[1].set_xlabel('Time (ms)')
        plt.show()
        plt.close()
    return r


def get_2D_rate_map(xoff, yoff, lam, theta_off, dim=100., dx=None, plot=True):
    """

    :param xoff: float: x offset for this grid cell (cm)
    :param yoff: float: y offset for this grid cell (cm)
    :param lam: float: grid spacing for this grid cell (cm)
    :param theta_off: float: orientation offset for this grid cell (rads)
    :param dim: float: length of side of square 2D space (cm)
    :param dx: float: spatial resolution to compute rate map (cm)
    :param plot: bool
    :return: array: rate
    """
    if dx is None:
        dx = space_res
    G = lambda x, y: g(np.sum([np.cos(4. * np.pi / np.sqrt(3.) / lam * np.dot(u(theta - theta_off),
                                                                              (x - xoff, y - yoff))) for theta in ori]))
    G = np.vectorize(G)
    X, Y = np.meshgrid(np.arange(-dim, dim, dx), np.arange(-dim, dim, dx))
    peak = G(xoff, yoff)
    r = G(X, Y) / peak * grid_peak_rate
    if plot:
        plt.pcolor(X, Y, r)
        plt.ylabel('Y distance (cm)')
        plt.xlabel('X distance (cm)')
        plt.show()
        plt.close()
    return r


def compute_single_rate_map((id, grid_cell)):
    """

    :param grid_cell: dict {string: value}
    :return: dict: {string: array}
    """
    start_time = time.time()
    this_grid_waveforms = {}
    this_grid_waveforms['cont'] = get_1D_rate_map(grid_cell['x_off'], grid_cell['y_off'], grid_cell['lambda'],
                                             grid_cell['theta'], X, Y, plot=False)
    this_grid_waveforms['discrete'] = get_1D_rate_map_discrete(grid_cell['x_off'], grid_cell['y_off'],
                                                               grid_cell['lambda'], grid_cell['theta'], X_low, Y_low,
                                                               stepwise=False, plot=False)
    print 'Process:', os.getpid(), 'computed rate maps for grid cell', id, 'in', time.time() - start_time, 's'

    return id, this_grid_waveforms