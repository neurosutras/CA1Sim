TITLE minimal model of GABAB receptors (NET_RECEIVE version)

COMMENT
Minimal kinetic model for GABA-B receptors (Destexhe & Sejnowski 1995).
Rewritten to use NET_RECEIVE for compatibility with NetCon/VecStim,
replacing the old POINTER pre / voltage-threshold approach.

Kinetics:
  dR/dt = K1 * T * (1-R) - K2 * R
  dG/dt = K3 * R - K4 * G
  Open fraction: O = G^n / (G^n + KD)
  i = gmax * O * (v - Erev)

Written by Destexhe, Laval University, 1995.
NET_RECEIVE adaptation: compatible with Synapse() / VecStim.
ENDCOMMENT

NEURON {
    POINT_PROCESS GABAb
    RANGE R, G, g, gmax, lastrelease
    NONSPECIFIC_CURRENT i
    GLOBAL Cmax, Cdur, Deadtime
    GLOBAL K1, K2, K3, K4, KD, n, Erev
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM)  = (milli/liter)
}

PARAMETER {
    Cmax     = 1      (mM)   : max transmitter concentration
    Cdur     = 1      (ms)   : transmitter pulse duration
    Deadtime = 1      (ms)   : refractory period between events

    K1  = 0.09   (/ms mM)   : receptor binding rate
    K2  = 0.0012 (/ms)      : receptor unbinding rate
    K3  = 0.18   (/ms)      : G-protein production rate
    K4  = 0.034  (/ms)      : G-protein decay rate
    KD  = 100               : K+ channel dissociation constant
    n   = 4                 : G-protein binding cooperativity
    Erev = -95    (mV)      : reversal potential (E_K)
    gmax = 0      (umho)    : maximum conductance (set by user)
}

ASSIGNED {
    v           (mV)
    i           (nA)
    g           (umho)
    C           (mM)        : current transmitter concentration
    Gn                      : G^n
    lastrelease (ms)        : time of last NET_RECEIVE event
    on                      : 1 = transmitter being released
}

STATE {
    R                       : fraction of activated receptors
    G                       : fraction of activated G-protein
}

INITIAL {
    R           = 0
    G           = 0
    C           = 0
    on          = 0
    lastrelease = -1e9
}

BREAKPOINT {
    SOLVE bindkin METHOD cnexp
    Gn = G^n
    g  = gmax * Gn / (Gn + KD)
    i  = g * (v - Erev)
}

DERIVATIVE bindkin {
    : transmitter pulse: on for Cdur ms after each event
    if (on && (t - lastrelease) > Cdur) {
        C  = 0
        on = 0
    }
    R' = K1 * C * (1 - R) - K2 * R
    G' = K3 * R - K4 * G
}

NET_RECEIVE(weight) {
    : weight is ignored - gmax controls amplitude
    : only fire if outside the dead time
    if ((t - lastrelease) > (Cdur + Deadtime)) {
        C           = Cmax
        on          = 1
        lastrelease = t
    }
}
