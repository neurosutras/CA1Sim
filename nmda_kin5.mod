COMMENT
-----------------------------------------------------------------------------

GluN-R (NMDA-type Glutamate Receptors)

Postsynaptic current is additionally constrained by voltage-dependent channel
block by extracellular Mg.

-----------------------------------------------------------------------------
Synaptic mechanism based on a simplified model of transmitter binding to
postsynaptic receptors.

Modified by Aaron Milstein, 2015.

Modification of original code by:

A. Destexhe & Z. Mainen, The Salk Institute, March 12, 1993.
Last modif. Sept 8, 1993.

Reference:

Destexhe, A., Mainen, Z. and Sejnowski, T.J.  An efficient method for
computing synaptic conductances based on a kinetic model of receptor binding.
Neural Computation, 6: 14-18, 1994.

-----------------------------------------------------------------------------

Upon arrival of a presynaptic spike (a net_event), the concentration of the
neurotransmitter C in the synaptic cleft is briefly stepped to concentration
Cmax for duration Cdur.

C     _____ . . . . . . Cmax
      |   |
 _____|   |______ . . . 0
     t0   t0 + Cdur

The receptors then bind transmitter, change conformation, and open their
channel according to the following kinetic scheme:

       ---[C] * kon-->        -------CC----->        -----Beta----->
C + Ru <-----koff-----   Rb   <------CO------   Rc   <----Alpha-----   Ro

where Ru, Rb, Rc, and Ro are respectively the fraction of channels in the
unbound, closed bound, closed cleft, and open states of the postsynaptic
receptor. kon and koff are the binding and unbinding rate constants, CC and
CO are the closing and opening rates of the receptor ligand-binding cleft,
and Beta and Alpha are the opening and closing rates of the channel.

The maximal conductance of the channel can be set for each instance of the
synaptic mechanism by specifying gmax, and the relative weight of events
can be set through the weight parameter of the associated netcon object.

The postsynaptic current is given by:

i = weight * gmax * Ro * (V-Erev)

-----------------------------------------------------------------------------

ENDCOMMENT

NEURON {
	POINT_PROCESS NMDA_KIN5
	RANGE Cmax, Cdur, kon, koff, CC, CO, Beta, Alpha, Erev, Kd, gamma, mg, gmax, g, B, kin_scale, Ro
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax        = 1.        (mM)        : transmitter concentration during release event
	Cdur	    = 0.3       (ms)		: transmitter duration (rising phase)
	kon         = 86.89     (/mM/ms)    : unbound receptor ligand-binding rate
    koff        = 0.69      (/ms)       : bound receptor ligand-unbinding rate
    CC          = 9.64      (/ms)       : bound receptor cleft closing rate
    CO          = 2.60      (/ms)       : bound receptor cleft opening rate
    Beta	    = 0.68      (/ms)	    : channel opening rate
    Alpha       = 0.079     (/ms)       : open channel closing rate
	Erev	    = 0.        (mV)		: reversal potential
	Kd          = 9.98      (mM)        : modulate Mg concentration dependence
    gamma       = 0.101     (/mV)       : modulate slope of Mg sensitivity
    mg          = 1.0       (mM)        : extracellular Mg concentration
    gmax	    = 0.003026  (umho)	    : maximum conductance
    kin_scale   = 1.83      (1)         : scale voltage sensitivity of decay kinetics
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g * (v - Erev)
	g 		(umho)		: conductance
	C       (mM)        : unbound transmitter concentration
    B                   : fraction of channels not blocked by extracellular Mg
    kin_factor          : exponential scaling of decay rates with voltage
    scale               : allow netcon weight to scale conductance
}

STATE {
    Ru                  : fraction of receptors not bound to transmitter
    Rb                  : fraction of receptors bound to transmitter
    Rc                  : fraction of receptors in closed cleft state
    Ro                  : fraction of channels in open state
}

INITIAL {
	C = 0.
    Ru = 1.
    Rb = 0.
    Rc = 0.
    Ro = 0.
    scale = 1.
    mgblock(v)
}

BREAKPOINT {
    mgblock(v)
    SOLVE kstates METHOD sparse
	g = scale * gmax * B * Ro
    i = g * (v - Erev)
}

KINETIC kstates {
    ~ Ru <-> Rb     (C * kon, koff)
    ~ Rb <-> Rc     (CC, CO * kin_factor)
    ~ Rc <-> Ro     (Beta, Alpha * kin_factor)
}

NET_RECEIVE(weight) {
    if (flag == 0) { : a new spike received
        C = Cmax
        scale = weight
        net_send(Cdur, 1)
    } else {    : a self event
        C = 0.
    }
}

PROCEDURE mgblock(v(mV)) {
	: from Jahr & Stevens
    B = 1. / (1. + exp(gamma * (-v)) * (mg / Kd))
    kin_factor = (1 - kin_scale) * B + kin_scale
}