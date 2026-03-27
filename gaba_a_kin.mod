COMMENT
-----------------------------------------------------------------------------

GABA-A Receptors

-----------------------------------------------------------------------------
Synaptic mechanism based on a simplified model of transmitter binding to
postsynaptic receptors.

Written by Aaron Milstein, 2015.

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
	POINT_PROCESS GABA_A_KIN
	RANGE Cmax, Cdur, kon, koff, CC, CO, Beta, Alpha, Erev, gmax, g, Ro
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax    = 1.        (mM)    	    : transmitter concentration during release event
	Cdur	= 0.5       (ms)		    : transmitter duration (rising phase)
	kon     = 5.397     (/ms/mM)        : unbound receptor ligand-binding rate
    koff    = 4.433     (/ms)           : bound receptor ligand-unbinding rate
    CC      = 20.945    (/ms)           : bound receptor cleft closing rate
    CO      = 1.233     (/ms)           : bound receptor cleft opening rate
    Beta	= 283.090   (/ms)	        : channel opening rate
    Alpha   = 254.520   (/ms)           : open channel closing rate
	Erev	= -73.0     (mV)		    : reversal potential (-70.0)
	gmax	= 0.000590  (umho)	        : maximum conductance [0.000603, 0.000492, 0.000590]
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g * (v - Erev)
	g 		(umho)		: conductance
	C                   : unbound transmitter concentration
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
}

BREAKPOINT {
    SOLVE kstates METHOD sparse
	g = scale * gmax * Ro
    i = g * (v - Erev)
}

KINETIC kstates {
    ~ Ru <-> Rb     (C * kon, koff)
    ~ Rb <-> Rc     (CC, CO)
    ~ Rc <-> Ro     (Beta, Alpha)
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