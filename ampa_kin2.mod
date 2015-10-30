COMMENT
-----------------------------------------------------------------------------

GluA-R (AMPA-type Glutamate Receptors)

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

       ---[C] * kon-->        -----Beta----->
C + Ru <-----koff-----   Rb   <----Alpha-----   Ro

where Ru, Rb, and Ro are respectively the fraction of channels in the
unbound, closed bound, and open states of the postsynaptic
receptor. kon and koff are the binding and unbinding rate constants,
and Beta and Alpha are the opening and closing rates of the channel.

The maximal conductance of the channel can be set for each instance of the
synaptic mechanism by specifying gmax, and the relative weight of events
can be set through the weight parameter of the associated netcon object.

The postsynaptic current is given by:

i = weight * gmax * Ro * (V-Erev)

-----------------------------------------------------------------------------

ENDCOMMENT

NEURON {
	POINT_PROCESS AMPA_KIN2
	RANGE Cmax, Cdur, kon, koff, Beta, Alpha, Erev, gmax, g
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
	Cdur	= 0.3       (ms)		    : transmitter duration (rising phase)
	kon     = 62.88     (/ms/mM)        : unbound receptor ligand-binding rate
    koff    = 16.63     (/ms)           : bound receptor ligand-unbinding rate
    Beta	= 20.53     (/ms)	        : channel opening rate
    Alpha   = 0.71      (/ms)           : open channel closing rate
	Erev	= 0.        (mV)		    : reversal potential
	gmax	= 0.001     (umho)	        : maximum conductance
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
    Ro                  : fraction of channels in open state
}

INITIAL {
	C = 0.
    Ru = 1.
    Rb = 0.
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
    ~ Rb <-> Ro     (Beta, Alpha)
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