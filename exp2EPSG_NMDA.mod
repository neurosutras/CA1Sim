COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS Exp2EPSG_NMDA
	RANGE tau1, tau2, e, gmax, Kd, gamma, mg, M, i
	NONSPECIFIC_CURRENT i

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1    = .1    (ms) <1e-9,1e9>
	tau2    = 10    (ms) <1e-9,1e9>
    Kd      = 9.98  (mM)    : modulate Mg concentration dependence
    gamma   = 0.101 (/mV)   : modulate slope of Mg sensitivity
    mg      = 1.0   (mM)    : extracellular Mg concentration
    e       = 0	    (mV)    : reversal potential
    gmax    = 1     (1)     : peak conductance
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
    M       : fraction of channels not blocked by extracellular Mg
	factor
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
    mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
    mgblock(v)
	g = B - A
	i = g*gmax*M*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}

PROCEDURE mgblock(v(mV)) {
	: from Jahr & Stevens
    M = 1. / (1. + exp(gamma * (-v)) * (mg / Kd))
}
