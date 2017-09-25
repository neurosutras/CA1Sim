COMMENT

Adapted exp2syn to generate a current rather than a conductance (Milstein 2015):

--------------------------------------------------------------------------------

Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak current is 1.
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
	POINT_PROCESS EPSC
	RANGE tau1, tau2, i, inward, imax
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	tau1 = 0.2   (ms)    <1e-9,1e9>
	tau2 = 5.0    (ms)   <1e-9,1e9>
    inward = -1.0                     : for inward currents
    imax = 1.0           <0,1>        : to scale amplitude
}

ASSIGNED {
	v (mV)
	i (nA)
	factor
}

STATE {
	A (nA)
	B (nA)
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
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	i = inward * imax * (B - A)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (nA)) {
	A = A + weight*factor
    B = B + weight*factor
}
