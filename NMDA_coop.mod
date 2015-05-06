COMMENT

Milstein 2015:

Adapted exp2syn to generate a conductance with tunable saturation / cooperativity.

--------------------------------------------------------------------------------

Two state kinetic scheme synapse described by rise time tau1, and decay time
constant tau2. The normalized peak conductance for a single pulse is gmax, but
it can grow during repeated activation depending on the value of a cooperativity
variable.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The normalization factor required to normalize the peak conductance of a
single event can be calculated with a root finding algorithm, or by simply
measuring the peak value of the gating variable B. The value of norm can also
be set via hoc or python.


ENDCOMMENT

NEURON {
	POINT_PROCESS AMPA_coop
	RANGE tau1, tau2, gmax, g, norm, max_coop
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	tau1 = 0.05     (ms)    <1e-9,1e9>
	tau2 = 3.16     (ms)    <1e-9,1e9>
    gmax = 0.0005   (umho)  <0,1>       : peak conductance of single EPSP, depends on norm
    e = 0           (mV)
    norm = 0.0468                       : peak of gating variable B, must be recalculated for each new set of time
                                        : constants to normalize the peak conductance of a single EPSP to gmax
    max_coop = 0            <0, 1>      : max cooperativity (degree of post-synaptic facilitation during a train)
}

ASSIGNED {
	v (mV)
	i (nA)
    g (umho)
    coop
    scale
}

STATE {
	A
	B
}

INITIAL {
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
    coop = 0
    scale = 1
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gmax * scale * B / norm
    i = g * (v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2 + A
}

NET_RECEIVE(weight) {
    if (B / norm / (1 + max_coop) < 0.05) {
        coop = 0.
    } else {
        coop = B / norm * max_coop / (1 + max_coop)
    }
    A = 1 + coop - B / norm
    scale = weight
}