COMMENT

Milstein 2015:

Adapted exp2syn to generate a conductance with a level of (approximate)
saturation that can be modulated with a cooperativity variable.

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

The factor is evaluated in the initial block such that an event of weight 1
generates a peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS AMPA_coop
	RANGE tau1, tau2, gmax, max_coop
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	tau1 = 0.05      (ms)   <1e-9,1e9>
	tau2 = 3.16      (ms)   <1e-9,1e9>
    gmax = 0.0005   (umho)  <0,1>       : to scale amplitude
    e = 0           (mV)
    max_coop = 1            <1, 1e9>    : max cooperativity. not exact. should be tuned to match target level of
                                        : peak facilitation at a given stimulation frequency and set of rate constants
}

ASSIGNED {
	v (mV)
	i (nA)
    g (umho)
	tp
    factor
    inc
    current_coop
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
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
    inc = 0
    current_coop = 1
    scale = 1
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gmax * scale * (B - A)
    i = g * (v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight) {
    scale = weight
    if (((B - A) / current_coop) < 0.05) {
        current_coop = 1
    } else {
        current_coop = (max_coop - 1) * (B - A) / current_coop + 1
        if (current_coop > max_coop) {
            current_coop = max_coop
        }
    }
    inc = factor * (current_coop - (B * exp(-tp/tau2) - A * exp(-tp/tau1)))
    A = A + inc
    B = B + inc
}