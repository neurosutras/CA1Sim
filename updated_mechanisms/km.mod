: Updated for CVODE compatibility on 2026-04-27
: Modified for CVODE compatibility on 2026-04-27
COMMENT
km.mod
Potassium channel, Hodgkin-Huxley style kinetics
Based on I-M (muscarinic K channel)
Slow, noninactivating
Author: Zach Mainen, Salk Institute, 1995, zach@salk.edu
	
ENDCOMMENT


NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE n, gk, gbar
	RANGE ninf, ntau
	GLOBAL Ra, Rb
	GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	v 		(mV)
	gbar = 10   	(pS/um2)	: 0.03 mho/cm2
	tha  = -30	(mV)		: v 1/2 for inf
	qa   = 9	(mV)		: inf slope		
	Ra   = 0.001	(/ms)		: max act rate  (slow)
	Rb   = 0.001	(/ms)		: max deact rate  (slow)
	celsius		(degC)
	temp = 23	(degC)		: original temp 	
	q10  = 2.3			: temperature sensitivity
	vmin = -120	(mV)
	vmax = 100	(mV)
} 


ASSIGNED {
	a		(/ms)
	b		(/ms)
	ik 		(mA/cm2)
	gk		(pS/um2)
	ek		(mV)
	ninf
	ntau (ms)	
	tadj
}
 

STATE { n }

INITIAL { 
	trates(v)
	n = ninf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = tadj*gbar*n
	ik = (1e-4) * gk * (v - ek)
} 

DERIVATIVE states {
        trates(v)
        n' = (ninf - n)/ntau
}

PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                       :Call once from HOC to initialize inf at resting v.
        TABLE ninf, ntau
	DEPEND celsius, temp, Ra, Rb, tha, qa
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable_hh == 1
        tadj = q10^((celsius - temp)/10)  :temperature adjastment
}


PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

        a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
        b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
        ntau = 1/(a+b)
	ninf = a*ntau
}
