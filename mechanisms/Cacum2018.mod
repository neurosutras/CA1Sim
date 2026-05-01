: Ca channels (T,N,L-type)


NEURON {
	SUFFIX Cacum
	USEION ca READ ica, cai WRITE cao, cai
	RANGE depth, tau, cai0, cao0
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	F = (faraday) (coulomb)
}

PARAMETER {
	tau = 100		(ms)
	depth = 0.2 	(um)  : stuck with depth from Aradi (1999)
	cai0 = 7e-5 	(mM)
 	cao0 = 1.3 		(mM)
}

ASSIGNED {
	v			(mV)
	ica			(mA/cm2)
 	diam 		(um)
	VSR 		(um)
 	B 			(mM*cm2/mA)
}

STATE {
	cai (mM) <1e-5>
 	cao (mM)
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	cai' = -ica * B - (cai-cai0)/tau
 	cao' = 0.
}

INITIAL {
	cai = cai0
 	cao = cao0
	: From Beining (2017)
	if (2. * depth >= diam) {
		VSR = 0.25 * diam : if diameter gets less than double the depth,
						: the surface to volume ratio (here volume to surface ratio VSR)
						: cannot be less than 1/(0.25diam) (instead of diam/(d*(diam-d)) )
	} else{
		VSR = depth*(1-depth/diam)
	}
	B = (1e4)/(2*F*VSR)
}
