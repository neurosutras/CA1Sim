: Ca channels (T,N,L-type)


NEURON {
	SUFFIX Cacum
	USEION ca READ ica WRITE cai
	RANGE depth, tau, taucadiv, cai0
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	F = (faraday) (coulomb)
}

PARAMETER {
	tau = 9			(ms)
	taucadiv = 1
	depth = 200 	(nm)
	cai0 = 7e-5 	(mM)
}

ASSIGNED {
	v		(mV)
    eca     (mV)
	ica	(mA/cm2)
    cao     (mM)
	celsius	(degC)
}

STATE {
	cai (mM) <1e-5>
}

BREAKPOINT {
    SOLVE state METHOD derivimplicit
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	cai' = -(1e7)*ica/2./F/depth - taucadiv*(cai-cai0)/tau
}

INITIAL {
	cai = cai0
}
