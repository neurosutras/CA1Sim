TITLE Ca L-type channel with high treshold of activation
: inserted in distal dendrites to account for distally
: restricted initiation of Ca++ spikes
: uses channel conductance (not permeability)
: written by Yiota Poirazi, 1/8/00 poirazi@LNC.usc.edu

NEURON {
	SUFFIX calH
	USEION ca READ eca WRITE ica
        RANGE gcalbar, m, h
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {          : parameters that can be entered when function is called in cell-setup
        v               (mV)
        celsius = 34	(degC)
	dt              (ms)
        gcalbar = 0     (mho/cm2) : initialized conductance
	eca = 140       (mV)      : Ca++ reversal potential
        }

STATE {	m h }                     : unknown activation and inactivation parameters to be solved in the DEs  

ASSIGNED {                        : parameters needed to solve DE
	ica (mA/cm2)
        inf[2]
	fac[2]
	tau[2]
}

BREAKPOINT {
	SOLVE states
	ica = gcalbar*m*m*m*h*(v - eca)       
	}

INITIAL {
        m = 0    : initial activation parameter value
	h = 1    : initial inactivation parameter value
        states()
	ica = gcalbar*m*m*m*h*(v - eca) : initial Ca++ current value
     	}

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	}	

PROCEDURE states() {	: exact when v held constant
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}

FUNCTION varss(v, i) {
	if (i==0) { 
             varss = 1 / (1 + exp((v+37)/(-1)))  : Ca activation 
	}
	else if (i==1) { 
             varss = 1 / (1 + exp((v+41)/(0.5))) : Ca inactivation 
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
           vartau = 3.6  : activation variable time constant
        }
	else if (i==1) {
:           vartau = 25   : inactivation variable time constant
           vartau = 29   : inactivation variable time constant
        }
}	

PROCEDURE mhn(v) {LOCAL a, b :rest = -70
:      TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}















