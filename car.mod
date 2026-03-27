TITLE Ca R-type channel with medium threshold for activation
: used in distal dendritic regions, together with calH.mod, to help
: the generation of Ca++ spikes in these regions
: uses channel conductance (not permeability)
: written by Yiota Poirazi on 11/13/00 poirazi@LNC.usc.edu

NEURON {
	SUFFIX car
	USEION ca READ eca WRITE ica
        RANGE gcabar, m, h
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {              : parameters that can be entered when function is called in cell-setup
        v               (mV)
        celsius = 34	(degC)
	dt              (ms)
        gcabar = 0      (mho/cm2) : initialized conductance
	eca = 140       (mV)      : Ca++ reversal potential
        }  

STATE {	m h }            : unknown activation and inactivation parameters to be solved in the DEs  

ASSIGNED {               : parameters needed to solve DE
	ica (mA/cm2)
        inf[2]
	fac[2]
	tau[2]
}

BREAKPOINT {
	SOLVE states
	ica = gcabar*m*m*m*h*(v - eca)
	}

INITIAL {
        m = 0    : initial activation parameter value
	h = 1    : initial inactivation parameter value
	states()
	ica = gcabar*m*m*m*h*(v - eca) : initial Ca++ current value
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
	    varss = 1 / (1 + exp((v+48.5)/(-3))) : Ca activation
	}
	else if (i==1) {
             varss = 1/ (1 + exp((v+53)/(1)))    : Ca inactivation
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
           vartau = 50  : activation variable time constant
        }
	else if (i==1) {
           vartau = 5   : inactivation variable time constant
       }
	
}	

PROCEDURE mhn(v) {LOCAL a, b :rest = -70
:	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}















