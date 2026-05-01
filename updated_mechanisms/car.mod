: Updated for CVODE compatibility on 2026-04-27
TITLE Ca R-type channel with medium threshold for activation
: used in distal dendritic regions, together with calH.mod, to help
: the generation of Ca++ spikes in these regions
: uses channel conductance (not permeability)
: written by Yiota Poirazi on 11/13/00 poirazi@LNC.usc.edu

NEURON {
	SUFFIX car
	USEION ca READ eca WRITE ica
        RANGE gcabar, m, h
	RANGE inf, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {              : parameters that can be entered when function is called in cell-setup
        v               (mV)
        celsius = 34	(degC)
        gcabar = 0      (mho/cm2) : initialized conductance
	eca = 140       (mV)      : Ca++ reversal potential
        }  

STATE {	m h }            : unknown activation and inactivation parameters to be solved in the DEs  

ASSIGNED {               : parameters needed to solve DE
	ica (mA/cm2)
        inf[2]
	tau[2]
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar*m*m*m*h*(v - eca)
	}

INITIAL {
        mhn(v)
        m = inf[0]    : initial activation parameter value
	h = inf[1]    : initial inactivation parameter value
	ica = gcabar*m*m*m*h*(v - eca) : initial Ca++ current value
        }

DERIVATIVE states {
        mhn(v)
        m' = (inf[0] - m) / tau[0]
        h' = (inf[1] - h) / tau[1]
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

PROCEDURE mhn(v) {
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
	}
}
