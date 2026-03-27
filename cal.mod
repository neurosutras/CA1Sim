TITLE L-type calcium channel with low threshold for activation
: used in somatic and proximal dendritic regions 
: it calculates I_Ca using channel permeability instead of conductance

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {		:parameters that can be entered when function is called in cell-setup 
        dt              (ms)
	v               (mV)
	celsius = 34	(degC)
	gcalbar = 0     (mho/cm2) : initialized conductance
	ki  = 0.001     (mM)  
	cai = 5.e-5     (mM)      : initial internal Ca++ concentration
	cao = 2         (mM)      : initial external Ca++ concentration
        tfa = 5                   : time constant scaling factor
        eca = 140                 : Ca++ reversal potential
}

NEURON {
	SUFFIX cal
	USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf,taum
}

STATE {	m }                      : unknown parameter to be solved in the DEs 

ASSIGNED {                       : parameters needed to solve DE
	ica (mA/cm2)
        gcal  (mho/cm2) 
        minf
        taum
        }

INITIAL {                        : initialize the following parameter using rates()
        rates(v)
        m = minf
	gcal = gcalbar*m*h2(cai)
}

BREAKPOINT {
	SOLVE states
	gcal = gcalbar*m*h2(cai) : maximum channel permeability
	ica = gcal*ghk(v,cai,cao): calcium current induced by this channel
}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f
        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) { : temperature-dependent adjustment factor
        KTF = ((25./293.15)*(celsius + 273.15))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = 0.055*(-27.01 - v)/(exp((-27.01-v)/3.8) - 1)
}


FUNCTION betm(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
        betm =0.94*exp((-63.01-v)/17)
}

UNITSON
LOCAL facm
:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
:only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        m = m + facm*(minf - m)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a
        a = alpm(v)
        taum = 1/(tfa*(a+betm(v))) : estimation of activation tau
        minf = a/(a+betm(v))       : estimation of activation steady state value
        facm = (1 - exp(-dt/taum))
}



