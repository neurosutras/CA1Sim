TITLE t-type calcium channel with high threshold for activation
: used in somatic and dendritic regions 
: it calculates I_Ca using channel permeability instead of conductance

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {           :parameters that can be entered when function is called in cell-setup 
        dt            (ms)
	v             (mV)
        tBase = 23.5  (degC)
	celsius = 22  (degC)
	gcatbar = 0   (mho/cm2)  : initialized conductance
	ki = 0.001    (mM)
	cai = 5.e-5   (mM)       : initial internal Ca++ concentration
	cao = 2       (mM)       : initial external Ca++ concentration
        tfa = 1                  : activation time constant scaling factor
        tfi = 0.68               : inactivation time constant scaling factor
        eca = 140                : Ca++ reversal potential
}

NEURON {
	SUFFIX cat
	USEION ca READ cai,cao 
        USEION Ca WRITE iCa VALENCE 2
        : The T-current does not activate calcium-dependent currents.
        : The construction with dummy ion Ca prevents the updating of the 
        : internal calcium concentration. 
        RANGE gcatbar, hinf, minf, taum, tauh, iCa
}

STATE {	m h }  : unknown activation and inactivation parameters to be solved in the DEs 

ASSIGNED {     : parameters needed to solve DE
	iCa (mA/cm2)
        gcat  (mho/cm2) 
        minf
        hinf
        taum
        tauh
}

INITIAL {
:        tadj = 3^((celsius-tBase)/10)   : assume Q10 of 3
	rates(v)
        m = minf
        h = hinf
	gcat = gcatbar*m*m*h*h2(cai)
}

BREAKPOINT {
	SOLVE states
	gcat = gcatbar*m*m*h*h2(cai) : maximum channel permeability
	iCa = gcat*ghk(v,cai,cao)    : dummy calcium current induced by this channel

}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) { LOCAL nu,f
        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {   : temperature-dependent adjustment factor
        KTF = ((25./293.15)*(celsius + 273.15))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alph(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alph = 1.6e-4*exp(-(v+57)/19)
}

FUNCTION beth(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
	beth = 1/(exp((-v+15)/10)+1.0)
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	betm = 0.046*exp(-v/22.73)
}

UNITSON
LOCAL facm,fach

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE states() {     : exact when v held constant; integrates over dt step
        rates(v)
        m = m + facm*(minf - m)
        h = h + fach*(hinf - h)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a
        a = alpm(v)
        taum = 1/(tfa*(a + betm(v))) : estimation of activation tau
        minf =  a/(a+betm(v))        : estimation of activation steady state
        facm = (1 - exp(-dt/taum))
        a = alph(v)
        tauh = 1/(tfi*(a + beth(v))) : estimation of inactivation tau
        hinf = a/(a+beth(v))         : estimation of inactivation steady state
        fach = (1 - exp(-dt/tauh))
}
