TITLE CA1 KM channel from Mala Shah
: M. Migliore June 2006
: option to have faster activation than inactivation kinetics based on Chen & Johnston, 2004.

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkmbar 	= .0001 	(mho/cm2)
    vhalfl 	= -40   	(mV)
	kl 		= -7
    vhalft 	= -42   	(mV)
    a0t_f 	= 0.009     (/ms)
	a0t_s 	= 0.036
    zetat 	= 7    		(1)
    gmt 	= .4   		(1)
	q10 	= 5
	t0_f 	= 15
	t0_s 	= 60
	st 		= 1
}

NEURON {
	SUFFIX km2
	USEION k READ ek WRITE ik
    RANGE  gkmbar,ik,gk,inf,tau
}

STATE {
    m
}

ASSIGNED {
    v 	    (mV)
	ek
	celsius (degC)
	ik      (mA/cm2)
    gk 		(S/cm2)
    inf
	tau
    tau_f
	tau_s
}

INITIAL {
	rate(v)
	m=inf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gk = gkmbar*m^st
	ik = gk*(v-ek)
}

FUNCTION alpt(v(mV)) {
    alpt = exp(0.0378*zetat*(v-vhalft))
}

FUNCTION bett(v(mV)) {
    bett = exp(0.0378*zetat*gmt*(v-vhalft))
}

DERIVATIVE state {
    rate(v)
	if (m<inf) {tau=tau_s} else {tau=tau_f}
    :if (m<inf) {tau=tau_f} else {tau=tau_s}
	:tau=tau_s
	:tau=tau_f
	m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
    LOCAL a,qt
    qt=q10^((celsius-35)/10)
    inf = (1/(1 + exp((v-vhalfl)/kl)))
    a = alpt(v)
    tau_f = t0_f + bett(v)/(a0t_f*(1+a))
    tau_s = t0_s + bett(v)/(a0t_s*(1+a))
}
