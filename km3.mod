TITLE CA1 KM channel from Mala Shah
: M. Migliore June 2006
: faster activation kinetics based on Chen & Johnston, 2004.

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkmbar 	= .0001 	(mho/cm2)
    vhalfl 	= -42   	(mV)
	kl 		= -4
    vhalft 	= -42   	(mV)
    a0t 	= 0.009     (/ms)
    zetat 	= 7    		(1)
    gmt 	= .4   		(1)
	q10 	= 5
	t0 		= 15
	st 		= 1
}

NEURON {
	SUFFIX km3
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
}

INITIAL {
	rate(v)
	m=inf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gk = gkmbar*m^st
	ik = gk*m^st*(v-ek)
}

FUNCTION alpt(v(mV)) {
    alpt = exp(0.0378*zetat*(v-vhalft))
}

FUNCTION bett(v(mV)) {
    bett = exp(0.0378*zetat*gmt*(v-vhalft))
}

DERIVATIVE state {
    rate(v)
	m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
    LOCAL a,qt
    qt=q10^((celsius-35)/10)
    inf = (1/(1 + exp((v-vhalfl)/kl)))
    a = alpt(v)
    tau = t0 + bett(v)/(a0t*(1+a))
}
