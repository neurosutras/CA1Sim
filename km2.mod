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
    a0t_a 	= 0.009     (/ms)
	a0t_d 	= 0.036
    zetat 	= 7    		(1)
    gmt 	= .4   		(1)
	q10 	= 5
	t0_a 	= 15
	t0_d 	= 60
	st 		= 1
}

NEURON {
	SUFFIX km2
	USEION k READ ek WRITE ik
    RANGE  gkmbar,ik,inf,tau
}

STATE {
    m
}

ASSIGNED {
    v 	    (mV)
	ek
	celsius (degC)
	ik      (mA/cm2)
    inf
	tau
    tau_a
	tau_d
}

INITIAL {
	rate(v)
	m=inf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gkmbar*m^st*(v-ek)
}

FUNCTION alpt(v(mV)) {
    alpt = exp(0.0378*zetat*(v-vhalft))
}

FUNCTION bett(v(mV)) {
    bett = exp(0.0378*zetat*gmt*(v-vhalft))
}

DERIVATIVE state {
    rate(v)
	:if (m<inf) {tau=tau_d} else {tau=tau_a}
    :if (m<inf) {tau=tau_a} else {tau=tau_d}
	tau=tau_d
	m' = (inf - m)/tau
}

PROCEDURE rate(v (mV)) { :callable from hoc
    LOCAL a,qt
    qt=q10^((celsius-35)/10)
    inf = (1/(1 + exp((v-vhalfl)/kl)))
    a = alpt(v)
    :tau_a = t0_a + bett(v)/(a0t_a*(1+a))
    tau_d = t0_d + bett(v)/(a0t_d*(1+a))
}
