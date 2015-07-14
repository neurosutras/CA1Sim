TITLE I-h channel from Magee 1998 for distal dendrites

NEURON {
	SUFFIX h
	NONSPECIFIC_CURRENT i
    RANGE ghbar, vhalfl, eh
    RANGE linf,taul
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    eh=-30.		    (mV)
	ghbar=.0001 	(mho/cm2)
    vhalfl=-90   	(mV)
    vhalft=-75   	(mV)
    a0t=0.011      	(/ms)
    zetal=4    	    (1)
    zetat=2.2    	(1)
    gmt=.4   	    (1)
	q10=4.5
	qtl=1
}

STATE {
    l
}

ASSIGNED {
	v 	    (mV)
    i       (mA/cm2)
    linf
    taul
    g
    celsius (degC)
}

INITIAL {
	rate(v)
	l=linf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = ghbar*l
	i = g*(v-eh)

}

DERIVATIVE states {
    rate(v)
    l' =  (linf - l)/taul
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION alpt(v(mV)) {
  alpt = exp(1e-3*zetat*(v-vhalft)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bett(v(mV)) {
  bett = exp(1e-3*zetat*gmt*(v-vhalft)*9.648e4/(8.315*(273.16+celsius)))
}

PROCEDURE rate(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-33)/10)
        a = alpt(v)
        linf = 1/(1+ alpl(v))
        taul = bett(v)/(qtl*qt*a0t*(1+a))
}