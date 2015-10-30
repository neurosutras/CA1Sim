TITLE nax
: Na current for axon. No slow inact.
: M.Migliore Jul. 1997
: added sh to account for higher threshold M.Migliore, Apr.2002
: Aaron Milstein modified October 2015, added additional sha that applies only to activation threshold

NEURON {
	THREADSAFE
    SUFFIX nax
	USEION na READ ena WRITE ina
	RANGE  gbar, sh, sha, minf, hinf, mtau, htau
	GLOBAL thinf, qinf
}

PARAMETER {
	sh    = 0	        (mV)
    sha   = 0           (mV)
	gbar  = 0.010   	(mho/cm2)
								
	tha   =  -30	    (mV)		: act vhalf
	qa    = 7.2	        (/mV)		: act slope (4.5)
	Ra    = 0.4	        (/ms)		: open (v)
	Rb    = 0.124 	    (/ms)		: close (v)

	thi1  = -45	        (mV)		: v 1/2 for inact
	thi2  = -45 	    (mV)		: v 1/2 for inact
	qd    = 1.5	        (/mV)	    : inact tau slope
	qg    = 1.5         (/mV)
	mmin  = 0.02
	hmin  = 0.5
	q10   = 2
	Rg    = 0.01 	    (/ms)		: inact recov (v)
	Rd    = .03 	    (/ms)		: inact (v)

	thinf  = -50 	    (mV)		: inact inf vhalf
	qinf   = 4 	        (/mV)		: inact inf slope
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ena		        (mV)
	celsius
	v 		        (mV)
    ina 		    (mA/cm2)
	thegna		    (mho/cm2)
	minf
    hinf
	mtau            (ms)
    htau            (ms)
}
 

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
} 

INITIAL {
	trates(v,sh,sha)
	m=minf  
	h=hinf
    thegna = gbar*m*m*m*h
	ina = thegna * (v - ena)
}

DERIVATIVE states {   
    trates(v,sh,sha)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
}

PROCEDURE trates(vm,sh2,sha2) {
    LOCAL  a, b, qt
    qt=q10^((celsius-24)/10)
	a = trap0(vm,tha+sh2+sha2,Ra,qa)
	b = trap0(-vm,-tha-sh2-sha2,Rb,qa)
	mtau = 1/(a+b)/qt
    if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = trap0(vm,thi1+sh2,Rd,qd)
	b = trap0(-vm,-thi2-sh2,Rg,qg)
    :a = trap0(vm,thi1,Rd,qd)
	:b = trap0(-vm,-thi2,Rg,qg)
	htau =  1/(a+b)/qt
    if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf-sh2)/qinf))
    :hinf = 1/(1+exp((vm-thinf)/qinf))
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
        trap0 = a * q
 	}
}