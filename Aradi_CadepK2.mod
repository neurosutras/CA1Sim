: Ca-dependent K channels (BK and SK)


NEURON {
	SUFFIX CadepK
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gbkbar, gskbar, gbar, i, ask, bsk, gsk, gbk, isk, ibk, gcakmult
	GLOBAL stau, tauskdiv
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gbkbar = .0003	(S/cm2)	: maximum permeability
	gskbar = .0005	(S/cm2)	: maximum permeability
	gcakmult = 1.
	tauskdiv = 1
	ask = 1
	bsk = 1
	alphar = 7.5	(/ms)
	stau = 10		(ms)
}

ASSIGNED {
	v			(mV)
	ek		(mV)
	ik		(mA/cm2)
	isk		(mA/cm2)
	ibk		(mA/cm2)
	i 		(mA/cm2)
	cai		(mM)
	area	(microm2)
  	gbk		(S/cm2)
  	gsk		(S/cm2)
  	gbar  (S/cm2)
}

STATE {
	q 
	r 
	s 
}

BREAKPOINT {
    SOLVE state METHOD cnexp
	gbk = gbkbar*gcakmult*r*s*s
	gsk = gskbar*gcakmult*q*q
	isk = gsk*(v - ek)
	ibk = gbk*(v - ek)
	ik = isk + ibk
	i = ik
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	q' = tauskdiv*(ask*alphaq(cai)*(1-q)-bsk*betaq(cai)*q)
	r' = alphar*(1-r)-betar(v)*r
	s' = (sinf(cai)-s)/stau
}

INITIAL {
	q = alphaq(cai)/(alphaq(cai)+betaq(cai))
	r = alphar/(alphar+betar(v))
  	s = sinf(cai)
	gbk = gbkbar*gcakmult*r*s*s
	gsk = gskbar*gcakmult*q*q
	isk = gsk*(v - ek)
	ibk = gbk*(v - ek)
	ik = isk + ibk
	i = ik
  	gbar = gcakmult * (gbkbar + gskbar)
}

FUNCTION exp1(A (/ms), d, k, x (mM)) (/ms) {
	if (x > 1e-7) {
		exp1 = A/exp((12*log10(x)+d)/k)
	} else {
		exp1 = A/exp((12*(-7.)+d)/k)
	}
}

FUNCTION alphaq(x (mM)) (/ms) {
	alphaq = exp1(0.00246,28.48,-4.5,x)	:28
}

FUNCTION betaq(x (mM)) (/ms) {
	betaq = exp1(0.006,60.4,35,x)
}

FUNCTION betar(v (mV)) (/ms) {
	betar = 0.11/exp((v-35)/14.9)
}

FUNCTION sinf(x (mM)) {
	if (x > 1e-7) {
		sinf = 1/(1+4/(1000*x))
	} else {
		sinf = 1/(1+4/(1000*(1e-7)))
	}
}
