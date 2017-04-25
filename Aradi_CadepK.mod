: Ca-dependent K channels (BK and SK)


NEURON {
	SUFFIX CadepK
	USEION ca READ ica
	USEION k READ ek WRITE ik
	RANGE gbkbar, gskbar, gbar, i, ask, bsk, gsk, gbk, isk, ibk, gcakmult, ca_i
	GLOBAL ca0, tau, stau, taucadiv, tauskdiv
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	B = .26 (mM-cm2/mA-ms)
}

PARAMETER {
	gbkbar = .0003	(S/cm2)	: maximum permeability
	gskbar = .0005	(S/cm2)	: maximum permeability
	gcakmult = 1.
	ca0 = .00007	(mM)
	tau = 9		(ms)
	taucadiv = 1
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
	ica		(mA/cm2)
	area	(microm2)
  	gbk		(S/cm2)
  	gsk		(S/cm2)
  	gbar  (S/cm2)
}

STATE { 
	ca_i (mM)		<1e-5> 
	q 
	r 
	s 
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
	gbk = gbkbar*gcakmult*r*s*s
	gsk = gskbar*gcakmult*q*q
	isk = gsk*(v - ek)
	ibk = gbk*(v - ek)
	ik = isk + ibk
	i = ik
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	ca_i' = -B*ica - decay_rate(ca_i)
	q' = tauskdiv*(ask*alphaq(ca_i)*(1-q)-bsk*betaq(ca_i)*q)
	r' = alphar*(1-r)-betar(v)*r
	s' = (sinf(ca_i)-s)/stau
}

INITIAL {
	ca_i = ca0
	q = alphaq(ca_i)/(alphaq(ca_i)+betaq(ca_i))
	r = alphar/(alphar+betar(v))
  	s = sinf(ca_i)
  	gbar = gbkbar + gskbar
}

FUNCTION decay_rate(ca_i (mM)) {
	LOCAL rate, min_rate
	min_rate = 0.0015
	rate = taucadiv*(ca_i-ca0)/tau
	if (rate < min_rate) {
		decay_rate = min_rate
	} else {
		decay_rate = rate
	}
}


FUNCTION exp1(A (/ms), d, k, x (mM)) (/ms) {
	UNITSOFF
	exp1 = A/exptrap((12*log10(x)+d)/k)
	UNITSON
}

FUNCTION alphaq(x (mM)) (/ms) {
	alphaq = exp1(0.00246,28.48,-4.5,x)	:28
}

FUNCTION betaq(x (mM)) (/ms) {
	betaq = exp1(0.006,60.4,35,x)
}

FUNCTION betar(v (mV)) (/ms) {
	UNITSOFF
	betar = 0.11/exptrap((v-35)/14.9)
	UNITSON
}

FUNCTION sinf(x (mM)) {
	UNITSOFF
	sinf = 1/(1+4/(1000*x))
	UNITSON
}

FUNCTION exptrap(x) {
  if (x>=700.0) {
    exptrap = exp(700.0)
  } else {
    exptrap = exp(x)
  }
}
