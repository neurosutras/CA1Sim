: Ca channels (T,N,L-type)


NEURON {
	SUFFIX Ca
	USEION ca WRITE ica
	RANGE gtcabar, gncabar, glcabar, gtca, gnca, glca
	RANGE ainf, taua, binf, taub, e_ca, gbar, gcamult, i, ca_i
	GLOBAL ca0, cao, taucadiv, B, tauctdiv, Vshift
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	ca0 = .00007	(mM)		: initial calcium concentration inside
	cao = 1.6		  (mM)	  : calcium concentration outside
	tau = 9		    (ms)
	taucadiv = 1
	tauctdiv = 1
	gtcabar = .00015	(S/cm2)
	gncabar = .002	(S/cm2)
	glcabar = .01	(S/cm2)
    gcamult = 1.
	B = .26 (mM-cm2/mA-ms)
	Vshift = 0    (mV)
}

ASSIGNED {
	v		(mV)
	e_ca		(mV)
	ica			(mA/cm2)
	i 			(mA/cm2)
	gtca		(S/cm2)
	gnca		(S/cm2)
	glca		(S/cm2)
	gbar		(S/cm2)
	celsius	(degC)
	taua    (ms)
	ainf
	taub    (ms)
	binf
}

STATE { 
	ca_i (mM)  <1e-5>
	a 
	b 
	c 
	d 
	e
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit
	e_ca = (1000)*(celsius+273.15)*R/(2*F)*log(cao/ca_i)
	gtca = gtcabar*gcamult*a*a*b
	gnca = gncabar*gcamult*c*c*d
	glca = glcabar*gcamult*e*e
	ica = (gtca+gnca+glca)*(v - e_ca)
	i = ica
}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	ca_i' = -B*ica - taucadiv*(ca_i-ca0)/tau
	ainf = alphaa(v)/(alphaa(v) + betaa(v))
	taua = 1/(alphaa(v) + betaa(v))
	a' = (ainf - a)/taua
	binf = alphab(v)/(alphab(v) + betab(v))
	taub = 1/(alphab(v+Vshift) + betab(v+Vshift))
	b' = (binf - b)/taub
	c' = alphac(v)*(1-c)-betac(v)*c
	d' = alphad(v)*(1-d)-betad(v)*d
	e' = alphae(v)*(1-e)-betae(v)*e
}

INITIAL {
	ca_i = ca0
	a = alphaa(v)/(alphaa(v)+betaa(v))
	b = alphab(v)/(alphab(v)+betab(v))
	c = alphac(v)/(alphac(v)+betac(v))
	d = alphad(v)/(alphad(v)+betad(v))
	e = alphae(v)/(alphae(v)+betae(v))
	gbar = gtcabar + gncabar + glcabar 
}

FUNCTION alphaa(v (mV)) (/ms) {
	alphaa = f(2,0.1,v,19.26)
}

FUNCTION betaa(v (mV)) (/ms) {
	betaa = exponential(0.009,-0.045393,v,0)
}

FUNCTION alphab(v (mV)) (/ms) {
	alphab = exponential(1e-6,-0.061501,v,0)
}

FUNCTION betab(v (mV)) (/ms) {
	betab = logistic(1,-0.1,v,29.79)
}

FUNCTION alphac(v (mV)) (/ms) {
	alphac = f(1.9,0.1,v,19.88)
}

FUNCTION betac(v (mV)) (/ms) {
	betac = exponential(0.046,-0.048239,v,0)
}

FUNCTION alphad(v (mV)) (/ms) {
	alphad = exponential(1.6e-4,-0.020661,v,0)
}

FUNCTION betad(v (mV)) (/ms) {
	betad = logistic(1,-0.1,v,39)
}

FUNCTION alphae(v (mV)) (/ms) {
	alphae = f(156.9,0.1,v,81.5)
}

FUNCTION betae(v (mV)) (/ms) {
	betae = exponential(0.29,-0.092081,v,0)
}

FUNCTION f(A, k, v (mV), D) (/ms) {
	LOCAL x
	UNITSOFF
	x = k*(v-D)
	if (fabs(x) > 1e-6) {
		f = A*x/(1-exptrap(-x))
	}else{
		f = A/(1-0.5*x)
	}
	UNITSON
}

FUNCTION logistic(A, k, v (mV), D) (/ms) {
	UNITSOFF
	logistic = A/(1+exptrap(k*(v-D)))
	UNITSON
}

FUNCTION exponential(A, k, v (mV), D) (/ms) {
	UNITSOFF
	exponential = A*exptrap(k*(v-D))
	UNITSON
}

FUNCTION exptrap(x) {
  if (x>=700.0) {
    exptrap = exp(700.0)
  } else {
    exptrap = exp(x)
  }
}



