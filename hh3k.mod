TITLE HH K channel
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)
: Updated by K. Archie:
:    removed leak current
:    changed to derivatives rather than explicit step calculation to
:    support NEURON's spiffy, smarter integration
: BFB Cleaned (2007)
: Aaron Milstein modified 2015

NEURON {
    SUFFIX hh3k
    USEION k READ ek WRITE ik
    RANGE gbar, gbar2
    GLOBAL vmin, vmax
    GLOBAL taun, pown, taun2
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gbar = .12 (mho/cm2)
    gbar2 = .12 (mho/cm2)
    vmin = -120 (mV)
    vmax = 100 (mV)
    taun = 1 (ms)
    taun2 = 10 (ms)
    pown = 2
}

STATE {    
    n   <1e-1>
    n2  <1e-1>
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
}

INITIAL {
    n = ssn(v)
    n2 = ssn(v)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gbar*(v - ek)*n^pown + gbar2*(v - ek)*n2^pown
}

DERIVATIVE states {
    n' = (ssn(v) - n)/taun
    n2' = (ssn(v) - n2)/taun2
}

FUNCTION ssn(v(mV)) {  : K activation steady state
    TABLE FROM vmin TO vmax WITH 200
    ssn = 1/(1 + exp((v + 40 (mV))/(-3 (mV))))
}
