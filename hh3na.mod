TITLE HH na channel
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)
: Updated by K. Archie:
:    removed leak current
:    changed to derivatives rather than explicit step calculation to
:    support NEURON's spiffy, smarter integration
: BFB Cleaned (2007)
: Aaron Milstein modified 2015

NEURON {
    SUFFIX hh3na
    USEION na READ ena WRITE ina
    RANGE gbar, g
    GLOBAL vmin, vmax
    GLOBAL taum, tauh, taus, tausvh, tausvs, tausb
    GLOBAL powm, powh, pows
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    gbar = .20 (mho/cm2)
    vmin = -120 (mV)
    vmax = 100 (mV)
    tausvh = 30 (mV)
    tausvs = 1 (mV)
    taus = 50 (ms)
    tausb = .5 (ms)
    taum = .05 (ms)
    tauh = .5 (ms)
    pows = 1
    powm = 3
    powh = 1
}

STATE {    
    h  <1e-1> 
    m  <1e-1> 
    s  <1e-1> 
}

ASSIGNED {
    v (mV)
    ina (mA/cm2)
    ena (mV)
    g (mho/cm2)
}

INITIAL {
    m = ssm(v)
    h = ssh(v)
    s = sss(v)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gbar*(m^powm)*(h^powh)*(s^pows)
    ina = g*(v - ena)
}

DERIVATIVE states {
    m' = (ssm(v) - m)/taum
    h' = (ssh(v) - h)/tauh
    s' = (sss(v) - s)/tauss(v)
}

FUNCTION ssm(v (mV)) {  : Na activation steady state
    TABLE FROM vmin TO vmax WITH 200
    ssm = 1/(1 + exp((v + 40 (mV))/(-3 (mV))))
}

FUNCTION ssh(v (mV)) {  : Na inactivation steady state
    TABLE FROM vmin TO vmax WITH 200
    ssh = 1/(1 + exp((v + 45 (mV))/3 (mV)))
}

FUNCTION sss(v (mV)) {  : Na ... steady state
    TABLE FROM vmin TO vmax WITH 200
    sss = 1/(1 + exp((v + 44 (mV))/3 (mV)))
}

FUNCTION tauss(v (mV)) (ms) {  : Na ... tau
    TABLE FROM vmin TO vmax WITH 200
    tauss = tausb + taus/(1 + exp((v + tausvh)/tausvs))
}
