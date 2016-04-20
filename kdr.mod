TITLE K-DR channel
: from Klee Ficker and Heinemann
: modified to account for Dax et al.
: M.Migliore 1997
: Aaron Milstein modified in 2015:
: removed q10: was set to 1, so was unnecessary

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (mol) = (1)
}

NEURON {
        SUFFIX kdr
        USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar,ik
        RANGE ninf,taun
        GLOBAL nscale
}

PARAMETER {
        temp    = 24            (degC)
        gkdrbar = 0.003         (mho/cm2)
        vhalfn  = 13            (mV)
        a0n     = 0.02          (/ms)
        zetan   = -3            (1)
        gmn     = 0.7           (1)
        nmin    = 2             (ms)
        nscale  = 1
}

STATE {
        n
}

ASSIGNED {
        v                       (mV)
        ik                      (mA/cm2)
        ninf
        gkdr                    (mho/cm2)
        taun                    (ms)
        ek                      (mV)
        celsius                 (degC)
}

INITIAL {
        rates(v)
        n=ninf
}        

BREAKPOINT {
        SOLVE states METHOD cnexp
        gkdr = gkdrbar*n
        ik = gkdr*(v-ek)
}

DERIVATIVE states {
        rates(v)
        n' = (ninf-n)/taun
}

FUNCTION alpn(v(mV)) {
        alpn = exp(zetan*(v-vhalfn)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/degC/mol)*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) {
        betn = exp(zetan*gmn*(v-vhalfn)*1.e-3(V/mV)*9.648e4(coulomb/mol)/(8.315(joule/degC/mol)*(273.16(degC)+celsius))) 
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(a0n*(1+a))
        if (taun<nmin) {taun=nmin}
        taun=taun/nscale
}