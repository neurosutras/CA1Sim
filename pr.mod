TITLE stochastic release probability

COMMENT

Milstein 2015. When this point process receives a spike, it requests a random number from a random number generator and
compares it to an internal release probability variable 'P' to decide if a single vesicle should be released at this
synapse. The release probability is then updated to follow specified dynamics of facilitation and depression. In order
to make use of this event, additional synaptic mechanisms must be connected to this point process via a NetCon object.

Dynamics based on:

Implementation of the model of short-term facilitation and depression described in
  Varela, J.A., Sen, K., Gibson, J., Fost, J., Abbott, L.R., and Nelson, S.B.
  A quantitative description of short-term plasticity at excitatory synapses
  in layer 2/3 of rat primary visual cortex
  Journal of Neuroscience 17:7926-7940, 1997

ENDCOMMENT

NEURON {
	POINT_PROCESS Pr
	RANGE P, P0, random, f, tau_F, d1, tau_D1, F, D1, tlast
    THREADSAFE
    POINTER randObjPtr
}

PARAMETER {
    : the (1) is needed for the range limits to be effective
    P0 = 0.200          (1)     < 0, 1 >        : basal release probability
    f = 1.769           (1)     < 0, 1e9 >      : additive facilitation per spike
    tau_F = 67.351      (ms)    < 1e-9, 1e9 >   : rate of decay back to baseline following facilitation
    d1 = 0.878          (1)     < 0, 1 >        : multiplicative fast depression per spike
    tau_D1 = 92.918     (ms)    < 1e-9, 1e9 >   : rate of decay back to baseline following fast depression
    :d2 = 0.975         (1)     < 0, 1 >        : multiplicative slow depression per spike
    :tau_D2 = 9200      (ms)    < 1e-9, 1e9 >   : rate of decay back to baseline following slow depression
}

ASSIGNED {
	P				        : instantaneous release probability
    randObjPtr              : pointer to a hoc random number generator Random.uniform(0,1)
    random                  : individual instance of random number
    F                       : current level of facilitation
    D1                      : current level of fast depression
    :D2                      : current level of slow depression
    tlast (ms)              : time of last spike
}

INITIAL {
	P = P0
    random = 1
    F = 1
    D1 = 1
    :D2 = 1
}

NET_RECEIVE(weight) {
    INITIAL {
        tlast = t
    }
    F = 1 + (F-1)*exp(-(t - tlast)/tau_F)
    D1 = 1 - (1-D1)*exp(-(t - tlast)/tau_D1)
:    D2 = 1 - (1-D2)*exp(-(t - tlast)/tau_D2)
:    if (P0*F*D1*D2 > 1) {
    if (P0*F*D1 > 1) {
        P = 1
    } else {
:        P = P0*F*D1*D2
        P = P0*F*D1
    }
    random = randGen()
    if (random <= P) {
        net_event(t)
    }
    tlast = t
    F = F + f
    D1 = D1 * d1
:    D2 = D2 * d2
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION randGen() {
VERBATIM
   if (_p_randObjPtr) {
      /*
      :Supports separate independent but reproducible streams for
      : each instance. However, the corresponding hoc Random
      : distribution MUST be set to Random.uniform(0,1)
      */
      _lrandGen = nrn_random_pick(_p_randObjPtr);
   }else{
      hoc_execerror("Random object ref not set correctly for randObjPtr"," only via hoc Random");
   }
ENDVERBATIM
}

PROCEDURE setRandObjRef() {
VERBATIM
   void** pv4 = (void**)(&_p_randObjPtr);
   if (ifarg(1)) {
      *pv4 = nrn_random_arg(1);
   }else{
      *pv4 = (void*)0;
   }
ENDVERBATIM
}