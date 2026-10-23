
import argparse
import subprocess
import sys
import traceback
import numpy as np

MODULE = 'optimize_dendritic_excitability_102715'


def finite_scalar(value):
    assert value is not None, 'Function returned None'
    assert np.isscalar(value), f'Expected scalar, got {type(value)}'
    assert np.isfinite(value), f'Expected finite result, got {value}'
    return value


def run_case(name):
    model = __import__(MODULE)
    if name == 'offset_vm':
        for compartment in ('soma', 'trunk', 'tuft'):
            value = finite_scalar(model.offset_vm(compartment))
            assert -120 < value < 50, (compartment, value)
            print(f'  {compartment}: initial resting voltage {value:.3f} mV; holding current {model.i_holding[compartment]:.4f} nA')
    elif name == 'zero_h':
        model.zero_h()
        assert np.isclose(model.cell.mech_dict['soma']['h']['ghbar']['value'], 0)
        assert np.isclose(model.cell.mech_dict['trunk']['h']['ghbar']['value'], 0)
        assert np.isclose(model.cell.mech_dict['trunk']['h']['ghbar']['slope'], 0)
    elif name == 'zero_na':
        model.zero_na()
        for section, mech in [('soma', 'nas'), ('axon_hill', 'nax'), ('ais', 'nax')]:
            value = model.cell.mech_dict[section][mech]['gbar']['value']
            assert np.isclose(value, 0), (section, mech, value)
    elif name == 'get_spike_shape':
        # Synthetic action potential; tests numerical processing, not biological accuracy.
        t = np.arange(0, model.duration, model.dt)
        vm = np.full_like(t, -67., dtype=float)
        # Smooth rising spike, followed by hyperpolarization.
        vm += 107. * np.exp(-0.5 * ((t - (model.equilibrate + 5.)) / 0.45)**2)
        vm -= 4. * np.exp(-0.5 * ((t - (model.equilibrate + 9.)) / 1.2)**2)
        peak, threshold, adp, ahp = model.get_spike_shape(vm)
        for v in (peak, threshold, adp, ahp):
            finite_scalar(v)
        assert 25 < peak < 50, peak
        print(f'  peak={peak:.3f}, threshold={threshold:.3f}, ADP={adp:.3f}, AHP={ahp:.3f}')
    elif name == 'get_EPSP_dur':
        # Baseline -67 mV; single synthetic EPSP with exponential decay.
        t = np.arange(0, 300, 0.05)
        vm = np.full_like(t, -67., dtype=float)
        after = t >= 250
        s = t[after] - 250
        vm[after] += 5 * (1 - np.exp(-s / 1.0)) * np.exp(-s / 10.0)
        result = finite_scalar(model.get_EPSP_dur(t, vm, 250., 300.))
        assert 0 < result < 50, result
        print(f'  EPSP duration={result:.3f} ms')
    elif name == 'pas_error':
        result = finite_scalar(model.pas_error(np.asarray(model.x0['pas'], dtype=float)))
        assert result >= 0
        print(f'  error={result:.6g}')
    elif name == 'h_error':
        result = finite_scalar(model.h_error(np.asarray(model.x0['h'], dtype=float)))
        assert result >= 0
        print(f'  error={result:.6g}')
    elif name == 'v_rest_error':
        result = finite_scalar(model.v_rest_error())
        assert result >= 0
        print(f'  error={result:.6g}')
    elif name in ('na_ka_error', 'na_ka_error_soma', 'na_ka_error_dend'):
        x = np.asarray(model.x0['na_ka'], dtype=float)
        if name == 'na_ka_error_soma':
            x = x[:3]
        elif name == 'na_ka_error_dend':
            x = x[[0, 3]]
        result = finite_scalar(getattr(model, name)(x))
        assert result >= 0
        print(f'  error={result:.6g}')
    elif name == 'optimize_polish':
        result = model.optimize_polish('na_ka', np.asarray(model.x0['na_ka'], dtype=float), maxfev=2)
        assert 'na_ka' in result
        finite_scalar(result['na_ka']['Err'])
        assert np.all(np.isfinite(result['na_ka']['x']))
        print(f"  error={result['na_ka']['Err']:.6g}")
    elif name == 'optimize_explore':
        # Some versions of the original script pass a float as niter_success.
        # If that raises TypeError, fix the source to use maxfev // 2.
        result = model.optimize_explore('na_ka', np.asarray(model.x0['na_ka'], dtype=float),
                                        model.xmin['na_ka'], model.xmax['na_ka'], maxfev=2)
        assert 'na_ka' in result
        finite_scalar(result['na_ka']['Err'])
        assert np.all(np.isfinite(result['na_ka']['x']))
        print(f"  error={result['na_ka']['Err']:.6g}")
    else:
        raise ValueError(f'Unknown case: {name}')
    print(f'PASS: {name}')


CASES = ['offset_vm', 'zero_h', 'zero_na', 'get_spike_shape', 'get_EPSP_dur',
         'pas_error', 'h_error', 'v_rest_error', 'na_ka_error',
         'na_ka_error_soma', 'na_ka_error_dend']
OPTIMIZERS = ['optimize_polish', 'optimize_explore']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=CASES + OPTIMIZERS)
    parser.add_argument('--include-optimizers', action='store_true')
    parser.add_argument('--timeout', type=int, default=180)
    args = parser.parse_args()
    if args.case:
        try:
            run_case(args.case)
        except Exception:
            traceback.print_exc()
            sys.exit(1)
        return
    selected = CASES + (OPTIMIZERS if args.include_optimizers else [])
    failed = []
    for name in selected:
        print(f'\n=== {name} ===', flush=True)
        try:
            result = subprocess.run([sys.executable, __file__, '--case', name], timeout=args.timeout)
            if result.returncode != 0:
                failed.append(name)
        except subprocess.TimeoutExpired:
            print(f'TIMEOUT: {name} (>{args.timeout} seconds)')
            failed.append(name)
    print('\nSUMMARY:', len(selected)-len(failed), 'passed,', len(failed), 'failed/timed out')
    if failed:
        print('Needs investigation:', ', '.join(failed))
        sys.exit(1)


if __name__ == '__main__':
    main()
