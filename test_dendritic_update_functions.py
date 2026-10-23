
import argparse
import ast
import math
from pathlib import Path
import sys
import traceback

SOURCE_NAME = 'optimize_dendritic_excitability_102715.py'


def load_model(path):
    """Execute the original module only up to its bottom-level x1 setup."""
    source = path.read_text()
    tree = ast.parse(source, filename=str(path))
    # Keep imports, functions, CA1 cell initialization, and x0 definitions;
    # omit the original bottom-level run block (which currently overwrites x1).
    tree.body = [node for node in tree.body if node.lineno < 863]
    ast.fix_missing_locations(tree)
    namespace = {'__name__': '_ca1_update_test_', '__file__': str(path)}
    exec(compile(tree, str(path), 'exec'), namespace)
    return namespace


def close(a, b, *, label):
    if not math.isclose(float(a), float(b), rel_tol=1e-6, abs_tol=1e-12):
        raise AssertionError(f'{label}: expected {b:.9g}, got {a:.9g}')
    print(f'    PASS {label}: {float(a):.9g}')


def soma_value(m, mech, param):
    return getattr(getattr(m['cell'].tree.root.sec(0.5), mech), param)


def dend_samples(m, mech, param):
    """Return (distance, value) from all trunk segments that contain mech."""
    values = []
    for node in m['cell'].trunk:
        for seg in node.sec:
            if not hasattr(seg, mech):
                continue
            distance = m['cell'].get_distance_to_node(m['cell'].tree.root, node, seg.x)
            values.append((float(distance), float(getattr(getattr(seg, mech), param))))
    return sorted(values)


def assert_changed(before, after, label):
    if len(before) != len(after):
        raise AssertionError(f'{label}: segment count changed')
    changed = sum(not math.isclose(x[1], y[1], rel_tol=1e-6, abs_tol=1e-12)
                  for x, y in zip(before, after))
    if not changed:
        raise AssertionError(f'{label}: no sampled trunk values changed')
    print(f'    PASS {label}: {changed}/{len(before)} trunk segments changed')


def run_test(name, callback, results):
    print(f'\n=== {name} ===')
    try:
        callback()
        results.append((name, True))
    except Exception as exc:
        results.append((name, False))
        print(f'    FAIL: {type(exc).__name__}: {exc}')
        traceback.print_exc(limit=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', default=SOURCE_NAME, help='Original optimizer script path')
    parser.add_argument('--simulate', action='store_true', help='Also run slow passive/HCN/spike simulation checks')
    args = parser.parse_args()
    path = Path(args.source).resolve()
    if not path.is_file():
        sys.exit(f'Cannot find {path}. Put this test file in the CA1Sim folder, or use --source PATH.')

    print('Loading model setup (without starting optimization)...')
    m = load_model(path)
    x0 = m['x0']
    results = []

    def pas():
        base = list(x0['pas'])
        m['update_pas'](base)
        close(soma_value(m, 'pas', 'g'), base[0], label='soma pas.g')
        before = dend_samples(m, 'pas', 'g')
        alt = base.copy()
        alt[0] *= 1.2
        alt[1] *= 1.2
        m['update_pas'](alt)
        close(soma_value(m, 'pas', 'g'), alt[0], label='changed soma pas.g')
        assert_changed(before, dend_samples(m, 'pas', 'g'), 'trunk pas.g response')
        m['update_pas'](base)

    def hcn():
        base = list(x0['h'])
        m['update_h'](base)
        close(soma_value(m, 'h', 'ghbar'), base[0], label='soma h.ghbar')
        before = dend_samples(m, 'h', 'ghbar')
        alt = base.copy()
        alt[0] *= 1.2
        alt[1] *= 1.2
        m['update_h'](alt)
        close(soma_value(m, 'h', 'ghbar'), alt[0], label='changed soma h.ghbar')
        assert_changed(before, dend_samples(m, 'h', 'ghbar'), 'trunk HCN response')
        m['update_h'](base)

    def naka():
        base = list(x0['na_ka'])
        m['update_na_ka'](base)
        close(soma_value(m, 'nas', 'sh'), base[0], label='soma nas.sh')
        close(soma_value(m, 'kap', 'gkabar'), base[1], label='soma kap.gkabar')
        close(soma_value(m, 'kdr', 'gkdrbar'), base[2], label='soma kdr.gkdrbar')
        before = dend_samples(m, 'kad', 'gkabar')
        alt = base.copy()
        alt[3] *= 1.2
        m['update_na_ka'](alt)
        assert_changed(before, dend_samples(m, 'kad', 'gkabar'), 'dendritic kad gradient response')
        m['update_na_ka'](base)

    def naka_soma():
        base = list(x0['na_ka'][:3])
        m['update_na_ka_soma'](base)
        close(soma_value(m, 'nas', 'sh'), base[0], label='soma-only nas.sh')
        close(soma_value(m, 'kap', 'gkabar'), base[1], label='soma-only kap.gkabar')
        close(soma_value(m, 'kdr', 'gkdrbar'), base[2], label='soma-only kdr.gkdrbar')
        alt = base.copy()
        alt[0] += 0.2
        m['update_na_ka_soma'](alt)
        close(soma_value(m, 'nas', 'sh'), alt[0], label='changed soma-only nas.sh')
        m['update_na_ka_soma'](base)

    def naka_dend():
        base = [x0['na_ka'][0], x0['na_ka'][3]]
        m['update_na_ka_dend'](base)
        before = dend_samples(m, 'kad', 'gkabar')
        alt = [base[0], base[1] * 1.2]
        m['update_na_ka_dend'](alt)
        assert_changed(before, dend_samples(m, 'kad', 'gkabar'), 'dend-only kad gradient response')
        m['update_na_ka_dend'](base)

    def vrest():
        # There is no x0['v_rest'] in this script. The function currently
        # requires four entries and uses x[2], x[3] to set pas.e gradient.
        # We use HCN's four-value vector only as a structural smoke test,
        # NOT as validated physiological v_rest parameters.
        before = dend_samples(m, 'pas', 'e')
        test = list(x0['h'])
        m['update_v_rest'](test)
        after = dend_samples(m, 'pas', 'e')
        assert_changed(before, after, 'trunk pas.e response (structural only)')
        m['update_pas'](x0['pas'])
        print('    WARNING: update_v_rest has no defined parameter vector or dedicated optimizer target.')

    for name, fn in [('update_pas', pas), ('update_h', hcn),
                     ('update_na_ka', naka), ('update_na_ka_soma', naka_soma),
                     ('update_na_ka_dend', naka_dend), ('update_v_rest', vrest)]:
        run_test(name, fn, results)

    if args.simulate:
        print('\nSimulation checks take longer and test model output, not just parameter assignment.')

        def passive_sim():
            error = m['pas_error'](m['np'].array(x0['pas'], dtype=float), plot=0)
            if not math.isfinite(error):
                raise AssertionError('nonfinite passive error')
            print(f'    PASS passive simulation completed; objective error = {error:.6g}')

        def hcn_sim():
            m['update_pas'](x0['pas'])
            error = m['h_error'](m['np'].array(x0['h'], dtype=float), plot=0)
            if not math.isfinite(error):
                raise AssertionError('nonfinite HCN error')
            print(f'    PASS HCN simulation completed; objective error = {error:.6g}')

        def naka_sim():
            m['update_h'](x0['h'])
            error = m['na_ka_error'](m['np'].array(x0['na_ka'], dtype=float), plot=0)
            if not math.isfinite(error):
                raise AssertionError('nonfinite Na/K error')
            print(f'    PASS Na/K simulation completed; objective error = {error:.6g}')

        for name, fn in [('passive simulation', passive_sim),
                         ('HCN simulation', hcn_sim), ('Na/K simulation', naka_sim)]:
            run_test(name, fn, results)

    print('\n=== SUMMARY ===')
    for name, passed in results:
        print(f'{"PASS" if passed else "FAIL"}: {name}')
    failures = sum(not passed for _, passed in results)
    print(f'{len(results) - failures}/{len(results)} checks passed')
    if failures:
        sys.exit(1)


if __name__ == '__main__':
    main()
