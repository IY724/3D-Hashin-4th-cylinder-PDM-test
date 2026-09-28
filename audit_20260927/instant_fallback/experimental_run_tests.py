# -*- coding: utf-8 -*-
"""实际编译Fortran + 原版回归 + 张量参考 + 转换器测试。
所有构建与测试产物隔离到runs；不使用Abaqus许可证。
"""
from pathlib import Path
import ctypes as ct
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import unittest
import numpy as np
from build import build, check_fixed_format
from model import (MODE, NSTATE, make_props, elastic, rotation, tensor,
                   vector, rotate_strain, rotate_stress_back, smeared_stiffness, error_text)
from convert_wcm_inp import transform, blocks, material_groups

ROOT = Path(__file__).resolve().parents[1]
PLY = [141000., 11400., 11400., .28, .28, .28, 7100., 7100., 7100.]
STRENGTH = [2080., 1250., 60., 290., 110., 110., 110.]
BASE = {'A': PLY+STRENGTH+[.07, .14, .2, .4, .9, .5],
        'B': PLY+STRENGTH+[133., 10., .5, 1.6, 1., .9, .5, 1e-4, 5e-3]}
DLL_HANDLES = []


def ptr(a):
    return a.ctypes.data_as(ct.POINTER(ct.c_double))


class Core:
    def __init__(self, directory):
        compiler = shutil.which('gfortran')
        if not compiler:
            raise RuntimeError('未找到gfortran；不能将未运行测试标为通过')
        if os.name == 'nt':
            DLL_HANDLES.append(os.add_dll_directory(str(Path(compiler).parent)))
        self.run = directory
        directory.mkdir(parents=True)
        sources = []
        (directory/'ABA_PARAM.INC').write_text('      IMPLICIT REAL*8(A-H,O-Z)\n', encoding='ascii')
        for mode, name in (('A', 'hashin_constant.for'), ('B', 'hashin.for')):
            text = (ROOT/'baseline'/name).read_text(encoding='utf-8')
            text = text.replace('SUBROUTINE UMAT(', 'SUBROUTINE LEGACY_'+mode+'(', 1)
            target = directory/('legacy_'+mode+'.for')
            target.write_text(text, encoding='utf-8')
            sources.append(str(target))
        # 旧A代码有长计算行，仅测试原版时按完整行编译。
        legacy = directory/'legacy.dll'
        cmds = [
            [compiler, '-shared', '-O1', '-frecursive', '-fcheck=all',
             '-ffixed-line-length-none', '-I'+str(directory), *sources,
             str(ROOT/'tests/baseline_driver.for'), '-o', str(legacy)],
            [compiler, '-shared', '-O1', '-frecursive', '-fcheck=all',
             '-Wall', '-Wextra', '-ffixed-line-length-72',
             str(ROOT/'src/wcm_core.for'), '-o', str(directory/'core.dll')]]
        for i, cmd in enumerate(cmds):
            result = subprocess.run(cmd, cwd=directory, capture_output=True, text=True)
            (directory/f'compile_{i}.txt').write_text(result.stdout+result.stderr, encoding='utf-8')
            if result.returncode:
                raise RuntimeError(result.stdout+result.stderr)
        self.library = ct.CDLL(str(directory/'core.dll'))
        self.legacy = ct.CDLL(str(legacy))
        self.library.wcm_update_.restype = None
        self.legacy.legacy_eval_.restype = None

    def update(self, variant, props, old, e0, de, dt=.01, length=.1, advance=1, allow_error=False):
        props, old, e0, de = [np.array(a, dtype=np.float64, order='F', copy=True)
                              for a in (props, old, e0, de)]
        new = old.copy()
        s, c = np.zeros(6), np.zeros((6, 6), order='F')
        mode, np_, adv, err = map(ct.c_int, (MODE[variant], len(props), advance, 0))
        dt_, lc = ct.c_double(dt), ct.c_double(length)
        self.library.wcm_update_(ct.byref(mode), ptr(props), ct.byref(np_), ptr(old),
            ptr(e0), ptr(de), ct.byref(dt_), ct.byref(lc), ct.byref(adv),
            ptr(new), ptr(s), ptr(c), ct.byref(err))
        if err.value and not allow_error:
            raise RuntimeError(error_text(err.value))
        return new, s, c, err.value

    def original(self, variant, old, e0, de, dt=.01, length=.1):
        props, old, e0, de = [np.array(a, dtype=np.float64, order='F', copy=True)
                              for a in (BASE[variant], old, e0, de)]
        new, s, c = old.copy(), np.zeros(6), np.zeros((6, 6), order='F')
        mode, n = ct.c_int(MODE[variant]), ct.c_int(len(props))
        dt_, lc = ct.c_double(dt), ct.c_double(length)
        self.legacy.legacy_eval_(ct.byref(mode), ptr(props), ct.byref(n), ptr(old),
            ptr(e0), ptr(de), ct.byref(dt_), ct.byref(lc), ptr(new), ptr(s), ptr(c))
        return new, s, c


def stamped(variant, props):
    s = np.zeros(84)
    it = 22 if variant == 'A' else 25
    s[80:] = [props[it], props[it+1], MODE[variant], 202609]
    return s


class MaterialTests(unittest.TestCase):
    core = None
    metrics = {}

    def assertClose(self, a, b, tolerance=1e-9):
        np.testing.assert_allclose(a, b, rtol=tolerance, atol=tolerance)

    def test_rotation_and_power(self):
        rng = np.random.default_rng(109)
        for angle in (0, 15, 37, 45, 79, 90):
            r = rotation(angle)
            for _ in range(15):
                e, s = rng.normal(size=6), rng.normal(size=6)
                el = rotate_strain(e, angle)
                sl = vector(r @ tensor(s) @ r.T)
                self.assertClose(np.dot(s, e), np.dot(sl, el))
                self.assertClose(rotate_strain(el, -angle), e)
                self.assertClose(rotate_stress_back(sl, angle), s)

    def test_actual_fortran_all_elastic_angles_and_weights(self):
        rng = np.random.default_rng(29)
        for variant in ('A', 'B'):
            base = BASE[variant].copy()
            base[6:9] = [4900., 6200., 3800.]
            for angle in (0, 15, 37, 45, 65, 90):
                for w in (0., .3, .5, 1.):
                    p = make_props(variant, base, angle, weight=w)
                    expected = smeared_stiffness(base[:9], angle, w)
                    for e in np.vstack([np.eye(6)*1e-5, rng.normal(size=(2, 6))*1e-5]):
                        state, s, c, _ = self.core.update(variant, p, np.zeros(84), np.zeros(6), e)
                        self.assertClose(c, expected, 2e-10)
                        self.assertClose(s, expected @ e)
                        if w:
                            self.assertClose(state[34:40], rotate_strain(e, angle))

    def test_original_formula_and_history_regression(self):
        # 真实旧源码对照：0度两族相同，应力、切线、旧SDV逐增量相等。
        for variant in ('A', 'B'):
            for component, peak in ((0, .025), (0, -.02), (1, .012),
                                    (1, -.04), (3, .04), (4, .04), (5, .04)):
                p = make_props(variant, BASE[variant], 0, capf=.95, capm=.85)
                new, old, e0 = np.zeros(84), np.zeros(24), np.zeros(6)
                load = np.r_[np.linspace(0, peak, 31), np.linspace(peak, 0, 15),
                             np.linspace(0, peak*1.1, 15)]
                direction = np.linalg.solve(elastic(PLY), np.eye(6)[component])
                direction /= direction[component]
                for value in load:
                    # 自由泊松收缩的单轴应力路径，避免把横向零应变
                    # 却有横向起始应力的等效位移退化奇点当作有效基准。
                    e = direction*value
                    old, so, co = self.core.original(variant, old, e0, e-e0)
                    new, sn, cn, _ = self.core.update(variant, p, new, e0, e-e0)
                    # B intentionally differs after the historical .95/.85 cap.
                    if variant == 'B' and np.any(new[:4] >= [.95,.95,.85,.85]):
                        break
                    count = 12 if variant == 'A' else 24
                    self.assertClose(new[:count], old[:count], 2e-10)
                    self.assertClose(new[40:40+count], old[:count], 2e-10)
                    self.assertClose(sn, so, 2e-10)
                    self.assertClose(cn, co, 2e-10)
                    e0 = e

    def test_90_degree_fiber_not_matrix(self):
        # 纯环向100MPa应当恢复纤维拉伸，而非以YT=60判断基体失效。
        for variant in ('A', 'B'):
            p = make_props(variant, BASE[variant], 90)
            c = smeared_stiffness(PLY, 90)
            e = np.linalg.solve(c, [0, 100, 0, 0, 0, 0])
            state, _, _, _ = self.core.update(variant, p, np.zeros(84), np.zeros(6), e)
            self.assertClose(state[24], (100/2080)**2)
            self.assertLess(abs(state[26]), 1e-15)
            self.assertClose(state[28], 100.)

    def test_independent_families_and_coupling(self):
        for variant in ('A', 'B'):
            p = make_props(variant, BASE[variant], 35)
            state = stamped(variant, p)
            if variant == 'A':
                state[6] = 1  # 仅+族纤维拉伸失效。
            else:
                state[20] = .6  # 仅+族纤维黏性损伤。
            out, stress, tangent, _ = self.core.update(variant, p, state,
                np.zeros(6), np.ones(6)*1e-5, advance=0)
            self.assertGreater(abs(tangent[0, 3]), 1.)
            self.assertClose(out, state)  # 纯刚度请求不推进历史。
            self.assertClose(tangent, tangent.T)
            self.assertGreater(np.linalg.eigvalsh(tangent)[0], 0)

    def test_initiation_asymmetry(self):
        p = make_props('A', BASE['A'], 45)
        e = np.array([0., 0., 0., .04, 0., 0.])
        state, _, _, _ = self.core.update('A', p, np.zeros(84), np.zeros(6), e)
        self.assertNotEqual(state[6], state[46])
        self.assertNotEqual(state[7], state[47])

    def test_current_update_A(self):
        p = make_props('A', BASE['A'], 0, update=1)
        state, s, c, _ = self.core.update('A', p, np.zeros(84), np.zeros(6), [.02, 0, 0, 0, 0, 0])
        self.assertEqual(state[6], 1)
        self.assertClose(c[0, 0], .93*elastic(PLY)[0, 0])
        self.assertClose(s, c @ np.array([.02, 0, 0, 0, 0, 0]))

    def test_current_update_B_tangent(self):
        base = BASE['B'].copy()
        base[-2:] = [0., 0.]
        p = make_props('B', base, 27, update=1)
        state, _, _, _ = self.core.update('B', p, np.zeros(84), np.zeros(6), [.025, .001, 0, .003, 0, 0])
        e0 = np.array([.025, .001, 0, .003, 0, 0])
        de = np.array([.002, .0001, 0, .0002, 0, 0])
        _, s, c, _ = self.core.update('B', p, state, e0, de)
        fd = np.zeros((6, 6))
        for j in range(6):
            h = np.eye(6)[j]*2e-8
            sp = self.core.update('B', p, state, e0, de+h)[1]
            sm = self.core.update('B', p, state, e0, de-h)[1]
            fd[:, j] = (sp-sm)/4e-8
        self.assertClose(c, fd, 5e-5)

    def test_caps_and_no_healing(self):
        for caps in ((.95, .85), (.97, .97)):
            base = BASE['B'].copy()
            base[-2:] = [0, 0]
            p = make_props('B', base, 0, capf=caps[0], capm=caps[1])
            state, e0 = np.zeros(84), np.zeros(6)
            for value in np.r_[np.linspace(0, .4, 500), np.linspace(.4, 0, 50)]:
                e = np.array([value, value, 0, value, 0, 0])
                previous = state.copy()
                state, _, _, _ = self.core.update('B', p, state, e0, e-e0, length=.1)
                self.assertTrue(np.all(state[:4] >= previous[:4]-1e-14))
                self.assertTrue(np.all(state[:2] <= 1.0))
                self.assertTrue(np.all(state[2:4] <= 1.0))
                e0 = e

    def test_input_validation_and_energy_condition(self):
        for variant in ('A', 'B'):
            good = make_props(variant, BASE[variant], 15)
            theta = 22 if variant == 'A' else 25
            invalid = [(0, -1), (1, float('nan')), (-1, 9),
                       (theta, 91), (theta+1, -1),
                       (20 if variant == 'A' else 27, -1)]
            for index, value in invalid:
                p = good.copy()
                p[index] = value
                self.assertNotEqual(self.core.update(variant, p, np.zeros(84),
                    np.zeros(6), np.zeros(6), allow_error=True)[3], 0)
            bad = np.zeros(84)
            bad[0] = 1
            self.assertEqual(self.core.update(variant, good, bad, np.zeros(6),
                np.zeros(6), allow_error=True)[3], 5)
        p = make_props('B', BASE['B'], 0)
        state, _, _, error = self.core.update('B', p, np.zeros(84),
            np.zeros(6), [0, .01, 0, 0, 0, 0], length=100)
        self.assertEqual(error, 0)
        self.assertEqual(state[2], 1.0)
        self.assertGreater(state[22], 0.0)
        self.assertLess(state[22], 1.0)
        self.assertLess(2*p[18]/state[18], state[14])

    def test_full_input_and_actual_fortran_all_bins(self):
        config = json.loads((ROOT/'config/project.json').read_text(encoding='utf-8'))
        maximum = 0.
        for variant, filename in (('A', 'Job-hashinnewconstant.inp'), ('B', 'Job-hashinnewcdm.inp')):
            source = (ROOT.parent/filename).read_bytes().decode('latin-1')
            result, records = transform(source, variant, config)
            self.assertEqual(len(records), 138)
            old_blocks, new_blocks = blocks(source), blocks(result)
            # 逐块检查节点/单元/取向/分布/截面/载荷等没有改变。
            protected = {'node', 'element', 'orientation', 'distribution', 'distribution table',
                         'solid section', 'boundary', 'dsload', 'tie', 'step', 'static'}
            for key in protected:
                self.assertEqual([b.text for b in old_blocks if b.key == key],
                                 [b.text for b in new_blocks if b.key == key])
            for rec in records:
                p = np.array(rec['props'])
                _, _, c, _ = self.core.update(variant, p, np.zeros(84), np.zeros(6), np.ones(6)*1e-6)
                self.assertClose(c, smeared_stiffness(PLY, rec['angle']), 1e-9)
                maximum = max(maximum, rec['relative_elastic_error'])
            self.assertEqual(len(material_groups(new_blocks)), len(material_groups(old_blocks)))
        self.metrics['bin_count_each'] = 138
        self.metrics['max_wcm_elastic_error'] = maximum

    def test_converter_rejects_missing_angle_and_repeat(self):
        config = json.loads((ROOT/'config/project.json').read_text(encoding='utf-8'))
        p = ', '.join(map(str, PLY))
        snippet = ('*Solid Section, elset=x, material=WCM_Tank1_Mat1_Bin0, orientation=o\n1.,\n'
                   '*Material, name=t-700\n*Elastic, type=ENGINEERING CONSTANTS\n'+p+'\n'
                   '*Material, name=WCM_Tank1_Mat1_Bin0\n*Elastic, type=ENGINEERING CONSTANTS\n'+p+'\n')
        with self.assertRaisesRegex(ValueError, 'Beta'):
            transform(snippet, 'A', config)
        with self.assertRaisesRegex(ValueError, 'Include'):
            transform('*Include, input=missing.inp\n', 'A', config)
        config['angles'] = {'WCM_Tank1_Mat1_Bin0': 0}
        converted, _ = transform(snippet, 'A', config)
        with self.assertRaisesRegex(ValueError, '重复转换'):
            transform(converted, 'A', config)

    def test_original_files_unchanged(self):
        manifest = json.loads((ROOT/'baseline/source_manifest.json').read_text(encoding='utf-8'))
        for name, rec in manifest.items():
            self.assertEqual(hashlib.sha256((ROOT.parent/name).read_bytes()).hexdigest(), rec['sha256'])


def main():
    build()
    directory = ROOT/'runs'/('core_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    MaterialTests.core = Core(directory)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MaterialTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    summary = {'tests_run': result.testsRun, 'passed': result.wasSuccessful(),
               'failures': [str(t) for t, _ in result.failures],
               'errors': [str(t) for t, _ in result.errors],
               'build_directory': str(directory), 'metrics': MaterialTests.metrics,
               'scope': '实际Fortran核心及原版回归，非Abaqus求解验证'}
    (ROOT/'validation/core_tests.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
