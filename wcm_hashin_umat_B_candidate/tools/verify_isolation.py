# -*- coding: utf-8 -*-
"""候选工程隔离核查：正式文件哈希保护、写路径逃逸测试、环境摘要。

用法：
  python tools/verify_isolation.py --freeze   # 首次冻结正式文件哈希基线
  python tools/verify_isolation.py            # 复核哈希零变化 + 逃逸测试 + 环境摘要

逃逸测试的临时沙箱只建在候选 runs/ 下并自动清理；除该沙箱与
validation/ 报告外本工具不写任何文件。正式文件（PROJECT 下的
wcm_hashin_umat 源码、test 顶层 INP/Fortran、rebuild 几何输入）
在本候选工作期间必须保持字节零变化。
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from isolation import CAND_ROOT, PROJECT_ROOT, IsolationError, inside, require_inside

FORMAL_FILES = [
    'wcm_hashin_umat/src/wcm_core.for',
    'wcm_hashin_umat/src/umat_entry.for.in',
    'wcm_hashin_umat/tools/bootstrap.py',
    'wcm_hashin_umat/tools/build.py',
    'wcm_hashin_umat/tools/abaqus_verify.py',
    'wcm_hashin_umat/tools/convert_wcm_inp.py',
    'wcm_hashin_umat/tools/extract_failure.py',
    'wcm_hashin_umat/tools/model.py',
    'wcm_hashin_umat/tools/prepare_test_cases.py',
    'wcm_hashin_umat/tools/run_job.py',
    'wcm_hashin_umat/tools/run_tests.py',
    'wcm_hashin_umat/tests/baseline_driver.for',
    'wcm_hashin_umat/tests/extract_odb.py',
    'wcm_hashin_umat/config/project.json',
    'wcm_hashin_umat/config/failure_thresholds_A.json',
    'wcm_hashin_umat/config/failure_thresholds_B.json',
    'wcm_hashin_umat/baseline/hashin.for',
    'wcm_hashin_umat/baseline/hashin_constant.for',
    'wcm_hashin_umat/baseline/source_manifest.json',
    'test/A_stiffness_reduction/Rebuild_PDM_A.inp',
    'test/A_stiffness_reduction/hashin_constant_wcm.for',
    'test/B_energy_softening/Rebuild_PDM_B.inp',
    'test/B_energy_softening/hashin_energy_wcm.for',
    'geometry_zhuning_joint_rebuild/WCM_Job.inp',
    'geometry_zhuning_joint_rebuild/WCM_Zhuning_Joint_Rebuilt_WindAngles.ang',
]
CANDIDATE_COPIED = [
    ('src/wcm_core.for', 'wcm_hashin_umat/src/wcm_core.for'),
    ('src/umat_entry.for.in', 'wcm_hashin_umat/src/umat_entry.for.in'),
    ('tools/bootstrap.py', 'wcm_hashin_umat/tools/bootstrap.py'),
    ('tools/convert_wcm_inp.py', 'wcm_hashin_umat/tools/convert_wcm_inp.py'),
    ('tests/baseline_driver.for', 'wcm_hashin_umat/tests/baseline_driver.for'),
    ('config/project.json', 'wcm_hashin_umat/config/project.json'),
    ('config/failure_thresholds_A.json', 'wcm_hashin_umat/config/failure_thresholds_A.json'),
    ('config/failure_thresholds_B.json', 'wcm_hashin_umat/config/failure_thresholds_B.json'),
    ('baseline/hashin.for', 'wcm_hashin_umat/baseline/hashin.for'),
    ('baseline/hashin_constant.for', 'wcm_hashin_umat/baseline/hashin_constant.for'),
    ('baseline/source_manifest.json', 'wcm_hashin_umat/baseline/source_manifest.json'),
]
MANIFEST = CAND_ROOT/'validation'/'source_manifest.json'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def formal_hashes():
    table = {}
    for rel in FORMAL_FILES:
        path = PROJECT_ROOT/rel
        if not path.is_file():
            raise FileNotFoundError('受保护的正式文件缺失：'+str(path))
        table[rel] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
    return table


def freeze():
    copied = {}
    for cand_rel, formal_rel in CANDIDATE_COPIED:
        cand = CAND_ROOT/cand_rel
        formal = sha256(PROJECT_ROOT/formal_rel)
        if cand.is_file():
            if sha256(cand) != formal:
                raise RuntimeError('候选副本与正式父版本不一致：'+cand_rel)
            copied[cand_rel] = {'parent_sha256': formal, 'parent': formal_rel}
        else:
            copied[cand_rel] = {'parent_sha256': formal, 'parent': formal_rel,
                                'note': '已按计划派生修改，不再与父版本一致'}
    manifest = {
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
        'project': str(PROJECT_ROOT),
        'candidate': str(CAND_ROOT),
        'parent_source': 'wcm_hashin_umat',
        'candidate_parent_version_note': '本候选由 PROJECT/wcm_hashin_umat 于 2026-09-27 复制隔离；'
                                        '此后仅候选内文件允许修改。',
        'formal_protected': formal_hashes(),
        'candidate_copied_provenance': copied,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n',
                        encoding='utf-8')
    print('已冻结正式文件哈希基线：', MANIFEST)
    return manifest


def check_formal(baseline):
    changed, missing = [], []
    for rel, rec in baseline['formal_protected'].items():
        path = PROJECT_ROOT/rel
        if not path.is_file():
            missing.append(rel)
        elif sha256(path) != rec['sha256']:
            changed.append(rel)
    return changed, missing


def escape_tests():
    """写路径防逃逸测试：全部必须被 require_inside 拒绝/正确放行。"""
    results = []
    sandbox = CAND_ROOT/'runs'/('_isolation_tmp_'+datetime.datetime.now()
                                .strftime('%Y%m%d_%H%M%S_%f'))

    def record(name, raised, expect_reject, detail=''):
        ok = raised if expect_reject else not raised
        results.append({'case': name, 'rejected': raised,
                        'expected_reject': expect_reject, 'pass': ok, 'detail': detail})

    try:
        sandbox.mkdir(parents=True)
        # 1. .. 上跳逃出候选根（runs/沙箱→runs→CAND→PROJECT）：必须拒绝。
        try:
            require_inside(sandbox/'..'/'..'/'..'/'escaped.txt')
            record('dotdot_escape', False, True)
        except IsolationError:
            record('dotdot_escape', True, True)
        # 2. .. 上跳但仍留在候选内：必须放行。
        try:
            require_inside(sandbox/'..'/'ok.txt')
            record('dotdot_inside', False, False)
        except IsolationError:
            record('dotdot_inside', True, False)
        # 3. 目录联接指向候选外：真实路径重解析后必须拒绝。
        link = sandbox/'junction_out'
        target = str(PROJECT_ROOT/'test')
        made = False
        try:
            subprocess.run(['cmd', '/c', 'mklink', '/J', str(link), target],
                           capture_output=True, check=True)
            made = True
        except (subprocess.CalledProcessError, FileNotFoundError):
            try:
                import _winapi
                _winapi.CreateJunction(target, str(link))
                made = True
            except Exception:
                record('junction_escape', False, True,
                       'SKIPPED：无法创建目录联接')
        if made:
            try:
                require_inside(link/'leaked.txt')
                record('junction_escape', False, True)
            except IsolationError:
                record('junction_escape', True, True)
        # 4. 字符串前缀陷阱：root=a 时路径 ab 不能靠前缀放行。
        try:
            require_inside(CAND_ROOT/'ab_cand_prefix_trap', root=CAND_ROOT/'a')
            record('prefix_trap', False, True)
        except IsolationError:
            record('prefix_trap', True, True)
        # 5. 候选外绝对路径：必须拒绝。
        try:
            require_inside(PROJECT_ROOT/'test'/'outside.txt')
            record('absolute_outside', False, True)
        except IsolationError:
            record('absolute_outside', True, True)
        # 6. 候选内普通路径：必须放行并返回重解析路径。
        resolved = require_inside(sandbox/'normal.txt')
        record('normal_inside', False, False, str(resolved))
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)
    return results


def environment_summary():
    info = {'python': sys.executable, 'python_version': sys.version,
            'platform': sys.platform}
    try:
        import numpy
        info['numpy_version'] = numpy.__version__
    except ImportError:
        info['numpy_version'] = None
    gfortran = shutil.which('gfortran')
    info['gfortran'] = gfortran
    if gfortran:
        try:
            out = subprocess.run([gfortran, '--version'], capture_output=True,
                                 text=True, timeout=60)
            info['gfortran_version'] = out.stdout.splitlines()[0] if out.stdout else None
        except Exception as exc:
            info['gfortran_version'] = 'ERROR: '+str(exc)
    abaqus = Path(r'E:\ABAQUS2025\Commands\abaqus.BAT')
    info['abaqus_bat'] = str(abaqus) if abaqus.is_file() else None
    info['cand_root'] = str(CAND_ROOT)
    info['project_root'] = str(PROJECT_ROOT)
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze', action='store_true',
                        help='冻结/刷新正式文件哈希基线')
    args = parser.parse_args()
    if args.freeze:
        manifest = freeze()
    else:
        if not MANIFEST.is_file():
            raise SystemExit('尚未冻结基线：先运行 --freeze')
        manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    changed, missing = check_formal(manifest)
    escapes = escape_tests()
    env = environment_summary()
    passed = not changed and not missing and all(e['pass'] for e in escapes)
    report = {
        'time': datetime.datetime.now().isoformat(timespec='seconds'),
        'formal_hash_zero_change': not changed and not missing,
        'formal_changed': changed, 'formal_missing': missing,
        'escape_tests': escapes,
        'escape_all_rejected': all(e['pass'] for e in escapes),
        'environment': env,
        'gate_G0_pass': passed,
    }
    out = CAND_ROOT/'validation'/'isolation_check.json'
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n',
                   encoding='utf-8')
    print(json.dumps({'gate_G0_pass': passed,
                      'formal_zero_change': report['formal_hash_zero_change'],
                      'escape_tests': {e['case']: e['pass'] for e in escapes},
                      'environment': {k: env.get(k) for k in
                                      ('python_version', 'numpy_version',
                                       'gfortran_version', 'abaqus_bat')},
                      'report': str(out)}, ensure_ascii=False, indent=2))
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
