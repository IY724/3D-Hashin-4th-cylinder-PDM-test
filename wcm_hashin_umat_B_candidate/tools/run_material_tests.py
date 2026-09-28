# -*- coding: utf-8 -*-
"""材料点测试编排器（规划4.3接口）。

用法：
  python tools/run_material_tests.py --suite reproduce|onset|matrix|energy|all

每次调用创建唯一运行目录 runs/material_<suite>_<stamp>/：gfortran 编译
产物与全套逐点历史CSV/JSON都落在其中；证据摘要写 validation/ 下对应
JSON（reproduction/onset/matrix/energy）。构建走 build.py 的互斥锁。

套件与阶段对应：
  reproduce  阶段1  严格版1033重编译复现+旧版静默跳过对照+解析基准
  onset      阶段2  起始定位修正的独立参考、A回归、状态回滚
  matrix     阶段3A 材料点验证矩阵
  energy     阶段3B 能量闭合账本与解析负对照
"""
import argparse
import datetime
import hashlib
import importlib
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
TOOLS = Path(__file__).resolve().parent
TESTS = TOOLS.parent/'tests'
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TESTS))
from isolation import ensure_writable_dir, require_inside  # noqa: E402

ROOT = TOOLS.parent
SUITES = {
    'reproduce': ('reproduce_suite', 'validation/reproduction.json'),
    'onset': ('test_onset', 'validation/onset.json'),
    'matrix': ('test_matrix', 'validation/matrix.json'),
    'energy': ('test_energy', 'validation/energy.json'),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=tuple(SUITES)+('all',), required=True)
    parser.add_argument('--keep-going', action='store_true',
                        help='all模式下单个套件失败不中断其余套件')
    args = parser.parse_args()
    names = list(SUITES) if args.suite == 'all' else [args.suite]
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    results = {}
    for name in names:
        module_name, report_rel = SUITES[name]
        run_dir = ensure_writable_dir(
            ROOT/'runs'/f'material_{name}_{stamp}', what='材料点运行目录')
        try:
            module = importlib.import_module(module_name)
        except ImportError as exc:
            results[name] = {'status': 'BLOCKED_ENV', 'error': f'套件模块未就绪：{exc}'}
            continue
        from run_tests import Core
        core = Core(run_dir/'build')
        evidence = module.run(core, run_dir)
        evidence.setdefault('suite', name)
        evidence['run_dir'] = str(run_dir)
        evidence['core_src_sha256'] = hashlib.sha256(
            (ROOT/'src/wcm_core.for').read_bytes()).hexdigest()
        evidence['dist_umat_sha256'] = hashlib.sha256(
            (ROOT/'dist/hashin_energy_wcm.for').read_bytes()).hexdigest()
        report = require_inside(ROOT/report_rel, what=f'{name}证据')
        report.write_text(json.dumps(evidence, ensure_ascii=False, indent=2)+'\n',
                          encoding='utf-8')
        results[name] = {'status': evidence.get('status', 'UNKNOWN'),
                         'report': str(report), 'run_dir': str(run_dir)}
        print(json.dumps({name: results[name]}, ensure_ascii=False))
        if results[name]['status'] not in ('PASS', 'LENGTH_LIMITED',
                                           'EXPECTED_INCOMPATIBLE'):
            if not args.keep_going and args.suite == 'all':
                break
    print(json.dumps({'suite': args.suite, 'results': results},
                     ensure_ascii=False, indent=2))
    ok = all(v['status'] in ('PASS', 'LENGTH_LIMITED', 'EXPECTED_INCOMPATIBLE')
             for v in results.values())
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
