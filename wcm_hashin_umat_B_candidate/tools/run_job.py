# -*- coding: utf-8 -*-
"""隔离运行已转换的气瓶模型；每次生成独立作业目录。

在wcm_hashin_umat中使用：
  python tools/run_job.py --variant A --prepare-only
  python tools/run_job.py --variant B --datacheck
  python tools/run_job.py --variant B --cpus 4

最后一条会实际提交原250MPa工况，可能耗时很长；不保证爆破预测有效。
默认不自动改变材料、时间步、载荷、内衬属性或稳定化。
Git回退时仅回退源代码/配置；生成模型应重新转换，不重用不兼容SDV。
原始基线标签：baseline-original；回看用git show，撤销提交优先git revert。
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
from build import build
from model import validate_props
from convert_wcm_inp import blocks

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', choices=('A', 'B'), required=True)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--cpus', type=int, default=2)
    parser.add_argument('--abaqus', default=shutil.which('abaqus'))
    parser.add_argument('--datacheck', action='store_true')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.cpus < 1:
        parser.error('cpus必须为正整数')
    source = (args.input or ROOT/'generated'/f'Job_wcm_{args.variant}.inp').resolve()
    from isolation import require_inside, IsolationError
    try:
        require_inside(source, what='run_job输入')
    except IsolationError:
        parser.error('只允许运行候选内已转换的inp，防止误用旧材料卡或逃逸写路径')
    content = source.read_bytes()
    materials = [b for b in blocks(content.decode('latin-1')) if b.key == 'user material']
    if not materials:
        parser.error('未找到用户材料')
    for block in materials:
        validate_props(args.variant, block.numbers())
    build()
    name = 'hashin_constant_wcm.for' if args.variant == 'A' else 'hashin_energy_wcm.for'
    run = require_inside(ROOT/'runs'/('tank_'+args.variant+'_'+datetime.datetime.now()
                         .strftime('%Y%m%d_%H%M%S_%f')), what='run_job作业目录')
    run.mkdir(parents=True)
    (run/'scratch').mkdir()
    (run/'tank.inp').write_bytes(content)
    shutil.copyfile(ROOT/'dist'/name, run/name)
    command = [args.abaqus or 'abaqus', 'job=tank', 'input=tank.inp', 'user='+name,
               'cpus='+str(args.cpus), 'scratch='+str(run/'scratch'), 'interactive']
    if args.datacheck:
        command.append('datacheck')
    manifest = {'input_source': str(source), 'variant': args.variant, 'command': command,
                'cwd': str(run), 'input_sha256': hashlib.sha256(content).hexdigest(),
                'umat_sha256': hashlib.sha256((run/name).read_bytes()).hexdigest()}
    (run/'run_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print('作业目录：', run, flush=True)
    print('命令：', subprocess.list2cmdline(command), flush=True)
    if args.prepare_only:
        return 0
    if not args.abaqus:
        raise RuntimeError('没有找到Abaqus，可通过--abaqus指定启动器')
    with (run/'launcher.txt').open('wb') as log:
        result = subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT)
    # Windows启动批处理可能掩盖子进程退出码，必须同时检查完成标志。
    text = (run/'launcher.txt').read_text(errors='replace')
    if result.returncode or 'Abaqus JOB tank COMPLETED' not in text:
        print('作业未正常完成，请检查launcher.txt及tank.msg/dat/sta。')
        return 1
    print('Abaqus作业完成；求解完成不等于材料模型已通过实验验证。')
    return 0


if __name__ == '__main__':
    sys.exit(main())
