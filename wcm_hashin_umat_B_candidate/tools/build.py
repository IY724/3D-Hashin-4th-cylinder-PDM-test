# -*- coding: utf-8 -*-
"""候选工程内生成可单独提交的两版 UMAT；不读取或修改旧源码。

dist 只写候选目录；构建用互斥锁串行化，避免并发构建互相覆盖
同名摘要（规划4.2：构建/写同名摘要的任务必须串行）。
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

ROOT = Path(__file__).resolve().parents[1]


def props_layout(mode):
    """逐项列出输入；实际命名赋值见同文件 WCM_CHECK。"""
    common_names = (
        'E11', 'E22', 'E33', 'PR12', 'PR13', 'PR23',
        'G12', 'G13', 'G23', 'XT', 'XC', 'YT', 'YC',
        'S12', 'S13', 'S23',
    )
    if mode == 1:
        specific_names = (
            'DFT0', 'DFC0', 'DMT0', 'DMC0', 'SMT', 'SMC',
            'THETA', 'WPLUS', 'UPDATE',
        )
    else:
        specific_names = (
            'GFT', 'GFC', 'GMT', 'GMC', 'ALPHA', 'SMT', 'SMC',
            'ETA_F', 'ETA_M', 'THETA', 'WPLUS', 'CAPF', 'CAPM',
            'UPDATE',
        )
    names = common_names + specific_names
    return '\n'.join(
        ['C 完整材料参数表；下方 WCM_CHECK 包含逐项可执行赋值：']
        + [f'C {name}=PROPS({index})' for index, name in enumerate(names, 1)]
        + ['C Units: modulus/strength MPa; GFT/GFC/GMT/GMC N/mm;',
           'C THETA deg; ETA_F/ETA_M step-time units.',
           'C S13=PROPS(15) 保留；当前 tau13 判据使用 S12。']
        + (['C B: GFT/GFC/GMT/GMC 为峰后软化功，非全程总断裂能。',
            'C B: deltaf=delta0+2G/sigma0；禁止混用旧 restart。']
           if mode == 2 else [])
    )


def check_fixed_format(text):
    for number, line in enumerate(text.splitlines(), 1):
        if not line or line[0] in 'Cc*!':
            continue
        if '\t' in line or len(line) > 72 or not line.isascii():
            raise ValueError(f'Fortran固定格式不合法：第{number}行 {line!r}')


class build_lock:
    """候选内互斥锁：独占创建式，最多等10分钟。"""

    def __init__(self, path):
        self.path = Path(path)
        self.handle = None

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        deadline = time.time() + 600
        while True:
            try:
                self.handle = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(self.handle, str(os.getpid()).encode())
                return self
            except FileExistsError:
                if time.time() > deadline:
                    raise RuntimeError('构建锁超时：'+str(self.path))
                time.sleep(0.2)

    def __exit__(self, *exc):
        if self.handle is not None:
            os.close(self.handle)
            try:
                self.path.unlink()
            except OSError:
                pass
        return False


def render(mode):
    """构建和交付共用同一渲染入口，避免 .for 与源码脱节。"""
    template = (ROOT / 'src/umat_entry.for.in').read_text(encoding='utf-8')
    core = (ROOT / 'src/wcm_core.for').read_text(encoding='utf-8')
    variant = 'A（直接刚度折减）' if mode == 1 else 'B（峰后跨度软化）'
    text = (template.replace('@MODE@', str(mode))
            .replace('@VERSION@', variant)
            .replace('C @PROPS_LAYOUT@', props_layout(mode)) + '\n' + core)
    check_fixed_format(text)
    return text


def build(variants=('A', 'B')):
    from isolation import require_inside
    dist = require_inside(ROOT / 'dist', what='build输出')
    dist.mkdir(exist_ok=True)
    records = {}
    with build_lock(ROOT/'validation'/'.build.lock'):
        for mode, variant, name in ((1, 'A', 'hashin_constant_wcm.for'),
                                    (2, 'B', 'hashin_energy_wcm.for')):
            if variant not in variants:
                continue
            text = render(mode)
            path = require_inside(ROOT / 'dist' / name, what='build产物')
            path.write_text(text, encoding='utf-8', newline='\n')
            records[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest = ROOT / 'dist/build_manifest.json'
        previous = json.loads(manifest.read_text(encoding='utf-8')) if manifest.exists() else {}
        previous.update(records)
        manifest.write_text(json.dumps(previous, indent=2) + '\n', encoding='utf-8')
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', choices=('A', 'B', 'all'), default='all')
    args = parser.parse_args()
    variants = ('A', 'B') if args.variant == 'all' else (args.variant,)
    print(json.dumps(build(variants), indent=2))
