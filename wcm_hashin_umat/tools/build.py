# -*- coding: utf-8 -*-
"""在新仓库中生成可单独提交的两版 UMAT；不读取或修改旧源码。"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def props_layout(mode):
    """Short, variant-specific guide placed at the top of each .for file."""
    common = [
        'C Elastic: E11/E22/E33=1:3, PR12/PR13/PR23=4:6,',
        'C          G12/G13/G23=7:9 (MPa; Poisson ratios unitless).',
        'C Strength: XT/XC/YT/YC=10:13, S12/S13/S23=14:16 (MPa).',
    ]
    if mode == 1:
        specific = [
            'C A: DFT0/DFC0/DMT0/DMC0=17:20 (fixed loss fractions).',
            'C    SMT/SMC=21:22 (matrix shear loss factors).',
            'C    THETA/WPLUS/UPDATE=23:25 (deg, fraction, 0 or 1).',
        ]
    else:
        specific = [
            'C B: GFT/GFC/GMT/GMC=17:20 (fracture energy, N/mm).',
            'C    ALPHA=21 (reserved), SMT/SMC=22:23 (shear factors).',
            'C    ETA_F/ETA_M=24:25 (viscosity in step-time units).',
            'C    THETA/WPLUS=26:27 (deg, fraction).',
            'C    CAPF/CAPM=28:29 (stiffness loss limits, 0..1).',
            'C    UPDATE=30 (0: next increment, 1: current feedback).',
        ]
    return '\n'.join(common + specific)


def check_fixed_format(text):
    for number, line in enumerate(text.splitlines(), 1):
        if not line or line[0] in 'Cc*!':
            continue
        if '\t' in line or len(line) > 72 or not line.isascii():
            raise ValueError(f'Fortran固定格式不合法：第{number}行 {line!r}')


def build():
    template = (ROOT / 'src/umat_entry.for.in').read_text(encoding='utf-8')
    core = (ROOT / 'src/wcm_core.for').read_text(encoding='utf-8')
    (ROOT / 'dist').mkdir(exist_ok=True)
    records = {}
    for mode, variant, name in ((1, 'A（直接刚度折减）', 'hashin_constant_wcm.for'),
                                (2, 'B（能量线性软化）', 'hashin_energy_wcm.for')):
        text = (template.replace('@MODE@', str(mode))
                .replace('@VERSION@', variant)
                .replace('C @PROPS_LAYOUT@', props_layout(mode)) + '\n' + core)
        check_fixed_format(text)
        path = ROOT / 'dist' / name
        path.write_text(text, encoding='utf-8', newline='\n')
        records[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    (ROOT / 'dist/build_manifest.json').write_text(
        json.dumps(records, indent=2) + '\n', encoding='utf-8')
    return records


if __name__ == '__main__':
    print(json.dumps(build(), indent=2))
