# -*- coding: utf-8 -*-
"""仅首次执行：归档旧源码和输入哈希，不修改原项目。"""
from pathlib import Path
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent


def main():
    for name in ('baseline', 'src', 'dist', 'tests', 'config', 'generated', 'runs', 'validation'):
        (ROOT / name).mkdir(exist_ok=True)
    manifest = {}
    for relative in ('hashin_self/hashin.for', 'hashin_self/hashin_constant.for',
                     'Job-hashinnewcdm.inp', 'Job-hashinnewconstant.inp'):
        source = SOURCE / relative
        content = source.read_bytes()
        manifest[relative] = {'sha256': hashlib.sha256(content).hexdigest(), 'bytes': len(content)}
        if source.suffix == '.for':
            target = ROOT / 'baseline' / source.name
            if target.exists() and target.read_bytes() != content:
                raise RuntimeError('基线已存在且内容不同，拒绝覆盖：' + str(target))
            if not target.exists():
                shutil.copyfile(source, target)
    target = ROOT / 'baseline' / 'source_manifest.json'
    text = json.dumps(manifest, ensure_ascii=False, indent=2) + '\n'
    if target.exists() and target.read_text(encoding='utf-8') != text:
        raise RuntimeError('原项目已变化，拒绝更新原始基线。')
    target.write_text(text, encoding='utf-8')
    print('基线归档完成：', ROOT)


if __name__ == '__main__':
    main()
