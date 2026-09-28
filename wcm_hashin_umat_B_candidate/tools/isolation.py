# -*- coding: utf-8 -*-
"""候选工程写路径隔离：规范化目录包含检查。

所有写文件的入口必须经 require_inside 校验。检查基于 os.path.realpath
重解析后的真实路径（解析符号链接/目录联接/..），不做字符串前缀比较。
"""
import os
from pathlib import Path


class IsolationError(PermissionError):
    """写路径逃出允许的根目录或经过重解析点。"""


def _real(path):
    return Path(os.path.realpath(str(path)))


CAND_ROOT = _real(Path(__file__).resolve().parents[1])
PROJECT_ROOT = _real(CAND_ROOT.parent)


def inside(path, root=None):
    """path 重解析后是否位于 root 重解析路径之内（含根本身）。"""
    root = CAND_ROOT if root is None else _real(root)
    target = _real(path)
    try:
        target.relative_to(root)
    except ValueError:
        return False
    return True


def require_inside(path, root=None, what=''):
    """拒绝逃出 root 的写路径；通过时返回重解析后的 Path。"""
    root = CAND_ROOT if root is None else _real(root)
    if not inside(path, root):
        raise IsolationError(
            f'写路径逃出候选根目录{root}：{path}{("（" + what + "）") if what else ""}')
    return _real(path)


def ensure_writable_dir(path, root=None, what=''):
    """校验目录在 root 内并确保存在；返回重解析 Path。"""
    path = require_inside(path, root, what)
    path.mkdir(parents=True, exist_ok=True)
    return path
