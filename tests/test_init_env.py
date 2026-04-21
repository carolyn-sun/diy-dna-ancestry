"""
tests/test_init_env.py — 基础环境检查测试
"""

import sys
from unittest.mock import patch, MagicMock
import pytest

from dna.init_env import _check_python, _check_conda, _check_tool


def test_check_python_ok():
    """当前 Python 应满足 3.11+ 要求（conda 环境已配置）。"""
    status = _check_python()
    assert status.ok, f"Python 版本不足 3.11：{status.version}"
    assert status.name == "Python"
    assert len(status.version) > 0


def test_check_python_version_string():
    """验证版本字符串格式。"""
    status = _check_python()
    parts = status.version.split(".")
    assert len(parts) >= 2
    major, minor = int(parts[0]), int(parts[1])
    assert (major, minor) >= (3, 11)


def test_check_conda_in_env(monkeypatch):
    """模拟已激活 conda 环境。"""
    monkeypatch.setenv("CONDA_DEFAULT_ENV", "dna-ancestry")
    monkeypatch.setenv("CONDA_PREFIX", "/opt/conda/envs/dna-ancestry")
    status = _check_conda()
    assert status.ok
    assert status.version == "dna-ancestry"


def test_check_conda_not_activated(monkeypatch):
    """模拟未激活 conda 环境。"""
    monkeypatch.delenv("CONDA_DEFAULT_ENV", raising=False)
    monkeypatch.delenv("CONDA_PREFIX", raising=False)
    status = _check_conda()
    assert not status.ok
    assert "conda activate" in status.note


def test_check_tool_not_found():
    """检测不存在的命令。"""
    status = _check_tool("nonexistent_tool_xyz")
    assert not status.ok
    assert "PATH" in status.note


def test_check_tool_found(tmp_path):
    """模拟命令存在。"""
    with patch("shutil.which", return_value="/usr/bin/plink"):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                returncode=0, stdout="PLINK v1.90b6.21", stderr=""
            )
            status = _check_tool("plink")
    assert status.ok
    assert status.path == "/usr/bin/plink"
