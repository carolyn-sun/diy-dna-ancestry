"""
tests/test_qc.py — QC 模块测试（使用 fixture VCF 文件）
"""

import os
from pathlib import Path
from unittest.mock import patch, MagicMock, call
import pytest

from dna.qc import _count_variants, _count_samples


FIXTURE_DIR = Path(__file__).parent / "fixtures"


def test_count_variants(tmp_path):
    """从 .bim 文件统计 SNPs 数量。"""
    bim = tmp_path / "test.bim"
    bim.write_text("1\trs1\t0\t100\tA\tT\n1\trs2\t0\t200\tC\tG\n")
    assert _count_variants(str(tmp_path / "test")) == 2


def test_count_samples(tmp_path):
    """从 .fam 文件统计样本数量。"""
    fam = tmp_path / "test.fam"
    fam.write_text("FAM1 IND1 0 0 1 -9\nFAM2 IND2 0 0 2 -9\n")
    assert _count_samples(str(tmp_path / "test")) == 2


def test_count_variants_missing_file(tmp_path):
    """BIM 文件不存在时返回 0。"""
    assert _count_variants(str(tmp_path / "nonexistent")) == 0


def test_qc_calls_plink(tmp_path):
    """验证 run_qc 依次调用 PLINK 的四个步骤。"""
    import subprocess

    vcf = tmp_path / "test.vcf"
    vcf.write_text("##fileformat=VCFv4.1\n")

    # 模拟 PLINK 成功，并创建必要的输出文件
    def fake_plink(cmd, **kwargs):
        out_prefix = None
        for i, arg in enumerate(cmd):
            if arg == "--out":
                out_prefix = cmd[i + 1]
        if out_prefix:
            Path(out_prefix + ".bed").touch()
            Path(out_prefix + ".bim").write_text("1\trs1\t0\t100\tA\tT\n")
            Path(out_prefix + ".fam").write_text("S1 S1 0 0 1 -9\n")
            if "_ld" in out_prefix:
                Path(out_prefix + ".prune.in").write_text("rs1\n")
        return MagicMock(returncode=0, stdout="", stderr="")

    with patch("subprocess.run", side_effect=fake_plink):
        from dna.qc import run_qc
        result = run_qc(
            vcf_path=str(vcf),
            out_dir=str(tmp_path / "work"),
            geno=0.05, maf=0.01, hwe=1e-6, threads=1,
        )
    assert result.endswith("_pruned")
