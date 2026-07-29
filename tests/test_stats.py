"""Stats 模块单元测试（使用模拟数据）"""

import os

import pytest

from local_frequency.stats.dnafusion_stat import stat_dnafusion
from local_frequency.stats.rnafusion_stat import stat_rnafusion
from local_frequency.stats.snv_stat import stat_snv


def _write_xls(dir_path: str, filename: str, data: list[dict]) -> str:
    import pandas as pd

    df = pd.DataFrame(data)
    path = os.path.join(dir_path, filename)
    df.to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def snv_data_dir(tmp_path):
    """创建模拟 SNV 数据"""
    data_dir = tmp_path / "snv"
    data_dir.mkdir()

    data = [
        {
            "Chr": "chr1",
            "Pos": 100,
            "End_pos": 100,
            "Ref": "A",
            "Alt": "T",
            "Gene.refGeneWithVer": "GENE1:transcript1",
            "Func.refGeneWithVer": "exonic",
            "ExonicFunc.refGeneWithVer": "nonsynonymous SNV",
            "AAChange.refGeneWithVer": "GENE1:transcript1:exon1:c.100A>T:p.Lys34Ter",
            "TAF": "0.5",
            "new_gene": None,
            "new_transcript": None,
            "new_exon": None,
            "new_cdna": None,
        },
    ]

    _write_xls(data_dir, "sample1.hg38_multianno.filter.xls", data)
    _write_xls(data_dir, "sample2.hg38_multianno.filter.xls", data)
    return str(data_dir)


@pytest.fixture
def fusion_data_dir(tmp_path):
    """创建模拟 Fusion 数据"""
    data_dir = tmp_path / "fusion"
    data_dir.mkdir()

    data = [
        {
            "gene1_chr": "chr1",
            "gene1_pos": 1000,
            "gene1": "GENE1",
            "gene2_chr": "chr2",
            "gene2_pos": 2000,
            "gene2": "GENE2",
            "exon": "1-3",
            "freq": "0.3",
        },
    ]

    _write_xls(data_dir, "sample1.total.fusion.xls", data)
    _write_xls(data_dir, "sample2.total.fusion.xls", data)
    return str(data_dir)


@pytest.fixture
def rna_fusion_data_dir(tmp_path):
    """创建模拟 RNA Fusion 数据"""
    data_dir = tmp_path / "rnafusion"
    data_dir.mkdir()

    data = [
        {
            "gene1_chr": "chr1",
            "gene1_pos": 1000,
            "gene1": "GENE1",
            "gene2_chr": "chr2",
            "gene2_pos": 2000,
            "gene2": "GENE2",
            "exon": "1-3",
            "freq": "0.3",
        },
    ]

    _write_xls(data_dir, "sample1.tsv.redup.xls", data)
    return str(data_dir)


class TestStatSNV:
    def test_stat_snv_basic(self, snv_data_dir, tmp_path):
        out_dir = str(tmp_path / "out")
        result = stat_snv(snv_data_dir, out_dir, sample_prefix="84")

        assert result["total"] == 2
        assert result["output_file"]
        assert os.path.exists(result["output_file"])
        assert "84_snv_somatic" in result["output_file"]

        # 检查输出内容
        with open(result["output_file"]) as f:
            lines = f.readlines()
        assert lines[0].startswith("Chr\t")
        # 2 个样本，每个变异出现 2 次
        assert "2\t1.0" in lines[1]

    def test_stat_snv_no_files(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        result = stat_snv(str(empty_dir), str(tmp_path / "out"), sample_prefix="84")
        assert result["total"] == 0
        assert result["output_file"] == ""


class TestStatDnaFusion:
    def test_stat_dnafusion_basic(self, fusion_data_dir, tmp_path):
        out_dir = str(tmp_path / "out")
        result = stat_dnafusion(fusion_data_dir, out_dir, sample_prefix="84")

        assert result["total"] == 2
        assert result["output_file"]
        assert os.path.exists(result["output_file"])
        assert "84_dnafusion" in result["output_file"]

        with open(result["output_file"]) as f:
            lines = f.readlines()
        assert lines[0].startswith("gene1_chr\t")

    def test_stat_dnafusion_no_files(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        result = stat_dnafusion(str(empty_dir), str(tmp_path / "out"), sample_prefix="84")
        assert result["total"] == 0


class TestStatRnaFusion:
    def test_stat_rnafusion_basic(self, rna_fusion_data_dir, tmp_path):
        out_dir = str(tmp_path / "out")
        result = stat_rnafusion(rna_fusion_data_dir, out_dir, sample_prefix="84")

        assert result["total"] == 1
        assert result["output_file"]
        assert os.path.exists(result["output_file"])
        assert "84_rnafusion" in result["output_file"]

    def test_stat_rnafusion_no_files(self, tmp_path):
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        result = stat_rnafusion(str(empty_dir), str(tmp_path / "out"), sample_prefix="84")
        assert result["total"] == 0
