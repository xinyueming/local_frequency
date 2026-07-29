"""DB Builder 模块单元测试"""

import json
import os

import pytest

from local_frequency.db_builder import (
    FILENAME_TO_KEY,
    _parse_taf_values,
    avinput_to_txt,
    build_frequency_db,
    vcf_to_avinput,
    xls_to_vcf,
)


@pytest.fixture
def mock_xls(tmp_path):
    """创建模拟 mutation_frequency.xls"""
    xls = tmp_path / "test.mutation_frequency.xls"
    xls.write_text(
        "Chr\tPos\tEnd_pos\tRef\tAlt\tGene.refGeneWithVer\tFunc.refGeneWithVer\t"
        "ExonicFunc.refGeneWithVer\tAAChange.refGeneWithVer\tpopulation\t"
        "mutation_frequency\tTAF\n"
        "chr1\t100\t100\tA\tT\tGENE1:NM_001\texonic\tnonsynonymous SNV\t"
        "GENE1:NM_001:exon1:c.100A>T\t2\t1.0\t50%,50%\n"
    )
    return str(xls)


class TestParseTafValues:
    def test_percentage_values(self):
        assert _parse_taf_values("50%,30%,20%") == pytest.approx([0.5, 0.3, 0.2])

    def test_plain_floats(self):
        assert _parse_taf_values("0.5,0.3") == pytest.approx([0.5, 0.3])

    def test_empty(self):
        assert _parse_taf_values("") == []

    def test_dot_value(self):
        assert _parse_taf_values(".,50%") == pytest.approx([0.5])


class TestXlsToVcf:
    def test_xls_to_vcf(self, mock_xls, tmp_path):
        vcf = xls_to_vcf(mock_xls, "test_version", 10)
        assert vcf.endswith(".vcf")
        assert os.path.exists(vcf)

        with open(vcf) as f:
            lines = f.readlines()
        assert lines[0].startswith("#CHROM\t")
        assert "LOCAL_FREQUENCY=" in lines[1]
        assert "AAChange=" in lines[1]


class TestVcfToAvinput:
    def test_vcf_to_avinput(self, mock_xls, tmp_path):
        vcf = xls_to_vcf(mock_xls, "test_version", 10)
        avinput = vcf_to_avinput(vcf)

        assert avinput.endswith(".avinput")
        assert os.path.exists(avinput)
        with open(avinput) as f:
            content = f.read()
        assert "chr1" in content


class TestAvinputToTxt:
    def test_avinput_to_txt(self, mock_xls, tmp_path):
        vcf = xls_to_vcf(mock_xls, "test_version", 10)
        avinput = vcf_to_avinput(vcf)
        txt = avinput_to_txt(avinput)

        assert txt.endswith(".raw.txt")
        assert os.path.exists(txt)
        with open(txt) as f:
            header = f.readline()
        assert header.startswith("#Chr\t")
        assert "LOCAL_FREQUENCY" in header


class TestBuildFrequencyDb:
    def test_build_db(self, tmp_path):
        # 创建输入目录
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        # 使用已知的文件名
        xls = input_dir / "84.mutation_frequency.xls"
        xls.write_text(
            "Chr\tPos\tEnd_pos\tRef\tAlt\tGene\tFunc\tExonicFunc\t"
            "AAChange\tpopulation\tmutation_frequency\tTAF\n"
            "chr1\t100\t100\tA\tT\tGENE1\texonic\tcoding\t"
            "GENE1:c.100A>T\t1\t1.0\t0.5\n"
        )

        # 创建 stat.json
        stat = {"84": 1}
        stat_path = tmp_path / "stat.json"
        stat_path.write_text(json.dumps(stat))

        output_dir = tmp_path / "output"

        results = build_frequency_db(
            str(input_dir), str(output_dir), str(stat_path), processes=1
        )

        # 应有 .txt + .idx 两个文件
        assert len(results) >= 2

    def test_build_db_no_stat(self, tmp_path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        xls = input_dir / "84.mutation_frequency.xls"
        xls.write_text(
            "Chr\tPos\tEnd_pos\tRef\tAlt\tGene\tFunc\tExonicFunc\t"
            "AAChange\tpopulation\tmutation_frequency\tTAF\n"
            "chr1\t100\t100\tA\tT\tG\tE\tC\tA\t1\t1.0\t0.5\n"
        )

        output_dir = tmp_path / "output"
        results = build_frequency_db(str(input_dir), str(output_dir), processes=1)
        assert len(results) >= 2

    def test_build_db_empty_dir(self, tmp_path):
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        output_dir = tmp_path / "output"

        results = build_frequency_db(str(input_dir), str(output_dir), processes=1)
        assert results == []

    def test_filename_mapping_complete(self):
        # 确保所有常用模式都有映射
        assert "84" in FILENAME_TO_KEY.values()
        assert "624" in FILENAME_TO_KEY.values()
        assert "wes" in FILENAME_TO_KEY.values()
