"""FileCopier 模块单元测试"""

import os
from unittest.mock import MagicMock, patch

import pytest

from local_frequency.file_copier import FileCopier


@pytest.fixture
def mock_sftp():
    return MagicMock()


@pytest.fixture
def source_dict():
    return {
        "proj001": {
            ".filter.xls": ["/data/a/b/proj001/sample1.filter.xls"],
            ".cnv.vcf": ["/data/a/b/proj001/sample1.cnv.vcf"],
        },
        "proj002": {
            ".total.fusion.xls": ["/data/a/b/proj002/sample2.total.fusion.xls"],
        },
    }


class TestFileCopier:
    @patch("os.makedirs")
    def test_copy_single(self, mock_makedirs, mock_sftp, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))
        local = str(tmp_path / "proj001" / "data" / "sample.filter.xls")

        result = copier.copy_single("/remote/sample.filter.xls", local)

        assert result is True
        mock_sftp.get.assert_called_once()

    @patch("os.makedirs")
    def test_copy_single_failure(self, mock_makedirs, mock_sftp, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))
        mock_sftp.get.side_effect = Exception("Permission denied")
        local = str(tmp_path / "proj001" / "data" / "fail.xls")

        result = copier.copy_single("/remote/fail.xls", local)

        assert result is False

    def test_copy_files_all(self, mock_sftp, source_dict, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))

        result = copier.copy_files(source_dict, "/data")

        assert len(result) == 3
        # 每个文件都应有对应的本地路径
        for local in result.values():
            assert "proj001" in local or "proj002" in local

    def test_copy_files_project_filter(self, mock_sftp, source_dict, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))

        result = copier.copy_files(source_dict, "/data", project_ids=["proj002"])

        assert len(result) == 1
        assert "proj002" in next(iter(result.values()))

    def test_skip_existing(self, mock_sftp, source_dict, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))
        # 预创建文件模拟已存在
        existing = tmp_path / "proj001" / "data" / "sample1.filter.xls"
        existing.parent.mkdir(parents=True, exist_ok=True)
        existing.touch()

        result = copier.copy_files(source_dict, "/data")

        # 该文件应被跳过
        assert "/data/a/b/proj001/sample1.filter.xls" not in result

    def test_resolve_local_path(self, mock_sftp, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))
        path = copier._resolve_local_path(
            "/data/a/b/proj001/sample.filter.xls", "/data", "proj001"
        )
        assert path.endswith(os.path.join("proj001", "data", "sample.filter.xls"))

    def test_skip_copy_existing(self, mock_sftp, tmp_path):
        copier = FileCopier(mock_sftp, str(tmp_path))
        existing = tmp_path / "exists.xls"
        existing.touch()

        assert copier._skip_copy(str(existing)) is True
        assert copier._skip_copy(str(tmp_path / "no_exist.xls")) is False
