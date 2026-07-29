"""Pipeline 模块单元测试"""

from unittest.mock import MagicMock, patch

import pytest

from local_frequency.pipeline import Pipeline, PipelineConfig, create_pipeline


@pytest.fixture
def pipeline_config():
    return PipelineConfig(
        ssh={"host": "test.example.com", "port": 22, "user": "testuser"},
        remote={"base_path": "/data"},
        local={"base_path": "/tmp/local_data"},
        stats={"max_workers": 2},
        db_builder={"processes": 1},
    )


class TestPipelineDryRun:
    def test_dry_run_success(self, pipeline_config):
        p = Pipeline(pipeline_config)
        assert p.run(dry_run=True) is True

    def test_dry_run_missing_host(self, pipeline_config):
        pipeline_config.ssh = {}
        p = Pipeline(pipeline_config)
        assert p.run(dry_run=True) is False

    def test_dry_run_missing_remote(self, pipeline_config):
        pipeline_config.remote = {}
        p = Pipeline(pipeline_config)
        assert p.run(dry_run=True) is False

    def test_dry_run_missing_local(self, pipeline_config):
        pipeline_config.local = {}
        p = Pipeline(pipeline_config)
        assert p.run(dry_run=True) is False


class TestPipelineCreate:
    def test_create_pipeline_from_config(self, tmp_path):
        cfg_path = tmp_path / "config.yaml"
        cfg_path.write_text("ssh:\n  host: h\n  user: u\nremote:\n  base_path: /remote\nlocal:\n  base_path: /local\n")

        p = create_pipeline(str(cfg_path))
        assert isinstance(p, Pipeline)
        assert p.config.ssh["host"] == "h"


class TestPipelineCollectOnly:
    @patch("local_frequency.pipeline.SSHClient")
    def test_collect_only(self, MockClient, pipeline_config, tmp_path):
        mock_instance = MagicMock()
        mock_instance.__enter__ = MagicMock(return_value=mock_instance)
        mock_instance.__exit__ = MagicMock(return_value=False)
        mock_instance.walk_remote_dir.return_value = [
            "/data/a/b/c/d/e/proj001/sample.filter.xls",
        ]
        MockClient.return_value = mock_instance

        p = Pipeline(pipeline_config)
        # 模拟 collect_only: 只执行前两步
        # 但我们的 run(collect_only=True) 只执行 step 1 和 2
        # step 2 需要 sftp，mock 它
        mock_instance.sftp = MagicMock()

        # 由于 collect_only=True 只做收集+拷贝，
        # 而拷贝会因无实际连接失败，我们测试 step_collect 单独调用
        result = p._step_collect()
        assert result is True
        assert "proj001" in p.source_dict


class TestPipelineParallelStats:
    @patch("local_frequency.pipeline.stat_snv")
    @patch("local_frequency.pipeline.stat_dnafusion")
    @patch("local_frequency.pipeline.stat_rnafusion")
    def test_parallel_stats(self, mock_rna, mock_dna, mock_snv, pipeline_config):
        mock_snv.return_value = {"total": 5}
        mock_dna.return_value = {"total": 3}
        mock_rna.return_value = {"total": 2}

        p = Pipeline(pipeline_config)
        p.source_dict = {"proj001": {".filter.xls": ["/x"]}}
        result = p._step_parallel_stats(project_ids=["proj001"])
        assert result is True
        assert mock_snv.call_count == 1
        assert mock_dna.call_count == 1
        assert mock_rna.call_count == 1


class TestPipelineBuildDb:
    @patch("local_frequency.pipeline.build_frequency_db")
    def test_build_db(self, mock_build, pipeline_config):
        p = Pipeline(pipeline_config)
        p.source_dict = {"proj001": {".filter.xls": ["/x"]}}
        result = p._step_build_db(project_ids=["proj001"])
        assert result is True
        mock_build.assert_called_once()


class TestPipelineCnv:
    @patch("local_frequency.pipeline.stat_cnv")
    def test_cnv_step(self, mock_cnv, pipeline_config):
        mock_cnv.return_value = {"total": 1, "output_file": "/out.gz"}

        p = Pipeline(pipeline_config)
        p.source_dict = {"proj001": {".cnv.vcf": ["/x"]}}
        result = p._step_cnv(project_ids=["proj001"])
        assert result is True
        mock_cnv.assert_called_once()
