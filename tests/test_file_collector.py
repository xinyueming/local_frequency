"""FileCollector 模块单元测试"""


from local_frequency.file_collector import FileCollector


class TestFileCollector:
    def test_classify_snv_files(self):
        collector = FileCollector(base_path="/data")
        files = [
            "/data/a/b/c/d/proj001/sample.filter.xls",
            "/data/a/b/c/d/proj001/sample.filter.germline.xls",
        ]
        result = collector.classify_files(files)

        assert "proj001" in result
        assert ".filter.xls" in result["proj001"]
        assert ".filter.germline.xls" in result["proj001"]

    def test_classify_fusion_files(self):
        collector = FileCollector(base_path="/data")
        files = [
            "/data/a/b/c/d/proj002/sample.total.fusion.xls",
            "/data/a/b/c/d/proj002/sample.tsv.redup.xls",
        ]
        result = collector.classify_files(files)

        assert "proj002" in result
        assert ".total.fusion.xls" in result["proj002"]
        assert ".tsv.redup.xls" in result["proj002"]

    def test_classify_cnv_files(self):
        collector = FileCollector(base_path="/data")
        files = ["/data/a/b/c/d/proj003/sample.cnv.vcf"]
        result = collector.classify_files(files)

        assert "proj003" in result
        assert ".cnv.vcf" in result["proj003"]

    def test_classify_multiple_projects(self):
        collector = FileCollector(base_path="/data")
        files = [
            "/data/a/b/c/d/proj001/sample.filter.xls",
            "/data/a/b/c/d/proj002/sample.filter.xls",
        ]
        result = collector.classify_files(files)

        assert "proj001" in result
        assert "proj002" in result
        assert len(result["proj001"][".filter.xls"]) == 1

    def test_skip_short_paths(self):
        collector = FileCollector(base_path="/data")
        # 相对路径: a/b （不足 project_level=5 层）
        files = ["/data/a/b/sample.filter.xls"]
        result = collector.classify_files(files)
        assert result == {}

    def test_unknown_suffix_ignored(self):
        collector = FileCollector(base_path="/data")
        files = ["/data/a/b/c/d/proj001/readme.txt"]
        result = collector.classify_files(files)
        assert result == {}

    def test_custom_suffixes(self):
        collector = FileCollector(base_path="/data", suffixes=[".bam", ".bai"])
        files = [
            "/data/a/b/c/d/proj001/sample.bam",
            "/data/a/b/c/d/proj001/sample.bai",
        ]
        result = collector.classify_files(files)
        assert "proj001" in result
        assert ".bam" in result["proj001"]

    def test_match_suffix_precedence(self):
        """更具体的后缀优先匹配"""
        collector = FileCollector(base_path="/data")
        files = [
            "/data/a/b/c/d/proj001/sample.filter.germline.xls",
        ]
        result = collector.classify_files(files)
        # 应该匹配 .filter.germline.xls 而不是 .filter.xls
        assert ".filter.germline.xls" in result["proj001"]
        assert ".filter.xls" not in result["proj001"]

    def test_custom_project_level(self):
        collector = FileCollector(base_path="/data", project_level=2)
        # 相对路径: b/c/file.xls → project_level=2 取 'c'
        files = ["/data/a/b/c/sample.filter.xls"]
        result = collector.classify_files(files)

        assert "c" in result
