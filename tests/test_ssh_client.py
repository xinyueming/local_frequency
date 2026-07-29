"""SSH 模块单元测试"""

import os
from unittest.mock import MagicMock, patch

import pytest

from local_frequency.ssh_client import SSHClient


@pytest.fixture
def ssh_config():
    return {
        "host": "test-server.example.com",
        "port": 22,
        "user": "testuser",
        "key_path": "~/.ssh/id_ed25519",
    }


@pytest.fixture
def mock_sftp():
    sftp = MagicMock()
    sftp.listdir.return_value = ["file1.txt", "file2.filter.xls", "subdir"]
    return sftp


@pytest.fixture
def mock_client(mock_sftp):
    client = MagicMock()
    client.get_transport().is_active.return_value = True
    client.open_sftp.return_value = mock_sftp
    return client


class TestSSHClientConnect:
    """SSH 连接测试"""

    @patch("paramiko.SSHClient")
    def test_connect_with_key(self, MockSSHClient, ssh_config):
        client = SSHClient(**ssh_config)
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_instance.open_sftp.return_value = MagicMock()
        MockSSHClient.return_value = mock_instance

        assert client.connect() is True
        mock_instance.connect.assert_called_once()
        call_kwargs = mock_instance.connect.call_args[1]
        assert call_kwargs["key_filename"] == os.path.expanduser(ssh_config["key_path"])

    @patch("paramiko.SSHClient")
    def test_connect_with_password(self, MockSSHClient):
        client = SSHClient(
            host="test.example.com", user="testuser", password="secret123"
        )
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_instance.open_sftp.return_value = MagicMock()
        MockSSHClient.return_value = mock_instance

        assert client.connect() is True
        call_kwargs = mock_instance.connect.call_args[1]
        assert call_kwargs["password"] == "secret123"
        assert "key_filename" not in call_kwargs

    @patch("paramiko.SSHClient")
    def test_connect_failure(self, MockSSHClient):
        MockSSHClient.return_value.connect.side_effect = Exception("Connection refused")

        client = SSHClient(host="bad-host", user="testuser")
        assert client.connect() is False

    @patch("paramiko.SSHClient")
    def test_context_manager(self, MockSSHClient):
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_instance.open_sftp.return_value = MagicMock()
        MockSSHClient.return_value = mock_instance

        with SSHClient(host="test.example.com", user="testuser") as client:
            assert client._client is not None
        mock_instance.close.assert_called()


class TestWalkRemoteDir:
    """远程目录遍历测试"""

    def _make_attr(self, name, is_dir=True):
        attr = MagicMock()
        attr.filename = name
        import stat

        if is_dir:
            attr.st_mode = stat.S_IFDIR
        else:
            attr.st_mode = stat.S_IFREG
        return attr

    @patch("paramiko.SSHClient")
    def test_walk_flat_directory(self, MockSSHClient):
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_sftp = MagicMock()
        mock_instance.open_sftp.return_value = mock_sftp
        MockSSHClient.return_value = mock_instance

        # /data 下有 3 个文件
        mock_sftp.listdir_attr.return_value = [
            self._make_attr("a.xls", False),
            self._make_attr("b.filter.xls", False),
            self._make_attr("c.txt", False),
        ]

        client = SSHClient(host="test", user="test")
        client.connect()
        files = client.walk_remote_dir("/data")

        assert len(files) == 3
        assert "/data/a.xls" in files


class TestGetFilesBySuffix:
    """按后缀过滤测试"""

    def _make_attr(self, name, is_dir=True):
        attr = MagicMock()
        attr.filename = name
        import stat

        if is_dir:
            attr.st_mode = stat.S_IFDIR
        else:
            attr.st_mode = stat.S_IFREG
        return attr

    @patch("paramiko.SSHClient")
    def test_filter_by_suffix(self, MockSSHClient):
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_sftp = MagicMock()
        mock_instance.open_sftp.return_value = mock_sftp
        MockSSHClient.return_value = mock_instance

        # remote_path = /base
        # 第 6 层: /base/a/b/c/d/proj123/file.filter.xls
        files_flat = [self._make_attr("file.filter.xls", False)]

        def walk_side_effect(path):
            if path == "/base":
                return [self._make_attr("a", True)]
            elif path == "/base/a":
                return [self._make_attr("b", True)]
            elif path == "/base/a/b":
                return [self._make_attr("c", True)]
            elif path == "/base/a/b/c":
                return [self._make_attr("d", True)]
            elif path == "/base/a/b/c/d":
                return [self._make_attr("proj123", True)]
            elif path == "/base/a/b/c/d/proj123":
                return files_flat
            return []

        mock_sftp.listdir_attr.side_effect = walk_side_effect

        client = SSHClient(host="test", user="test")
        client.connect()
        result = client.get_files_by_suffix("/base", [".filter.xls"])

        assert "proj123" in result
        assert ".filter.xls" in result["proj123"]
        assert len(result["proj123"][".filter.xls"]) == 1


class TestEnsureConnected:
    """连接检查测试"""

    @patch("paramiko.SSHClient")
    def test_ensure_connected_active(self, MockSSHClient):
        """活跃连接应直接返回"""
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_instance.open_sftp.return_value = MagicMock()
        MockSSHClient.return_value = mock_instance

        client = SSHClient(host="test", user="test")
        client.connect()
        assert client._ensure_connected() is True
        # Should not have called connect again
        assert mock_instance.connect.call_count == 1


class TestClose:
    """关闭连接测试"""

    @patch("paramiko.SSHClient")
    def test_close_sftp_and_client(self, MockSSHClient):
        mock_instance = MagicMock()
        mock_instance.get_transport.return_value = MagicMock()
        mock_instance.get_transport().is_active.return_value = True
        mock_sftp = MagicMock()
        mock_instance.open_sftp.return_value = mock_sftp
        MockSSHClient.return_value = mock_instance

        client = SSHClient(host="test", user="test")
        client.connect()
        client.close()

        mock_sftp.close.assert_called_once()
        mock_instance.close.assert_called_once()
