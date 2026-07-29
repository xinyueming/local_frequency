"""SSH 连接与远程文件管理模块"""

from __future__ import annotations

import logging
import os
from collections import defaultdict

import paramiko

logger = logging.getLogger(__name__)


class SSHClient:
    """SSH 连接管理器，支持密钥/密码认证和自动重连"""

    def __init__(
        self,
        host: str,
        port: int = 22,
        user: str = "root",
        key_path: str | None = None,
        password: str | None = None,
    ):
        self.host = host
        self.port = port
        self.user = user
        self.key_path = key_path
        self.password = password
        self._client: paramiko.SSHClient | None = None
        self._sftp = None

    def connect(self) -> bool:
        """建立 SSH 连接

        Returns:
            连接是否成功
        """
        self._client = paramiko.SSHClient()
        self._client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        connect_kwargs = {
            "hostname": self.host,
            "port": self.port,
            "username": self.user,
            "timeout": 30,
            "allow_agent": True,
            "look_for_keys": True,
        }

        if self.key_path:
            key_path = os.path.expanduser(self.key_path)
            connect_kwargs["key_filename"] = key_path
            logger.info("使用密钥认证: %s", key_path)
        elif self.password:
            connect_kwargs["password"] = self.password
            logger.info("使用密码认证")
        else:
            logger.info("使用默认 SSH 认证（agent/密钥）")

        try:
            self._client.connect(**connect_kwargs)
            self._sftp = self._client.open_sftp()
            logger.info("SSH 连接成功: %s@%s:%d", self.user, self.host, self.port)
            return True
        except Exception as e:  # noqa: BLE001 — 连接层需捕获所有异常
            logger.error("SSH 连接失败: %s", e)
            self._client = None
            return False

    def _ensure_connected(self) -> bool:
        """确保连接存活，断线自动重连"""
        if self._client is None:
            return self.connect()

        try:
            transport = self._client.get_transport()
            if transport is None or not transport.is_active():
                logger.warning("连接已断开，尝试重连")
                self._client = None
                return self.connect()
            return True
        except Exception:  # noqa: BLE001 — 连接层需捕获所有异常
            self._client = None
            return self.connect()

    def walk_remote_dir(self, remote_path: str) -> list[str]:
        """递归遍历远程目录

        Args:
            remote_path: 远程目录路径

        Returns:
            所有文件的绝对路径列表
        """
        if not self._ensure_connected():
            raise ConnectionError("SSH 连接不可用")

        files = []
        self._walk_recursive(remote_path, files)
        logger.info("遍历完成: %s → %d 个文件", remote_path, len(files))
        return files

    def _walk_recursive(self, path: str, files: list[str]) -> None:
        """递归遍历辅助"""
        try:
            items = self._sftp.listdir_attr(path)
        except OSError as e:
            logger.warning("无法读取目录 %s: %s", path, e)
            return

        for item in items:
            full_path = os.path.join(path, item.filename)
            if self._is_dir(item):
                self._walk_recursive(full_path, files)
            else:
                files.append(full_path)

    def _is_dir(self, attr) -> bool:
        """判断是否为目录"""
        import stat

        return stat.S_ISDIR(attr.st_mode)

    def get_files_by_suffix(
        self,
        remote_path: str,
        suffixes: list[str],
    ) -> dict[str, dict[str, list[str]]]:
        """按后缀过滤文件，并按目录层级提取项目编号

        Args:
            remote_path: 远程根目录
            suffixes: 文件后缀列表，如 ['.filter.xls', '.total.fusion.xls']

        Returns:
            {项目编号: {后缀: [文件路径列表]}}
        """
        all_files = self.walk_remote_dir(remote_path)
        result: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))

        suffixes_set = set(suffixes)
        # 规范化 remote_path 确保结尾无斜杠
        norm_base = remote_path.rstrip("/")

        for fpath in all_files:
            rel_path = os.path.relpath(fpath, norm_base)
            parts = rel_path.split(os.sep)

            # 第 5 层目录（索引 4）为项目编号:
            # a/b/c/d/proj123/file.filter.xls → parts[4] = "proj123"
            if len(parts) < 6:
                continue
            project_id = parts[4]

            # 匹配后缀
            fname = os.path.basename(fpath)
            for suffix in suffixes_set:
                if fname.endswith(suffix):
                    result[project_id][suffix].append(fpath)
                    break

        return dict(result)

    def listdir(self, remote_path: str) -> list[str]:
        """列出远程目录下的文件名

        Args:
            remote_path: 远程目录路径

        Returns:
            文件名列表
        """
        if not self._ensure_connected():
            raise ConnectionError("SSH 连接不可用")
        return self._sftp.listdir(remote_path)

    def stat(self, remote_path: str):
        """获取远程文件/目录属性"""
        if not self._ensure_connected():
            raise ConnectionError("SSH 连接不可用")
        return self._sftp.stat(remote_path)

    @property
    def sftp(self):
        """获取 SFTP 会话"""
        if not self._ensure_connected():
            raise ConnectionError("SSH 连接不可用")
        return self._sftp

    def close(self) -> None:
        """关闭连接"""
        if self._sftp:
            self._sftp.close()
            self._sftp = None
        if self._client:
            self._client.close()
            self._client = None
        logger.info("SSH 连接已关闭")

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
