"""配置加载模块"""

import yaml


def load_config(config_path: str) -> dict:
    """加载 YAML 配置文件

    Args:
        config_path: 配置文件路径

    Returns:
        配置字典

    Raises:
        FileNotFoundError: 配置文件不存在
        yaml.YAMLError: YAML 格式错误
    """
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)
