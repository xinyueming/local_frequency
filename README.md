# Local Frequency
从远程服务器收集突变数据，统计频率并构建本地annovar频率库。

## 简介
Local Frequency 是一个用于分析本地突变频率的工具。它通过 SSH 连接远程服务器，
收集 SNV、DNA fusion、RNA fusion 和 CNV 数据，进行统计分析，
最终生成可用于 annovar 注释的频率数据库。

## 安装
```bash
pip install -e .
```

开发依赖：
```bash
pip install -e ".[dev]"
```

## 配置
复制并编辑配置文件：
```bash
cp config/config.yaml.example config/config.yaml
```

## 使用
```bash
# 完整运行
python -m local_frequency --config config/config.yaml

# 仅收集文件
python -m local_frequency --config config/config.yaml --collect-only

# 指定项目
python -m local_frequency --config config/config.yaml --project 项目编号

# dry-run 模式
python -m local_frequency --config config/config.yaml --dry-run
```

## 目录结构
```
local_frequency/
├── src/local_frequency/    # 源代码
│   ├── config.py           # 配置加载
│   ├── ssh_client.py       # SSH 连接管理
│   ├── file_collector.py   # 文件收集与分类
│   ├── file_copier.py      # SFTP 文件拷贝
│   ├── stats/              # 统计模块
│   ├── db_builder.py       # 频率库构建
│   └── pipeline.py         # 流程编排
├── config/                 # 配置文件
├── tests/                  # 测试
├── scripts/                # 参考脚本（不纳入包）
└── pyproject.toml
```

## 开发
```bash
# 代码检查
ruff check src/ tests/

# 运行测试
pytest

# 测试覆盖率
pytest --cov=local_frequency --cov-report=term-missing
```

## License
MIT
