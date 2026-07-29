# Local Frequency

从远程服务器收集突变数据，统计频率并构建本地 annovar 频率库。

## 简介

Local Frequency 是一个用于分析本地突变频率的工具。它通过 SSH 连接远程服务器，
收集 SNV、DNA fusion、RNA fusion 和 CNV 数据，进行统计分析，
最终生成可用于 annovar 注释的频率数据库。

### 功能特性

- **SSH 远程文件遍历** — 支持密钥/密码认证，自动重连
- **文件分类与拷贝** — 按后缀自动分类（SNV/Fusion/CNV），增量拷贝
- **并行统计** — SNV、DNA fusion、RNA fusion 使用线程池并行统计
- **CNV 统计** — VCF 解析 + bgzip 压缩 + tabix 索引
- **频率库构建** — 将统计结果转换为 annovar 数据库格式
- **CLI 入口** — 支持 dry-run、项目过滤、仅收集等模式

## 安装

### 方式一：pip 安装

```bash
pip install -e .
```

开发依赖：
```bash
pip install -e ".[dev]"
```

### 方式二：使用 uv（推荐）

```bash
# 安装 uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# 创建虚拟环境并安装
uv venv .venv
source .venv/bin/activate
uv pip install -e ".[dev]"
```

## 配置

复制并编辑配置文件：

```bash
cp config/config.yaml.example config/config.yaml
```

配置项说明：

```yaml
ssh:
  host: "your-server.example.com"   # SSH 服务器地址
  port: 22                          # SSH 端口
  user: "your_username"             # 用户名
  key_path: "~/.ssh/id_rsa"         # 密钥路径（可选）
  # password: "your_password"       # 密码认证（二选一）

remote:
  base_path: "/path/to/remote/data" # 远程数据根目录

local:
  base_path: "/path/to/local/data"  # 本地存储目录

stats:
  max_workers: 4                    # 统计并行线程数

db_builder:
  processes: 2                      # 频率库构建进程数
```

## 使用

```bash
# 完整运行
python -m local_frequency --config config/config.yaml

# 仅收集文件（不统计）
python -m local_frequency --config config/config.yaml --collect-only

# 指定项目
python -m local_frequency --config config/config.yaml --project 项目编号1 项目编号2

# dry-run 模式（预检查配置）
python -m local_frequency --config config/config.yaml --dry-run
```

### 输出目录结构

```
本地存储目录/
├── 项目编号/
│   ├── data/                           # 拷贝的原始数据
│   │   ├── sample.filter.xls
│   │   ├── sample.total.fusion.xls
│   │   ├── sample.tsv.redup.xls
│   │   └── sample.cnv.vcf
│   ├── mutation_frequency_result/      # 统计结果
│   │   ├── mutation_frequency.xls      # SNV 频率
│   │   └── stat.json                   # 样本数统计
│   └── frequency_db/                   # annovar 频率库
│       ├── local_freq.84.txt
│       ├── local_freq.84.txt.idx
│       └── ...
```

## 开发

### 代码检查

```bash
ruff check src/ tests/
```

### 运行测试

```bash
pytest                    # 运行所有测试
pytest -v                 # 详细输出
pytest --cov=local_frequency  # 测试覆盖率
```

### 目录结构

```
local_frequency/
├── src/local_frequency/    # 源代码
│   ├── __init__.py         # 包入口
│   ├── config.py           # 配置加载
│   ├── ssh_client.py       # SSH 连接管理
│   ├── file_collector.py   # 文件收集与分类
│   ├── file_copier.py      # SFTP 文件拷贝
│   ├── stats/              # 统计模块
│   │   ├── snv_stat.py     # SNV 统计
│   │   ├── dnafusion_stat.py  # DNA Fusion 统计
│   │   ├── rnafusion_stat.py  # RNA Fusion 统计
│   │   └── cnv_stat.py     # CNV 统计
│   ├── db_builder.py       # 频率库构建
│   └── pipeline.py         # 流程编排 + CLI
├── config/                 # 配置文件
├── tests/                  # 测试
├── scripts/                # 参考脚本（不纳入包）
└── pyproject.toml
```

## License

MIT
