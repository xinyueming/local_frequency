#!/usr/bin/env bash
# Local Frequency 安装脚本
set -e

echo "=== Local Frequency 安装 ==="

# 检查 Python
if ! command -v python3 &>/dev/null; then
    echo "错误: 未找到 python3，请先安装 Python 3.10+"
    exit 1
fi

python_version=$(python3 --version 2>&1 | grep -oP '\d+\.\d+')
if [[ "$python_version" < "3.10" ]]; then
    echo "错误: Python 版本需 >= 3.10，当前: $python_version"
    exit 1
fi

echo "✓ Python $(python3 --version)"

# 检查 uv（如果已安装）
if command -v uv &>/dev/null; then
    echo "✓ uv 已安装: $(uv --version)"
    echo "使用 uv 创建环境..."
    uv venv .venv
    source .venv/bin/activate
    uv pip install -e ".[dev]"
else
    echo "提示: 未找到 uv，使用 pip 安装..."
    echo "  推荐安装: curl -LsSf https://astral.sh/uv/install.sh | sh"
    python3 -m venv .venv
    source .venv/bin/activate
    pip install -e ".[dev]"
fi

echo "✓ 依赖安装完成"

# 检查 ruff
if command -v ruff &>/dev/null; then
    echo "✓ ruff: $(ruff check src/ tests/ 2>&1 | tail -1)"
fi

# 运行测试
echo ""
echo "运行测试..."
python -m pytest tests/ -q

echo ""
echo "=== 安装完成 ==="
echo ""
echo "下一步:"
echo "  1. cp config/config.yaml.example config/config.yaml"
echo "  2. 编辑 config/config.yaml 填写 SSH 和目录配置"
echo "  3. python -m local_frequency --config config/config.yaml --dry-run"
