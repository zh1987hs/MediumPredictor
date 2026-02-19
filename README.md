# Genome Nutrition Predictor

一个面向**在线环境**（可联网下载依赖与数据库）的细菌基因组营养需求与底物偏好预测工具（Python 3.10+）。

## 功能
- 基于 KO/EC 规则库评估通路完整性（Present/Partial/Absent + 0-100 分）
- 预测可能外源需求（维生素/氨基酸）
- 预测碳源与氮源可利用能力与偏好分数（0-1）
- 输出结构化 JSON、可读报告、证据表、培养基建议
- 支持自定义 YAML 规则追加覆盖

## 安装（在线）
```bash
python -m pip install -U pip
pip install -e .
```

## 快速开始
```bash
genome_nutrition_predictor run \
  --annotations examples/toy_ko.tsv \
  --annotation_mode provided_ko \
  --out examples/output_toy
```

## CLI
- `genome_nutrition_predictor run --genome xxx.fna --out outdir`
- `genome_nutrition_predictor run --proteins xxx.faa --annotations ann.tsv --out outdir`
- `genome_nutrition_predictor build-rules-template --out rules_template.yaml`
- `genome_nutrition_predictor explain --factor biotin_B7`
- `genome_nutrition_predictor batch --input_list examples/batch_input.tsv --out outdir`

## 输入格式
- `provided_ko`: TSV 至少包含 `gene_id`、`ko`
- `eggnog_tsv`/`dram`/`prokka_gff`/`bakta`: 通用 TSV 解析（自动识别 gene/KO/EC 列）
- `--genome`: 若无蛋白文件，可调用 prodigal 预测 CDS

## 输出文件
- `metadata.json`: 版本、参数、时间戳、数据库路径
- `summary.json`: 聚合结果
- `report.md`: 报告
- `pathways.tsv`: 通路评分明细
- `evidence.tsv`: 基因证据
- `medium_suggestion.tsv`: 建议培养基

## 自定义规则
通过 `--rules custom.yaml` 覆盖默认规则。顶层键：
- `pathways`
- `auxotrophy_factors`
- `substrate_rules`
- `metal_rules`
- `preference_model`

## 说明
- `kofam` 模式预留接口，当前推荐先使用 `provided_ko` 或来自 eggNOG/DRAM 的 KO 注释结果。
- 规则来源为通用代谢框架经验规则，适合可解释初筛，不替代实验验证。


## 安装故障排查（Windows 常见）
如果你看到：`does not appear to be a Python project: neither setup.py nor pyproject.toml found`，通常是以下原因之一：

1. 当前目录不是项目根目录（请先确认同级能看到 `pyproject.toml`）。
2. 你在 GitHub 上下载/切换到了不完整分支（缺少仓库根文件）。
3. 本地目录名正确但文件未同步（例如只下载了子目录）。

建议执行：
```powershell
# 1) 确认在项目根目录
Get-ChildItem

# 2) 确认关键文件存在
Test-Path .\pyproject.toml
Test-Path .\setup.py

# 3) 然后再安装
python -m pip install -U pip
pip install -e .
```
