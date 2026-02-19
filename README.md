# Genome Nutrition Predictor

一个离线优先的细菌基因组营养需求与底物偏好预测工具（Python 3.10+）。

## 功能
- 基于 KO/EC 规则库评估通路完整性（Present/Partial/Absent + 0-100 分）
- 预测可能外源需求（维生素/氨基酸）
- 预测碳源与氮源可利用能力与偏好分数（0-1）
- 输出结构化 JSON、可读报告、证据表、培养基建议
- 支持自定义 YAML 规则追加覆盖

## 安装
```bash
pip install -e .
```

## 快速开始
使用示例 KO 表：
```bash
genome_nutrition_predictor run \
  --annotations examples/toy_ko.tsv \
  --annotation_mode provided_ko \
  --out examples/output_toy
```

解释某个因子规则：
```bash
genome_nutrition_predictor explain --factor biotin_B7
```

生成自定义规则模板：
```bash
genome_nutrition_predictor build-rules-template --out my_rules.yaml
```

批处理：
```bash
genome_nutrition_predictor batch --input_list examples/batch_input.tsv --out examples/batch_out
```

## 输入格式
- `provided_ko`（推荐起步）: TSV 至少包含 `gene_id`、`ko`
- `eggnog_tsv`/`dram`/`prokka_gff`/`bakta`: 当前为通用 TSV 解析模式（自动寻找 gene/KO 列）
- `--genome`：若无蛋白文件，可调用 prodigal 预测 CDS

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
- KOfamScan/HMM 自动注释预留了接口，建议当前先用 `provided_ko` 快速验证。
- 规则来源是通用代谢框架经验规则，适合可解释初筛，不替代实验验证。
