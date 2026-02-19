from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_rules(rules_path: str | Path | None = None) -> dict[str, Any]:
    """Load built-in and optional user rules, merging by top-level keys."""
    built_in = yaml.safe_load(Path(__file__).resolve().parent.parent.parent.joinpath("rules", "default_rules.yaml").read_text())
    if rules_path is None:
        return built_in
    custom = yaml.safe_load(Path(rules_path).read_text())
    return merge_rules(built_in, custom)


def merge_rules(base: dict[str, Any], custom: dict[str, Any]) -> dict[str, Any]:
    """Deep merge nested dictionaries/lists for extensible rule customization."""
    merged = dict(base)
    for key, value in custom.items():
        if key not in merged:
            merged[key] = value
            continue
        if isinstance(value, dict) and isinstance(merged[key], dict):
            merged[key] = merge_rules(merged[key], value)
        elif isinstance(value, list) and isinstance(merged[key], list):
            existing_ids = {x.get("id") for x in merged[key] if isinstance(x, dict)}
            combined = list(merged[key])
            for item in value:
                if isinstance(item, dict) and item.get("id") in existing_ids:
                    combined = [item if (isinstance(x, dict) and x.get("id") == item.get("id")) else x for x in combined]
                else:
                    combined.append(item)
            merged[key] = combined
        else:
            merged[key] = value
    return merged


def build_template() -> str:
    """Return a minimal template for user-defined rules."""
    return """pathways:\n  - id: custom_example\n    name: Custom pathway\n    category: custom\n    essential_weight: 0.8\n    optional_weight: 0.2\n    essential_steps:\n      - id: key_step\n        any_of: [K00001, K00002]\n    optional_steps:\n      - id: aux_step\n        any_of: [K01000]\n"""
