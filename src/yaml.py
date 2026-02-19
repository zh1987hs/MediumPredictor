"""Tiny yaml compatibility layer for offline environments.
Supports JSON-formatted YAML files.
"""
import json
from typing import Any


def safe_load(text: str) -> Any:
    return json.loads(text)
