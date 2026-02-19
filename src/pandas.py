"""Very small subset of pandas API used in this project."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any


class Series(list):
    def value_counts(self) -> "ValueCounts":
        return ValueCounts(Counter(self))


class ValueCounts(dict):
    def to_dict(self) -> dict[str, int]:
        return dict(self)


class DataFrame:
    def __init__(self, rows: list[dict[str, Any]]):
        self.rows = rows
        self.columns = list(rows[0].keys()) if rows else []
        self.empty = len(rows) == 0

    def fillna(self, value: str) -> "DataFrame":
        out = []
        for row in self.rows:
            out.append({k: (value if v is None else v) for k, v in row.items()})
        return DataFrame(out)

    def iterrows(self):
        for i, row in enumerate(self.rows):
            yield i, row

    def to_csv(self, path: str | Path, sep: str = "\t", index: bool = False) -> None:
        _ = index
        with Path(path).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.columns, delimiter=sep)
            writer.writeheader()
            for row in self.rows:
                writer.writerow(row)

    def __getitem__(self, key: str) -> Series:
        return Series([row.get(key) for row in self.rows])


def read_csv(path: str | Path, sep: str = "\t", comment: str | None = None, dtype=str) -> DataFrame:
    _ = dtype
    rows: list[dict[str, Any]] = []
    with Path(path).open() as handle:
        lines = [ln for ln in handle.read().splitlines() if not (comment and ln.startswith(comment))]
    if not lines:
        return DataFrame([])
    reader = csv.DictReader(lines, delimiter=sep)
    for row in reader:
        rows.append(dict(row))
    return DataFrame(rows)
