from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Record:
    id: str
    seq: str


def parse(path: str, fmt: str):
    if fmt != "fasta":
        raise ValueError("Only fasta supported")
    current_id = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    yield Record(id=current_id, seq="".join(seq))
                current_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if current_id is not None:
        yield Record(id=current_id, seq="".join(seq))
