"""Small Typer compatibility shim for offline use."""
from __future__ import annotations

from typing import Any, Callable


class Typer:
    def __init__(self, help: str = ""):
        self.help = help

    def command(self, name: str | None = None):
        def deco(func: Callable[..., Any]) -> Callable[..., Any]:
            return func
        return deco

    def __call__(self):
        raise SystemExit("CLI shim active; install typer for full command-line parsing.")


def Option(default: Any = None, **kwargs: Any) -> Any:
    _ = kwargs
    return default


def echo(msg: str) -> None:
    print(msg)
