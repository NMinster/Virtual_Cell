"""Configuration loading, dotted overrides, output paths, provenance.

The pattern is Mouse_Compensatory's: one ``config/params.yaml``, side files
for the things that grow (the corpus manifest, the crosswalk, the
benchmark), ``--set a.b.c=value`` overrides that refuse unknown keys, the
environment winning for the roots that move between machines, and a
provenance json written by every stage before it writes anything else.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PARAMS = REPO_ROOT / "config" / "params.yaml"
ENV_PREFIX = "OVARY_FM"


def _read_yaml(path) -> dict:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"config file not found: {path}")
    try:
        return yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ValueError(f"could not parse {path}\n\n{exc}") from None


def _coerce(text: str):
    """Turn a CLI --set value into a Python scalar, JSON first."""
    try:
        return json.loads(text)
    except (ValueError, TypeError):
        pass
    lowered = text.strip().lower()
    if lowered in {"none", "null", "~"}:
        return None
    if lowered in {"true", "yes"}:
        return True
    if lowered in {"false", "no"}:
        return False
    return text


def apply_override(cfg: dict, assignment: str) -> dict:
    """Apply one ``a.b.c=value`` override in place; an unknown key is an error."""
    if "=" not in assignment:
        raise ValueError(f"--set expects key=value, got {assignment!r}")
    key, _, raw = assignment.partition("=")
    parts = key.strip().split(".")
    node = cfg
    for part in parts[:-1]:
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"--set {key}: no such config section {part!r}")
        node = node[part]
    if not isinstance(node, dict) or parts[-1] not in node:
        raise KeyError(f"--set {key}: no such config key {parts[-1]!r}")
    node[parts[-1]] = _coerce(raw)
    return cfg


def resolve_path(path, prefer: str = "repo") -> Path:
    """Anchor a relative path to the repo root first (or the cwd, when asked)."""
    text = str(path)
    if Path(text).is_absolute():
        return Path(text)
    roots = [REPO_ROOT, Path.cwd()]
    if prefer == "cwd":
        roots.reverse()
    for root in roots:
        if (root / text).exists():
            return root / text
    return roots[0] / text


def load_config(params=None, overrides=()) -> dict:
    """params.yaml, then the overrides, then the environment, then the side files."""
    cfg = _read_yaml(params or DEFAULT_PARAMS)
    for assignment in overrides:
        apply_override(cfg, assignment)
    for key, env in (("data", f"{ENV_PREFIX}_DATA"), ("results", f"{ENV_PREFIX}_RESULTS")):
        value = os.environ.get(env)
        if value:
            cfg["paths"][key] = value
    for key in ("data", "results"):
        cfg["paths"][key] = str(resolve_path(cfg["paths"][key]))
    for key, side in (("corpus_manifest", "corpus"), ("crosswalk", "crosswalk"), ("benchmark", "benchmark")):
        cfg["paths"][key] = str(resolve_path(cfg["paths"][key], prefer="cwd"))
        cfg[side] = _read_yaml(cfg["paths"][key])
    return cfg


def data_dir(cfg, *parts) -> Path:
    out = Path(cfg["paths"]["data"]).joinpath(*parts)
    out.mkdir(parents=True, exist_ok=True)
    return out


def results_dir(cfg, *parts) -> Path:
    out = Path(cfg["paths"]["results"]).joinpath(*parts)
    out.mkdir(parents=True, exist_ok=True)
    return out


def tables_dir(cfg, *parts) -> Path:
    return results_dir(cfg, "tables", *parts)


def git_revision() -> str | None:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return None


def package_versions(*names) -> dict:
    import importlib

    out = {}
    for name in names:
        try:
            out[name] = getattr(importlib.import_module(name), "__version__", "?")
        except Exception:
            out[name] = None
    return out


def write_provenance(cfg, stage: str, extra: dict | None = None) -> Path:
    """results/provenance/<stage>_<utc>.json: git revision, the live config, versions."""
    from . import __version__

    folder = results_dir(cfg, "provenance")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    record = {
        "stage": stage,
        "utc": stamp,
        "ovary_fm": __version__,
        "git_revision": git_revision(),
        "python": sys.version.split()[0],
        "versions": package_versions("numpy", "scipy", "pandas", "anndata", "scanpy", "sklearn"),
        "config": {k: v for k, v in cfg.items() if not k.startswith("_")},
        "extra": extra or {},
    }
    path = folder / f"{stage}_{stamp}.json"
    path.write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
    return path


_SAFE = re.compile(r"[^A-Za-z0-9_.-]+")


def slug(text: str) -> str:
    return _SAFE.sub("_", str(text)).strip("_")
