"""Command line entry point: ``ovary-fm <stage>``."""

from __future__ import annotations

import argparse
import sys
import time

from . import __version__
from .config import load_config, resolve_path

STAGES = ("inspect", "fetch", "qc", "harmonise", "bench", "predict", "score", "t1", "fm", "fm-compare",
          "train", "report", "all")
# `all` on the plumbing: qc -> harmonise -> bench -> the three baseline predictions scored -> t1 -> report
CHAIN = ("qc", "harmonise", "bench", "predict", "score", "t1", "report")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="ovary-fm",
        description="An ovary-specific foundation model: corpus, held-out-context benchmark, model arms.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="stages:\n"
               "  inspect     the manifest: what is on disk, what is QC'd, what is still fill_in\n"
               "  fetch       download each enabled dataset (--dataset NAME for one)\n"
               "  qc          per dataset: modality, filters, MAD bounds, the obs schema -> data/corpus/qc/\n"
               "  harmonise   the crosswalk, the gene space, corpus.h5ad and the held-out shards\n"
               "  bench       build the T3 reference bundles: b0, b1, r per (context, condition, compartment)\n"
               "              and the ceilings table\n"
               "  predict     a prediction object (--model b0|b1|stub) for every bundle\n"
               "  score       score a prediction (--pred PATH --name NAME [--scorer lite|cell-eval2])\n"
               "  t1          leave-one-atlas-out label transfer (--rep an embedding.parquet to score instead)\n"
               "  fm          one generalist's zero-shot embedding of the corpus (--tool stub|scprint2|uce|transcriptformer)\n"
               "  fm-compare  t1 on every results/fm/*/embedding.parquet\n"
               "  train       --arm a|b writes the training recipe; --arm c runs scVI/scANVI (the baseline)\n"
               "  report      results/report/first_results.md from every table above\n"
               "  all         qc harmonise bench predict score t1 report, with the stub model\n")
    parser.add_argument("stage", choices=STAGES)
    parser.add_argument("--params", default=None, help="alternative config/params.yaml")
    parser.add_argument("--set", dest="overrides", action="append", default=[],
                        metavar="key.path=value", help="override one config value; repeatable")
    parser.add_argument("--data", default=None, help="data root (downloads and the corpus)")
    parser.add_argument("--results", default=None, help="output root")
    parser.add_argument("--dataset", default=None, help="fetch / qc: one dataset of the manifest")
    parser.add_argument("--model", default="stub", help="predict: b0 | b1 | stub")
    parser.add_argument("--pred", default=None, help="score: the prediction .h5ad")
    parser.add_argument("--name", default=None, help="score / t1: a name for the output tables")
    parser.add_argument("--scorer", default="lite", choices=["lite", "cell-eval2"])
    parser.add_argument("--rep", default=None, help="t1: an embedding.parquet to transfer labels inside")
    parser.add_argument("--tool", default=None, help="fm: which model")
    parser.add_argument("--arm", default=None, choices=["a", "b", "c"], help="train: which arm")
    parser.add_argument("--adata", default=None, help="fm / train: an object instead of the corpus")
    parser.add_argument("--version", action="version", version=f"ovary-fm {__version__}")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    cfg = load_config(args.params, args.overrides)
    if args.data:
        cfg["paths"]["data"] = str(resolve_path(args.data, prefer="cwd"))
    if args.results:
        cfg["paths"]["results"] = str(resolve_path(args.results, prefer="cwd"))

    print(f"ovary-fm {__version__} :: stage {args.stage}")
    print(f"  data      = {cfg['paths']['data']}")
    print(f"  results   = {cfg['paths']['results']}")
    print(f"  manifest  = {cfg['paths']['corpus_manifest']}")
    print(f"  benchmark = {cfg['paths']['benchmark']}")
    print()

    if args.stage == "inspect":
        from .io import inspect
        inspect(cfg)
        return 0

    from . import fetch, fm, harmonise, qc, report
    from .bench import runner, t1
    from .models import recipes, stub

    def _predict():
        if args.stage == "all":
            for model in ("b0", "b1", "stub"):
                stub.predict(cfg, model)
        else:
            stub.predict(cfg, args.model)

    def _score():
        if args.stage == "all":
            for model in ("b0", "b1", "stub"):
                runner.score(cfg, stub.prediction_path(cfg, model), model, scorer="lite")
        else:
            if not args.pred:
                raise SystemExit("score needs --pred PATH (a prediction .h5ad from `predict`)")
            runner.score(cfg, resolve_path(args.pred, prefer="cwd"), args.name or "prediction", scorer=args.scorer)

    def _train():
        if args.arm is None:
            raise SystemExit("train needs --arm a|b|c")
        if args.arm == "c":
            from .models import scvi_arm
            scvi_arm.train(cfg, adata_path=args.adata)
        else:
            recipes.write(cfg, args.arm)

    runners = {
        "fetch": lambda: fetch.run(cfg, args.dataset),
        "qc": lambda: qc.run(cfg, args.dataset),
        "harmonise": lambda: harmonise.run(cfg),
        "bench": lambda: runner.build(cfg),
        "predict": _predict,
        "score": _score,
        "t1": lambda: t1.run(cfg, rep=args.rep, name=args.name),
        "fm": lambda: fm.embed(cfg, args.tool, adata_path=args.adata),
        "fm-compare": lambda: fm.compare(cfg),
        "train": _train,
        "report": lambda: report.run(cfg),
    }
    order = list(CHAIN) if args.stage == "all" else [args.stage]
    failures = []
    for name in order:
        start = time.time()
        print(f"== {name} " + "=" * (60 - len(name)))
        try:
            runners[name]()
        except Exception as exc:
            if len(order) == 1:
                raise
            failures.append((name, exc))
            print(f"!! {name} failed: {exc}")
            break
        print(f"   {name} done in {time.time() - start:.0f} s\n")
    if failures:
        print(f"stopped at {failures[0][0]}: {failures[0][1]}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
