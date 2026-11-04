"""``ovary-fm train --arm a|b``: write the training recipe, never launch it.

Arm A is scPRINT-2 medium continued-pretrained on the corpus (GPL-3; the
package ships its pretraining documentation and a SLURM recipe). Arm B is
a <= 30 M-parameter model from scratch in the Tahoe-x1 stack (Apache-2.0;
Composer / FSDP). Both recipes are filled from config/train/*.yml templates
and written with an sbatch script and a README under results/train/<arm>/.
The templates carry every value this plan fixes; the keys of the tools'
own configs change between releases, so the README says to transcribe,
not to run blind.
"""

from __future__ import annotations

from pathlib import Path

from ..config import REPO_ROOT, resolve_path, results_dir, write_provenance
from ..io import corpus_path

SBATCH = """#!/bin/bash
#SBATCH --job-name=ovary_fm_{arm}
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:{gpus}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output={out}/slurm_%j.out
set -euo pipefail
cd {repo}
export HF_HOME=${{HF_HOME:-/data/$USER/hf}}
export http_proxy=${{http_proxy:-http://dtn03-e0:3128}} https_proxy=${{https_proxy:-http://dtn03-e0:3128}}
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate {env}
echo ">> $(date) on $(hostname): {arm}"
{command}
echo ">> finished $(date)"
"""

README = """# Training recipe: arm {arm_upper}

Written by `ovary-fm train --arm {arm}` on {stamp}. Nothing has been run.

* corpus: `{corpus}`
* template: `{template}` (copied here as `config.yml` with the paths filled in)
* environment: `{env}` (`bash biowulf/setup_env.sh --train {setup}`)
* submit: `sbatch {out}/train.sbatch` once `config.yml` has been checked
  against the tool's current documentation ({docs})

Before submitting (docs/first_results_checklist.md, step 10):

1. open the tool's own config documentation and transcribe every key of
   `config.yml` into the tool's format; the keys below are this plan's
   values, not a file the tool loads as is;
2. run a few hundred steps on ONE GPU with `--set models.{arm}.epochs=1`
   and a corpus subset before asking for a node;
3. record the checkpoint hash and the wall time in the provenance json.
"""


def write(cfg, arm: str) -> Path:
    from datetime import datetime, timezone

    write_provenance(cfg, "train", extra={"arm": arm})
    spec = cfg["models"][f"arm_{arm}"]
    template = resolve_path(spec["template"])
    if not template.exists():
        raise FileNotFoundError(template)
    corpus = corpus_path(cfg)
    out = results_dir(cfg, "train", f"arm_{arm}")
    text = template.read_text(encoding="utf-8")
    text = text.replace("{corpus}", str(corpus)).replace("{out}", str(out)).replace("{epochs}", str(spec.get("epochs", 1)))
    text = text.replace("{base_checkpoint}", str(spec.get("base_checkpoint", "")))
    (out / "config.yml").write_text(text, encoding="utf-8")
    if arm == "a":
        command = (f"scprint2 fit --config {out}/config.yml   # transcribe config.yml into scPRINT-2's "
                   f"pretrain_medium.yml format first (docs/pretrain.md in their repo)")
        docs, setup, gpus, cpus, mem, time = "cantinilab/scPRINT-2 docs/pretrain.md", "scprint2", 4, 16, "360g", "3-00:00:00"
    else:
        command = (f"composer scripts/train.py {out}/config.yml   # from a clone of tahoebio/tahoe-x1; "
                   f"transcribe config.yml into their configs/ format first")
        docs, setup, gpus, cpus, mem, time = "tahoebio/tahoe-x1 README and configs/", "tahoe", 4, 16, "240g", "2-00:00:00"
    (out / "train.sbatch").write_text(SBATCH.format(arm=f"arm_{arm}", gpus=gpus, cpus=cpus, mem=mem, time=time, out=out,
                                                    repo=REPO_ROOT, env=spec["env"], command=command), encoding="utf-8")
    (out / "README.md").write_text(README.format(arm=arm, arm_upper=arm.upper(),
                                                 stamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
                                                 corpus=corpus, template=template, env=spec["env"], setup=setup, out=out,
                                                 docs=docs), encoding="utf-8")
    print(f"  arm {arm.upper()} recipe -> {out}/config.yml, train.sbatch, README.md (nothing launched)")
    return out
