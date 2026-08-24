"""Command line interface."""
from __future__ import annotations

import argparse
import logging
import sys

from .config import Config
from . import pipeline


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="tart-transient",
        description="Transient search for the TART radio telescope array.")
    p.add_argument("stage", choices=["download", "build", "search", "run"],
                   help="download = fetch HDF; build = MS + calibration; "
                        "search = fit, peel, image, search; run = all three")
    p.add_argument("--config", required=True, help="path to a run YAML")
    p.add_argument("--runs-dir", default="runs", help="output root (default: runs)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(name)-28s %(message)s",
        datefmt="%H:%M:%S")

    cfg = Config.load(args.config, args.runs_dir)
    logging.info("run '%s' -> %s", cfg.name, cfg.root)

    try:
        if args.stage in ("download", "run"):
            pipeline.stage_download(cfg)
        if args.stage in ("build", "run"):
            pipeline.stage_build(cfg)
        if args.stage in ("search", "run"):
            rec = pipeline.stage_search(cfg)
            s = rec["search"]
            print("\n" + "=" * 62)
            print("  power removed        {:.2%}".format(rec["fit"]["power_removed"]))
            print("  image RMS            {:.4f} -> {:.4f}".format(
                rec["images"]["rms_before"], rec["images"]["rms_after"]))
            print("  residual peaks       {}".format(s["n_peaks"]))
            print("  passed single look   {}".format(s["n_passed_single_look"]))
            print("  TRANSIENT CANDIDATES {}".format(s["n_confirmed"]))
            print("  elapsed              {:.0f} s".format(rec["elapsed_s"]))
            print("=" * 62)
            print("  results: {}".format(cfg.results_dir))
    except Exception as exc:            # noqa: BLE001
        logging.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
