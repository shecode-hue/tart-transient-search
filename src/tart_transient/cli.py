"""Command line interface."""
from __future__ import annotations

import argparse
import logging
import sys

from .config import Config
from . import pipeline, burst, timeseries


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="tart-transient",
        description="Transient search for the TART radio telescope array.")
    p.add_argument("stage",
                   choices=["download", "build", "search", "run", "transient",
                            "burst", "compare", "timeseries"],
                   help="download = fetch HDF; build = MS + calibration; "
                        "search = fit, peel, image, search; run = all three; "
                        "transient = light curve at a target direction, from HDF "
                        "(burst is a synonym); "
                        "compare = difference two completed runs; "
                        "timeseries = per-interval fit of every catalogued source")
    p.add_argument("--config", help="path to a run YAML (not needed for burst)")
    p.add_argument("--sites", default="na-unam,za-rhodes,ghana",
                   help="burst: comma-separated telescope names")
    p.add_argument("--peak", default="2025-11-11T10:00:00+00:00",
                   help="burst: reported peak time, ISO UTC")
    p.add_argument("--window", type=float, default=30.0,
                   help="burst: minutes of data centred on the peak")
    p.add_argument("--out", default="runs/burst", help="burst/compare: output directory")
    p.add_argument("--quiet-run", help="compare: the reference run directory")
    p.add_argument("--burst-run", help="compare: the run to test against it")
    p.add_argument("--target", help="compare: 'RA,DEC' in degrees to circle, "
                                    "or 'sun@ISOTIME'")
    p.add_argument("--interval", type=float, default=30.0,
                   help="timeseries: seconds per interval")
    p.add_argument("--min-elevation", type=float, default=10.0,
                   help="timeseries: ignore sources below this")
    p.add_argument("--runs-dir", default="runs", help="output root (default: runs)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(name)-28s %(message)s",
        datefmt="%H:%M:%S")

    if args.stage in ("transient", "burst"):
        return _burst(args)
    if args.stage == "compare":
        return _compare(args)
    if args.stage == "timeseries":
        return _timeseries(args)

    if not args.config:
        p.error("--config is required for this stage")
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
            exp = s.get("expected_false_positives", 0.0)
            n = s["n_confirmed"]
            print("  TRANSIENT CANDIDATES {}".format(n))
            print("  expected false alarms {:.2f}  ({} looks at alpha={})".format(
                exp, s["n_peaks"],
                cfg.section("significance").get("false_alarm_rate", 0.01)))
            if n and n <= max(1.0, 2.0 * exp):
                print("  NOTE: {} candidate(s) is consistent with the false-alarm "
                      "rate.".format(n))
                print("        Do not treat this as a detection without an "
                      "independent epoch or site.")
            print("  elapsed              {:.0f} s".format(rec["elapsed_s"]))
            print("=" * 62)
            print("  results: {}".format(cfg.results_dir))
    except Exception as exc:
        logging.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())


def _burst(args) -> int:
    import json
    from datetime import datetime, timedelta
    from pathlib import Path

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    peak = datetime.fromisoformat(args.peak)
    start = (peak - timedelta(minutes=args.window / 2.0)).isoformat()

    curves, record = {}, {}
    for site in [s.strip() for s in args.sites.split(",") if s.strip()]:
        try:
            files = burst.fetch(site, out / "data" / site, start, args.window)
            if not files:
                logging.warning("%s: no data in that window", site)
                continue
            lc = burst.light_curve(files, args.peak)
            st = burst.summarise(lc)
            curves[site] = {"lc": lc, "summary": st}
            record[site] = st
            logging.info("%s: sun %.0f deg, %d integrations, peak x%.1f at %+.2f min, "
                         "%.1f sigma over controls", site, st["elev_deg"],
                         st["n_integrations"], st["ratio"], st["peak_t_min"],
                         st["sigma_over_controls"])
        except Exception as exc:
            logging.error("%s failed: %s", site, exc)

    if not curves:
        logging.error("no site produced a light curve")
        return 1

    fig = burst.figure(curves, args.peak, out / "burst_lightcurves.png")
    with open(out / "burst_summary.json", "w") as f:
        json.dump({"peak_utc": args.peak, "sites": record}, f, indent=2)

    print("\n" + "=" * 66)
    print("  site           sun    peak/base   peak time   vs controls")
    print("  " + "-" * 62)
    for site, st in record.items():
        print("  %-14s %3.0f deg   %6.1fx    %+6.2f min   %5.1f sigma"
              % (site, st["elev_deg"], st["ratio"], st["peak_t_min"],
                 st["sigma_over_controls"]))
    peaks = [st["peak_t_min"] for st in record.values()]
    if len(peaks) > 1:
        print("\n  peak times agree to %.2f min across %d sites"
              % (max(peaks) - min(peaks), len(peaks)))
    print("=" * 66)
    print("  figure:  %s" % fig)
    print("  summary: %s" % (out / "burst_summary.json"))
    return 0


def _compare(args) -> int:
    from pathlib import Path

    if not (args.quiet_run and args.burst_run):
        logging.error("compare needs --quiet-run and --burst-run")
        return 1

    ra = dec = None
    if args.target:
        if args.target.lower().startswith("sun@"):
            from datetime import datetime
            from astropy.coordinates import get_body
            from astropy.time import Time
            s = get_body("sun", Time(datetime.fromisoformat(args.target[4:])))
            ra, dec = float(s.ra.deg), float(s.dec.deg)
            logging.info("target: Sun at RA %.3f dec %.3f", ra, dec)
        else:
            ra, dec = [float(v) for v in args.target.split(",")]

    rec = pipeline.stage_compare(Path(args.quiet_run), Path(args.burst_run),
                                 Path(args.out), ra, dec)

    print("\n" + "=" * 66)
    print("  quiet : %s  %s" % (rec["quiet"]["run"], rec["quiet"]["t_start"][11:19]))
    print("  burst : %s  %s" % (rec["burst"]["run"], rec["burst"]["t_start"][11:19]))
    t = rec.get("target")
    if t:
        print("\n  AT THE TARGET")
        print("    quiet peak   %.5f" % t["quiet"])
        print("    burst peak   %.5f" % t["burst"])
        print("    ratio        %.1fx" % t["ratio"])
        print("    change       %.1f rms" % t["change_rms"])
    print("\n  positions that changed by >5 sigma: %d" % rec["n_changed"])
    for c in rec["changed"][:5]:
        print("    RA %7.2f dec %7.2f   %.5f -> %.5f  (%.1fx, %.1f rms)"
              % (c["ra_deg"], c["dec_deg"], c["quiet"], c["burst"],
                 c["ratio"], c["change_rms"]))
    print("=" * 66)
    print("  fits   : %s/fits/04_epoch_difference.fits" % args.out)
    print("  figure : %s/plots/07_epoch_comparison.png" % args.out)
    print("  json   : %s/results/comparison.json" % args.out)
    return 0


def _timeseries(args) -> int:
    import csv
    import json
    from datetime import datetime, timedelta
    from pathlib import Path

    out = Path(args.out)
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    site = [s.strip() for s in args.sites.split(",") if s.strip()][0]
    peak = datetime.fromisoformat(args.peak)
    start = (peak - timedelta(minutes=args.window / 2.0)).isoformat()
    files = burst.fetch(site, out / "data", start, args.window)
    if not files:
        logging.error("no data for %s in that window", site)
        return 1

    res = timeseries.run(files, site, interval_s=args.interval,
                         min_elevation_deg=args.min_elevation)
    ranked = timeseries.summarise(res)
    fig = timeseries.figure(res, ranked, out / "figures" / "timeseries.png")

    with open(out / "results" / "timeseries.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "unix_time", "amplitude", "snr", "elevation_deg"])
        for name, rec in res["series"].items():
            for i in range(len(rec["t"])):
                w.writerow([name, "%.3f" % rec["t"][i], "%.6f" % rec["amp"][i],
                            "%.3f" % rec["snr"][i], "%.2f" % rec["el"][i]])
    with open(out / "results" / "ranked_variability.json", "w") as f:
        json.dump(ranked, f, indent=2)

    events = []
    for r in ranked[:5]:
        ev = timeseries.onset_offset(res, r["name"])
        if ev and ev["duration_s"] > 0:
            events.append(ev)
            timeseries.event_figure(
                res, ev, out / "figures" /
                ("event_%s.png" % r["name"].split()[0].replace("/", "_")))
    with open(out / "results" / "events.json", "w") as f:
        json.dump(events, f, indent=2)

    from astropy.time import Time
    print("\n" + "=" * 74)
    print("  %s   %s .. %s UTC   %d intervals of %.0f s"
          % (site, Time(res["t_start"], format="unix").iso[:19],
             Time(res["t_end"], format="unix").iso[11:19],
             res["n_intervals"], res["interval_s"]))
    print("  " + "-" * 70)
    print("   sigma   peak/med   elev   peak at            source")
    for r in ranked[:10]:
        print("  %6.1f   %7.2f   %4.0f   %s   %s"
              % (r["sigma"], r["ratio"], r["median_el"],
                 Time(r["peak_t"], format="unix").iso[11:19], r["name"][:30]))
    if events:
        print("\n  EVENTS  (rise above 3 sigma of own baseline)")
        print("   onset      peak       end        dur    rise   decay  source")
        for ev in events:
            print("   %s  %s  %s  %5.0fs %5.0fs %5.0fs  %s"
                  % (Time(ev["onset_t"], format="unix").iso[11:19],
                     Time(ev["peak_t"], format="unix").iso[11:19],
                     Time(ev["offset_t"], format="unix").iso[11:19],
                     ev["duration_s"], ev["rise_s"], ev["decay_s"],
                     ev["name"][:26]))
    print("=" * 74)
    print("  figure : %s" % fig)
    print("  series : %s" % (out / "results" / "timeseries.csv"))
    print("  ranked : %s" % (out / "results" / "ranked_variability.json"))
    return 0
