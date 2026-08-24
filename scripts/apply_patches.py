#!/usr/bin/env python3
"""Apply the upstream patches described in PATCHES.md."""
from __future__ import annotations

import sys
from pathlib import Path


def find_package(name: str) -> Path:
    import importlib.util
    spec = importlib.util.find_spec(name)
    if spec is None or not spec.submodule_search_locations:
        raise SystemExit(f"{name} is not installed")
    return Path(list(spec.submodule_search_locations)[0])


def patch_tart2ms() -> None:
    root = find_package("tart2ms")
    core, cli = root / "tart2ms.py", root / "scripts" / "tart2ms.py"

    c = cli.read_text()
    if "--filter-elevation" not in c:
        anchor = '    parser.add_argument(\n        "--uncalibrated",'
        if anchor not in c:
            raise SystemExit("tart2ms CLI layout changed; patch by hand")
        c = c.replace(anchor,
            '    parser.add_argument(\n'
            '        "--filter-elevation", dest="filter_elevation", type=float,\n'
            '        default=45.0, required=False,\n'
            '        help="Minimum satellite elevation (deg) for the GNSS catalogue. "\n'
            '             "45 (the historical default) discards most satellites in a "\n'
            '             "170deg TART field of view; use ~5-10 for the whole sky.",\n'
            '    )\n' + anchor, 1)
        c = c.replace("            cat_name_prefix=ARGS.model_catalog_name_prefix,",
                      "            cat_name_prefix=ARGS.model_catalog_name_prefix,\n"
                      "            filter_elevation=ARGS.filter_elevation,")
        cli.write_text(c)
        print("  patched scripts/tart2ms.py")
    else:
        print("  scripts/tart2ms.py already patched")

    s = core.read_text()
    changed = False
    for fn in ("ms_from_hdf5", "ms_from_json", "ms_create"):
        i = s.index(f"def {fn}(")
        j = s.index("):", i)
        sig = s[i:j]
        if "filter_elevation" in sig:
            continue
        for anchor in ('cat_name_prefix=None', 'cat_name_prefix="model_sources_",'):
            if anchor in sig:
                s = s[:i] + sig.replace(
                    anchor, anchor.rstrip(",") + ",\n    filter_elevation=45.0") + s[j:]
                changed = True
                break
    for old, new in [
        ("""                online_sources, online_sources_timestamps = __fetch_sources(
                    timestamps=timestamps,
                    observer_lat=lat,
                    observer_lon=lon,
                    force_recache=catalog_recache,
                )""",
         """                online_sources, online_sources_timestamps = __fetch_sources(
                    timestamps=timestamps,
                    observer_lat=lat,
                    observer_lon=lon,
                    force_recache=catalog_recache,
                    filter_elevation=filter_elevation,
                )"""),
        ("                    filter_elevation=20.0,",
         "                    filter_elevation=filter_elevation,"),
    ]:
        if old in s:
            s = s.replace(old, new)
            changed = True

    i = s.index("            model_data = predict_model(")
    j = s.index("            )", i)
    if "filter_elevation" not in s[i:j]:
        s = s[:i] + s[i:j].replace(
            "                cat_name_prefix=cat_name_prefix,",
            "                filter_elevation=filter_elevation,\n"
            "                cat_name_prefix=cat_name_prefix,") + s[j:]
        changed = True
    start = 0
    while True:
        try:
            k = s.index("    ms_create(", start)
        except ValueError:
            break
        end = s.index("\n    )", k)
        call = s[k:end]
        if "filter_elevation" not in call and "cat_name_prefix=cat_name_prefix," in call:
            new = call.replace("        cat_name_prefix=cat_name_prefix,",
                               "        cat_name_prefix=cat_name_prefix,\n"
                               "        filter_elevation=filter_elevation,")
            s = s[:k] + new + s[end:]
            changed = True
            start = k + len(new)
        else:
            start = end
    if changed:
        core.write_text(s)
        print("  patched tart2ms.py")
    else:
        print("  tart2ms.py already patched")


def patch_tartcargo() -> None:
    try:
        root = find_package("tartcargo")
    except SystemExit:
        print("  tartcargo not installed - skipping (only needed for Stimela recipes)")
        return
    cab = root / "tart2ms.yml"
    s = cab.read_text()
    if "filter-elevation" in s:
        print("  tartcargo/tart2ms.yml already patched")
        return
    old = """            write-model-catalog:
                dtype: bool
                default: false"""
    s = s.replace(old,
        """            filter-elevation:
                dtype: float
                default: 45.0
                required: false
                info: Minimum satellite elevation in degrees for the GNSS catalogue.
""" + old, 1)
    cab.write_text(s)
    print("  patched tartcargo/tart2ms.yml")


if __name__ == "__main__":
    print("applying upstream patches (see PATCHES.md)")
    patch_tart2ms()
    patch_tartcargo()
    print("done")
