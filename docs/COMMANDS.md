# Commands

## Start

```bash
open -a "Docker Desktop"
docker start shain_tart_jammy
docker exec -it shain_tart_jammy bash
```

Inside the container:

```bash
source /workspace/gro_demo/.venv/bin/activate
cd /workspace/tart
git checkout solar-burst-detection
git pull
```

## Check

```bash
python3 scripts/apply_patches.py
PYTHONPATH=src python3 -m pytest tests/ -q
```

Expect three "already patched" lines, then `7 passed`.

## Transient light curves — 4 sites, 15 min

```bash
PYTHONPATH=src python3 -m tart_transient burst \
  --sites na-unam,za-rhodes,ghana,mu-udm \
  --peak "2025-11-11T10:00:00+00:00" \
  --window 30 --out runs/burst
```

Out: `runs/burst/burst_lightcurves.png`, `runs/burst/burst_summary.json`

## Before/after images — 1 site, 1 hour

```bash
bash scripts/run_all_sites.sh mu-udm
```

Out: `runs/mu-udm-1111q/`, `runs/mu-udm-1111b/`, `runs/burst_before_after.png`

More sites, 1 hour each:

```bash
bash scripts/run_all_sites.sh na-unam
bash scripts/run_all_sites.sh ghana
bash scripts/run_all_sites.sh za-rhodes
```

All four at once (4 hours):

```bash
bash scripts/run_all_sites.sh
```

## Look

```bash
exit
open ~/Desktop/tart-transient-search/runs/burst/burst_lightcurves.png
open ~/Desktop/tart-transient-search/runs/burst_before_after.png
```

FITS for CARTA:

```
runs/mu-udm-1111q/fits/01_before.fits
runs/mu-udm-1111b/fits/01_before.fits
```

## Another burst

```bash
PYTHONPATH=src python3 -m tart_transient burst \
  --sites za-rhodes,ghana,mu-udm \
  --peak "2025-11-09T07:12:00+00:00" \
  --window 30 --out runs/burst-1109
```

## If a run fails

```bash
cat runs/<name>.log
```

## If the container is gone

```bash
docker run -dit --name shain_tart_jammy \
  -v ~/Desktop/gro_demo:/workspace/gro_demo \
  -v ~/Desktop/tart-transient-search:/workspace/tart \
  shain_tart:snapshot-20260824 bash
```
