FITS for the decay epoch are on disk here but not committed — 16 MB each.

Regenerate:

    PYTHONPATH=src python3 -m tart_transient run --config config/unam-ev-decay.yaml

Produces 01_before.fits, 02_after.fits, 03_removed.fits in runs/unam-ev-decay/fits/.
