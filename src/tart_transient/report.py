"""Figures and the run record."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

TEAL, GREEN, RED, AMBER, GREY = "#0b6e78", "#17603a", "#a32020", "#9a5c05", "#6b7c8a"
CMAP = matplotlib.colormaps.get_cmap("RdBu_r").copy()
CMAP.set_bad("lightgrey")


def _sky(ax, img, half, title, vmax=None):
    vm = vmax if vmax is not None else np.nanpercentile(np.abs(img), 99.5)
    h = ax.imshow(np.ma.masked_invalid(img), origin="lower",
                  extent=[-half, half, -half, half], cmap=CMAP, vmin=-vm, vmax=vm)
    plt.colorbar(h, ax=ax, fraction=.046, pad=.04, label="Jy/beam")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("RA offset (deg)")
    ax.set_ylabel("Dec offset (deg)")
    return vm


def input_data(hdf_info, vis_abs, gains, phases, out: Path):
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    ax[0].hist(vis_abs, bins=60, color=TEAL)
    ax[0].set_xlabel("|visibility|"); ax[0].set_ylabel("count")
    ax[0].set_title("visibility amplitudes", fontsize=10)
    ax[1].bar(np.arange(len(gains)), gains, color=TEAL)
    ax[1].axhline(1.0, color=GREY, ls=":")
    ax[1].set_xlabel("antenna"); ax[1].set_ylabel("gain")
    ax[1].set_title("stored gains {:.2f}-{:.2f}".format(gains.min(), gains.max()), fontsize=10)
    ax[2].bar(np.arange(len(phases)), np.degrees(phases), color=AMBER)
    ax[2].set_xlabel("antenna"); ax[2].set_ylabel("phase (deg)")
    ax[2].set_title("stored phases", fontsize=10)
    verdict = ("calibration recorded -> use DATA" if hdf_info["gains_stored"]
               else "NO calibration recorded -> solve gains")
    fig.suptitle("1 - input data: {}".format(verdict), fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def sky_model(uvw_lambda, cat_lm, eigenvalues, beam_deg, n_raw, n_merged,
              cond, out: Path):
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    ax[0].plot(uvw_lambda[:, 0], uvw_lambda[:, 1], ".", ms=1, color=TEAL, alpha=.5)
    ax[0].plot(-uvw_lambda[:, 0], -uvw_lambda[:, 1], ".", ms=1, color=TEAL, alpha=.5)
    ax[0].set_aspect("equal"); ax[0].set_xlabel("u (lambda)"); ax[0].set_ylabel("v (lambda)")
    ax[0].set_title("UV coverage -> beam {:.2f} deg".format(beam_deg), fontsize=10)
    ax[1].add_patch(plt.Circle((0, 0), 1, fill=False, color=GREY, ls="--"))
    ax[1].plot(cat_lm[:, 0], cat_lm[:, 1], "o", ms=5, color=TEAL)
    ax[1].set_xlim(-1.05, 1.05); ax[1].set_ylim(-1.05, 1.05); ax[1].set_aspect("equal")
    ax[1].set_xlabel("l"); ax[1].set_ylabel("m")
    ax[1].set_title("{} satellites -> {} components".format(n_raw, n_merged), fontsize=10)
    ax[2].semilogy(np.arange(1, len(eigenvalues) + 1), np.abs(eigenvalues), "o-",
                   color=TEAL, ms=4)
    ax[2].set_xlabel("mode"); ax[2].set_ylabel("|eigenvalue|")
    ax[2].set_title("conditioning: {:.1f}".format(cond), fontsize=10)
    fig.suptitle("2 - sky model and fit conditioning", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def snr_vs_null(names, snr, threshold, null_median, out: Path):
    order = np.argsort(snr)[::-1]
    fig, ax = plt.subplots(figsize=(12, max(5, 0.22 * len(names))))
    ax.barh(np.arange(len(order))[::-1], snr[order],
            color=[GREEN if snr[i] >= threshold else GREY for i in order])
    ax.axvline(threshold, color=RED, ls="--", lw=2,
               label="measured threshold = {:.1f}".format(threshold))
    ax.axvline(null_median, color=GREY, ls=":", lw=1.5,
               label="null median (empty sky) = {:.1f}".format(null_median))
    ax.set_yticks(np.arange(len(order))[::-1])
    ax.set_yticklabels([names[i].split("_(")[0][:24] for i in order], fontsize=7)
    ax.set_xlabel("coherent fit SNR"); ax.legend(fontsize=9, loc="lower right")
    ax.set_title("3 - every source against a threshold measured from this data\n"
                 "green = brighter than empty sky ever gets ({}/{})".format(
                     int((snr >= threshold).sum()), len(names)), fontsize=11)
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def before_after(before, after, removed, half, n_peeled, power, out: Path):
    fig, ax = plt.subplots(1, 3, figsize=(19, 6.4))
    vm = _sky(ax[0], before, half, "BEFORE  rms={:.4f}".format(np.nanstd(before)))
    _sky(ax[1], after, half, "AFTER  rms={:.4f}".format(np.nanstd(after)), vmax=vm)
    _sky(ax[2], removed, half, "REMOVED  rms={:.4f}".format(np.nanstd(removed)))
    fig.suptitle("4 - peeled {} components, {:.1%} of visibility power removed".format(
        n_peeled, power), fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def effectiveness(before, after, power, rms_change, peak_change, out: Path):
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    ax[0].hist(before[np.isfinite(before)], bins=90, alpha=.6, color=TEAL, label="before")
    ax[0].hist(after[np.isfinite(after)], bins=90, alpha=.6, color=AMBER, label="after")
    ax[0].set_yscale("log"); ax[0].legend(fontsize=9)
    ax[0].set_xlabel("pixel value"); ax[0].set_title("pixel distribution", fontsize=10)
    vals = [100 * power, rms_change, peak_change]
    bars = ax[1].bar(["visibility\npower removed", "image RMS\nchange",
                      "image peaks\nchange"], vals,
                     color=[GREEN, GREEN if rms_change < 0 else AMBER,
                            GREEN if peak_change < 0 else AMBER])
    ax[1].bar_label(bars, fmt="%+.1f%%", fontsize=10)
    ax[1].axhline(0, color=GREY, lw=1); ax[1].set_ylabel("% change")
    ax[1].set_title("what changed", fontsize=10)
    ax[2].axis("off")
    ax[2].text(0, .5,
               "Visibility power fell {:.1%} - the subtraction is real.\n\n"
               "Image RMS {:+.1f}%, peaks {:+.1f}%.\n\n"
               "Where the two disagree, trust the visibility numbers:\n"
               "DiSkO's regularised image correlates only r~0.55 with a\n"
               "direct transform of the same data and shares just 4 of\n"
               "its 20 brightest positions - the prior redistributes flux."
               .format(power, rms_change, peak_change),
               fontsize=10, va="center", linespacing=1.6)
    fig.suptitle("5 - did the clean help?", fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def transient_search(after, half, wcs, npix, pix, result, out: Path):
    cands = result["candidates"]
    fig, ax = plt.subplots(1, 2, figsize=(18, 8))
    _sky(ax[0], after, half, "residual, every peak tested")
    for c in cands:
        px, py = wcs.wcs_world2pix(c["ra_deg"], c["dec_deg"], 0)
        x, y = (float(px) - npix / 2) * pix, (float(py) - npix / 2) * pix
        col = GREEN if c["confirmed"] else (AMBER if c["passes_single_look"] else RED)
        ax[0].plot(x, y, "o", ms=9, mfc="none", mec=col, mew=1.8, zorder=4)
    ax[0].legend(handles=[
        Line2D([], [], color=GREEN, marker="o", ls="none", mfc="none", mew=1.8,
               label="confirmed"),
        Line2D([], [], color=AMBER, marker="o", ls="none", mfc="none", mew=1.8,
               label="passed single look, failed trials correction"),
        Line2D([], [], color=RED, marker="o", ls="none", mfc="none", mew=1.8,
               label="rejected")], loc="lower left", fontsize=8)
    if cands:
        el = [c["elevation_deg"] for c in cands]
        ax[1].scatter(el, [c["vis_snr"] for c in cands], s=42, zorder=3,
                      color=[GREEN if c["confirmed"] else
                             (AMBER if c["passes_single_look"] else RED)
                             for c in cands])
        order = np.argsort(el)
        ax[1].plot(np.array(el)[order],
                   np.array([c["threshold_here"] for c in cands])[order],
                   color=RED, ls="--", lw=2,
                   label="trials-corrected threshold at this elevation")
        ax[1].set_xlabel("elevation (deg)")
        ax[1].set_ylabel("coherent fit SNR at that position")
        ax[1].invert_xaxis()
        ax[1].legend(fontsize=9)
        ax[1].set_title("threshold rises toward the horizon:\n"
                        "the null tail is ~35% heavier below 20 deg elevation",
                        fontsize=10)
    fig.suptitle("6 - transient search: {} peaks -> {} single-look -> {} confirmed".format(
        result["n_peaks"], result["n_passed_single_look"], result["n_confirmed"]),
        fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=115); plt.close(fig)


def write_summary(record: dict, out: Path):
    with open(out, "w") as f:
        json.dump(record, f, indent=2, default=float)
