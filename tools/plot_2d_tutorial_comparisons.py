#!/usr/bin/env python3
"""
Generate 2D-vs-3D comparison plots for the MOCASSIN 2D mode tutorial
(TUTORIAL_2D_MODE.md).

Parses MOCASSIN spatial grid outputs (grid0.out, grid1.out, grid2.out,
grid3.out) for single-chemistry (e.g., PN150 with `Hdensity`) and
multi-chemistry (e.g., test_Problem with `densityFile` and `Ndust file`)
simulations, and produces comparison figures of:
  1. Electron temperature (Te), electron density (Ne), and hydrogen density (nH)
  2. Ionic fractions (H, He, O, N, S) across ionization stages

Both radial/directional cuts (with residuals) and meridional (R, Z) spatial
maps are generated at 300 DPI.

Usage example:
  python3 tools/plot_2d_tutorial_comparisons.py \\
      --pn150-3d /path/to/PN150_3D \\
      --pn150-2d /path/to/PN150_2D_cylinder \\
      --pn150-2d-res2x /path/to/PN150_2D_cylinder_res2x \\
      --test-3d /path/to/test_Problem_3D \\
      --test-2d /path/to/test_Problem_2D_cylinder \\
      --out-dir docs/tutorial_2d_figures
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Arc


# Publication-clean styling matching paper_2d aesthetics
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9.5,
    "axes.titlesize": 10.0,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 7.8,
    "axes.linewidth": 0.85,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.06,
})

# Colorblind-safe Okabe-Ito palette
C_3D = "#1f1f1f"
C_CYL = "#d55e00"
C_RES2X = "#0072b2"
C_SPEC = "#009e73"
C_BUF = "#cc79a7"

ELEM_SYMBOLS = {
    1: "H", 2: "He", 6: "C", 7: "N", 8: "O",
    10: "Ne", 12: "Mg", 14: "Si", 16: "S", 18: "Ar", 26: "Fe",
}


def _resolve_abundance_path(model_dir, ab_rel, repo_root):
    """Resolve abundance file path referenced in grid3.out."""
    candidates = [
        ab_rel if os.path.isabs(ab_rel) else os.path.join(model_dir, ab_rel),
        os.path.join(model_dir, "input", os.path.basename(ab_rel)),
        os.path.join(repo_root, "benchmarks", "gas", "PN150", os.path.basename(ab_rel)),
        os.path.join(repo_root, "benchmarks", "test_Problem", "input", os.path.basename(ab_rel)),
        os.path.join(repo_root, ab_rel),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    raise FileNotFoundError(
        f"Cannot find abundance file '{ab_rel}' referenced by '{model_dir}/output/grid3.out'."
    )


def load_spatial_grid(model_dir, repo_root=None):
    """
    Parse MOCASSIN grid0.out, grid1.out, grid2.out, grid3.out into a dict.

    Supports both single-chemistry and multiChemistry (nAbComponents > 1) runs,
    as well as 3D Cartesian and 2D (planar or cylinder) grids.
    """
    if repo_root is None:
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_dir = os.path.join(model_dir, "output")

    g3_path = os.path.join(out_dir, "grid3.out")
    if not os.path.isfile(g3_path):
        raise FileNotFoundError(f"Missing grid3.out in {out_dir}")
    with open(g3_path) as f:
        g3_lines = [line.strip() for line in f]

    ab_files = []
    nstages = 7
    lg_sym = True
    for idx, line in enumerate(g3_lines):
        if "lgSymmetricXYZ" in line:
            lg_sym = line.split()[0].upper().startswith("T")
        if "nAbComponents" in line:
            n_ab = int(line.split()[0])
            for k in range(1, n_ab + 1):
                ab_rel = g3_lines[idx + k].split()[0].strip("\"'")
                ab_files.append(_resolve_abundance_path(model_dir, ab_rel, repo_root))
        if "emittingGrid" in line and idx >= 30:
            try:
                nstages = int(line.split()[0])
            except ValueError:
                pass

    # Active element union mirrors lgElementOn(1:nElements) in source/grid_mod.f90
    elem_set = set()
    for af in ab_files:
        with open(af) as f:
            for elem_idx in range(1, 31):
                raw = f.readline()
                if not raw:
                    break
                toks = raw.split()
                if toks and float(toks[0]) > 1e-12:
                    elem_set.add(elem_idx)
    active_elems = sorted(elem_set)

    with open(os.path.join(out_dir, "grid0.out")) as f0:
        _ = int(f0.readline().split()[0])
        hdr = f0.readline().split()
        nx, ny, nz, ncells = int(hdr[0]), int(hdr[1]), int(hdr[2]), int(hdr[3])
        r_out = float(hdr[5])
        x_ax = np.array([float(f0.readline().strip()) for _ in range(nx)])
        y_ax = np.array([float(f0.readline().strip()) for _ in range(ny)])
        z_ax = np.array([float(f0.readline().strip()) for _ in range(nz)])
        active = np.zeros((nx, ny, nz), dtype=int)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    active[i, j, k] = int(f0.readline().split()[0])

    te = np.zeros((nx, ny, nz), dtype=float)
    ne = np.zeros((nx, ny, nz), dtype=float)
    hden = np.zeros((nx, ny, nz), dtype=float)
    ab_idx = np.ones((nx, ny, nz), dtype=int)
    with open(os.path.join(out_dir, "grid1.out")) as f1:
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    parts = f1.readline().split()
                    te[i, j, k] = float(parts[0])
                    ne[i, j, k] = float(parts[1])
                    hden[i, j, k] = float(parts[2])
                    if len(parts) >= 4:
                        ab_idx[i, j, k] = int(parts[3])

    ion_fracs = {
        el: np.zeros((nx, ny, nz, min(el + 1, nstages)), dtype=float)
        for el in active_elems
    }
    with open(os.path.join(out_dir, "grid2.out")) as f2:
        lines2 = f2.readlines()
    expected = nx * ny * nz * len(active_elems)
    if len(lines2) != expected:
        raise ValueError(
            f"Line count mismatch in {out_dir}/grid2.out: found {len(lines2)}, "
            f"expected {expected} (nx*ny*nz*n_active_elems = {nx}*{ny}*{nz}*{len(active_elems)})."
        )

    ptr = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for el in active_elems:
                    vals = [float(v) for v in lines2[ptr].split()]
                    n_st = min(el + 1, nstages)
                    ion_fracs[el][i, j, k, :n_st] = vals[:n_st]
                    ptr += 1

    X, Y, Z = np.meshgrid(x_ax, y_ax, z_ax, indexing="ij")
    R_sph = np.sqrt(X**2 + Y**2 + Z**2)

    return {
        "dir": model_dir,
        "nx": nx, "ny": ny, "nz": nz, "ncells": ncells,
        "r_out": r_out, "lg_sym": lg_sym,
        "x_ax": x_ax, "y_ax": y_ax, "z_ax": z_ax,
        "X": X, "Y": Y, "Z": Z, "R_sph": R_sph,
        "active": active, "te": te, "ne": ne, "hden": hden,
        "ab_idx": ab_idx, "ions": ion_fracs, "active_elems": active_elems,
    }


def cell_edges_from_centers(centers, symmetric=True):
    """
    Construct MOCASSIN cell boundary edges from node coordinates.
    Under symmetricXYZ in MOCASSIN, boundary nodes sit on outer/inner walls
    (half-width cells), matching getVolume() in source/grid_mod.f90.
    """
    n = len(centers)
    if n == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5])
    edges = np.zeros(n + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    if symmetric and centers[0] == 0.0:
        edges[0] = 0.0
    else:
        edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] if symmetric else centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def radial_bin_profile(r_vals, data_vals, mask, bin_edges):
    """Compute binned mean profile and bin centers for active masked cells."""
    r = r_vals[mask]
    v = data_vals[mask]
    nbins = len(bin_edges) - 1
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    means = np.full(nbins, np.nan)
    for b in range(nbins):
        if b < nbins - 1:
            sel = (r >= bin_edges[b]) & (r < bin_edges[b + 1])
        else:
            sel = (r >= bin_edges[b]) & (r <= bin_edges[b + 1])
        if np.any(sel):
            means[b] = np.mean(v[sel])
    valid = ~np.isnan(means)
    return centers[valid], means[valid]


def active_gas_mask(g):
    """Mask selecting active cells containing gas (Hden > 0)."""
    return (g["active"] > 0) & (g["hden"] > 0.0)


def mask_slice_y0(g):
    """Extract 2D (nx, nz) slice at j = 0 (y = 0) with NaN on inactive cells."""
    m = (g["active"][:, 0, :] > 0) & (g["hden"][:, 0, :] > 0.0)
    return m


def _plot_pn150_panel_map(ax, x_ed, z_ed, arr2d, mask2d, cmap, norm=None, vmin=None, vmax=None,
                          r_in=1.0e17, r_out=4.07e17, bad_color="#e6e6e6"):
    """Helper to plot an (R, Z) / (X, Z) meridional heatmap on clean MOCASSIN cell edges."""
    masked = np.ma.masked_where(~mask2d, arr2d)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(color=bad_color)
    kw = {"cmap": cm, "shading": "flat"}
    if norm is not None:
        kw["norm"] = norm
    else:
        kw["vmin"] = vmin
        kw["vmax"] = vmax
    mesh = ax.pcolormesh(x_ed / 1e17, z_ed / 1e17, masked.T, **kw)
    for r_circle, ls in [(r_in / 1e17, ":"), (r_out / 1e17, "--")]:
        arc = Arc((0.0, 0.0), 2 * r_circle, 2 * r_circle, theta1=0, theta2=90,
                  color="#111111", lw=0.85, ls=ls, alpha=0.85)
        ax.add_patch(arc)
    ax.set_aspect("equal")
    ax.set_xlim(x_ed[0] / 1e17, x_ed[-1] / 1e17)
    ax.set_ylim(z_ed[0] / 1e17, z_ed[-1] / 1e17)
    return mesh


def plot_pn150_te_density(models, out_path, r_in=1.0e17, r_out=4.07e17):
    """
    Generate Case 1 (PN150, `Hdensity 3000.`) Te, Ne, and nH comparison figure.
    """
    g3d = models["3d"]
    g2d = models["2d_cyl"]
    g2dx = models.get("2d_cyl_2x")
    gspec = models.get("2d_spec")
    gbuf = models.get("2d_buf")

    m3d = active_gas_mask(g3d)
    m2d = active_gas_mask(g2d)
    bin_edges = np.linspace(r_in, r_out, 11)

    r_3d, te_3d = radial_bin_profile(g3d["R_sph"], g3d["te"], m3d, bin_edges)
    r_2d, te_2d = radial_bin_profile(g2d["R_sph"], g2d["te"], m2d, bin_edges)
    _, ne_3d = radial_bin_profile(g3d["R_sph"], g3d["ne"], m3d, bin_edges)
    _, ne_2d = radial_bin_profile(g2d["R_sph"], g2d["ne"], m2d, bin_edges)
    _, nh_3d = radial_bin_profile(g3d["R_sph"], g3d["hden"], m3d, bin_edges)
    _, nh_2d = radial_bin_profile(g2d["R_sph"], g2d["hden"], m2d, bin_edges)

    fig = plt.figure(figsize=(12.2, 8.3), dpi=300)
    outer = GridSpec(2, 1, height_ratios=[1.18, 0.92], hspace=0.30, figure=fig)

    # Top row: 2 radial profile stacks (main + % diff)
    top_gs = outer[0].subgridspec(2, 2, height_ratios=[3.0, 1.15], hspace=0.06, wspace=0.22)
    ax_te = fig.add_subplot(top_gs[0, 0])
    ax_te_res = fig.add_subplot(top_gs[1, 0], sharex=ax_te)
    ax_ne = fig.add_subplot(top_gs[0, 1])
    ax_ne_res = fig.add_subplot(top_gs[1, 1], sharex=ax_ne)

    # Scatter cloud of 3D grid nodes shows intrinsic Cartesian spread
    ax_te.scatter(g3d["R_sph"][m3d] / 1e17, g3d["te"][m3d] / 1e3, s=6, color="#bbbbbb",
                  alpha=0.35, edgecolors="none", zorder=1, label="3D active cells ($13^3$)")
    ax_te.plot(r_3d / 1e17, te_3d / 1e3, color=C_3D, lw=2.2, zorder=4, label="3D shell mean ($13^3$)")
    ax_te.plot(r_2d / 1e17, te_2d / 1e3, color=C_CYL, lw=1.9, ls="--", marker="o", ms=4.5,
               mfc="white", mew=1.3, zorder=5, label="2D Cylinder ($13\\times 13$)")
    if g2dx is not None:
        m2dx = active_gas_mask(g2dx)
        r_2dx, te_2dx = radial_bin_profile(g2dx["R_sph"], g2dx["te"], m2dx, np.linspace(r_in, r_out, 17))
        ax_te.plot(r_2dx / 1e17, te_2dx / 1e3, color=C_RES2X, lw=1.6, ls="-.", zorder=4,
                   label="2D Cylinder ($25\\times 25$)")
    if gspec is not None:
        ms = active_gas_mask(gspec)
        rs, tes = radial_bin_profile(gspec["R_sph"], gspec["te"], ms, bin_edges)
        ax_te.plot(rs / 1e17, tes / 1e3, color=C_SPEC, lw=1.2, ls=":", alpha=0.85,
                   label="2D XZ Specular ($13\\times 13$)")
    if gbuf is not None:
        mb = active_gas_mask(gbuf)
        rb, teb = radial_bin_profile(gbuf["R_sph"], gbuf["te"], mb, bin_edges)
        ax_te.plot(rb / 1e17, teb / 1e3, color=C_BUF, lw=1.2, ls=":", alpha=0.85,
                   label="2D XZ Buffer ($13\\times 13$)")

    ax_te.set_ylabel(r"Electron Temperature $T_e$ [$10^3\ \mathrm{K}$]")
    ax_te.set_title(r"(a) PN150 Radial Electron Temperature $T_e(r)$ ($n_{\mathrm{H}} = 3000\ \mathrm{cm}^{-3}$)", fontweight="bold")
    ax_te.set_ylim(9.3, 18.4)
    ax_te.grid(alpha=0.22, ls="--")
    ax_te.legend(loc="upper right", frameon=True, framealpha=0.92, ncol=1)
    plt.setp(ax_te.get_xticklabels(), visible=False)

    # Te residuals vs 3D shell mean
    ax_te_res.axhline(0.0, color=C_3D, lw=1.0, ls="-")
    common_te = np.interp(r_2d, r_3d, te_3d)
    ax_te_res.plot(r_2d / 1e17, (te_2d - common_te) / common_te * 100.0,
                   color=C_CYL, lw=1.5, marker="o", ms=3.8, label="Cyl $13\\times 13$")
    if g2dx is not None:
        _, te_2dx_b = radial_bin_profile(g2dx["R_sph"], g2dx["te"], m2dx, bin_edges)
        ax_te_res.plot(r_2d / 1e17, (te_2dx_b - common_te) / common_te * 100.0,
                       color=C_RES2X, lw=1.3, ls="-.", marker="^", ms=3.2, label="Cyl $25\\times 25$")
    if gspec is not None:
        _, tes_b = radial_bin_profile(gspec["R_sph"], gspec["te"], ms, bin_edges)
        ax_te_res.plot(r_2d / 1e17, (tes_b - common_te) / common_te * 100.0,
                       color=C_SPEC, lw=1.1, ls=":", label="Specular")
    if gbuf is not None:
        _, teb_b = radial_bin_profile(gbuf["R_sph"], gbuf["te"], mb, bin_edges)
        ax_te_res.plot(r_2d / 1e17, (teb_b - common_te) / common_te * 100.0,
                       color=C_BUF, lw=1.1, ls=":", label="Buffer")
    ax_te_res.set_ylabel(r"$\Delta T_e / T_{e,\mathrm{3D}}$ [%]")
    ax_te_res.set_xlabel(r"Spherical Radius $r$ [$10^{17}\ \mathrm{cm}$]")
    ax_te_res.set_ylim(-8.0, 12.5)
    ax_te_res.grid(alpha=0.22, ls="--")

    # (b) Ne and nH radial profile
    ax_ne.scatter(g3d["R_sph"][m3d] / 1e17, g3d["ne"][m3d], s=6, color="#bbbbbb",
                  alpha=0.35, edgecolors="none", zorder=1, label="3D active cells ($N_e$)")
    ax_ne.plot(r_3d / 1e17, ne_3d, color=C_3D, lw=2.2, zorder=4, label="3D $N_e$ mean ($13^3$)")
    ax_ne.plot(r_2d / 1e17, ne_2d, color=C_CYL, lw=1.9, ls="--", marker="o", ms=4.5,
               mfc="white", mew=1.3, zorder=5, label="2D Cylinder $N_e$ ($13\\times 13$)")
    if g2dx is not None:
        r_2dx, ne_2dx = radial_bin_profile(g2dx["R_sph"], g2dx["ne"], m2dx, np.linspace(r_in, r_out, 17))
        ax_ne.plot(r_2dx / 1e17, ne_2dx, color=C_RES2X, lw=1.6, ls="-.", zorder=4,
                   label="2D Cylinder $N_e$ ($25\\times 25$)")
    ax_ne.plot(r_3d / 1e17, nh_3d, color="#666666", lw=1.5, ls="-", label=r"$n_{\mathrm{H}}$ = 3000 $\mathrm{cm}^{-3}$ (3D & 2D)")
    ax_ne.set_ylabel(r"Density [$N_e,\ n_{\mathrm{H}}$ in $\mathrm{cm}^{-3}$]")
    ax_ne.set_title(r"(b) PN150 Electron & Hydrogen Density ($N_e$, $n_{\mathrm{H}}$)", fontweight="bold")
    ax_ne.set_ylim(1600, 3800)
    ax_ne.grid(alpha=0.22, ls="--")
    ax_ne.legend(loc="lower left", frameon=True, framealpha=0.92)
    plt.setp(ax_ne.get_xticklabels(), visible=False)

    ax_ne_res.axhline(0.0, color=C_3D, lw=1.0, ls="-")
    common_ne = np.interp(r_2d, r_3d, ne_3d)
    ax_ne_res.plot(r_2d / 1e17, (ne_2d - common_ne) / common_ne * 100.0,
                   color=C_CYL, lw=1.5, marker="o", ms=3.8, label="Cyl $13\\times 13$")
    if g2dx is not None:
        _, ne_2dx_b = radial_bin_profile(g2dx["R_sph"], g2dx["ne"], m2dx, bin_edges)
        ax_ne_res.plot(r_2d / 1e17, (ne_2dx_b - common_ne) / common_ne * 100.0,
                       color=C_RES2X, lw=1.3, ls="-.", marker="^", ms=3.2, label="Cyl $25\\times 25$")
    if gspec is not None:
        _, nes_b = radial_bin_profile(gspec["R_sph"], gspec["ne"], ms, bin_edges)
        ax_ne_res.plot(r_2d / 1e17, (nes_b - common_ne) / common_ne * 100.0,
                       color=C_SPEC, lw=1.1, ls=":", label="Specular")
    ax_ne_res.set_ylabel(r"$\Delta N_e / N_{e,\mathrm{3D}}$ [%]")
    ax_ne_res.set_xlabel(r"Spherical Radius $r$ [$10^{17}\ \mathrm{cm}$]")
    ax_ne_res.set_ylim(-6.0, 6.0)
    ax_ne_res.grid(alpha=0.22, ls="--")

    # Bottom row: 2D (X/R, Z) meridional maps at 13x13
    bot_gs = outer[1].subgridspec(1, 3, wspace=0.22)
    ax_m1 = fig.add_subplot(bot_gs[0, 0])
    ax_m2 = fig.add_subplot(bot_gs[0, 1], sharex=ax_m1, sharey=ax_m1)
    ax_m3 = fig.add_subplot(bot_gs[0, 2], sharex=ax_m1, sharey=ax_m1)

    x_ed = cell_edges_from_centers(g3d["x_ax"], g3d["lg_sym"])
    z_ed = cell_edges_from_centers(g3d["z_ax"], g3d["lg_sym"])
    m_3d_y0 = mask_slice_y0(g3d)
    m_2d_y0 = mask_slice_y0(g2d)
    m_both = m_3d_y0 & m_2d_y0

    te_3d_y0 = g3d["te"][:, 0, :] / 1e3
    te_2d_y0 = g2d["te"][:, 0, :] / 1e3
    rel_diff_te = np.where(m_both, (te_2d_y0 - te_3d_y0) / np.maximum(te_3d_y0, 1e-6) * 100.0, 0.0)

    vmin_t = min(te_3d_y0[m_3d_y0].min(), te_2d_y0[m_2d_y0].min())
    vmax_t = max(te_3d_y0[m_3d_y0].max(), te_2d_y0[m_2d_y0].max())

    im1 = _plot_pn150_panel_map(ax_m1, x_ed, z_ed, te_3d_y0, m_3d_y0, "inferno", vmin=vmin_t, vmax=vmax_t)
    im2 = _plot_pn150_panel_map(ax_m2, x_ed, z_ed, te_2d_y0, m_2d_y0, "inferno", vmin=vmin_t, vmax=vmax_t)
    norm_d = TwoSlopeNorm(vmin=-3.0, vcenter=0.0, vmax=3.0)
    im3 = _plot_pn150_panel_map(ax_m3, x_ed, z_ed, rel_diff_te, m_both, "coolwarm", norm=norm_d)

    ax_m1.set_title(r"(c) 3D ($y=0$ Slice): $T_e$ [$10^3\ \mathrm{K}$]", fontweight="bold")
    ax_m2.set_title(r"(d) 2D Cylinder: $T_e$ [$10^3\ \mathrm{K}$]", fontweight="bold")
    ax_m3.set_title(r"(e) Relative Diff: $(T_{e,\mathrm{2D}}-T_{e,\mathrm{3D}})/T_{e,\mathrm{3D}}$", fontweight="bold")
    ax_m1.set_ylabel(r"$Z$ [$10^{17}\ \mathrm{cm}$]")
    for axm in (ax_m1, ax_m2, ax_m3):
        axm.set_xlabel(r"$R$ or $X$ [$10^{17}\ \mathrm{cm}$]")

    fig.colorbar(im1, ax=[ax_m1, ax_m2], fraction=0.025, pad=0.02, label=r"$T_e$ [$10^3\ \mathrm{K}$]")
    fig.colorbar(im3, ax=ax_m3, fraction=0.048, pad=0.03, label=r"$\Delta T_e / T_{e,\mathrm{3D}}$ [%]")

    fig.suptitle(
        r"Tutorial Case 1 (Meudon PN150, Hdensity = 3000 $\mathrm{cm}^{-3}$): Temperature & Density Structure (2D vs 3D)",
        fontsize=11.2, fontweight="bold", y=0.995
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_pn150_ionic_fractions(models, out_path, r_in=1.0e17, r_out=4.07e17):
    """
    Generate Case 1 (PN150, `Hdensity 3000.`) Ionic Fractions comparison figure.
    """
    g3d = models["3d"]
    g2d = models["2d_cyl"]
    g2dx = models.get("2d_cyl_2x")

    m3d = active_gas_mask(g3d)
    m2d = active_gas_mask(g2d)
    m2dx = active_gas_mask(g2dx) if g2dx is not None else None
    bins_13 = np.linspace(r_in, r_out, 11)
    bins_25 = np.linspace(r_in, r_out, 17)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.6), dpi=300)
    plt.subplots_adjust(hspace=0.28, wspace=0.22)

    def plot_species_group(ax, species_specs, yscale="linear", ylim=(0, 1.05)):
        for label, el, st, color in species_specs:
            r3, v3 = radial_bin_profile(g3d["R_sph"], g3d["ions"][el][..., st], m3d, bins_13)
            r2, v2 = radial_bin_profile(g2d["R_sph"], g2d["ions"][el][..., st], m2d, bins_13)
            ax.plot(r3 / 1e17, v3, color=color, lw=2.0, ls="-", label=label)
            ax.plot(r2 / 1e17, v2, color=color, lw=1.2, ls="--", marker="o", ms=3.8, mfc="white", mew=1.1)
            if g2dx is not None:
                rx, vx = radial_bin_profile(g2dx["R_sph"], g2dx["ions"][el][..., st], m2dx, bins_25)
                ax.plot(rx / 1e17, vx, color=color, lw=1.1, ls=":", alpha=0.85)
        ax.set_yscale(yscale)
        ax.set_ylim(ylim)
        ax.set_xlim(r_in / 1e17 - 0.05, r_out / 1e17 + 0.05)
        ax.grid(alpha=0.22, ls="--", which="both")
        # Dummy handles for line-style legend key
        h_3d = plt.Line2D([0], [0], color="#444444", lw=2.0, ls="-", label="3D ($13^3$)")
        h_2d = plt.Line2D([0], [0], color="#444444", lw=1.2, ls="--", marker="o", ms=3.8, mfc="white", label="2D Cyl ($13\\times 13$)")
        handles, labels = ax.get_legend_handles_labels()
        if g2dx is not None:
            h_25 = plt.Line2D([0], [0], color="#444444", lw=1.2, ls=":", label="2D Cyl ($25\\times 25$)")
            ax.legend(handles + [h_3d, h_2d, h_25], labels + ["3D ($13^3$)", "2D Cyl ($13\\times 13$)", "2D Cyl ($25\\times 25$)"],
                      loc="best", frameon=True, framealpha=0.90, ncol=2)
        else:
            ax.legend(handles + [h_3d, h_2d], labels + ["3D ($13^3$)", "2D Cyl ($13\\times 13$)"],
                      loc="best", frameon=True, framealpha=0.90, ncol=2)

    # (a) H & He (log scale makes H0 trace transition visible alongside H+ and He ions)
    plot_species_group(
        axes[0, 0],
        [
            (r"$\mathrm{H}^+$", 1, 1, "#1b9e77"),
            (r"$\mathrm{H}^0$", 1, 0, "#d95f02"),
            (r"$\mathrm{He}^{2+}$", 2, 2, "#7570b3"),
            (r"$\mathrm{He}^+$", 2, 1, "#e7298a"),
            (r"$\mathrm{He}^0$", 2, 0, "#66a61e"),
        ],
        yscale="log", ylim=(1e-4, 1.4),
    )
    axes[0, 0].set_title(r"(a) Hydrogen & Helium Ionic Fractions $X(\mathrm{ion})$", fontweight="bold")
    axes[0, 0].set_ylabel(r"Ionic Fraction $n(\mathrm{X}^{i+}) / n(\mathrm{X})$")
    axes[0, 0].set_xlabel(r"Spherical Radius $r$ [$10^{17}\ \mathrm{cm}$]")

    # (b) Oxygen ions
    plot_species_group(
        axes[0, 1],
        [
            (r"$\mathrm{O}^{4+}$", 8, 4, "#e7298a"),
            (r"$\mathrm{O}^{3+}$", 8, 3, "#7570b3"),
            (r"$\mathrm{O}^{2+}$", 8, 2, "#0072b2"),
            (r"$\mathrm{O}^+$", 8, 1, "#d55e00"),
            (r"$\mathrm{O}^0$", 8, 0, "#111111"),
        ],
        yscale="log", ylim=(1e-4, 1.4),
    )
    axes[0, 1].set_title(r"(b) Oxygen Ionization Stratification ($\mathrm{O}^0 \to \mathrm{O}^{4+}$)", fontweight="bold")
    axes[0, 1].set_ylabel(r"Ionic Fraction $n(\mathrm{O}^{i+}) / n(\mathrm{O})$")
    axes[0, 1].set_xlabel(r"Spherical Radius $r$ [$10^{17}\ \mathrm{cm}$]")

    # (c) Nitrogen & Sulfur ions
    plot_species_group(
        axes[1, 0],
        [
            (r"$\mathrm{N}^{3+}$", 7, 3, "#7570b3"),
            (r"$\mathrm{N}^{2+}$", 7, 2, "#0072b2"),
            (r"$\mathrm{N}^+$", 7, 1, "#d55e00"),
            (r"$\mathrm{S}^{3+}$", 16, 3, "#66a61e"),
            (r"$\mathrm{S}^{2+}$", 16, 2, "#e69f00"),
            (r"$\mathrm{S}^+$", 16, 1, "#a65628"),
        ],
        yscale="log", ylim=(1e-4, 1.4),
    )
    axes[1, 0].set_title(r"(c) Nitrogen & Sulfur Ionic Fractions ($\mathrm{N}^{i+}$, $\mathrm{S}^{i+}$)", fontweight="bold")
    axes[1, 0].set_ylabel(r"Ionic Fraction $n(\mathrm{X}^{i+}) / n(\mathrm{X})$")
    axes[1, 0].set_xlabel(r"Spherical Radius $r$ [$10^{17}\ \mathrm{cm}$]")

    # (d) Side-by-side 2D meridional O2+ map inside subplot (3D above diagonal / 2D below, or split sub-axes)
    axes[1, 1].axis("off")
    sub = axes[1, 1].get_subplotspec().subgridspec(1, 2, wspace=0.08)
    ax_l = fig.add_subplot(sub[0, 0])
    ax_r = fig.add_subplot(sub[0, 1], sharex=ax_l, sharey=ax_l)

    x_ed = cell_edges_from_centers(g3d["x_ax"], g3d["lg_sym"])
    z_ed = cell_edges_from_centers(g3d["z_ax"], g3d["lg_sym"])
    m_3d_y0 = mask_slice_y0(g3d)
    map_src = g2dx if g2dx is not None else g2d
    x_ed_r = cell_edges_from_centers(map_src["x_ax"], map_src["lg_sym"])
    z_ed_r = cell_edges_from_centers(map_src["z_ax"], map_src["lg_sym"])
    m_2d_y0 = mask_slice_y0(map_src)

    im_o2_l = _plot_pn150_panel_map(ax_l, x_ed, z_ed, g3d["ions"][8][:, 0, :, 2], m_3d_y0, "viridis", vmin=0.0, vmax=1.0)
    im_o2_r = _plot_pn150_panel_map(ax_r, x_ed_r, z_ed_r, map_src["ions"][8][:, 0, :, 2], m_2d_y0, "viridis", vmin=0.0, vmax=1.0)
    # Overlay X(He2+) = 0.5 contour
    ax_l.contour(g3d["x_ax"] / 1e17, g3d["z_ax"] / 1e17,
                 np.where(m_3d_y0, g3d["ions"][2][:, 0, :, 2], 0.0).T,
                 levels=[0.5], colors=["white"], linewidths=1.3, linestyles="-")
    ax_r.contour(map_src["x_ax"] / 1e17, map_src["z_ax"] / 1e17,
                 np.where(m_2d_y0, map_src["ions"][2][:, 0, :, 2], 0.0).T,
                 levels=[0.5], colors=["white"], linewidths=1.3, linestyles="-")

    ax_l.set_title(r"(d) 3D ($y=0$): $X(\mathrm{O}^{2+})$", fontweight="bold")
    right_lbl = r"2D Cyl $25^2$: $X(\mathrm{O}^{2+})$" if g2dx is not None else r"2D Cyl $13^2$: $X(\mathrm{O}^{2+})$"
    ax_r.set_title(right_lbl, fontweight="bold")
    ax_l.set_ylabel(r"$Z$ [$10^{17}\ \mathrm{cm}$]")
    ax_l.set_xlabel(r"$X$ [$10^{17}\ \mathrm{cm}$]")
    ax_r.set_xlabel(r"$R$ [$10^{17}\ \mathrm{cm}$]")
    plt.setp(ax_r.get_yticklabels(), visible=False)
    fig.colorbar(im_o2_r, ax=[ax_l, ax_r], fraction=0.036, pad=0.025,
                 label=r"$X(\mathrm{O}^{2+})$ (white contour: $X(\mathrm{He}^{2+})=0.5$)")

    fig.suptitle(
        r"Tutorial Case 1 (Meudon PN150, Hdensity = 3000 $\mathrm{cm}^{-3}$): Ionic Fraction Stratification (2D vs 3D)",
        fontsize=11.2, fontweight="bold", y=0.985
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def _plot_test_map(ax, x_ed, z_ed, arr2d, mask2d, cmap, norm=None, vmin=None, vmax=None,
                   r_out=8.0e17, bad_color="#e6e6e6"):
    """Plot (R/X, Z) heatmap for test_Problem with torus and outer sphere overlays."""
    masked = np.ma.masked_where(~mask2d, arr2d)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(color=bad_color)
    kw = {"cmap": cm, "shading": "flat"}
    if norm is not None:
        kw["norm"] = norm
    else:
        kw["vmin"] = vmin
        kw["vmax"] = vmax
    mesh = ax.pcolormesh(x_ed / 1e17, z_ed / 1e17, masked.T, **kw)
    arc = Arc((0.0, 0.0), 2 * r_out / 1e17, 2 * r_out / 1e17, theta1=0, theta2=90,
              color="#222222", lw=0.8, ls="--", alpha=0.8)
    ax.add_patch(arc)
    # Torus outline: 1.0..3.8 x 1e17 in R (nodes 1.2..3.6), 0..2.2 x 1e17 in Z (nodes 0..2.0)
    ax.plot([1.0, 3.8, 3.8, 1.0, 1.0], [0.0, 0.0, 2.2, 2.2, 0.0],
            color="white", lw=0.95, ls=":", alpha=0.9)
    ax.set_aspect("equal")
    ax.set_xlim(x_ed[0] / 1e17, x_ed[-1] / 1e17)
    ax.set_ylim(z_ed[0] / 1e17, z_ed[-1] / 1e17)
    return mesh


def plot_test_problem_te_density(g3d, g2d, out_path):
    """
    Generate Case 2 (test_Problem, `densityFile` + `Ndust file`)
    nH, Ne, and Te structure comparison figure (2D Cylinder vs 3D).
    """
    fig = plt.figure(figsize=(12.6, 10.6), dpi=300)
    gs = GridSpec(3, 3, height_ratios=[1.0, 1.0, 1.08], hspace=0.33, wspace=0.25, figure=fig)

    x_ed = cell_edges_from_centers(g3d["x_ax"], g3d["lg_sym"])
    z_ed = cell_edges_from_centers(g3d["z_ax"], g3d["lg_sym"])
    m3 = mask_slice_y0(g3d)
    m2 = mask_slice_y0(g2d)
    mb = m3 & m2

    # Row 1: log10 nH (input density), 3D log10 Ne (y=0), 2D log10 Ne
    ax_nh = fig.add_subplot(gs[0, 0])
    ax_ne3 = fig.add_subplot(gs[0, 1], sharex=ax_nh, sharey=ax_nh)
    ax_ne2 = fig.add_subplot(gs[0, 2], sharex=ax_nh, sharey=ax_nh)

    log_nh = np.log10(np.maximum(g2d["hden"][:, 0, :], 1.0))
    log_ne3 = np.log10(np.maximum(g3d["ne"][:, 0, :], 0.1))
    log_ne2 = np.log10(np.maximum(g2d["ne"][:, 0, :], 0.1))

    im_nh = _plot_test_map(ax_nh, x_ed, z_ed, log_nh, m2, "cividis", vmin=2.8, vmax=5.0)
    im_ne3 = _plot_test_map(ax_ne3, x_ed, z_ed, log_ne3, m3, "magma", vmin=-0.5, vmax=4.5)
    im_ne2 = _plot_test_map(ax_ne2, x_ed, z_ed, log_ne2, m2, "magma", vmin=-0.5, vmax=4.5)

    ax_nh.text(2.4, 1.0, r"$8\times 10^4$" + "\n" + r"$\mathrm{cm}^{-3}$", color="black",
               ha="center", va="center", fontsize=7.8, fontweight="bold")
    ax_nh.text(5.7, 0.9, r"$10^4$", color="white", ha="center", va="center", fontsize=7.8, fontweight="bold")
    ax_nh.text(2.0, 5.2, r"$10^3\ \mathrm{cm}^{-3}$", color="white", ha="center", va="center", fontsize=7.8, fontweight="bold")

    ax_nh.set_title(r"(a) Predefined $\log_{10} n_{\mathrm{H}}$ (densityFile)", fontweight="bold")
    ax_ne3.set_title(r"(b) 3D ($y=0$): $\log_{10} N_e$ [$\mathrm{cm}^{-3}$]", fontweight="bold")
    ax_ne2.set_title(r"(c) 2D Cylinder: $\log_{10} N_e$ [$\mathrm{cm}^{-3}$]", fontweight="bold")
    ax_nh.set_ylabel(r"$Z$ [$10^{17}\ \mathrm{cm}$]")
    fig.colorbar(im_nh, ax=ax_nh, fraction=0.046, pad=0.03, label=r"$\log_{10}(n_{\mathrm{H}}/\mathrm{cm}^{-3})$")
    fig.colorbar(im_ne2, ax=[ax_ne3, ax_ne2], fraction=0.024, pad=0.02, label=r"$\log_{10}(N_e/\mathrm{cm}^{-3})$")

    # Row 2: 3D Te (y=0), 2D Cylinder Te, and Delta Te
    ax_t3 = fig.add_subplot(gs[1, 0], sharex=ax_nh, sharey=ax_nh)
    ax_t2 = fig.add_subplot(gs[1, 1], sharex=ax_nh, sharey=ax_nh)
    ax_td = fig.add_subplot(gs[1, 2], sharex=ax_nh, sharey=ax_nh)

    te3 = g3d["te"][:, 0, :] / 1e3
    te2 = g2d["te"][:, 0, :] / 1e3
    dte = np.where(mb, (g2d["te"][:, 0, :] - g3d["te"][:, 0, :]) / 1e3, 0.0)

    im_t3 = _plot_test_map(ax_t3, x_ed, z_ed, te3, m3, "inferno", vmin=0.0, vmax=30.0)
    im_t2 = _plot_test_map(ax_t2, x_ed, z_ed, te2, m2, "inferno", vmin=0.0, vmax=30.0)
    im_td = _plot_test_map(ax_td, x_ed, z_ed, dte, mb, "coolwarm",
                           norm=TwoSlopeNorm(vmin=-1.5, vcenter=0.0, vmax=1.5))

    ax_t3.set_title(r"(d) 3D ($y=0$): $T_e$ [$10^3\ \mathrm{K}$]", fontweight="bold")
    ax_t2.set_title(r"(e) 2D Cylinder: $T_e$ [$10^3\ \mathrm{K}$]", fontweight="bold")
    ax_td.set_title(r"(f) Difference: $T_{e,\mathrm{2D}} - T_{e,\mathrm{3D},y=0}$", fontweight="bold")
    ax_t3.set_ylabel(r"$Z$ [$10^{17}\ \mathrm{cm}$]")
    for axm in (ax_t3, ax_t2, ax_td):
        axm.set_xlabel(r"$R$ or $X$ [$10^{17}\ \mathrm{cm}$]")
    fig.colorbar(im_t2, ax=[ax_t3, ax_t2], fraction=0.024, pad=0.02, label=r"$T_e$ [$10^3\ \mathrm{K}$]")
    fig.colorbar(im_td, ax=ax_td, fraction=0.046, pad=0.03, label=r"$\Delta T_e$ [$10^3\ \mathrm{K}$]")

    # Row 3: Two wider directional profile plots spanning 1.5 columns each via nested gridspec
    row3_gs = gs[2, :].subgridspec(1, 2, wspace=0.22)
    ax_pe = fig.add_subplot(row3_gs[0, 0])
    ax_pt = fig.add_subplot(row3_gs[0, 1])

    x_ax = g3d["x_ax"] / 1e17
    z_ax = g3d["z_ax"] / 1e17
    k_rim = int(np.argmin(np.abs(g3d["z_ax"] - 2.4e17)))

    # Equatorial (Z=0, exclude origin cavity cell i=0 where active==0) and Polar (R=0, k>=1) masks
    eq_3d_m = (g3d["active"][:, 0, 0] > 0) & (g3d["hden"][:, 0, 0] > 0)
    eq_2d_m = (g2d["active"][:, 0, 0] > 0) & (g2d["hden"][:, 0, 0] > 0)
    pol_3d_m = (g3d["active"][0, 0, :] > 0) & (g3d["hden"][0, 0, :] > 0)
    pol_2d_m = (g2d["active"][0, 0, :] > 0) & (g2d["hden"][0, 0, :] > 0)
    rim_3d_m = (g3d["active"][:, 0, k_rim] > 0) & (g3d["hden"][:, 0, k_rim] > 0)
    rim_2d_m = (g2d["active"][:, 0, k_rim] > 0) & (g2d["hden"][:, 0, k_rim] > 0)

    # (g) Equatorial & Polar nH / Ne profiles
    ax_pe.plot(x_ax[eq_3d_m], np.log10(g3d["hden"][eq_3d_m, 0, 0]), color="#666666", lw=1.4, ls="-",
               label=r"$n_{\mathrm{H}}$ Equator ($Z=0$)")
    ax_pe.plot(x_ax[eq_3d_m], np.log10(np.maximum(g3d["ne"][eq_3d_m, 0, 0], 0.1)), color=C_3D, lw=2.0,
               label=r"3D $N_e$ Equator ($Z=0$)")
    ax_pe.plot(x_ax[eq_2d_m], np.log10(np.maximum(g2d["ne"][eq_2d_m, 0, 0], 0.1)), color=C_CYL, lw=1.6, ls="--",
               marker="o", ms=3.8, mfc="white", label=r"2D Cyl $N_e$ Equator ($Z=0$)")
    ax_pe.plot(z_ax[pol_3d_m], np.log10(np.maximum(g3d["ne"][0, 0, pol_3d_m], 0.1)), color=C_RES2X, lw=2.0,
               label=r"3D $N_e$ Polar ($R=0$)")
    ax_pe.plot(z_ax[pol_2d_m], np.log10(np.maximum(g2d["ne"][0, 0, pol_2d_m], 0.1)), color="#009e73", lw=1.5, ls=":",
               marker="^", ms=3.8, mfc="white", label=r"2D Cyl $N_e$ Polar ($R=0$)")
    ax_pe.set_title(r"(g) Equatorial ($Z=0$) & Polar ($R=0$) Density Cuts", fontweight="bold")
    ax_pe.set_xlabel(r"Distance Along Cut ($R$ or $Z$) [$10^{17}\ \mathrm{cm}$]")
    ax_pe.set_ylabel(r"$\log_{10}(\mathrm{Density} / \mathrm{cm}^{-3})$")
    ax_pe.set_ylim(0.2, 5.2)
    ax_pe.grid(alpha=0.22, ls="--")
    ax_pe.legend(loc="upper right", frameon=True, framealpha=0.92, ncol=1)

    # (h) Te directional cuts: Polar (R=0), Rim (Z=2.4e17), Equatorial (Z=0)
    ax_pt.plot(z_ax[pol_3d_m], g3d["te"][0, 0, pol_3d_m] / 1e3, color=C_RES2X, lw=2.0,
               label=r"3D Polar ($R=0$)")
    ax_pt.plot(z_ax[pol_2d_m], g2d["te"][0, 0, pol_2d_m] / 1e3, color="#56b4e9", lw=1.5, ls="--",
               marker="^", ms=3.8, mfc="white", label=r"2D Cyl Polar ($R=0$)")

    ax_pt.plot(x_ax[rim_3d_m], g3d["te"][rim_3d_m, 0, k_rim] / 1e3, color=C_SPEC, lw=2.0,
               label=r"3D Above-Torus ($Z=2.4\times 10^{17}$)")
    ax_pt.plot(x_ax[rim_2d_m], g2d["te"][rim_2d_m, 0, k_rim] / 1e3, color="#a6d854", lw=1.5, ls="--",
               marker="s", ms=3.5, mfc="white", label=r"2D Cyl Above-Torus ($Z=2.4\times 10^{17}$)")

    ax_pt.plot(x_ax[eq_3d_m], g3d["te"][eq_3d_m, 0, 0] / 1e3, color=C_3D, lw=2.0,
               label=r"3D Equator ($Z=0$)")
    ax_pt.plot(x_ax[eq_2d_m], g2d["te"][eq_2d_m, 0, 0] / 1e3, color=C_CYL, lw=1.5, ls="--",
               marker="o", ms=3.8, mfc="white", label=r"2D Cyl Equator ($Z=0$)")

    ax_pt.set_title(r"(h) Electron Temperature Cuts Across Torus, Rim & Polar Axis", fontweight="bold")
    ax_pt.set_xlabel(r"Distance Along Cut ($R$ or $Z$) [$10^{17}\ \mathrm{cm}$]")
    ax_pt.set_ylabel(r"Electron Temperature $T_e$ [$10^3\ \mathrm{K}$]")
    ax_pt.set_ylim(-1.0, 31.5)
    ax_pt.grid(alpha=0.22, ls="--")
    ax_pt.legend(loc="upper right", frameon=True, framealpha=0.92, ncol=2)

    fig.suptitle(
        r"Tutorial Case 2 (test_Problem, densityFile + Ndust file): Density & Temperature Structure (2D vs 3D)",
        fontsize=11.4, fontweight="bold", y=0.992
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_test_problem_ionic_fractions(g3d, g2d, out_path):
    """
    Generate Case 2 (test_Problem, `densityFile` + `Ndust file`)
    Ionic Fractions comparison figure (2D Cylinder vs 3D).
    """
    fig = plt.figure(figsize=(12.6, 8.6), dpi=300)
    outer = GridSpec(2, 1, height_ratios=[1.12, 0.92], hspace=0.32, figure=fig)

    top = outer[0].subgridspec(1, 3, wspace=0.22)
    ax_eq = fig.add_subplot(top[0, 0])
    ax_rim = fig.add_subplot(top[0, 1])
    ax_pol = fig.add_subplot(top[0, 2])

    x_ax = g3d["x_ax"] / 1e17
    z_ax = g3d["z_ax"] / 1e17
    k_rim = int(np.argmin(np.abs(g3d["z_ax"] - 2.4e17)))

    eq_3m = (g3d["active"][:, 0, 0] > 0) & (g3d["hden"][:, 0, 0] > 0)
    eq_2m = (g2d["active"][:, 0, 0] > 0) & (g2d["hden"][:, 0, 0] > 0)
    rim_3m = (g3d["active"][:, 0, k_rim] > 0) & (g3d["hden"][:, 0, k_rim] > 0)
    rim_2m = (g2d["active"][:, 0, k_rim] > 0) & (g2d["hden"][:, 0, k_rim] > 0)
    pol_3m = (g3d["active"][0, 0, :] > 0) & (g3d["hden"][0, 0, :] > 0)
    pol_2m = (g2d["active"][0, 0, :] > 0) & (g2d["hden"][0, 0, :] > 0)

    def plot_cut(ax, coord, m3, m2, slice_fn_3, slice_fn_2, specs, legend_loc="best", bbox=None):
        for lbl, el, st, col in specs:
            y3 = slice_fn_3(g3d["ions"][el])[m3, st]
            y2 = slice_fn_2(g2d["ions"][el])[m2, st]
            ax.plot(coord[m3], y3, color=col, lw=2.0, ls="-", label=lbl)
            ax.plot(coord[m2], y2, color=col, lw=1.2, ls="--", marker="o", ms=3.5, mfc="white", mew=1.0)
        h3 = plt.Line2D([0], [0], color="#444444", lw=2.0, ls="-", label="3D ($y=0$)")
        h2 = plt.Line2D([0], [0], color="#444444", lw=1.2, ls="--", marker="o", ms=3.5, mfc="white", label="2D Cyl")
        handles, labels = ax.get_legend_handles_labels()
        leg_kw = dict(loc=legend_loc, frameon=True, framealpha=0.92, ncol=2)
        if bbox is not None:
            leg_kw["bbox_to_anchor"] = bbox
        ax.legend(handles + [h3, h2], labels + ["3D ($y=0$)", "2D Cyl"], **leg_kw)
        ax.set_ylim(-0.04, 1.06)
        ax.grid(alpha=0.22, ls="--")

    # (a) Equatorial Torus ionization front (Z = 0, zoom R <= 4e17 where front occurs)
    plot_cut(
        ax_eq, x_ax, eq_3m, eq_2m,
        lambda arr: arr[:, 0, 0], lambda arr: arr[:, 0, 0],
        [
            (r"$\mathrm{H}^+$", 1, 1, "#1b9e77"),
            (r"$\mathrm{He}^+$", 2, 1, "#e7298a"),
            (r"$\mathrm{O}^+$", 8, 1, "#d55e00"),
            (r"$\mathrm{O}^0$", 8, 0, "#111111"),
            (r"$\mathrm{N}^+$", 7, 1, "#0072b2"),
            (r"$\mathrm{S}^+$", 16, 1, "#a65628"),
        ],
        legend_loc="center right",
    )
    ax_eq.set_xlim(0.2, 4.0)
    ax_eq.axvspan(1.0, 4.0, color="#fff2cc", alpha=0.35, zorder=0)
    ax_eq.set_title(r"(a) Equatorial Torus Front ($Z=0$)", fontweight="bold")
    ax_eq.set_xlabel(r"$R$ or $X$ [$10^{17}\ \mathrm{cm}$]")
    ax_eq.set_ylabel(r"Ionic Fraction $n(\mathrm{X}^{i+}) / n(\mathrm{X})$")

    # (b) Above-Torus Rim / Shadow Cut (Z = 2.4e17)
    plot_cut(
        ax_rim, x_ax, rim_3m, rim_2m,
        lambda arr: arr[:, 0, k_rim], lambda arr: arr[:, 0, k_rim],
        [
            (r"$\mathrm{H}^+$", 1, 1, "#1b9e77"),
            (r"$\mathrm{O}^{3+}$", 8, 3, "#7570b3"),
            (r"$\mathrm{O}^{2+}$", 8, 2, "#0072b2"),
            (r"$\mathrm{O}^+$", 8, 1, "#d55e00"),
            (r"$\mathrm{N}^+$", 7, 1, "#e69f00"),
            (r"$\mathrm{S}^+$", 16, 1, "#a65628"),
        ],
        legend_loc="center right", bbox=(0.98, 0.50),
    )
    ax_rim.set_xlim(0.0, 5.2)
    ax_rim.set_title(r"(b) Torus Shadow Rim ($Z=2.4\times 10^{17}\ \mathrm{cm}$)", fontweight="bold")
    ax_rim.set_xlabel(r"$R$ or $X$ [$10^{17}\ \mathrm{cm}$]")

    # (c) Polar Axis Lobe (R = 0, Z = 0.4..8.0 x 1e17)
    plot_cut(
        ax_pol, z_ax, pol_3m, pol_2m,
        lambda arr: arr[0, 0, :], lambda arr: arr[0, 0, :],
        [
            (r"$\mathrm{He}^{2+}$", 2, 2, "#7570b3"),
            (r"$\mathrm{He}^+$", 2, 1, "#e7298a"),
            (r"$\mathrm{O}^{4+}$", 8, 4, "#e69f00"),
            (r"$\mathrm{O}^{3+}$", 8, 3, "#111111"),
            (r"$\mathrm{O}^{2+}$", 8, 2, "#0072b2"),
            (r"$\mathrm{N}^{2+}$", 7, 2, "#1b9e77"),
        ],
        legend_loc="center left", bbox=(0.02, 0.63),
    )
    ax_pol.set_xlim(0.2, 8.1)
    ax_pol.set_title(r"(c) Polar Axis High-Ionization Lobe ($R=0$)", fontweight="bold")
    ax_pol.set_xlabel(r"$Z$ [$10^{17}\ \mathrm{cm}$]")

    # Bottom Row: 4 Spatial Maps (H+ in 3D/2D and O2+ in 3D/2D) with dedicated cax & spacer column
    bot = outer[1].subgridspec(
        1, 7,
        width_ratios=[1.0, 1.0, 0.055, 0.38, 1.0, 1.0, 0.055],
        wspace=0.08,
    )
    ax_h3 = fig.add_subplot(bot[0, 0])
    ax_h2 = fig.add_subplot(bot[0, 1], sharex=ax_h3, sharey=ax_h3)
    cax_h = fig.add_subplot(bot[0, 2])
    ax_o3 = fig.add_subplot(bot[0, 4], sharex=ax_h3, sharey=ax_h3)
    ax_o2 = fig.add_subplot(bot[0, 5], sharex=ax_h3, sharey=ax_h3)
    cax_o = fig.add_subplot(bot[0, 6])

    x_ed = cell_edges_from_centers(g3d["x_ax"], g3d["lg_sym"])
    z_ed = cell_edges_from_centers(g3d["z_ax"], g3d["lg_sym"])
    m3 = mask_slice_y0(g3d)
    m2 = mask_slice_y0(g2d)

    im_h3 = _plot_test_map(ax_h3, x_ed, z_ed, g3d["ions"][1][:, 0, :, 1], m3, "viridis", vmin=0.0, vmax=1.0)
    im_h2 = _plot_test_map(ax_h2, x_ed, z_ed, g2d["ions"][1][:, 0, :, 1], m2, "viridis", vmin=0.0, vmax=1.0)
    im_o3 = _plot_test_map(ax_o3, x_ed, z_ed, g3d["ions"][8][:, 0, :, 2], m3, "plasma", vmin=0.0, vmax=1.0)
    im_o2 = _plot_test_map(ax_o2, x_ed, z_ed, g2d["ions"][8][:, 0, :, 2], m2, "plasma", vmin=0.0, vmax=1.0)

    ax_h3.set_title(r"(d) 3D ($y=0$): $X(\mathrm{H}^+)$", fontweight="bold")
    ax_h2.set_title(r"(e) 2D Cyl: $X(\mathrm{H}^+)$", fontweight="bold")
    ax_o3.set_title(r"(f) 3D ($y=0$): $X(\mathrm{O}^{2+})$", fontweight="bold")
    ax_o2.set_title(r"(g) 2D Cyl: $X(\mathrm{O}^{2+})$", fontweight="bold")

    ax_h3.set_ylabel(r"$Z$ [$10^{17}\ \mathrm{cm}$]")
    plt.setp(ax_h2.get_yticklabels(), visible=False)
    plt.setp(ax_o3.get_yticklabels(), visible=False)
    plt.setp(ax_o2.get_yticklabels(), visible=False)
    for axm in (ax_h3, ax_h2, ax_o3, ax_o2):
        axm.set_xlabel(r"$R$ or $X$ [$10^{17}\ \mathrm{cm}$]")

    fig.colorbar(im_h2, cax=cax_h, label=r"$X(\mathrm{H}^+)$")
    fig.colorbar(im_o2, cax=cax_o, label=r"$X(\mathrm{O}^{2+})$")

    fig.suptitle(
        r"Tutorial Case 2 (test_Problem, densityFile + Ndust file): Ionic Fractions (2D vs 3D)",
        fontsize=11.4, fontweight="bold", y=0.992
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def build_arg_parser():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(
        description="Plot 2D-vs-3D MOCASSIN temperature, density, and ionic fraction comparisons for TUTORIAL_2D_MODE.md."
    )
    parser.add_argument("--case", choices=["all", "pn150", "test_problem"], default="all",
                        help="Which tutorial case(s) to plot (default: all).")
    parser.add_argument("--out-dir", default=os.path.join(repo_root, "docs", "tutorial_2d_figures"),
                        help="Directory where PNG comparison plots will be written.")
    parser.add_argument("--pn150-3d", default=os.path.join(repo_root, "tutorial_runs", "PN150_3D"),
                        help="Path to PN150 3D model directory (containing input/ and output/).")
    parser.add_argument("--pn150-2d", default=os.path.join(repo_root, "tutorial_runs", "PN150_2D_cylinder"),
                        help="Path to PN150 2D XZ Z cylinder model directory.")
    parser.add_argument("--pn150-2d-res2x", default=None,
                        help="Optional path to PN150 2D XZ Z cylinder 25x25 double-resolution run.")
    parser.add_argument("--pn150-2d-spec", default=None,
                        help="Optional path to PN150 2D XZ Z specular run.")
    parser.add_argument("--pn150-2d-buffer", default=None,
                        help="Optional path to PN150 2D XZ Z buffer run.")
    parser.add_argument("--test-3d", default=os.path.join(repo_root, "tutorial_runs", "test_Problem_3D"),
                        help="Path to test_Problem 3D model directory.")
    parser.add_argument("--test-2d", default=os.path.join(repo_root, "tutorial_runs", "test_Problem_2D_cylinder"),
                        help="Path to test_Problem 2D XZ Z cylinder model directory.")
    return parser


def main():
    args = build_arg_parser().parse_args()
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(args.out_dir, exist_ok=True)

    if args.case in ("all", "pn150"):
        pn_models = {
            "3d": load_spatial_grid(args.pn150_3d, repo_root),
            "2d_cyl": load_spatial_grid(args.pn150_2d, repo_root),
        }
        if args.pn150_2d_res2x and os.path.isdir(args.pn150_2d_res2x):
            pn_models["2d_cyl_2x"] = load_spatial_grid(args.pn150_2d_res2x, repo_root)
        if args.pn150_2d_spec and os.path.isdir(args.pn150_2d_spec):
            pn_models["2d_spec"] = load_spatial_grid(args.pn150_2d_spec, repo_root)
        if args.pn150_2d_buffer and os.path.isdir(args.pn150_2d_buffer):
            pn_models["2d_buf"] = load_spatial_grid(args.pn150_2d_buffer, repo_root)

        f1 = os.path.join(args.out_dir, "pn150_te_density_2d_vs_3d.png")
        f2 = os.path.join(args.out_dir, "pn150_ionic_fractions_2d_vs_3d.png")
        plot_pn150_te_density(pn_models, f1)
        plot_pn150_ionic_fractions(pn_models, f2)
        print(f"Wrote {f1}")
        print(f"Wrote {f2}")

    if args.case in ("all", "test_problem"):
        t3d = load_spatial_grid(args.test_3d, repo_root)
        t2d = load_spatial_grid(args.test_2d, repo_root)
        f3 = os.path.join(args.out_dir, "test_problem_te_density_2d_vs_3d.png")
        f4 = os.path.join(args.out_dir, "test_problem_ionic_fractions_2d_vs_3d.png")
        plot_test_problem_te_density(t3d, t2d, f3)
        plot_test_problem_ionic_fractions(t3d, t2d, f4)
        print(f"Wrote {f3}")
        print(f"Wrote {f4}")


if __name__ == "__main__":
    main()
