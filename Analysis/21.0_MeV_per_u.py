# ============================================================
# CR-39 TWO-VECTOR SCATTERING PIPELINE
# DOT-PRODUCT ARCCOS METHOD
# SINGLE-ENERGY VERSION: 21 MeV/n
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator

import os
import time
from pathlib import Path

# ============================================================
# USER PARAMETERS (ONLY 21 MeV/n)
# ============================================================

energy_label = "21 MeV/n"

cfg = {
    "front_A": r"C:\Users\mahfuj\Learning Code\21 MeV\N212_Front.csv",
    "back_A":  r"C:\Users\mahfuj\Learning Code\21 MeV\N212_Back.csv",
    "front_B": r"C:\Users\mahfuj\Learning Code\21 MeV\N213_Front.csv",
    "back_B":  r"C:\Users\mahfuj\Learning Code\21 MeV\N213_Back.csv",

    "tag_A": "N212",
    "tag_B": "N213",

    "thick_A_post": 0.343,
    "thick_B_post": 0.354,

    "out_folder": (
        r"C:\Users\mahfuj\Learning Code\21 MeV"
        r"\CR39_report_opt"
    )
}
BASE_DIR = Path(__file__).resolve().parent.parent

cfg = {
    "front_A": BASE_DIR / "data/processed/21.0_MeV_per_u/N212_Front.csv",
    "back_A": BASE_DIR / "data/processed/21.0_MeV_per_u/N212_Back.csv",
    "front_B": BASE_DIR / "data/processed/21.0_MeV_per_u/N213_Front.csv",
    "back_B": BASE_DIR / "data/processed/21.0_MeV_per_u/N213_Back.csv",

    "tag_A": "N212",
    "tag_B": "N213",

    "thick_A_post": 0.384,
    "thick_B_post": 0.387,

    "out_folder": BASE_DIR / "results/21.0_MeV_per_u"
}

al_thickness_mm = 0.10
coords_unit = "micron"

max_kde_samples = 20000
alignment_radius_frac = 0.25

cross_match_radius_px = None
use_mutual_nn = True

rng = np.random.default_rng(12345)


# ============================================================
# OUTPUTS
# ============================================================

out_folder = cfg["out_folder"]

os.makedirs(
    out_folder,
    exist_ok=True
)

tag_A = cfg["tag_A"]
tag_B = cfg["tag_B"]


out_csv = os.path.join(
    out_folder,
    f"{tag_A}_{tag_B}_angles_dotrule.csv"
)


out_pdf = os.path.join(
    out_folder,
    f"{tag_A}_{tag_B}_report_dotrule.pdf"
)


out_plot_compare = os.path.join(
    out_folder,
    "Energy_vs_mu_vs_theory.png"
)


# ============================================================
# HELPERS
# ============================================================

def parse_energy_MeV_u(label):

    try:

        return float(
            label.split()[0]
        )

    except Exception:

        return np.nan


def load_csv(path):

    try:

        df = pd.read_csv(
            path
        )

    except Exception:

        df = pd.read_csv(
            path,
            sep=None,
            engine="python"
        )


    df = df.dropna(
        subset=[
            "GravX",
            "GravY"
        ]
    ).reset_index(
        names="row_id"
    )


    print(
        f"Loaded "
        f"{os.path.basename(path)} "
        f"-> {len(df)} rows"
    )


    return df


# ============================================================
# KDE ALIGNMENT
# ============================================================

def kd_alignment_sample(
    front,
    back,
    tag="DET",
    fig_idx_start=1,
    out_folder="."
):

    fpts = front[
        ["GravX", "GravY"]
    ].to_numpy()


    bpts = back[
        ["GravX", "GravY"]
    ].to_numpy()


    tree = cKDTree(
        bpts
    )


    span_x = max(
        np.ptp(
            fpts[:, 0]
        ),
        np.ptp(
            bpts[:, 0]
        )
    )


    span_y = max(
        np.ptp(
            fpts[:, 1]
        ),
        np.ptp(
            bpts[:, 1]
        )
    )


    search_radius = (
        max(
            span_x,
            span_y
        )
        * alignment_radius_frac
    )


    dists, idxs = tree.query(
        fpts,
        k=1,
        distance_upper_bound=search_radius
    )


    mask = np.isfinite(
        dists
    )


    dx = (
        bpts[
            idxs[mask],
            0
        ]
        - fpts[
            mask,
            0
        ]
    )


    dy = (
        bpts[
            idxs[mask],
            1
        ]
        - fpts[
            mask,
            1
        ]
    )


    # --------------------------------------------------------
    # SAMPLE FOR KDE IF NECESSARY
    # --------------------------------------------------------

    if len(dx) > max_kde_samples:

        sel = rng.choice(
            len(dx),
            size=max_kde_samples,
            replace=False
        )

        dx = dx[sel]
        dy = dy[sel]


    # --------------------------------------------------------
    # KDE
    # --------------------------------------------------------

    kde = gaussian_kde(
        np.vstack(
            [
                dx,
                dy
            ]
        )
    )


    cx = np.median(
        dx
    )

    cy = np.median(
        dy
    )


    sx = np.std(
        dx
    )

    sy = np.std(
        dy
    )


    gx = np.linspace(
        cx - 6 * sx,
        cx + 6 * sx,
        200
    )


    gy = np.linspace(
        cy - 6 * sy,
        cy + 6 * sy,
        200
    )


    X, Y = np.meshgrid(
        gx,
        gy
    )


    Z = kde(
        np.vstack(
            [
                X.ravel(),
                Y.ravel()
            ]
        )
    ).reshape(
        200,
        200
    )


    iy, ix = np.unravel_index(
        np.argmax(Z),
        Z.shape
    )


    peak_dx = gx[ix]
    peak_dy = gy[iy]


    sigma_x = np.std(
        dx
    )

    sigma_y = np.std(
        dy
    )


    sigma_mean = np.mean(
        [
            sigma_x,
            sigma_y
        ]
    )


    # ========================================================
    # KDE PLOT
    # ========================================================

    plt.figure(
        figsize=(6, 5)
    )


    plt.imshow(
        Z,
        extent=[
            gx.min(),
            gx.max(),
            gy.min(),
            gy.max()
        ],
        origin="lower",
        aspect="auto"
    )


    plt.scatter(
        peak_dx,
        peak_dy,
        color="red",
        label="peak"
    )


    plt.title(
        f"{tag} dx–dy KDE (peak)"
    )

    plt.xlabel("dx")
    plt.ylabel("dy")

    plt.legend()

    plt.colorbar(
        label="density"
    )


    png_path = os.path.join(
        out_folder,
        f"{fig_idx_start:02d}_{tag}_KDE.png"
    )


    plt.savefig(
        png_path,
        dpi=300
    )


    plt.show()

    plt.close()


    # ========================================================
    # dx/dy HISTOGRAMS
    # ========================================================

    plt.figure(
        figsize=(12, 4)
    )


    plt.subplot(
        1,
        2,
        1
    )


    plt.hist(
        dx,
        bins=120,
        histtype="step"
    )


    plt.axvline(
        peak_dx,
        color="r"
    )


    plt.title(
        f"{tag} dx "
        f"(μ≈{peak_dx:.2f}, "
        f"σ={sigma_x:.2f})"
    )


    plt.subplot(
        1,
        2,
        2
    )


    plt.hist(
        dy,
        bins=120,
        histtype="step"
    )


    plt.axvline(
        peak_dy,
        color="r"
    )


    plt.title(
        f"{tag} dy "
        f"(μ≈{peak_dy:.2f}, "
        f"σ={sigma_y:.2f})"
    )


    png_path = os.path.join(
        out_folder,
        f"{fig_idx_start + 1:02d}_{tag}_Hist.png"
    )


    plt.savefig(
        png_path,
        dpi=300
    )


    plt.show()

    plt.close()


    return (
        peak_dx,
        peak_dy,
        sigma_mean,
        dx,
        dy
    )


# ============================================================
# FRONT-BACK MATCHING
# ============================================================

def match_unique_front_back(
    df_f,
    df_b,
    radius
):

    tree = cKDTree(
        df_b[
            ["GravX", "GravY"]
        ].to_numpy()
    )


    pts_f = df_f[
        ["GravX", "GravY"]
    ].to_numpy()


    dists, idxs = tree.query(
        pts_f,
        k=2,
        distance_upper_bound=radius
    )


    matches = []

    used_b = set()


    for i, (
        dlist,
        idlist
    ) in enumerate(
        zip(
            dists,
            idxs
        )
    ):

        if (
            dlist[0] != np.inf
            and idlist[0] not in used_b
        ):

            if (
                dlist[1] != np.inf
                and idlist[1] != idlist[0]
            ):

                continue


            f = df_f.iloc[i]

            b = df_b.iloc[
                idlist[0]
            ]


            matches.append(
                {
                    "Fx": f["GravX"],
                    "Fy": f["GravY"],
                    "Bx": b["GravX"],
                    "By": b["GravY"]
                }
            )


            used_b.add(
                idlist[0]
            )


    return pd.DataFrame(
        matches
    )


# ============================================================
# GAUSSIAN
# ============================================================

def gaussian(
    x,
    A,
    mu,
    sigma,
    C
):

    return (
        A
        * np.exp(
            -0.5
            * (
                (
                    x - mu
                )
                / sigma
            ) ** 2
        )
        + C
    )


# ============================================================
# THEORY FUNCTIONS
# ============================================================

def beta_p_12C(
    En_MeV_u
):

    z = 6
    A = 12

    m_u = 931.5

    mc2 = A * m_u

    T = (
        En_MeV_u
        * A
    )

    E_tot = (
        mc2
        + T
    )

    p = np.sqrt(
        E_tot ** 2
        - mc2 ** 2
    )

    beta = (
        p
        / E_tot
    )

    return (
        beta,
        p
    )


def highland_theta0_deg(
    En_MeV_u
):

    z = 6

    x_cm = 0.01
    X0_cm = 8.9

    beta, p = beta_p_12C(
        En_MeV_u
    )

    term = (
        13.6
        / (
            beta
            * p
        )
        * z
        * np.sqrt(
            x_cm
            / X0_cm
        )
        * (
            1
            + 0.038
            * np.log(
                x_cm
                / X0_cm
            )
        )
    )

    return np.degrees(
        term
    )


def lynch_dahl_theta0_deg(
    En_MeV_u
):

    z = 6

    x_cm = 0.01
    X0_cm = 8.9

    beta, p = beta_p_12C(
        En_MeV_u
    )

    arg = (
        x_cm
        * z ** 2
    ) / (
        X0_cm
        * beta ** 2
    )

    term = (
        13.6
        / (
            beta
            * p
        )
        * z
        * np.sqrt(
            x_cm
            / X0_cm
        )
        * (
            1
            + 0.038
            * np.log(
                arg
            )
        )
    )

    return np.degrees(
        term
    )


def highland_theta0_deg_carbon(
    En_MeV_u,
    al_thickness_mm
):

    return highland_theta0_deg(
        En_MeV_u
    )


# ============================================================
# RUN
# ============================================================

start_total = time.time()


print(
    "=== CR-39 scattering "
    "(two-vector arccos) "
    "SINGLE-ENERGY ==="
)


print(
    f"====================== "
    f"{energy_label} "
    f"======================"
)


results_summary = []


# ============================================================
# GEOMETRY
# ============================================================

mm2u = (
    1000.0
    if "micro"
    in coords_unit.lower()
    else 1.0
)


dz_A = (
    cfg["thick_A_post"]
    * mm2u
)


dz_B = (
    cfg["thick_B_post"]
    * mm2u
)


sep_front = (
    cfg["thick_A_post"]
    + al_thickness_mm
) * mm2u


zA_f = 0.0

zA_b = dz_A

zB_f = sep_front

zB_b = (
    sep_front
    + dz_B
)


# ============================================================
# LOAD CSV FILES
# ============================================================

dfA_f, dfA_b, dfB_f, dfB_b = [

    load_csv(
        p
    )

    for p in [

        cfg["front_A"],

        cfg["back_A"],

        cfg["front_B"],

        cfg["back_B"]

    ]

]


# ============================================================
# PDF REPORT
# ============================================================

with PdfPages(
    out_pdf
) as pdf_one:


    # ========================================================
    # ALIGNMENT
    # ========================================================

    print(
        f"\n-- Alignment sampling "
        f"{tag_A} --"
    )


    (
        peak_dxA,
        peak_dyA,
        sigA,
        dxA_full,
        dyA_full
    ) = kd_alignment_sample(
        dfA_f,
        dfA_b,
        tag_A,
        1,
        out_folder
    )


    print(
        f"\n-- Alignment sampling "
        f"{tag_B} --"
    )


    (
        peak_dxB,
        peak_dyB,
        sigB,
        dxB_full,
        dyB_full
    ) = kd_alignment_sample(
        dfB_f,
        dfB_b,
        tag_B,
        3,
        out_folder
    )


    # ========================================================
    # APPLY ALIGNMENT SHIFTS
    # ========================================================

    dfA_b["GravX"] -= peak_dxA

    dfA_b["GravY"] -= peak_dyA


    dfB_b["GravX"] -= peak_dxB

    dfB_b["GravY"] -= peak_dyB


    # ========================================================
    # FRONT-BACK MATCHING
    # ========================================================

    rA = (
        4
        * sigA
    )

    rB = (
        4
        * sigB
    )


    print(
        f"\nFront–back radii: "
        f"rA={rA:.2f}, "
        f"rB={rB:.2f}"
    )


    A_pairs = match_unique_front_back(
        dfA_f,
        dfA_b,
        rA
    )


    B_pairs = match_unique_front_back(
        dfB_f,
        dfB_b,
        rB
    )


    print(
        f"{tag_A} pairs={len(A_pairs)}, "
        f"{tag_B} pairs={len(B_pairs)}"
    )


    # ========================================================
    # TRACK VECTORS
    # ========================================================

    Av = np.vstack(
        [
            [
                r.Bx - r.Fx,
                r.By - r.Fy,
                dz_A
            ]

            for r in A_pairs.itertuples()
        ]
    )


    Bv = np.vstack(
        [
            [
                r.Bx - r.Fx,
                r.By - r.Fy,
                dz_B
            ]

            for r in B_pairs.itertuples()
        ]
    )


    Ahat = (
        Av
        / np.linalg.norm(
            Av,
            axis=1,
            keepdims=True
        )
    )


    Bhat = (
        Bv
        / np.linalg.norm(
            Bv,
            axis=1,
            keepdims=True
        )
    )


    # ========================================================
    # CROSS-MATCHING
    # ========================================================

    dz_AB = (
        zB_f
        - zA_b
    )


    x_proj = (
        A_pairs["Bx"]
        + (
            A_pairs["Bx"]
            - A_pairs["Fx"]
        )
        / dz_A
        * dz_AB
    )


    y_proj = (
        A_pairs["By"]
        + (
            A_pairs["By"]
            - A_pairs["Fy"]
        )
        / dz_A
        * dz_AB
    )


    Bfront = B_pairs[
        ["Fx", "Fy"]
    ].to_numpy()


    treeB = cKDTree(
        Bfront
    )


    cross_r = (
        cross_match_radius_px
    )


    if cross_r is None:

        cross_r = (
            4
            * np.mean(
                [
                    sigA,
                    sigB
                ]
            )
        )


    dists, idxs = treeB.query(
        np.vstack(
            [
                x_proj,
                y_proj
            ]
        ).T,
        k=1,
        distance_upper_bound=cross_r
    )


    matched = [

        (
            i,
            j,
            d
        )

        for i, (
            j,
            d
        ) in enumerate(
            zip(
                idxs,
                dists
            )
        )

        if np.isfinite(d)

    ]


    print(
        f"Cross matches: "
        f"{len(matched)}"
    )


    # ========================================================
    # SCATTERING ANGLES
    # ========================================================

    angles = []


    for iA, jB, d in matched:

        cos_th = np.clip(
            np.dot(
                Ahat[iA],
                Bhat[jB]
            ),
            -1,
            1
        )


        angles.append(
            np.degrees(
                np.arccos(
                    cos_th
                )
            )
        )


    angles = np.array(
        angles
    )


    mean_angle = np.nanmean(
        angles
    )


    std_angle = np.nanstd(
        angles
    )


    # ========================================================
    # GAUSSIAN FIT
    # ========================================================

    hist_counts, hist_bins = np.histogram(
        angles,
        bins=120
    )


    bin_centers = (
        0.5
        * (
            hist_bins[:-1]
            + hist_bins[1:]
        )
    )


    try:

        popt, pcov = curve_fit(
            gaussian,
            bin_centers,
            hist_counts,
            p0=[
                hist_counts.max(),
                mean_angle,
                std_angle,
                hist_counts.min()
            ]
        )


        A_fit, mu_fit, sig_fit, C_fit = popt


        perr = np.sqrt(
            np.diag(
                pcov
            )
        )


        mu_err = perr[1]

        sig_err = perr[2]


        fwhm = (
            2.354820045
            * sig_fit
        )


        fit_ok = True


    except Exception as e:

        print(
            "Fit fail:",
            e
        )


        fit_ok = False

        mu_fit = np.nan

        sig_fit = np.nan

        fwhm = np.nan

        mu_err = np.nan

        sig_err = np.nan


    # ========================================================
    # SCATTERING ANGLE PLOT
    # ========================================================

    plt.figure(
        figsize=(8.0, 6.0)
    )


    ax = plt.gca()


    ax.hist(
        angles,
        bins=120,
        histtype="step",
        linewidth=2.5,
        color="navy",
        label="Experiment"
    )


    if fit_ok:

        xs = np.linspace(
            bin_centers.min(),
            bin_centers.max(),
            400
        )


        ax.plot(
            xs,
            gaussian(
                xs,
                *popt
            ),
            linewidth=2.2,
            color="teal",
            label=(
                f"Fit "
                f"(μ = {mu_fit:.3f}°)"
            )
        )


    bin_w = (
        hist_bins[1]
        - hist_bins[0]
    )


    ax.set_xlim(
        -3.0 * bin_w,
        9.3
    )


    ax.set_xlabel(
        "Scattering angle θ (deg)",
        fontsize=20,
        labelpad=10
    )


    ax.set_ylabel(
        "Counts",
        fontsize=20,
        labelpad=10
    )


    ax.set_title(
        f"Scattering Angle Distribution "
        f"(E = {energy_label})",
        fontsize=22,
        pad=14
    )


    for spine in ax.spines.values():

        spine.set_linewidth(
            1.8
        )

        spine.set_color(
            "black"
        )


    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=7,
        width=1.6,
        labelsize=18,
        top=True,
        right=True
    )


    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=4,
        width=1.2,
        top=True,
        right=True
    )


    ax.xaxis.set_minor_locator(
        AutoMinorLocator(5)
    )


    ax.yaxis.set_minor_locator(
        AutoMinorLocator(5)
    )


    ax.grid(
        False
    )


    leg = ax.legend(
        loc="upper right",
        fontsize=16,
        frameon=True
    )


    leg.get_frame().set_edgecolor(
        "black"
    )


    leg.get_frame().set_linewidth(
        1.0
    )


    leg.get_frame().set_alpha(
        1.0
    )


    plt.tight_layout()


    plt.savefig(
        os.path.join(
            out_folder,
            "05_ScatteringAngle.png"
        ),
        dpi=300
    )


    plt.show()


    pdf_one.savefig(
        plt.gcf()
    )


    plt.close()


    # ========================================================
    # N212 dx-dy PLOT
    # ========================================================

    plt.figure(
        figsize=(6, 5)
    )


    plt.hist2d(
        dxA_full,
        dyA_full,
        bins=150
    )


    plt.title(
        f"{tag_A} dx–dy "
        f"(full matched)"
    )


    plt.xlabel("dx")

    plt.ylabel("dy")


    plt.colorbar(
        label="counts"
    )


    plt.savefig(
        os.path.join(
            out_folder,
            f"06_{tag_A}_dxdy.png"
        ),
        dpi=300
    )


    plt.show()


    pdf_one.savefig(
        plt.gcf()
    )


    plt.close()


    # ========================================================
    # N213 dx-dy PLOT
    # ========================================================

    plt.figure(
        figsize=(6, 5)
    )


    plt.hist2d(
        dxB_full,
        dyB_full,
        bins=150
    )


    plt.title(
        f"{tag_B} dx–dy "
        f"(full matched)"
    )


    plt.xlabel("dx")

    plt.ylabel("dy")


    plt.colorbar(
        label="counts"
    )


    plt.savefig(
        os.path.join(
            out_folder,
            f"07_{tag_B}_dxdy.png"
        ),
        dpi=300
    )


    plt.show()


    pdf_one.savefig(
        plt.gcf()
    )


    plt.close()


    # ========================================================
    # SUMMARY
    # ========================================================

    sigma_theta_deg = np.degrees(
        np.sqrt(
            (
                sigA
                / dz_A
            ) ** 2
            +
            (
                sigB
                / dz_B
            ) ** 2
        )
    )


    E_MeV_u = parse_energy_MeV_u(
        energy_label
    )


    theta0_deg = highland_theta0_deg_carbon(
        E_MeV_u,
        al_thickness_mm
    )


    sigma_meas_deg = sig_fit

    sigma_instr_deg = sigma_theta_deg


    if (
        np.isfinite(
            sigma_meas_deg
        )
        and
        np.isfinite(
            sigma_instr_deg
        )
        and
        sigma_meas_deg
        > sigma_instr_deg
    ):

        sigma_phys_deg = np.sqrt(
            sigma_meas_deg ** 2
            -
            sigma_instr_deg ** 2
        )

    else:

        sigma_phys_deg = np.nan


    plt.figure(
        figsize=(
            8.27,
            11.69
        )
    )


    plt.axis(
        "off"
    )


    lines = [

        "CR-39 Two-Vector Scattering Report "
        "(dot-product arccos + Gaussian fit)",

        f"Energy: {energy_label}",

        f"Matched events: {len(angles)}",

        f"Mean angle = "
        f"{mean_angle:.6f}°",

        f"Std angle = "
        f"{std_angle:.6f}°",

        f"Gaussian "
        f"μ={mu_fit:.6f}±{mu_err:.6f}°  "
        f"σ={sig_fit:.6f}±{sig_err:.6f}°  "
        f"FWHM={fwhm:.6f}°",

        f"{tag_A} dx,dy="
        f"({peak_dxA:.3f}, "
        f"{peak_dyA:.3f}) "
        f"σ={sigA:.3f}",

        f"{tag_B} dx,dy="
        f"({peak_dxB:.3f}, "
        f"{peak_dyB:.3f}) "
        f"σ={sigB:.3f}",

        f"rA={rA:.2f} "
        f"rB={rB:.2f} "
        f"Cross radius={cross_r:.2f}",

        f"Global σθ≈"
        f"{sigma_theta_deg:.6f}°"

    ]


    y = 0.95


    for L in lines:

        plt.text(
            0.05,
            y,
            L,
            fontsize=10
        )

        y -= 0.04


    plt.savefig(
        os.path.join(
            out_folder,
            "08_Summary.png"
        ),
        dpi=300
    )


    plt.show()


    pdf_one.savefig(
        plt.gcf()
    )


    plt.close()


    # ========================================================
    # SAVE ANGLES CSV
    # ========================================================

    pd.DataFrame(
        {
            "theta_deg": angles
        }
    ).to_csv(
        out_csv,
        index=False
    )


    # ========================================================
    # STORE SUMMARY
    # ========================================================

    results_summary.append(
        {

            "Energy": energy_label,

            "E_MeV_u": E_MeV_u,

            "Events": int(
                len(angles)
            ),

            "Mean (deg)": float(
                mean_angle
            ),

            "Std (deg)": float(
                std_angle
            ),

            "mu_core_deg": float(
                mu_fit
            ),

            "mu_err_deg": float(
                mu_err
            ),

            "sigma (deg)": float(
                sig_fit
            ),

            "sigma_err": float(
                sig_err
            ),

            "FWHM (deg)": float(
                fwhm
            ),

            "σ_instr (deg)": float(
                sigma_theta_deg
            ),

            "σ_phys (deg)": float(
                sigma_phys_deg
            ),

            "θ0 (deg)": float(
                theta0_deg
            )

        }
    )


# ============================================================
# FINAL THEORY COMPARISON PLOT
# SINGLE-ENERGY OUTPUT
# ============================================================

df_results = pd.DataFrame(
    results_summary
).sort_values(
    "E_MeV_u"
)


print(
    "\nSummary table:"
)


print(
    df_results
)


energies_u = df_results[
    "E_MeV_u"
].values


mu_exp = df_results[
    "mu_core_deg"
].values


mu_err_arr = df_results[
    "mu_err_deg"
].values


high_vals = np.array(
    [

        highland_theta0_deg(
            E
        )

        for E in energies_u

    ]
)


ld_vals = np.array(
    [

        lynch_dahl_theta0_deg(
            E
        )

        for E in energies_u

    ]
)


plt.figure(
    figsize=(9, 6)
)


plt.errorbar(
    energies_u,
    mu_exp,
    yerr=mu_err_arr,
    fmt="d-",
    capsize=4,
    label="Experimental μ (core)",
    linewidth=2,
    markersize=8
)


plt.plot(
    energies_u,
    high_vals,
    "o--",
    label="Highland θ₀",
    linewidth=2,
    markersize=7
)


plt.plot(
    energies_u,
    ld_vals,
    "s--",
    label="Lynch–Dahl θ₀",
    linewidth=2,
    markersize=7
)


plt.xlabel(
    "Energy (MeV per nucleon)",
    fontsize=12
)


plt.ylabel(
    "Scattering angle θ (degrees)",
    fontsize=12
)


plt.title(
    "Experimental μ vs. Highland vs. Lynch–Dahl "
    "(Al foil)",
    fontsize=14
)


plt.grid(
    True,
    linestyle="--",
    alpha=0.5
)


plt.legend(
    fontsize=11
)


plt.tight_layout()


plt.savefig(
    out_plot_compare,
    dpi=300
)


plt.show()


plt.close()


# ============================================================
# FINISHED
# ============================================================

print(
    "\nSaved comparison plot:",
    out_plot_compare
)


print(
    "Saved PNGs, PDF, and CSV in:",
    out_folder
)


print(
    "Total run time: %.1f s"
    % (
        time.time()
        - start_total
    )
)
