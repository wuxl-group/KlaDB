from pathlib import Path
import math
from datetime import datetime

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from scipy.stats import fisher_exact


# ============================================================
# Basic settings
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]

POS_FILE = PROJECT_DIR / "Motif" / "KlaDB_motif_21aa.csv.gz"
BG_DIR = PROJECT_DIR / "background_proteomes"
RESULTS_DIR = PROJECT_DIR / "generated_figures" / "motif_analysis"

ENRICH_DIR = RESULTS_DIR / "background_enrichment"
ENRICH_SPECIES_DIR = ENRICH_DIR / "species"
ENRICH_GROUP_DIR = ENRICH_DIR / "group"
ENRICH_TABLE_DIR = ENRICH_DIR / "tables"
ENRICH_MAIN_DIR = ENRICH_DIR / "main_figures"
ENRICH_SUPP_DIR = ENRICH_DIR / "supplementary_figures"
CACHE_DIR = ENRICH_DIR / "background_cache"

for d in [
    ENRICH_SPECIES_DIR,
    ENRICH_GROUP_DIR,
    ENRICH_TABLE_DIR,
    ENRICH_MAIN_DIR,
    ENRICH_SUPP_DIR,
    CACHE_DIR,
]:
    d.mkdir(parents=True, exist_ok=True)


# 输出设置
SAVE_PDF = True
FIG_DPI = 600

# 热图颜色范围。固定为 -2 到 +2，更适合论文展示
HEATMAP_VMIN = -2
HEATMAP_VMAX = 2

# 星号标注阈值：减少星号密度
STAR_QVALUE_CUTOFF = 0.05
STAR_LOG2_ENRICHMENT_CUTOFF = 0.3

AAS = list("ACDEFGHIKLMNPQRSTVWY")
POSITIONS = list(range(-10, 11))
STANDARD_AA = set(AAS)
AA_TO_IDX = {aa: i for i, aa in enumerate(AAS)}

GROUP_MAP = {
    "Cavia_porcellus": "Mammals",
    "Homo_sapiens": "Mammals",
    "Mus_musculus": "Mammals",
    "Rattus_norvegicus": "Mammals",
    "Sus_scrofa": "Mammals",

    "Glycine_max": "Plants",
    "Oryza_sativa": "Plants",
    "Triticum_aestivum": "Plants",

    "Frankliniella_occidentalis": "Insects",
}


# ============================================================
# Utility functions
# ============================================================

def pretty_name(name):
    """
    Convert file-style names into publication-style labels.
    """
    name_map = {
        "All_species": "All species",
        "All species": "All species",
        "Cavia_porcellus": "Cavia porcellus",
        "Frankliniella_occidentalis": "Frankliniella occidentalis",
        "Glycine_max": "Glycine max",
        "Homo_sapiens": "Homo sapiens",
        "Mus_musculus": "Mus musculus",
        "Oryza_sativa": "Oryza sativa",
        "Rattus_norvegicus": "Rattus norvegicus",
        "Sus_scrofa": "Sus scrofa",
        "Triticum_aestivum": "Triticum aestivum",
        "Mammals": "Mammals",
        "Plants": "Plants",
        "Insects": "Insects",
    }
    return name_map.get(name, name.replace("_", " "))


def parse_fasta(path):
    header = None
    seq_lines = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines).upper()
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            yield header, "".join(seq_lines).upper()


def save_figure(fig, out_prefix):
    """
    Save figure as PNG and optionally PDF.
    """
    png_file = f"{out_prefix}.png"
    fig.savefig(png_file, dpi=FIG_DPI, bbox_inches="tight")
    print(f"[SAVED] {png_file}", flush=True)

    if SAVE_PDF:
        pdf_file = f"{out_prefix}.pdf"
        fig.savefig(pdf_file, bbox_inches="tight")
        print(f"[SAVED] {pdf_file}", flush=True)


def safe_save_dataframe(df, csv_file):
    """
    Save dataframe. If target files are occupied by Excel, save with timestamp.
    """
    try:
        df.to_csv(csv_file, index=False, encoding="utf-8-sig")
        print(f"[SAVED] {csv_file}", flush=True)

    except PermissionError:
        time_tag = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_file = csv_file.with_name(f"{csv_file.stem}_{time_tag}{csv_file.suffix}")
        df.to_csv(csv_file, index=False, encoding="utf-8-sig")

        print("[WARN] 原文件可能被 Excel 占用，已另存为新文件：", flush=True)
        print(f"[SAVED] {csv_file}", flush=True)


# ============================================================
# Data loading
# ============================================================

def load_positive_sequences():
    if not POS_FILE.exists():
        raise FileNotFoundError(f"Missing motif table: {POS_FILE}")
    frame = pd.read_csv(POS_FILE)
    data = {}
    for display_name, group in frame.groupby("species", sort=True):
        species = display_name.replace(" ", "_")
        seqs = [
            str(seq).upper() for seq in group["window_21"]
            if len(str(seq)) == 21 and str(seq)[10].upper() == "K"
            and all(aa in STANDARD_AA for aa in str(seq).upper())
        ]
        data[species] = seqs
        print(f"[POS] {species}: {len(seqs)} Kla 21-aa peptides", flush=True)

    return data


def count_21aa_sequences(seqs):
    counts = np.zeros((21, 20), dtype=int)

    for seq in seqs:
        for i, aa in enumerate(seq):
            if aa in AA_TO_IDX:
                counts[i, AA_TO_IDX[aa]] += 1

    return counts, len(seqs)


def find_background_file(species):
    candidates = [
        BG_DIR / f"{species}.fasta",
        BG_DIR / f"{species}.fa",
        BG_DIR / f"{species}.faa",
        BG_DIR / f"{species}.txt",
    ]

    for p in candidates:
        if p.exists():
            return p

    for p in BG_DIR.iterdir():
        if p.is_file() and species in p.stem:
            return p

    return None


def count_background_proteome(species):
    """
    Extract all K-centered 21-aa windows from the species proteome.
    Cache the background count matrix to speed up repeated runs.
    """
    cache_file = CACHE_DIR / f"{species}_background_counts.npz"

    if cache_file.exists():
        cache = np.load(cache_file)
        counts = cache["counts"]
        n_windows = int(cache["n_windows"])
        n_proteins = int(cache["n_proteins"])
        print(
            f"[CACHE] {species}: proteins={n_proteins}, background_K_windows={n_windows}",
            flush=True
        )
        return counts, n_windows, n_proteins

    bg_file = find_background_file(species)
    if bg_file is None:
        print(f"[WARN] 没有找到 {species} 的背景蛋白组 FASTA，跳过。", flush=True)
        return None, 0, 0

    print(f"[BG] 正在处理 {species} 背景蛋白组: {bg_file.name}", flush=True)

    counts = np.zeros((21, 20), dtype=int)
    n_windows = 0
    n_proteins = 0

    for _, protein in parse_fasta(bg_file):
        n_proteins += 1
        protein = protein.replace("*", "")

        for i, aa in enumerate(protein):
            if aa != "K":
                continue

            if i < 10 or i > len(protein) - 11:
                continue

            window = protein[i - 10:i + 11]

            if len(window) != 21:
                continue
            if window[10] != "K":
                continue
            if any(x not in STANDARD_AA for x in window):
                continue

            for j, x in enumerate(window):
                counts[j, AA_TO_IDX[x]] += 1

            n_windows += 1

    np.savez_compressed(
        cache_file,
        counts=counts,
        n_windows=n_windows,
        n_proteins=n_proteins
    )

    print(
        f"[BG OK] {species}: proteins={n_proteins}, background_K_windows={n_windows}",
        flush=True
    )

    return counts, n_windows, n_proteins


# ============================================================
# Statistics
# ============================================================

def bh_fdr(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    m = len(pvalues)

    order = np.argsort(pvalues)
    ranked = pvalues[order]

    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)

    result = np.empty(m, dtype=float)
    result[order] = q

    return result


def enrichment_analysis(pos_counts, pos_n, bg_counts, bg_n, dataset_name):
    """
    For each amino acid at each flanking position:
    Positive set: experimentally identified Kla-centered 21-aa peptides
    Background: all K-centered 21-aa windows from the corresponding proteome

    Central position 0 is fixed as K and excluded from enrichment interpretation.
    """
    rows = []
    pvalue_items = []

    for i, pos in enumerate(POSITIONS):
        for aa_idx, aa in enumerate(AAS):
            positive_count = int(pos_counts[i, aa_idx])
            background_count = int(bg_counts[i, aa_idx])

            positive_other = int(pos_n - positive_count)
            background_other = int(bg_n - background_count)

            # position 0 是固定中心乳酸化 K，不作为背景富集解释对象
            if pos == 0:
                pvalue = 1.0
                positive_freq = positive_count / pos_n if pos_n > 0 else 0
                background_freq = background_count / bg_n if bg_n > 0 else 0
                log2_enrichment = 0.0

            else:
                if pos_n == 0 or bg_n == 0:
                    pvalue = 1.0
                else:
                    _, pvalue = fisher_exact(
                        [
                            [positive_count, positive_other],
                            [background_count, background_other],
                        ],
                        alternative="greater"
                    )

                # 加伪计数，避免除以 0
                positive_freq = (positive_count + 0.5) / (pos_n + 0.5 * 20)
                background_freq = (background_count + 0.5) / (bg_n + 0.5 * 20)
                log2_enrichment = math.log2(positive_freq / background_freq)

                pvalue_items.append((len(rows), pvalue))

            rows.append({
                "dataset": dataset_name,
                "position": pos,
                "aa": aa,
                "positive_count": positive_count,
                "positive_n": pos_n,
                "positive_freq": positive_freq,
                "background_count": background_count,
                "background_n": bg_n,
                "background_freq": background_freq,
                "log2_enrichment": log2_enrichment,
                "pvalue": pvalue,
                "qvalue": 1.0,
                "significant_enriched": False,
            })

    # 只对 position != 0 的检验做 BH-FDR 校正
    if pvalue_items:
        indices = [x[0] for x in pvalue_items]
        pvalues = [x[1] for x in pvalue_items]
        qvalues = bh_fdr(pvalues)

        for row_index, qvalue in zip(indices, qvalues):
            rows[row_index]["qvalue"] = qvalue
            rows[row_index]["significant_enriched"] = (
                rows[row_index]["qvalue"] < 0.05
                and rows[row_index]["log2_enrichment"] > 0
            )

    return pd.DataFrame(rows)


def table_to_matrix(df, value_col):
    mat = pd.DataFrame(0.0, index=POSITIONS, columns=AAS)

    for _, row in df.iterrows():
        mat.loc[int(row["position"]), row["aa"]] = float(row[value_col])

    return mat


def get_star_df(df):
    """
    Publication-style star annotation:
    only show residues with adjusted P < 0.05 and log2 enrichment > 0.3.
    This avoids excessive stars in the figure.
    """
    return df[
        (df["position"] != 0)
        & (df["significant_enriched"] == True)
        & (df["qvalue"] < STAR_QVALUE_CUTOFF)
        & (df["log2_enrichment"] > STAR_LOG2_ENRICHMENT_CUTOFF)
    ].copy()


# ============================================================
# Output tables
# ============================================================

def save_enrichment_table(df, name):
    csv_file = ENRICH_TABLE_DIR / f"{name}_background_enrichment.csv"
    safe_save_dataframe(df, csv_file)


# ============================================================
# Plotting functions
# ============================================================

def save_enrichment_heatmap(df, title, out_prefix):
    """
    Single heatmap.
    """
    mat = table_to_matrix(df, "log2_enrichment").T

    fig, ax = plt.subplots(figsize=(10, 5))

    im = ax.imshow(
        mat.values,
        aspect="auto",
        interpolation="nearest",
        vmin=HEATMAP_VMIN,
        vmax=HEATMAP_VMAX,
        cmap="coolwarm"
    )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Position relative to lactylated lysine", fontsize=10)
    ax.set_ylabel("Amino acid", fontsize=10)

    ax.set_xticks(range(len(POSITIONS)))
    ax.set_xticklabels(POSITIONS, fontsize=8)
    ax.set_yticks(range(len(AAS)))
    ax.set_yticklabels(AAS, fontsize=8)

    sig_df = get_star_df(df)
    for _, row in sig_df.iterrows():
        x = POSITIONS.index(int(row["position"]))
        y = AAS.index(row["aa"])
        ax.text(x, y, "*", ha="center", va="center", fontsize=7, color="black")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log2 enrichment vs background", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    save_figure(fig, out_prefix)
    plt.close(fig)


def save_enrichment_logo(df, title, out_prefix):
    """
    Significant enrichment logo.
    Only residues passing publication-style star thresholds are shown.
    """
    logo_mat = pd.DataFrame(0.0, index=POSITIONS, columns=AAS)

    sig_df = get_star_df(df)

    for _, row in sig_df.iterrows():
        logo_mat.loc[int(row["position"]), row["aa"]] = max(
            0.0,
            float(row["log2_enrichment"])
        )

    if logo_mat.values.sum() == 0:
        print(f"[SKIP LOGO] {title}: no significant enriched residue", flush=True)
        return

    fig, ax = plt.subplots(figsize=(10, 3))

    logo = logomaker.Logo(logo_mat, ax=ax)
    logo.style_spines(visible=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Position relative to lactylated lysine", fontsize=10)
    ax.set_ylabel("log2 enrichment", fontsize=10)
    ax.set_xticks(POSITIONS)
    ax.set_xticklabels(POSITIONS, fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.set_xlim(-10.5, 10.5)

    plt.tight_layout()
    save_figure(fig, out_prefix)
    plt.close(fig)


def save_enrichment_heatmap_montage(df_dict, out_prefix, ncols=2):
    """
    Publication-style montage heatmap.
    No large figure title.
    Panel titles only.
    Fixed color scale: -2 to +2.
    Star annotation: q < 0.05 and log2 enrichment > 0.3.
    """
    valid_items = [
        (name, df)
        for name, df in df_dict.items()
        if df is not None and len(df) > 0
    ]

    if not valid_items:
        print(f"[SKIP] no valid datasets for montage", flush=True)
        return

    n = len(valid_items)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.6 * ncols, 4.2 * nrows)
    )

    axes = np.array(axes).reshape(-1)

    for ax in axes:
        ax.set_visible(False)

    last_im = None

    for ax, (name, df) in zip(axes, valid_items):
        ax.set_visible(True)

        mat = table_to_matrix(df, "log2_enrichment").T

        last_im = ax.imshow(
            mat.values,
            aspect="auto",
            interpolation="nearest",
            vmin=HEATMAP_VMIN,
            vmax=HEATMAP_VMAX,
            cmap="coolwarm"
        )

        ax.set_title(pretty_name(name), fontsize=11)

        ax.set_xticks([0, 5, 10, 15, 20])
        ax.set_xticklabels([-10, -5, 0, 5, 10], fontsize=8)

        ax.set_yticks(range(len(AAS)))
        ax.set_yticklabels(AAS, fontsize=7)

        sig_df = get_star_df(df)
        for _, row in sig_df.iterrows():
            x = POSITIONS.index(int(row["position"]))
            y = AAS.index(row["aa"])
            ax.text(x, y, "*", ha="center", va="center", fontsize=5.5, color="black")

        ax.tick_params(axis="both", length=3)

    # colorbar
    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes.tolist(), shrink=0.62)
        cbar.set_label("log2 enrichment vs background", fontsize=10)
        cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    save_figure(fig, out_prefix)
    plt.close(fig)


# ============================================================
# Main workflow
# ============================================================

def main():
    positive_data = load_positive_sequences()

    species_pos_counts = {}
    species_pos_n = {}
    species_bg_counts = {}
    species_bg_n = {}
    species_bg_proteins = {}
    species_results = {}

    summary_rows = []

    # --------------------------------------------------------
    # 1. Species-level analysis
    # --------------------------------------------------------
    for species, seqs in positive_data.items():
        print(f"\n========== {species} ==========", flush=True)

        pos_counts, pos_n = count_21aa_sequences(seqs)
        bg_counts, bg_n, bg_proteins = count_background_proteome(species)

        if bg_counts is None or bg_n == 0:
            print(f"[SKIP] {species}: background not available", flush=True)
            continue

        species_pos_counts[species] = pos_counts
        species_pos_n[species] = pos_n
        species_bg_counts[species] = bg_counts
        species_bg_n[species] = bg_n
        species_bg_proteins[species] = bg_proteins

        df = enrichment_analysis(
            pos_counts=pos_counts,
            pos_n=pos_n,
            bg_counts=bg_counts,
            bg_n=bg_n,
            dataset_name=species
        )

        species_results[species] = df

        save_enrichment_table(df, species)

        save_enrichment_heatmap(
            df,
            f"{pretty_name(species)} background-corrected Kla enrichment",
            ENRICH_SPECIES_DIR / f"{species}_background_enrichment_heatmap"
        )

        save_enrichment_logo(
            df,
            f"{pretty_name(species)} significant enrichment logo",
            ENRICH_SPECIES_DIR / f"{species}_background_enrichment_logo"
        )

        summary_rows.append({
            "dataset": species,
            "level": "species",
            "positive_21aa_windows": pos_n,
            "background_proteins": bg_proteins,
            "background_K_windows": bg_n,
            "significant_enriched_residues": int(df["significant_enriched"].sum()),
            "starred_residues_in_figures": int(len(get_star_df(df))),
        })

    # --------------------------------------------------------
    # 2. Group-level analysis
    # --------------------------------------------------------
    group_results = {}

    for group in ["Mammals", "Plants", "Insects"]:
        species_list = [
            species
            for species, mapped_group in GROUP_MAP.items()
            if mapped_group == group and species in species_pos_counts
        ]

        if not species_list:
            continue

        print(f"\n========== Group: {group} ==========", flush=True)

        pos_counts_sum = np.zeros((21, 20), dtype=int)
        bg_counts_sum = np.zeros((21, 20), dtype=int)
        pos_n_sum = 0
        bg_n_sum = 0
        bg_proteins_sum = 0

        for species in species_list:
            pos_counts_sum += species_pos_counts[species]
            bg_counts_sum += species_bg_counts[species]
            pos_n_sum += species_pos_n[species]
            bg_n_sum += species_bg_n[species]
            bg_proteins_sum += species_bg_proteins[species]

        df = enrichment_analysis(
            pos_counts=pos_counts_sum,
            pos_n=pos_n_sum,
            bg_counts=bg_counts_sum,
            bg_n=bg_n_sum,
            dataset_name=group
        )

        group_results[group] = df

        save_enrichment_table(df, group)

        save_enrichment_heatmap(
            df,
            f"{group} background-corrected Kla enrichment",
            ENRICH_GROUP_DIR / f"{group}_background_enrichment_heatmap"
        )

        save_enrichment_logo(
            df,
            f"{group} significant enrichment logo",
            ENRICH_GROUP_DIR / f"{group}_background_enrichment_logo"
        )

        summary_rows.append({
            "dataset": group,
            "level": "group",
            "positive_21aa_windows": pos_n_sum,
            "background_proteins": bg_proteins_sum,
            "background_K_windows": bg_n_sum,
            "significant_enriched_residues": int(df["significant_enriched"].sum()),
            "starred_residues_in_figures": int(len(get_star_df(df))),
        })

    # --------------------------------------------------------
    # 3. Overall analysis
    # --------------------------------------------------------
    if species_pos_counts:
        print(f"\n========== All species ==========", flush=True)

        all_pos_counts = np.zeros((21, 20), dtype=int)
        all_bg_counts = np.zeros((21, 20), dtype=int)
        all_pos_n = 0
        all_bg_n = 0
        all_bg_proteins = 0

        for species in species_pos_counts:
            all_pos_counts += species_pos_counts[species]
            all_bg_counts += species_bg_counts[species]
            all_pos_n += species_pos_n[species]
            all_bg_n += species_bg_n[species]
            all_bg_proteins += species_bg_proteins[species]

        all_df = enrichment_analysis(
            pos_counts=all_pos_counts,
            pos_n=all_pos_n,
            bg_counts=all_bg_counts,
            bg_n=all_bg_n,
            dataset_name="All_species"
        )

        group_results["All species"] = all_df

        save_enrichment_table(all_df, "All_species")

        save_enrichment_heatmap(
            all_df,
            "All species background-corrected Kla enrichment",
            ENRICH_GROUP_DIR / "All_species_background_enrichment_heatmap"
        )

        save_enrichment_logo(
            all_df,
            "All species significant enrichment logo",
            ENRICH_GROUP_DIR / "All_species_background_enrichment_logo"
        )

        summary_rows.append({
            "dataset": "All_species",
            "level": "overall",
            "positive_21aa_windows": all_pos_n,
            "background_proteins": all_bg_proteins,
            "background_K_windows": all_bg_n,
            "significant_enriched_residues": int(all_df["significant_enriched"].sum()),
            "starred_residues_in_figures": int(len(get_star_df(all_df))),
        })

    # --------------------------------------------------------
    # 4. Summary table
    # --------------------------------------------------------
    summary_df = pd.DataFrame(summary_rows)

    summary_csv = ENRICH_TABLE_DIR / "background_enrichment_summary.csv"
    safe_save_dataframe(summary_df, summary_csv)

    # --------------------------------------------------------
    # 5. Supplementary figure: species-level montage
    # --------------------------------------------------------
    print("\n[DRAW] Supplementary species enrichment heatmap montage", flush=True)

    save_enrichment_heatmap_montage(
        species_results,
        ENRICH_SUPP_DIR / "Supplementary_Figure_species_background_enrichment_heatmaps",
        ncols=3
    )

    # --------------------------------------------------------
    # 6. Main figure: overall and group-level montage
    # --------------------------------------------------------
    main_dict = {}

    if "All species" in group_results:
        main_dict["All species"] = group_results["All species"]

    for group in ["Mammals", "Plants", "Insects"]:
        if group in group_results:
            main_dict[group] = group_results[group]

    print("[DRAW] Main overall and group enrichment heatmap montage", flush=True)

    save_enrichment_heatmap_montage(
        main_dict,
        ENRICH_MAIN_DIR / "Main_Figure_overall_and_group_background_enrichment_heatmaps",
        ncols=2
    )

    print("\n背景校正 motif 富集分析完成。输出目录：", flush=True)
    print(f"物种背景富集图: {ENRICH_SPECIES_DIR}", flush=True)
    print(f"类群背景富集图: {ENRICH_GROUP_DIR}", flush=True)
    print(f"统计表: {ENRICH_TABLE_DIR}", flush=True)
    print(f"主图: {ENRICH_MAIN_DIR}", flush=True)
    print(f"补充图: {ENRICH_SUPP_DIR}", flush=True)


if __name__ == "__main__":
    main()
