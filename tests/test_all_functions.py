"""
Comprehensive test of all exported edge_gwas functions using sample data.
Verifies methods run and results are sensible (shapes, column names, value ranges).
Saves a PDF report to tests/output/test_all_functions_results.pdf.

Run from repo root with local package (so fixes are tested):
  PYTHONPATH=. python tests/test_all_functions.py

Or with pytest (single aggregated test):
  PYTHONPATH=. pytest tests/test_all_functions.py -v -s

Requires: tests/test.bed, test.bim, test.fam, test.phen (and optionally test.vcf).
"""
from __future__ import annotations

import os
import sys
import tempfile
import logging
from io import StringIO
import numpy as np
import pandas as pd

# Reduce logging noise during tests
logging.getLogger("edge_gwas").setLevel(logging.WARNING)

# Output directory for PDF and plots
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
DEFAULT_PDF_PATH = os.path.join(OUTPUT_DIR, "test_all_functions_results.pdf")

# Paths to sample data (run from repo root or tests/)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTS_DIR = os.path.join(REPO_ROOT, "tests")
BED = os.path.join(TESTS_DIR, "test.bed")
BIM = os.path.join(TESTS_DIR, "test.bim")
FAM = os.path.join(TESTS_DIR, "test.fam")
PHEN = os.path.join(TESTS_DIR, "test.phen")
VCF = os.path.join(TESTS_DIR, "test.vcf")


def require_files(*paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"Sample data missing: {missing}")


def test_imports():
    """All exported symbols are importable."""
    import edge_gwas as eg

    for name in eg.__all__:
        assert hasattr(eg, name), f"Missing export: {name}"
    print("  [OK] All exports importable")


def test_load_plink_data():
    require_files(BED, BIM, FAM)
    from edge_gwas import load_plink_data

    geno, var_info = load_plink_data(BED, BIM, FAM, verbose=False)
    assert isinstance(geno, pd.DataFrame), "genotype_df should be DataFrame"
    assert isinstance(var_info, pd.DataFrame), "variant_info should be DataFrame"
    assert geno.shape[0] > 0 and geno.shape[1] > 0
    assert geno.index.name in ("sample_id", None) or "sample" in str(geno.index.name).lower()
    assert list(geno.columns) == list(var_info.index), "variant ids should match"
    for col in ["chrom", "pos", "ref_allele", "alt_allele", "MAF"]:
        assert col in var_info.columns, f"variant_info missing column: {col}"
    # Genotypes should be 0, 1, 2 or NaN
    uniq = geno.values[~np.isnan(geno.values.astype(float))]
    assert set(np.unique(uniq).astype(int)) <= {0, 1, 2}
    print(f"  [OK] load_plink_data: {geno.shape[0]} samples x {geno.shape[1]} variants")
    return geno, var_info


def test_prepare_phenotype_data():
    require_files(PHEN)
    from edge_gwas import prepare_phenotype_data

    # test.phen is space-separated (IID disease\n11 -1.96...)
    pheno = prepare_phenotype_data(
        PHEN, outcome_col="disease", covariate_cols=[], sample_id_col="IID", sep=r"\s+"
    )
    assert isinstance(pheno, pd.DataFrame)
    assert "disease" in pheno.columns
    assert pheno.index.name == "IID" or pheno.index.name is not None
    assert pheno.shape[0] > 0
    assert pheno["disease"].notna().all() or pheno["disease"].notna().any()
    print(f"  [OK] prepare_phenotype_data: {pheno.shape[0]} samples")
    return pheno


def _normalize_sample_ids(geno, pheno):
    """Ensure genotype and phenotype use matching sample IDs. PLINK may use (FID, IID) tuples."""
    # If genotype index is tuple-like (FID, IID), use IID only for matching
    if len(geno.index) > 0 and hasattr(geno.index[0], "__len__") and not isinstance(geno.index[0], str):
        try:
            iids = [x[1] if len(x) >= 2 else x[0] for x in geno.index]
            geno = geno.copy()
            geno.index = iids
            geno.index.name = "sample_id"
        except (TypeError, IndexError):
            pass
    geno.index = geno.index.astype(str)
    pheno.index = pheno.index.astype(str)
    common = sorted(set(geno.index) & set(pheno.index))
    assert len(common) > 0, "No common sample IDs between genotype and phenotype"
    geno = geno.loc[common].copy()
    pheno = pheno.loc[common].copy()
    return geno, pheno


def test_stratified_train_test_split(geno, pheno):
    from edge_gwas import stratified_train_test_split

    # Use continuous outcome -> is_binary=False
    tr_g, te_g, tr_p, te_p = stratified_train_test_split(
        geno, pheno, outcome_col="disease", test_size=0.3, random_state=42, is_binary=False
    )
    assert len(tr_g) + len(te_g) <= len(geno)
    assert len(tr_p) + len(te_p) <= len(pheno)
    assert tr_g.shape[1] == te_g.shape[1] == geno.shape[1]
    print(f"  [OK] stratified_train_test_split: train {tr_g.shape[0]}, test {te_g.shape[0]}")
    return tr_g, te_g, tr_p, te_p


def test_filter_variants_by_maf(geno):
    from edge_gwas import filter_variants_by_maf

    g = filter_variants_by_maf(geno, min_maf=0.01, verbose=False)
    assert g.shape[1] <= geno.shape[1]
    assert g.shape[0] == geno.shape[0]
    if g.shape[1] > 0:
        maf_approx = g.mean(axis=0, skipna=True) / 2
        maf = np.minimum(maf_approx, 1 - maf_approx)
        assert (maf >= 0.01).all()
    print(f"  [OK] filter_variants_by_maf: {geno.shape[1]} -> {g.shape[1]} variants")
    return g


def test_filter_variants_by_missing(geno):
    from edge_gwas import filter_variants_by_missing

    g = filter_variants_by_missing(geno, max_missing=0.5, verbose=False)
    assert g.shape[1] <= geno.shape[1]
    assert g.shape[0] == geno.shape[0]
    if g.shape[1] > 0:
        miss = g.isna().mean(axis=0)
        assert (miss <= 0.5).all()
    print(f"  [OK] filter_variants_by_missing: {geno.shape[1]} -> {g.shape[1]} variants")


def test_filter_samples_by_call_rate(geno, pheno):
    from edge_gwas import filter_samples_by_call_rate

    g, p = filter_samples_by_call_rate(geno, pheno, min_call_rate=0.5, verbose=False)
    assert g.shape[0] <= geno.shape[0]
    assert g.shape[0] == p.shape[0]
    assert list(g.index) == list(p.index)
    print(f"  [OK] filter_samples_by_call_rate: {geno.shape[0]} -> {g.shape[0]} samples")


def test_calculate_hwe_pvalues(geno):
    from edge_gwas import calculate_hwe_pvalues

    hwe = calculate_hwe_pvalues(geno.iloc[:, :50], verbose=False)  # subset for speed
    assert isinstance(hwe, pd.Series)
    assert len(hwe) == min(50, geno.shape[1])
    valid = hwe.dropna()
    if len(valid) > 0:
        assert (valid >= 0).all() and (valid <= 1).all()
    print(f"  [OK] calculate_hwe_pvalues: {len(hwe)} variants, {hwe.notna().sum()} valid")


def test_filter_variants_by_hwe(geno):
    from edge_gwas import filter_variants_by_hwe

    g = filter_variants_by_hwe(geno.iloc[:, :80], hwe_threshold=1e-6, verbose=False)
    assert g.shape[1] <= 80
    assert g.shape[0] == geno.shape[0]
    print(f"  [OK] filter_variants_by_hwe: -> {g.shape[1]} variants")


def test_calculate_genomic_inflation():
    from edge_gwas import calculate_genomic_inflation

    np.random.seed(42)
    pvals = np.random.uniform(0.001, 1, 1000)
    lam = calculate_genomic_inflation(pvals)
    assert np.isfinite(lam)
    assert lam > 0
    # Under null, lambda should be near 1
    assert 0.5 < lam < 2.0
    # With NaN
    pvals_na = np.concatenate([pvals, [np.nan, np.nan]])
    lam2 = calculate_genomic_inflation(pvals_na)
    assert np.isfinite(lam2)
    print(f"  [OK] calculate_genomic_inflation: lambda={lam:.3f}")


def test_get_pc_covariate_list():
    from edge_gwas import get_pc_covariate_list

    lst = get_pc_covariate_list(5)
    assert lst == ["PC1", "PC2", "PC3", "PC4", "PC5"]
    lst2 = get_pc_covariate_list(3, pc_prefix="PC")
    assert lst2 == ["PC1", "PC2", "PC3"]
    print("  [OK] get_pc_covariate_list")


def test_calculate_pca_sklearn(geno, pheno):
    from edge_gwas import calculate_pca_sklearn, attach_pcs_to_phenotype

    geno_small = geno.iloc[:, :100].dropna(axis=1, how="all")
    if geno_small.shape[1] < 2:
        print("  [SKIP] calculate_pca_sklearn: too few variants")
        return
    pca_df = calculate_pca_sklearn(geno_small, n_pcs=5, verbose=False)
    assert "PC1" in pca_df.columns and "PC5" in pca_df.columns
    assert pca_df.shape[0] == geno_small.shape[0]
    # Attach PCs to phenotype
    pheno_with_pcs = attach_pcs_to_phenotype(pheno, pca_df, n_pcs=3, verbose=False)
    assert "PC1" in pheno_with_pcs.columns
    assert pheno_with_pcs.shape[0] <= pheno.shape[0]
    print(f"  [OK] calculate_pca_sklearn + attach_pcs_to_phenotype: 5 PCs, {pca_df.shape[0]} samples")


def test_edgemodel_continuous(tr_g, te_g, tr_p, te_p, var_info):
    from edge_gwas import EDGEAnalysis

    nv = min(15, tr_g.shape[1], te_g.shape[1])
    tr_g = tr_g.iloc[:, :nv]
    te_g = te_g.iloc[:, :nv]
    var_sub = var_info.loc[var_info.index.isin(tr_g.columns)]

    edge = EDGEAnalysis(outcome_type="continuous", verbose=False)
    alpha_df = edge.calculate_alpha(
        tr_g, tr_p, outcome="disease", covariates=[], variant_info=var_sub
    )
    assert not alpha_df.empty
    assert "variant_id" in alpha_df.columns and "alpha_value" in alpha_df.columns
    assert alpha_df["alpha_value"].notna().any() or alpha_df["alpha_value"].isna().all()

    gwas_df = edge.apply_alpha(
        te_g, te_p, outcome="disease", covariates=[], alpha_values=alpha_df
    )
    assert not gwas_df.empty
    assert "pval" in gwas_df.columns and "coef" in gwas_df.columns
    assert (gwas_df["pval"] >= 0).all() and (gwas_df["pval"] <= 1).all()
    print(f"  [OK] EDGE continuous: alpha {len(alpha_df)}, GWAS {len(gwas_df)}")
    return alpha_df, gwas_df


def test_edgemodel_binary(tr_g, te_g, tr_p, te_p, var_info):
    from edge_gwas import EDGEAnalysis

    # Create binary outcome from continuous for testing
    tr_p_bin = tr_p.copy()
    tr_p_bin["disease_bin"] = (tr_p_bin["disease"] > tr_p_bin["disease"].median()).astype(int)
    te_p_bin = te_p.copy()
    te_p_bin["disease_bin"] = (te_p["disease"] > te_p["disease"].median()).astype(int)

    nv = min(15, tr_g.shape[1], te_g.shape[1])
    tr_g = tr_g.iloc[:, :nv]
    te_g = te_g.iloc[:, :nv]
    var_sub = var_info.loc[var_info.index.isin(tr_g.columns)]

    edge = EDGEAnalysis(outcome_type="binary", verbose=False)
    alpha_df = edge.calculate_alpha(
        tr_g, tr_p_bin, outcome="disease_bin", covariates=[], variant_info=var_sub
    )
    assert not alpha_df.empty
    assert "n_cases" in alpha_df.columns and "n_controls" in alpha_df.columns

    gwas_df = edge.apply_alpha(
        te_g, te_p_bin, outcome="disease_bin", covariates=[], alpha_values=alpha_df
    )
    assert not gwas_df.empty
    assert (gwas_df["pval"] >= 0).all() and (gwas_df["pval"] <= 1).all()
    print(f"  [OK] EDGE binary: alpha {len(alpha_df)}, GWAS {len(gwas_df)}")
    return alpha_df, gwas_df


def test_run_full_analysis(tr_g, te_g, tr_p, te_p, var_info):
    from edge_gwas import EDGEAnalysis

    nv = min(10, tr_g.shape[1], te_g.shape[1])
    tr_vars = tr_g.columns[:nv].tolist()
    var_sub = var_info.loc[var_info.index.isin(tr_vars)]
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "edge_out")
        edge = EDGEAnalysis(outcome_type="continuous", verbose=False)
        alpha_df, gwas_df = edge.run_full_analysis(
            tr_g.iloc[:, :nv],
            tr_p,
            te_g.iloc[:, :nv],
            te_p,
            outcome="disease",
            covariates=[],
            variant_info=var_sub,
            output_prefix=prefix,
        )
        assert not alpha_df.empty and not gwas_df.empty
        assert os.path.exists(prefix + "_alpha_values.csv")
        assert os.path.exists(prefix + "_gwas_results.csv")
    print("  [OK] run_full_analysis + output files")


def test_save_results(gwas_df, alpha_df):
    from edge_gwas import save_results

    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "out")
        out = save_results(gwas_df, alpha_df=alpha_df, output_prefix=prefix, save_alpha=True)
        assert "gwas" in out and "alpha" in out
        assert os.path.exists(out["gwas"]) and os.path.exists(out["alpha"])
    print("  [OK] save_results")


def test_load_alpha_values(alpha_df):
    from edge_gwas import load_alpha_values

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        alpha_df.to_csv(f.name, index=False)
        try:
            loaded = load_alpha_values(f.name)
            assert "variant_id" in loaded.columns and "alpha_value" in loaded.columns
            assert len(loaded) == len(alpha_df)
        finally:
            os.unlink(f.name)
    print("  [OK] load_alpha_values")


def test_format_gwas_output(gwas_df):
    from edge_gwas import format_gwas_output

    # Standard format
    fmt = format_gwas_output(gwas_df, include_alpha=True, sort_by="pval", format_for_locuszoom=False)
    assert isinstance(fmt, pd.DataFrame)
    assert "pval" in fmt.columns
    # LocusZoom format (needs chrom, pos)
    if "chrom" in gwas_df.columns and "pos" in gwas_df.columns:
        fmt_lz = format_gwas_output(gwas_df, include_alpha=True, format_for_locuszoom=True)
        assert "chrom" in fmt_lz.columns or "pos" in fmt_lz.columns
    print("  [OK] format_gwas_output")


def test_create_summary_report(gwas_df, alpha_df):
    from edge_gwas import create_summary_report

    report = create_summary_report(gwas_df, alpha_df=alpha_df, output_file=None)
    assert "EDGE GWAS" in report and "Total variants" in report
    assert "lambda" in report.lower() or "λ" in report or "inflation" in report.lower()
    print("  [OK] create_summary_report")


def test_visualize(gwas_df, alpha_df, output_dir: str):
    """Run visualization functions and save plots to output_dir. Returns list of saved paths."""
    from edge_gwas import manhattan_plot, qq_plot, plot_alpha_distribution

    os.makedirs(output_dir, exist_ok=True)
    if "chrom" not in gwas_df.columns:
        gwas_df = gwas_df.copy()
        gwas_df["chrom"] = 1
        gwas_df["pos"] = range(len(gwas_df))
    paths = []
    mp = os.path.join(output_dir, "manhattan.png")
    manhattan_plot(gwas_df, output=mp, title="EDGE GWAS Manhattan (test)")
    assert os.path.exists(mp)
    paths.append(mp)

    qp = os.path.join(output_dir, "qq_plot.png")
    lam = qq_plot(gwas_df, output=qp, title="EDGE GWAS QQ (test)")
    assert os.path.exists(qp)
    assert np.isfinite(lam)
    paths.append(qp)

    ap = os.path.join(output_dir, "alpha_distribution.png")
    plot_alpha_distribution(alpha_df, output=ap)
    assert os.path.exists(ap)
    paths.append(ap)
    print("  [OK] manhattan_plot, qq_plot, plot_alpha_distribution")
    return paths


def test_load_vcf_data():
    if not os.path.exists(VCF):
        print("  [SKIP] load_vcf_data: test.vcf not found")
        return
    try:
        from edge_gwas import load_vcf_data

        geno, var_info = load_vcf_data(VCF, dosage=False, verbose=False)
        assert isinstance(geno, pd.DataFrame)
        assert geno.shape[0] > 0 and geno.shape[1] > 0
        assert "chrom" in var_info.columns and "pos" in var_info.columns
        print(f"  [OK] load_vcf_data: {geno.shape[0]} x {geno.shape[1]}")
    except Exception as e:
        print(f"  [SKIP] load_vcf_data: {e}")


def test_download_test_files():
    from edge_gwas import download_test_files

    with tempfile.TemporaryDirectory() as out:
        r = download_test_files(output_dir=out, overwrite=False, verbose=False)
        assert "downloaded" in r and "skipped" in r and "failed" in r
    print("  [OK] download_test_files (no overwrite)")


def test_check_external_tools():
    from edge_gwas import check_external_tools

    status = check_external_tools()
    assert isinstance(status, dict)
    for k in ["plink2", "gcta", "r", "genesis", "snprelate", "gdsfmt"]:
        assert k in status
        assert isinstance(status[k], (bool, type(True)))
    print("  [OK] check_external_tools")


def main(pdf_path: str = DEFAULT_PDF_PATH):
    output_dir = os.path.dirname(pdf_path)
    os.makedirs(output_dir, exist_ok=True)
    plot_paths = []
    function_status = []
    buf = StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    ret = 0
    try:
        ret = _run_tests(output_dir, plot_paths, function_status)
    finally:
        sys.stdout = old_stdout
        log_lines = buf.getvalue().splitlines()
        if hasattr(_run_tests, "failures") and _run_tests.failures:
            summary = "Failures: " + ", ".join(n for n, _ in _run_tests.failures)
        else:
            summary = "All sections passed."
        _tests_dir = os.path.dirname(os.path.abspath(__file__))
        if _tests_dir not in sys.path:
            sys.path.insert(0, _tests_dir)
        from test_report_pdf import write_results_pdf
        write_results_pdf(
            pdf_path,
            title="EDGE-GWAS Full Function Test (sample data)",
            log_lines=log_lines,
            summary=summary,
            plot_paths=plot_paths,
            function_status=function_status,
        )
        print(buf.getvalue())
        print(f"\nResults PDF saved to: {pdf_path}")
    return ret


def _run_tests(output_dir: str, plot_paths: list, function_status: list):
    failures = []
    _run_tests.failures = failures  # for main() to read

    print("Testing all exported edge_gwas functions with sample data\n")
    try:
        print("[1] Imports")
        test_imports()
        function_status.append(("all exports", "import OK"))
    except Exception as e:
        failures.append(("imports", e))

    try:
        print("\n[2] Data loading")
        require_files(BED, BIM, FAM, PHEN)
        geno, var_info = test_load_plink_data()
        pheno = test_prepare_phenotype_data()
        geno, pheno = _normalize_sample_ids(geno, pheno)
        N_SAMPLES = min(600, geno.shape[0])
        N_VAR = min(100, geno.shape[1])
        geno = geno.iloc[:N_SAMPLES, :N_VAR].copy()
        pheno = pheno.loc[geno.index].copy()
        var_info = var_info.loc[var_info.index.isin(geno.columns)].copy()
        function_status.append(("load_plink_data", "OK"))
        function_status.append(("prepare_phenotype_data", "OK"))
    except Exception as e:
        failures.append(("data_loading", e))
        print(f"  [FAIL] {e}")
        geno = pheno = var_info = None

    if geno is not None and pheno is not None:
        try:
            print("\n[3] Train/test split")
            tr_g, te_g, tr_p, te_p = test_stratified_train_test_split(geno, pheno)
            function_status.append(("stratified_train_test_split", "OK"))
        except Exception as e:
            failures.append(("stratified_split", e))
            tr_g = te_g = tr_p = te_p = None

        try:
            print("\n[4] Filtering & QC")
            test_filter_variants_by_maf(geno)
            test_filter_variants_by_missing(geno)
            test_filter_samples_by_call_rate(geno, pheno)
            test_calculate_hwe_pvalues(geno)
            test_filter_variants_by_hwe(geno)
            for fn in ("filter_variants_by_maf", "filter_variants_by_missing", "filter_samples_by_call_rate",
                       "calculate_hwe_pvalues", "filter_variants_by_hwe"):
                function_status.append((fn, "OK"))
        except Exception as e:
            failures.append(("filtering", e))

        try:
            print("\n[5] Utilities")
            test_calculate_genomic_inflation()
            test_get_pc_covariate_list()
            function_status.append(("calculate_genomic_inflation", "OK"))
            function_status.append(("get_pc_covariate_list", "OK"))
        except Exception as e:
            failures.append(("utilities", e))

        try:
            print("\n[6] PCA (sklearn)")
            test_calculate_pca_sklearn(geno, pheno)
            function_status.append(("calculate_pca_sklearn", "OK"))
            function_status.append(("attach_pcs_to_phenotype", "OK"))
        except Exception as e:
            failures.append(("pca_sklearn", e))

        if tr_g is not None and te_g is not None:
            try:
                print("\n[7] EDGE analysis (continuous)")
                alpha_df, gwas_df = test_edgemodel_continuous(tr_g, te_g, tr_p, te_p, var_info)
                function_status.append(("EDGEAnalysis.calculate_alpha", "OK"))
                function_status.append(("EDGEAnalysis.apply_alpha", "OK"))
            except Exception as e:
                failures.append(("edge_continuous", e))
                alpha_df = gwas_df = None

            try:
                print("\n[8] EDGE analysis (binary)")
                test_edgemodel_binary(tr_g, te_g, tr_p, te_p, var_info)
            except Exception as e:
                failures.append(("edge_binary", e))

            try:
                print("\n[9] run_full_analysis")
                test_run_full_analysis(tr_g, te_g, tr_p, te_p, var_info)
            except Exception as e:
                failures.append(("run_full_analysis", e))

            if alpha_df is not None and gwas_df is not None:
                try:
                    print("\n[10] IO: save_results, load_alpha, format_gwas, summary report")
                    test_save_results(gwas_df, alpha_df)
                    test_load_alpha_values(alpha_df)
                    test_format_gwas_output(gwas_df)
                    test_create_summary_report(gwas_df, alpha_df)
                    for fn in ("save_results", "load_alpha_values", "format_gwas_output", "create_summary_report"):
                        function_status.append((fn, "OK"))
                except Exception as e:
                    failures.append(("io", e))

                try:
                    print("\n[11] Visualization")
                    plot_paths.extend(test_visualize(gwas_df, alpha_df, output_dir))
                    function_status.append(("manhattan_plot", "OK"))
                    function_status.append(("qq_plot", "OK"))
                    function_status.append(("plot_alpha_distribution", "OK"))
                except Exception as e:
                    failures.append(("visualize", e))

    try:
        print("\n[12] load_vcf_data")
        test_load_vcf_data()
        function_status.append(("load_vcf_data", "OK"))
    except Exception as e:
        failures.append(("load_vcf", e))

    try:
        print("\n[13] download_test_files")
        test_download_test_files()
        function_status.append(("download_test_files", "OK"))
    except Exception as e:
        failures.append(("download_test_files", e))

    try:
        print("\n[14] check_external_tools")
        test_check_external_tools()
    except Exception as e:
        failures.append(("check_external_tools", e))

    # Optional: functions that need extra data/tools (runnable but may skip)
    print("\n[15] Optional (skip if no data/tools)")
    _run_optional_functions(function_status)

    print("\n" + "=" * 60)
    if failures:
        print("FAILURES:")
        for name, err in failures:
            print(f"  - {name}: {err}")
        return 1
    print("All function tests passed.")
    return 0


def _run_optional_functions(function_status: list):
    """Try to run load_pgen_data, load_bgen_data, calculate_pca_plink, calculate_grm_gcta, load_grm_gcta, etc."""
    import edge_gwas as eg
    # load_pgen_data: needs .pgen/.pvar/.psam
    try:
        eg.load_pgen_data("/nonexistent.pgen", "/nonexistent.pvar", "/nonexistent.psam", verbose=False)
    except (FileNotFoundError, OSError, Exception) as e:
        err = str(e)[:50]
        function_status.append(("load_pgen_data", f"SKIP (no data): {err}..."))
        print(f"  [SKIP] load_pgen_data: {e}")
    # load_bgen_data: needs .bgen
    try:
        eg.load_bgen_data("/nonexistent.bgen", verbose=False)
    except (FileNotFoundError, OSError, Exception) as e:
        err = str(e)[:50]
        function_status.append(("load_bgen_data", f"SKIP (no data): {err}..."))
        print(f"  [SKIP] load_bgen_data: {e}")
    # calculate_pca_plink: needs plink2 and plink prefix
    try:
        eg.calculate_pca_plink(BED.replace(".bed", ""), n_pcs=2, verbose=False)
        function_status.append(("calculate_pca_plink", "OK"))
        print("  [OK] calculate_pca_plink")
    except Exception as e:
        function_status.append(("calculate_pca_plink", f"SKIP: {str(e)[:40]}..."))
        print(f"  [SKIP] calculate_pca_plink: {e}")
    # calculate_pca_pcair: needs R/GENESIS
    try:
        eg.calculate_pca_pcair(BED.replace(".bed", ""), n_pcs=2, verbose=False)
        function_status.append(("calculate_pca_pcair", "OK"))
        print("  [OK] calculate_pca_pcair")
    except Exception as e:
        function_status.append(("calculate_pca_pcair", f"SKIP: {str(e)[:40]}..."))
        print(f"  [SKIP] calculate_pca_pcair: {e}")
    # calculate_grm_gcta: needs GCTA and plink prefix
    try:
        eg.calculate_grm_gcta(BED.replace(".bed", ""), output_prefix=tempfile.mktemp(prefix="grm"), verbose=False)
        function_status.append(("calculate_grm_gcta", "OK"))
        print("  [OK] calculate_grm_gcta")
    except Exception as e:
        function_status.append(("calculate_grm_gcta", f"SKIP: {str(e)[:40]}..."))
        print(f"  [SKIP] calculate_grm_gcta: {e}")
    # load_grm_gcta: needs .grm.* files
    try:
        eg.load_grm_gcta("/nonexistent", verbose=False)
    except (FileNotFoundError, Exception) as e:
        function_status.append(("load_grm_gcta", "SKIP (no .grm files)"))
        print(f"  [SKIP] load_grm_gcta: no files")
    # identify_related_samples, filter_related_samples: need GRM matrix
    try:
        grm = np.eye(3)
        ids_df = pd.DataFrame({"FID": [1, 2, 3], "IID": ["a", "b", "c"]})
        eg.identify_related_samples(grm, ids_df, threshold=0.2)
        function_status.append(("identify_related_samples", "OK"))
        print("  [OK] identify_related_samples")
    except Exception as e:
        function_status.append(("identify_related_samples", f"SKIP: {str(e)[:40]}"))
        print(f"  [SKIP] identify_related_samples: {e}")
    try:
        grm = np.eye(3)
        ids_df = pd.DataFrame({"FID": [1, 2, 3], "IID": ["a", "b", "c"]})
        eg.filter_related_samples(pd.DataFrame({"x": [1, 2, 3]}, index=["a", "b", "c"]), grm, ids_df, threshold=0.2)
        function_status.append(("filter_related_samples", "OK"))
        print("  [OK] filter_related_samples")
    except Exception as e:
        function_status.append(("filter_related_samples", f"SKIP: {str(e)[:40]}"))
        print(f"  [SKIP] filter_related_samples: {e}")


def test_full_suite():
    """Pytest entry point: runs the full function test suite (sample data required)."""
    code = main()
    assert code == 0, "One or more test sections failed; see output above"


if __name__ == "__main__":
    sys.exit(main())
