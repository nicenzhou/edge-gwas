"""
Quick smoke test of exported edge_gwas functions (synthetic data only).
Saves a PDF report to tests/output/test_all_functions_quick_results.pdf.

Run from repo root:  PYTHONPATH=. python tests/test_all_functions_quick.py
"""
from __future__ import annotations

import os
import sys
import tempfile
from io import StringIO
import numpy as np
import pandas as pd

import logging
logging.getLogger("edge_gwas").setLevel(logging.WARNING)

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
DEFAULT_PDF_PATH = os.path.join(OUTPUT_DIR, "test_all_functions_quick_results.pdf")


def main(pdf_path: str = DEFAULT_PDF_PATH):
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo not in sys.path:
        sys.path.insert(0, repo)
    os.chdir(repo)

    output_dir = os.path.dirname(pdf_path)
    os.makedirs(output_dir, exist_ok=True)
    plot_paths = []
    function_status = []
    buf = StringIO()
    old_stdout = sys.stdout
    sys.stdout = buf
    try:
        _run_quick_tests(output_dir, plot_paths, function_status)
        return 0
    except Exception as e:
        print(f"\nFAIL: {e}")
        import traceback
        traceback.print_exc()
        return 1
    finally:
        sys.stdout = old_stdout
        log_lines = buf.getvalue().splitlines()
        summary = "All quick smoke tests passed." if not any("FAIL" in L for L in log_lines) else "One or more tests failed."
        _tests_dir = os.path.dirname(os.path.abspath(__file__))
        if _tests_dir not in sys.path:
            sys.path.insert(0, _tests_dir)
        from test_report_pdf import write_results_pdf
        write_results_pdf(
            pdf_path,
            title="EDGE-GWAS Quick Function Test (synthetic data)",
            log_lines=log_lines,
            summary=summary,
            plot_paths=plot_paths,
            function_status=function_status,
        )
        print(buf.getvalue())
        print(f"\nResults PDF saved to: {pdf_path}")


def _run_quick_tests(output_dir: str, plot_paths: list, function_status: list):
    print("Quick smoke tests (synthetic/small data)...")
    import edge_gwas as eg

    # Imports
    for name in eg.__all__:
        assert hasattr(eg, name), name
    print("  [OK] imports")
    function_status.append(("all exports", "OK"))

    # calculate_genomic_inflation
    lam = eg.calculate_genomic_inflation(np.random.uniform(0.01, 1, 500))
    assert 0.3 < lam < 2.0
    print("  [OK] calculate_genomic_inflation")
    function_status.append(("calculate_genomic_inflation", "OK"))

    # get_pc_covariate_list
    assert eg.get_pc_covariate_list(3) == ["PC1", "PC2", "PC3"]
    print("  [OK] get_pc_covariate_list")
    function_status.append(("get_pc_covariate_list", "OK"))

    # format_gwas_output
    gwas = pd.DataFrame({
        "chrom": [1, 1], "pos": [100, 200], "variant_id": ["a", "b"],
        "coef": [0.1, -0.2], "std_err": [0.05, 0.05], "stat": [2.0, -4.0],
        "pval": [0.04, 1e-5], "conf_int_low": [0.0, -0.3], "conf_int_high": [0.2, -0.1],
        "n_samples": [100, 100], "alpha_value": [0.5, 0.5],
        "ref_allele": ["A", "C"], "alt_allele": ["G", "T"], "eaf": [0.2, 0.3], "maf": [0.2, 0.3],
    })
    fmt = eg.format_gwas_output(gwas, include_alpha=True, format_for_locuszoom=False)
    assert "pval" in fmt.columns and len(fmt) == 2
    fmt_lz = eg.format_gwas_output(gwas, include_alpha=True, format_for_locuszoom=True)
    assert "chrom" in fmt_lz.columns or "pos" in fmt_lz.columns
    print("  [OK] format_gwas_output")
    function_status.append(("format_gwas_output", "OK"))

    # create_summary_report
    report = eg.create_summary_report(gwas, alpha_df=None, output_file=None)
    assert "EDGE" in report and "variant" in report.lower()
    print("  [OK] create_summary_report")
    function_status.append(("create_summary_report", "OK"))

    # save_results
    alpha = pd.DataFrame({"variant_id": ["a", "b"], "alpha_value": [0.5, 0.6]})
    with tempfile.TemporaryDirectory() as d:
        out = eg.save_results(gwas, alpha_df=alpha, output_prefix=os.path.join(d, "x"), save_alpha=True)
        assert os.path.exists(out["gwas"]) and os.path.exists(out["alpha"])
    print("  [OK] save_results")
    function_status.append(("save_results", "OK"))

    # load_alpha_values
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
        alpha.to_csv(f.name, index=False)
        try:
            ld = eg.load_alpha_values(f.name)
            assert "alpha_value" in ld.columns and len(ld) == 2
        finally:
            os.unlink(f.name)
    print("  [OK] load_alpha_values")
    function_status.append(("load_alpha_values", "OK"))

    # EDGEAnalysis continuous + binary
    np.random.seed(42)
    n, p = 80, 5
    X = np.random.choice([0, 1, 2], size=(n, p))
    geno = pd.DataFrame(X, index=[str(i) for i in range(n)], columns=[f"v{j}" for j in range(p)])
    y = 1 + 0.1 * X[:, 0] + np.random.randn(n) * 0.5
    pheno = pd.DataFrame({"disease": y}, index=geno.index)
    edge = eg.EDGEAnalysis(outcome_type="continuous", verbose=False)
    alpha_df = edge.calculate_alpha(geno, pheno, outcome="disease", covariates=[])
    assert not alpha_df.empty and "alpha_value" in alpha_df.columns
    gwas_df = edge.apply_alpha(geno, pheno, outcome="disease", covariates=[], alpha_values=alpha_df)
    assert not gwas_df.empty and (gwas_df["pval"] >= 0).all() and (gwas_df["pval"] <= 1).all()
    print("  [OK] EDGEAnalysis (continuous, synthetic)")
    function_status.append(("EDGEAnalysis (continuous)", "OK"))

    pheno_bin = pheno.copy()
    pheno_bin["disease_bin"] = (pheno_bin["disease"] > pheno_bin["disease"].median()).astype(int)
    edge2 = eg.EDGEAnalysis(outcome_type="binary", verbose=False)
    alpha_df2 = edge2.calculate_alpha(geno, pheno_bin, outcome="disease_bin", covariates=[])
    assert not alpha_df2.empty and "n_cases" in alpha_df2.columns
    print("  [OK] EDGEAnalysis (binary, synthetic)")
    function_status.append(("EDGEAnalysis (binary)", "OK"))

    # Visualization (synthetic gwas_df / alpha_df)
    if "chrom" not in gwas_df.columns:
        gwas_df = gwas_df.copy()
        gwas_df["chrom"] = 1
        gwas_df["pos"] = range(len(gwas_df))
    mp = os.path.join(output_dir, "quick_manhattan.png")
    eg.manhattan_plot(gwas_df, output=mp, title="Quick test Manhattan")
    if os.path.exists(mp):
        plot_paths.append(mp)
    qp = os.path.join(output_dir, "quick_qq.png")
    eg.qq_plot(gwas_df, output=qp, title="Quick test QQ")
    if os.path.exists(qp):
        plot_paths.append(qp)
    ap = os.path.join(output_dir, "quick_alpha_dist.png")
    eg.plot_alpha_distribution(alpha_df, output=ap)
    if os.path.exists(ap):
        plot_paths.append(ap)
    print("  [OK] manhattan_plot, qq_plot, plot_alpha_distribution")
    function_status.append(("manhattan_plot", "OK"))
    function_status.append(("qq_plot", "OK"))
    function_status.append(("plot_alpha_distribution", "OK"))

    # check_external_tools
    status = eg.check_external_tools()
    assert isinstance(status, dict) and "plink2" in status
    print("  [OK] check_external_tools")
    function_status.append(("check_external_tools", "OK"))

    # download_test_files (no overwrite)
    with tempfile.TemporaryDirectory() as out:
        r = eg.download_test_files(output_dir=out, overwrite=False, verbose=False)
        assert "downloaded" in r and "skipped" in r and "failed" in r
    print("  [OK] download_test_files")
    function_status.append(("download_test_files", "OK"))

    print("\nAll quick smoke tests passed.")


if __name__ == "__main__":
    sys.exit(main())
