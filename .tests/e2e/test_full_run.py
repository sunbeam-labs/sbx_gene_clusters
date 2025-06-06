import os
import pytest
import shutil
import subprocess as sp
from pathlib import Path


@pytest.fixture
def setup(tmp_path):
    reads_fp = Path(".tests/data/reads/").resolve()
    db_fp = Path(".tests/data/dbs/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()

    project_dir = tmp_path / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_gene_clusters: {{genes_fp: {db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"qc: {{host_fp: {hosts_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmp_path, project_dir

    shutil.rmtree(tmp_path)


@pytest.fixture
def run_sunbeam(setup):
    tmp_path, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_gene_clusters",
            "--directory",
            tmp_path,
            "--include",
            "sbx_gene_clusters",
            "--show-failed-logs",
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT: ", sbx_proc.stdout)
    print("STDERR: ", sbx_proc.stderr)

    if os.getenv("GITHUB_ACTIONS") == "true":
        try:
            shutil.copytree(log_fp, "logs/")
            shutil.copytree(stats_fp, "stats/")
        except FileNotFoundError:
            print("No logs or stats directory found.")

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp, sbx_proc


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp, proc = run_sunbeam

    assert proc.returncode == 0, f"Sunbeam run failed with error: {proc.stderr}"

    kegg_results_fp = (
        output_fp / "mapping" / "sbx_gene_clusters" / "tiny_kegg" / "from_kegg_1.txt"
    )
    assert (
        kegg_results_fp.exists()
    ), f"KEGG results file does not exist: {kegg_results_fp}"
    assert (
        kegg_results_fp.stat().st_size > 0
    ), f"KEGG results file is empty: {kegg_results_fp}"
