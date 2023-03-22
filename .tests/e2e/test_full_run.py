import os
import pytest
import shutil
import subprocess as sp
import sys
import tempfile


@pytest.fixture
def setup():
    temp_dir = tempfile.mkdtemp()

    reads_fp = os.path.abspath(".tests/data/reads/")
    hosts_fp = os.path.abspath(".tests/data/hosts/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = os.path.join(project_dir, "sunbeam_config.yml")

    config_str = f"qc: {{host_fp: {hosts_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    output_fp = os.path.join(project_dir, "sunbeam_output")

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--conda-frontend",
                "conda",
                "--profile",
                project_dir,
                "--target_list",
                "all_annotate",
                "all_assembly",
                "all_coverage",
                "--directory",
                temp_dir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
        sys.exit(e)

    try:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
    except FileExistsError:
        pass

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield output_fp, benchmarks_fp


def test_full_run_assembly(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    lfinal_contigs_fp = os.path.join(output_fp, "assembly/contigs/LONG-contigs.fa")
    sfinal_contigs_fp = os.path.join(output_fp, "assembly/contigs/SHORT-contigs.fa")
    genes_fp = os.path.join(output_fp, "annotation/genes/prodigal")

    # Check output
    assert os.path.exists(lfinal_contigs_fp)
    assert os.stat(lfinal_contigs_fp).st_size > 0
    assert os.path.exists(sfinal_contigs_fp)
    for ext in ["_nucl.fa", "_prot.fa", ".gff"]:
        assert os.path.exists(os.path.join(genes_fp, f"LONG_genes{ext}"))
        assert os.path.exists(os.path.join(genes_fp, f"SHORT_genes{ext}"))


def test_full_run_annotation(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    all_samples_fp = os.path.join(output_fp, "annotation/all_samples.tsv")
    blastn_fp = os.path.join(output_fp, "annotation/blastn/bacteria/contig/LONG.btf")
    blastp_fp = os.path.join(output_fp, "annotation/blastp/prot/prodigal/LONG.btf")
    blastx_fp = os.path.join(output_fp, "annotation/blastx/prot/prodigal/LONG.btf")

    print(str(os.listdir(os.path.join(output_fp, "annotation"))))
    print(str(os.listdir(os.path.join(output_fp, "annotation/blastn"))))
    print(str(os.listdir(os.path.join(output_fp, "annotation/blastn/bacteria/contig"))))

    # Check output
    assert os.path.exists(all_samples_fp)
    assert os.stat(all_samples_fp).st_size > 0
    assert os.path.exists(blastn_fp)
    #assert os.stat(blastn_fp).st_size > 0
    assert os.path.exists(blastp_fp)
    #assert os.stat(blastp_fp).st_size > 0
    assert os.path.exists(blastx_fp)
    #assert os.stat(blastx_fp).st_size > 0


def test_full_run_coverage(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    contigs_coverage_fp = os.path.join(output_fp, "assembly/contigs_coverage.txt")

    # Check output
    assert os.path.exists(contigs_coverage_fp)
    with open(contigs_coverage_fp) as f:
        f.readline()  # Headers
        assert f.readline().strip() != ""  # Make sure there's something here
