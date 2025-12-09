from pathlib import Path

import pandas as pd
import pytest
from cloudpathlib import CloudPath, implementation_registry
from cloudpathlib.local import (
    LocalS3Client,
    LocalS3Path,
    local_s3_implementation,
)
from pandas.testing import assert_frame_equal

from gopher.parsers.tabular import read_diann


@pytest.fixture
def cloud_asset_file(monkeypatch):
    """Patch CloudPath dispatch and seed LocalS3Client test assets."""
    monkeypatch.setitem(implementation_registry, "s3", local_s3_implementation)

    # Option 1: Use LocalS3Path to set up test assets directly
    local_cloud_path = LocalS3Path(
        "s3://cloudpathlib-test-bucket/diann_report.pg_mat.tsv"
    )
    # Simulated DIANN output
    mock_data = (
        "Protein.Group\tProtein.Ids\tProtein.Names\tGenes\t"
        "First.Protein.Description\tIntensity.Sample1\tIntensity.Sample2",
        "PG1\tP12345;P67890\tProtein A\tGENE1\tDescription A\t1000\t2000",
        "PG2\tP23456\tProtein B\tGENE2\tDescription B\t1500\t2500",
    )
    local_cloud_path.write_text("\n".join(mock_data))

    local_cloud_path_genes = LocalS3Path(
        "s3://cloudpathlib-test-bucket/diann_report.gg_mat.tsv"
    )
    # Simulated DIANN output
    mock_data = (
        "Genes\tFirst.Protein.Description\tIntensity.Sample1\t"
        "Intensity.Sample2",
        "GENE1\tDescription A\t1000\t2000",
        "GENE2\tDescription B\t1500\t2500",
    )
    local_cloud_path_genes.write_text("\n".join(mock_data))

    cloud_path_1 = CloudPath(
        "s3://cloudpathlib-test-bucket/diann_report.pg_mat.tsv"
    )
    assert cloud_path_1.exists()

    # Expected DataFrame
    expected = pd.DataFrame(
        {
            # The real diann data has float values in the intensities.
            "Intensity.Sample1": [1000.0, 1500.0],
            "Intensity.Sample2": [2000.0, 2500.0],
        },
        index=["P12345", "P23456"],
    )
    expected.index.name = "Protein"

    yield {"cloud_path": cloud_path_1, "expected": expected}

    # Clean up temp directory and replace with new one.
    LocalS3Client.reset_default_storage_dir()


def test_read_diann_removes_metadata_and_sets_index_cloud(cloud_asset_file):
    """Read DIANN from cloud path and set protein index."""
    result = read_diann(
        "s3://cloudpathlib-test-bucket/diann_report.pg_mat.tsv"
    )
    assert_frame_equal(result, cloud_asset_file["expected"])


def test_read_diann_removes_metadata_and_sets_index_local(
    cloud_asset_file, tmpdir
):
    """Read DIANN from local path and set protein index."""
    local_path = Path(tmpdir) / "diann_report.pg_mat.tsv"
    with open(local_path, "w") as f:
        f.write(cloud_asset_file["cloud_path"].read_text())

    result = read_diann(local_path)
    assert_frame_equal(result, cloud_asset_file["expected"])


def test_read_diann_faile_with_gg(cloud_asset_file):
    """Fail when DIANN input lacks protein group columns."""
    with pytest.raises(ValueError) as e:
        read_diann("s3://cloudpathlib-test-bucket/diann_report.gg_mat.tsv")

    assert "Expected columns" in str(e.value.args[0])
