import logging

import pytest

from multiqc import config, report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.haplogrep3 import MultiqcModule
from multiqc.types import ColumnKey, SampleGroup, SectionKey

HEADER = '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"'
EXTENDED_HEADER = HEADER + '\t"Not_Found_Polys"\t"Found_Polys"\t"Remaining_Polys"\t"AAC_In_Remainings"\t"Input_Sample"'

ONE_SAMPLE = f"""{HEADER}
"HG002"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"
"""

TWO_SAMPLES = f"""{HEADER}
"HG003"\t"K1a1b1a"\t"1"\t"0.9815"\t"1-16569"
"HG004"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"
"""

HITS = f"""{HEADER}
"HG002"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"
"HG002"\t"H5a"\t"2"\t"0.9162"\t"1-16569"
"HG003"\t"K1a1b1a"\t"1"\t"0.9815"\t"1-16569"
"HG003"\t"K1a1b1"\t"2"\t"0.9457"\t"1-16569"
"""

EXTENDED = f"""{EXTENDED_HEADER}
"HG002"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"\t""\t"263G 456T"\t"302.1CC (globalPrivateMutation) 6800G (globalPrivateMutation)"\t""\t"263G 302.1CC 456T 6800G"
"HG003"\t"K1a1b1a"\t"1"\t"0.9815"\t"1-16569"\t"16093C"\t"73G 114T"\t"3107d (hotspot) 16223T (localPrivateMutation) 16519C (hotspot)"\t""\t"73G 114T 3107d 16223T 16519C"
"""

# The file that `--write-qc` writes. It is not a classification report.
QC_FILE = """"SampleID"\t"Haplogroup"\t"Type"\t"Message"\t"Missing Mutations"\t"Global Private Mutations"
"HG002"\t"H5a7"\t"warning"\t"The sample contains 4 global private mutation(s)"\t""\t"302.1CC"
"""


@pytest.fixture(autouse=True)
def reset():
    config.reset()
    report.reset()
    yield
    config.reset()
    report.reset()


def run_module(tmp_path, files: dict[str, str]) -> MultiqcModule:
    paths = []
    for name, content in files.items():
        path = tmp_path / name
        path.write_text(content)
        paths.append(path)
    report.analysis_files = paths
    report.search_files(["haplogrep3"])
    config.preserve_module_raw_data = True
    return MultiqcModule()


def module_data(module: MultiqcModule) -> dict:
    assert module.saved_raw_data is not None
    return module.saved_raw_data["multiqc_haplogrep3"]


def test_one_sample(tmp_path):
    # The file name must not have an effect on the sample name.
    m = run_module(tmp_path, {"any_name.tsv": ONE_SAMPLE})
    assert module_data(m) == {"HG002": {"haplogroup": "H5a7", "quality": 1.0, "range": "1-16569"}}
    assert len(m.sections) == 1


def test_more_than_one_sample(tmp_path):
    data = module_data(run_module(tmp_path, {"trio.txt": TWO_SAMPLES}))
    assert list(data) == ["HG003", "HG004"]
    assert data["HG003"]["haplogroup"] == "K1a1b1a"
    assert data["HG003"]["quality"] == 0.9815


def test_header_only(tmp_path):
    with pytest.raises(ModuleNoSamplesFound):
        run_module(tmp_path, {"empty.txt": HEADER + "\n"})


def test_header_only_with_other_report(tmp_path):
    data = module_data(run_module(tmp_path, {"empty.txt": HEADER + "\n", "HG002.txt": ONE_SAMPLE}))
    assert list(data) == ["HG002"]


def test_hits_keeps_rank_1(tmp_path):
    data = module_data(run_module(tmp_path, {"hits.txt": HITS}))
    assert {s: d["haplogroup"] for s, d in data.items()} == {"HG002": "H5a7", "HG003": "K1a1b1a"}
    assert data["HG002"]["quality"] == 1.0


def test_extended_report(tmp_path):
    data = module_data(run_module(tmp_path, {"extended.txt": EXTENDED}))
    assert data["HG002"]["not_found_polys"] == 0
    assert data["HG002"]["remaining_polys"] == 2
    assert data["HG003"]["not_found_polys"] == 1
    assert data["HG003"]["remaining_polys"] == 3


def test_haplogroup_with_apostrophe(tmp_path):
    content = f"""{HEADER}\n"SIM_L3"\t"L3e'i'k'x"\t"1"\t"1.0000"\t"1-16569"\n"""
    data = module_data(run_module(tmp_path, {"sim.txt": content}))
    assert data["SIM_L3"]["haplogroup"] == "L3e'i'k'x"


def test_sample_names_are_not_changed(tmp_path):
    # Haplogrep 3 writes a quote character in a sample name as two quote characters.
    content = f"""{HEADER}
"HG002.sorted.bam"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"
"/data/HG004.cram"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"
\"\"\"quoted\"\"\"\t"H5a7"\t"1"\t"1.0000"\t"1-16569"
"""
    data = module_data(run_module(tmp_path, {"names.txt": content}))
    assert list(data) == ["HG002.sorted.bam", "/data/HG004.cram", '"quoted"']


def test_sample_in_two_reports_gives_warning(tmp_path, caplog):
    second = ONE_SAMPLE.replace("H5a7", "H5a").replace("1.0000", "0.9162")
    with caplog.at_level(logging.WARNING, logger="multiqc.modules.haplogrep3.haplogrep3"):
        data = module_data(run_module(tmp_path, {"a.txt": ONE_SAMPLE, "b.txt": second}))
    assert list(data) == ["HG002"]
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    assert "'HG002' is in more than one Haplogrep 3 report" in warnings[0]


def test_qc_file_is_ignored(tmp_path):
    with pytest.raises(ModuleNoSamplesFound):
        run_module(tmp_path, {"sample.qc.txt": QC_FILE})


def test_general_stats(tmp_path):
    run_module(tmp_path, {"HG002.txt": ONE_SAMPLE})
    rows = report.general_stats_data[SectionKey("haplogrep3")][SampleGroup("HG002")]
    assert rows[0].data[ColumnKey("haplogroup")] == "H5a7"
    headers = report.general_stats_headers[SectionKey("haplogrep3")]
    assert list(headers) == ["haplogroup", "quality"]
    assert not headers[ColumnKey("haplogroup")].get("hidden")
    assert headers[ColumnKey("quality")]["hidden"] is True
