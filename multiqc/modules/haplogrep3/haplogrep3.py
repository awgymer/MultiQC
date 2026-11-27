import csv
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import ColumnDict
from multiqc.types import ColumnKey, LoadedFileDict

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module reads the report that `haplogrep3 classify` writes:

    ```bash
    haplogrep3 classify --tree phylotree-rcrs@17.2 --in sample.vcf.gz --out sample.txt
    ```

    The module finds the report by its header line, so the report can have any file name.
    One report can contain more than one sample.

    The sample names come from the `SampleID` column of the report, not from the file name.
    Haplogrep 3 takes them from the sample names of a VCF file, from the first word of each
    header line of a FASTA file, or from the first column of an HSD file.
    The module uses these names as they are, as the Haplocheck module does.

    With `--hits N`, the report contains the best `N` haplogroups for each sample.
    The module shows only the best haplogroup (`Rank` 1).

    With `--extend-report`, the report contains the polymorphisms of each sample.
    The module then adds two columns to the table: the number of expected polymorphisms
    that the sample does not contain, and the number of polymorphisms that the haplogroup
    does not explain.

    The report does not contain the name of the tree or the version of Haplogrep 3.

    Haplogrep 3 version 3.2.2 is tested.
    """

    def __init__(self):
        super().__init__(
            name="Haplogrep 3",
            anchor="haplogrep3",
            href="https://github.com/genepi/haplogrep3",
            info="Classifies mitochondrial DNA sequences into haplogroups.",
            doi=["10.1093/nar/gkad284"],
            license="MIT License",
            license_url="https://github.com/genepi/haplogrep3/blob/main/LICENSE",
        )

        self.haplogrep3_data: dict[str, dict[str, str | float | int]] = {}

        for f in self.find_log_files("haplogrep3"):
            for s_name, data in self._parse_report(f).items():
                if s_name in self.haplogrep3_data:
                    log.warning(
                        f"Sample '{s_name}' is in more than one Haplogrep 3 report. "
                        f"The result from '{f['fn']}' replaces the earlier result."
                    )
                self.add_data_source(f, s_name=s_name)
                self.add_software_version(None, s_name)
                self.haplogrep3_data[s_name] = data

        self.haplogrep3_data = self.ignore_samples(self.haplogrep3_data)

        if len(self.haplogrep3_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.haplogrep3_data)} samples")

        self._add_general_stats()
        self._add_table()

        self.write_data_file(self.haplogrep3_data, "multiqc_haplogrep3")

    def _parse_report(self, f: LoadedFileDict[str]) -> dict[str, dict[str, str | float | int]]:
        parsed: dict[str, dict[str, str | float | int]] = {}
        # The csv module also reads a quote character in a field, which Haplogrep 3 writes as two quote characters.
        for row in csv.DictReader(f["f"].splitlines(), delimiter="\t"):
            if row["Rank"] != "1":
                continue
            data: dict[str, str | float | int] = {
                "haplogroup": row["Haplogroup"],
                "quality": float(row["Quality"]),
                "range": row["Range"],
            }
            if "Not_Found_Polys" in row:
                data["not_found_polys"] = len(row["Not_Found_Polys"].split())
                # Each polymorphism has a category in parentheses after it, for example "16519C (hotspot)".
                data["remaining_polys"] = sum(1 for p in row["Remaining_Polys"].split() if not p.startswith("("))
            parsed[row["SampleID"]] = data

        if not parsed:
            log.debug(f"No samples in the Haplogrep 3 report '{f['fn']}'")
        return parsed

    def _add_general_stats(self):
        headers: dict[str, ColumnDict] = {
            "haplogroup": {
                "title": "Haplogroup",
                "description": "Mitochondrial haplogroup with the best quality",
                "scale": False,
            },
            "quality": {
                "title": "Haplogroup quality",
                "description": "Quality of the haplogroup classification, from 0 to 1",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
                "format": "{:,.3f}",
                "hidden": True,
            },
        }
        self.general_stats_addcols(self.haplogrep3_data, self.get_general_stats_headers(all_headers=headers))

    def _add_table(self):
        headers: dict[str | ColumnKey, ColumnDict] = {
            "haplogroup": {
                "title": "Haplogroup",
                "description": "Mitochondrial haplogroup with the best quality (rank 1)",
                "scale": False,
            },
            "quality": {
                "title": "Quality",
                "description": "Quality of the haplogroup classification, from 0 (bad) to 1 (good)",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
                "format": "{:,.3f}",
            },
            "range": {
                "title": "Range",
                "description": "Positions of the mitochondrial genome that the input contains",
                "scale": False,
                "hidden": True,
            },
            "not_found_polys": {
                "title": "Not found",
                "description": "Number of polymorphisms that the haplogroup expects and that the sample does not contain",
                "min": 0,
                "scale": "Oranges",
                "format": "{:,d}",
            },
            "remaining_polys": {
                "title": "Remaining",
                "description": "Number of polymorphisms in the sample that the haplogroup does not explain",
                "min": 0,
                "scale": "Purples",
                "format": "{:,d}",
            },
        }
        self.add_section(
            name="Haplogroups",
            anchor="haplogrep3-haplogroups",
            description="Mitochondrial haplogroup of each sample, from <code>haplogrep3 classify</code>.",
            helptext="""
            Haplogrep 3 compares the polymorphisms of each sample with a phylogenetic tree of
            mitochondrial haplogroups. The table shows the haplogroup with the best quality.

            The quality is a value from 0 to 1. A low quality shows that the sample does not
            agree well with the haplogroup. Possible causes are low coverage, contamination or
            a sample that contains only a part of the mitochondrial genome.

            The _Not found_ and _Remaining_ columns are only available for a report that
            Haplogrep 3 wrote with `--extend-report`.
            """,
            plot=table.plot(
                self.haplogrep3_data,
                headers,
                pconfig={
                    "id": "haplogrep3-haplogroups-table",
                    "title": "Haplogrep 3: Haplogroups",
                },
            ),
        )
