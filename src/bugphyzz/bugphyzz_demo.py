"""Five-minute demonstration of the user-facing Python bugphyzz package.

This offline version never contacts Zenodo.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd

from bugphyzz import get_taxon_signatures, make_signatures


EXPORT_DIR = Path(r"C:\Users\delan\OneDrive\Desktop\original_exports")
EXPORT_FILES = (
    "bugphyzz_binary.csv",
    "bugphyzz_multistate.csv",
    "bugphyzz_numeric.csv",
)


DISPLAY_COLUMNS = [
    "Taxon_name",
    "Rank",
    "Attribute",
    "Attribute_value",
    "Evidence",
    "Frequency",
    "Score",
    "Validation",
]


def heading(number: int, title: str) -> None:
    """Print a visible divider between demo sections."""
    print(f"\n{number}. {title}")
    print("-" * (len(title) + 3))


def import_local_bugphyzz(
    export_dir: Path = EXPORT_DIR,
    *,
    exclude_rarely: bool = True,
) -> dict[str, pd.DataFrame]:
    """Load downloaded Bugphyzz CSV exports into trait-specific DataFrames."""
    missing = [name for name in EXPORT_FILES if not (export_dir / name).is_file()]
    if missing:
        formatted = "\n".join(f"  - {export_dir / name}" for name in missing)
        raise FileNotFoundError(f"Missing required Bugphyzz export files:\n{formatted}")

    traits: dict[str, pd.DataFrame] = {}
    for filename in EXPORT_FILES:
        # The first line is the Bugphyzz release header; row 2 has CSV columns.
        table = pd.read_csv(
            export_dir / filename,
            skiprows=1,
            dtype={"NCBI_ID": "string", "Taxon_name": "string"},
            low_memory=False,
        )
        table["Attribute"] = table["Attribute"].astype("string").str.lower()
        if exclude_rarely:
            table = table.loc[table["Frequency"] != "rarely"].copy()

        for attribute, group in table.groupby("Attribute", sort=False):
            traits[str(attribute)] = group.copy()

    return traits


def signature_ending_in(
    signatures: dict[str, list[str]] | None,
    attribute_value: str,
) -> tuple[str, list[str]]:
    """Find a discrete signature without depending on its full prefix."""
    if not signatures:
        raise RuntimeError("No signatures were created with these filters.")

    suffix = f"|{attribute_value.lower()}"
    for name, taxa in signatures.items():
        if name.lower().endswith(suffix):
            return name, taxa

    available = ", ".join(signatures)
    raise KeyError(f"Could not find {attribute_value!r}. Available: {available}")


def main() -> None:
    # : these are previously downloaded published outputs. Reading
    # them locally avoids rerunning ASR and removes the live Zenodo dependency.
    heading(1, "Import the downloaded bugphyzz annotations")
    bp = import_local_bugphyzz(exclude_rarely=True)
    print(f"Source: {EXPORT_DIR}")
    print(f"Imported {len(bp)} microbial traits.")
    print("A few available traits:")
    print(", ".join(sorted(bp)[:12]))

    # : every row connects a taxon with an attribute value and
    # preserves the evidence and support behind that association.
    heading(2, "Inspect one trait: aerophilicity")
    aerophilicity = bp["aerophilicity"]
    visible_columns = [c for c in DISPLAY_COLUMNS if c in aerophilicity.columns]
    print(aerophilicity[visible_columns].head(8).to_string(index=False))

    print("\nAnnotations by evidence type:")
    print(aerophilicity["Evidence"].value_counts().to_string())

    # : including ASR expands a signature to taxa without direct
    # observations. Users can also request a more conservative observed-only set.
    heading(3, "Create an aerobic microbial signature")
    all_evidence = make_signatures(
        aerophilicity,
        tax_id_type="Taxon_name",
        tax_level="species",
        min_size=1,
    )
    observed_only = make_signatures(
        aerophilicity,
        tax_id_type="Taxon_name",
        tax_level="species",
        evidence=("exp", "igc", "tas", "nas", "tax"),
        min_size=1,
    )

    aerobic_name, aerobic_taxa = signature_ending_in(all_evidence, "aerobic")
    _, observed_aerobic_taxa = signature_ending_in(observed_only, "aerobic")

    print(f"Signature: {aerobic_name}")
    print(f"With ASR annotations: {len(aerobic_taxa):,} taxa")
    print(f"With observed annotations only: {len(observed_aerobic_taxa):,} taxa")
    print("Example members:")
    print("\n".join(f"  - {taxon}" for taxon in aerobic_taxa[:8]))

    # : numeric traits become user-defined ranges rather than
    # discrete labels. Here the package finds taxa associated with 30-40 C.
    heading(4, "Create a numeric trait signature")
    temperature_signatures = make_signatures(
        bp["growth temperature"],
        tax_id_type="Taxon_name",
        tax_level="species",
        min_value=30,
        max_value=40,
        min_size=1,
    )
    temperature_name, temperature_taxa = next(iter(temperature_signatures.items()))
    print(f"Signature: {temperature_name}")
    print(f"Members: {len(temperature_taxa):,} taxa")
    print("Example members:")
    print("\n".join(f"  - {taxon}" for taxon in temperature_taxa[:8]))

    # : users can reverse the question and ask which signatures
    # contain a particular organism. Taxid 562 is Escherichia coli.
    heading(5, "Ask which signatures contain Escherichia coli")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        ecoli_signatures = get_taxon_signatures(
            "562",
            bp,
            tax_id_type="NCBI_ID",
            tax_level="species",
            min_size=1,
        )

    print(f"E. coli appears in {len(ecoli_signatures)} signatures.")
    print("A few examples:")
    print("\n".join(f"  - {name}" for name in ecoli_signatures[:10]))

    print("\nDemo complete: published annotations became inspectable records,")
    print("trait-based microbial sets, and a trait profile for one organism.")


if __name__ == "__main__":
    main()
