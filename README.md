# Bugphyzz

Bugphyzz loads published microbial trait annotations and creates trait signatures: groups of taxa that share a trait. It runs independently of `bugphyzz_reconstruction`, the separate package that prepares data and produces release files.

## Install

Use Python 3.10 or newer. Open a terminal in this project folder and run:

```powershell
python -m pip install -e .
```

## Load annotations

```python
from bugphyzz import import_bugphyzz

traits = import_bugphyzz()
print(sorted(traits))
print(traits["motility"].head())
```

The result is a dictionary of pandas DataFrames, one for each trait. The first run downloads the published dataset and requires an internet connection. Later runs reuse the local cache.

By default, the importer excludes annotations with frequency `rarely` and retains reconstructed annotations with validation values of at least 0.8. To change those settings:

```python
traits = import_bugphyzz(v=0.9, exclude_rarely=False)
```

## Create signatures

```python
from bugphyzz import make_signatures

signatures = make_signatures(
    traits["motility"],
    tax_id_type="NCBI_ID",
    tax_level="species",
    min_size=10,
)
print(signatures)
```

Each signature maps a trait label to a list of taxa. Use `Taxon_name` instead of `NCBI_ID` for taxon names, or `mixed` instead of `species` to include multiple ranks. A signature must contain at least `min_size` unique taxa; restrictive filters can produce no signatures.

Numeric signatures use the bundled `resources/thresholds.tsv`. Import filtering uses the bundled `resources/validation_summary.tsv`.

This package reads released annotations; running ancestral-state reconstruction and generating release CSV/GMT files is handled by `bugphyzz_reconstruction`.
