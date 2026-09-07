"""Validate an ingested study without changing it."""
import argparse
import anndata
from epilepsy_microglia.validation import STUDIES, validate_adata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path")
    parser.add_argument("--study", required=True, choices=sorted(STUDIES))
    args = parser.parse_args()
    errors = validate_adata(anndata.read_h5ad(args.path), study=args.study)
    if errors:
        parser.exit(1, "\n".join(errors) + "\n")
    print("Validation passed")


if __name__ == "__main__":
    main()
