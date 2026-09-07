"""Pappalardo ingestion entry point, pending curated sources and inspection.

Implement only this study's explicit file layout, matrix orientation and
metadata joins here. Use epilepsy_microglia.metadata and io.write_study.
Do not download, mutate raw files, infer clinical labels or transform counts.
"""
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    parser.exit(2, "pappalardo: awaiting curated URLs and file inspection; no files written.\n")


if __name__ == "__main__":
    main()
