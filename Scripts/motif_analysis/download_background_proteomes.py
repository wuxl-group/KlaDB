"""Download the UniProt reference proteomes used as motif backgrounds."""

import argparse
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


PROTEOMES = {
    "Cavia_porcellus": "UP000005447",
    "Frankliniella_occidentalis": "UP000504606",
    "Glycine_max": "UP000008827",
    "Homo_sapiens": "UP000005640",
    "Mus_musculus": "UP000000589",
    "Oryza_sativa": "UP000059680",
    "Rattus_norvegicus": "UP000002494",
    "Sus_scrofa": "UP000008227",
    "Triticum_aestivum": "UP000019116",
}


def count_records(path):
    with path.open("rb") as handle:
        return sum(line.startswith(b">") for line in handle)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("background_proteomes"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for species, proteome in PROTEOMES.items():
        output = args.output_dir / f"{species}.fasta"
        if output.exists() and output.stat().st_size > 1000:
            print(f"skip {species}: {count_records(output)} records")
            continue
        url = "https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=" + quote(f"proteome:{proteome}")
        request = Request(url, headers={"User-Agent": "KlaDB-motif-analysis/1.0"})
        with urlopen(request, timeout=300) as response, output.open("wb") as handle:
            while block := response.read(1024 * 1024):
                handle.write(block)
        print(f"downloaded {species}: {count_records(output)} records")
        time.sleep(2)


if __name__ == "__main__":
    main()
