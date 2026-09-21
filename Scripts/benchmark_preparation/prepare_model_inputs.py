from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'outputs'
DESTINATION = ROOT / 'model_inputs'
MODELS = ('deepkla', 'autokla', 'hybridkla', 'pcbert')


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def validate(data: pd.DataFrame, name: str) -> None:
    required = {'sample_id', 'species', 'uniprot_id', 'position', 'window_51', 'label'}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f'{name}: missing columns {sorted(missing)}')
    if data['sample_id'].duplicated().any():
        raise ValueError(f'{name}: duplicate sample_id values')
    if not data['window_51'].astype(str).str.len().eq(51).all():
        raise ValueError(f'{name}: non-51-residue window')
    if not data['window_51'].astype(str).str[25].eq('K').all():
        counts = data['window_51'].astype(str).str[25].value_counts().to_dict()
        raise ValueError(f'{name}: center residues are not all K: {counts}')
    if set(data['label'].astype(int)) - {0, 1}:
        raise ValueError(f'{name}: invalid labels')
    conflicts = data.groupby('window_51')['label'].nunique().gt(1).sum()
    if conflicts:
        raise ValueError(f'{name}: {conflicts} exact windows have conflicting labels')


def write_view(source: Path, destination: Path, manifest: list[dict]) -> None:
    data = pd.read_csv(source, sep='\t', dtype={'position': str})
    validate(data, str(source))
    data.insert(data.columns.get_loc('window_51') + 1, 'window', data['window_51'])
    destination.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(destination, sep='\t', index=False)
    manifest.append({
        'source': str(source.relative_to(ROOT)),
        'input': str(destination.relative_to(ROOT)),
        'sha256': checksum(destination),
        'rows': len(data),
        'positive': int(data['label'].sum()),
        'putative_negative': int(data['label'].eq(0).sum()),
        'terminal_padded': int(data['terminal_padded'].astype(str).str.lower().eq('true').sum()),
    })


def write_delta(common_source: Path, specific_source: Path, destination: Path,
                manifest: list[dict]) -> None:
    common_ids = set(pd.read_csv(common_source, sep='\t', usecols=['sample_id'])['sample_id'])
    data = pd.read_csv(specific_source, sep='\t', dtype={'position': str}, low_memory=False)
    data = data[~data['sample_id'].isin(common_ids)].copy()
    validate(data, str(specific_source) + ' delta')
    data.insert(data.columns.get_loc('window_51') + 1, 'window', data['window_51'])
    destination.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(destination, sep='\t', index=False)
    manifest.append({
        'source': str(specific_source.relative_to(ROOT)),
        'input': str(destination.relative_to(ROOT)),
        'kind': 'model_specific_minus_common_strict',
        'sha256': checksum(destination), 'rows': len(data),
        'positive': int(data['label'].sum()),
        'putative_negative': int(data['label'].eq(0).sum()),
    })


def main() -> None:
    manifest = []
    common = SOURCE / 'common_strict' / 'identity_70' / 'all_samples.tsv'
    for model in MODELS:
        write_view(common, DESTINATION / 'common_strict' / model / 'all.tsv', manifest)
        specific = SOURCE / 'model_specific' / model / 'identity_70' / 'all_samples.tsv'
        write_view(specific, DESTINATION / 'model_specific' / model / 'all.tsv', manifest)
        write_delta(
            common, specific,
            DESTINATION / 'model_specific_delta' / model / 'delta.tsv', manifest,
        )
        sensitivity_90 = SOURCE / 'common_strict' / 'identity_90' / 'all_samples.tsv'
        write_delta(
            common, sensitivity_90,
            DESTINATION / 'identity_sensitivity_delta' / model / 'delta_90.tsv',
            manifest,
        )
    (DESTINATION / 'input_manifest.json').write_text(
        json.dumps(manifest, indent=2), encoding='utf-8'
    )


if __name__ == '__main__':
    main()
