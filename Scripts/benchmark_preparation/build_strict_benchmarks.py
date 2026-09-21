from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd


MODELS = ('deepkla', 'autokla', 'hybridkla', 'pcbert')
THRESHOLDS = (0.70, 0.80, 0.90)
TRAINING_FILES = {
    'deepkla': 'deepkla_train.fa',
    'autokla': 'autokla_train.fa',
    'hybridkla': 'hybridkla_train.csv',
    'pcbert': 'pcbert_train.csv',
}
CDHIT_IMAGE = 'quay.io/biocontainers/cd-hit:4.8.1--h43eeafb_10'
OUTPUT_COLUMNS = [
    'sample_id', 'species', 'uniprot_id', 'position', 'window_51', 'label',
    'label_name', 'source_pmids', 'source_record', 'terminal_padded',
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def slug(value: str) -> str:
    return re.sub(r'[^a-z0-9]+', '_', value.lower()).strip('_')


def word_size(identity: float) -> int:
    return 5 if identity >= 0.70 else 4 if identity >= 0.60 else 3


def normalize_sequence(value: str) -> str:
    return re.sub(r'\s+', '', str(value)).upper().replace('*', 'X')


def run(command: list[str], cwd: Path, commands: list[dict]) -> None:
    commands.append({'cwd': str(cwd), 'command': command})
    result = subprocess.run(
        command, cwd=cwd, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            + '\n'.join(result.stdout.splitlines()[-50:])
        )


def docker_cdhit(tool: str, work: Path, args: list[str], commands: list[dict]) -> None:
    docker = shutil.which('docker')
    if not docker:
        raise RuntimeError('Docker is required for the pinned CD-HIT image')
    run([
        docker, 'run', '--rm', '-v', f'{work.resolve()}:/work', '-w', '/work',
        CDHIT_IMAGE, tool, *args,
    ], work, commands)


def read_fasta(path: Path) -> list[str]:
    sequences = []
    current = []
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            if current:
                sequences.append(normalize_sequence(''.join(current)))
                current = []
        elif line.strip():
            current.append(line.strip())
    if current:
        sequences.append(normalize_sequence(''.join(current)))
    return sequences


def read_training(path: Path, model: str) -> list[str]:
    if model in {'deepkla', 'autokla'}:
        values = read_fasta(path)
    elif model == 'hybridkla':
        frame = pd.read_csv(path)
        values = [normalize_sequence(value) for value in frame['Sequence']]
    else:
        lines = path.read_text(encoding='utf-8-sig').splitlines()[1:]
        values = [normalize_sequence(lines[index + 1])
                  for index in range(0, len(lines) - 1, 2)
                  if lines[index].startswith('Protein ')]
    return sorted({value for value in values if value})


def read_candidates(input_dir: Path) -> pd.DataFrame:
    files = sorted(input_dir.glob('*.xlsx'))
    frames = []
    for path in files:
        frame = pd.read_excel(path, dtype={'uniprot_id': str, 'position': str})
        required = {'uniprot_id', 'position', 'label', 'window_seq', 'species'}
        if not required.issubset(frame.columns):
            continue
        frame = frame[list(required)].copy()
        frame['source_record'] = path.name
        frames.append(frame)
    if not frames:
        raise RuntimeError('No candidate spreadsheets were found')
    data = pd.concat(frames, ignore_index=True)
    data = data.rename(columns={'window_seq': 'window_51'})
    for column in ('uniprot_id', 'position', 'species', 'window_51'):
        data[column] = data[column].astype(str).str.strip()
    data['window_51'] = data['window_51'].str.upper()
    data['label'] = pd.to_numeric(data['label'], errors='raise').astype(int)
    if set(data['label']) - {0, 1}:
        raise ValueError('Candidate labels must be 0 or 1')
    if not data['window_51'].str.len().eq(51).all():
        raise ValueError('Every candidate window must contain exactly 51 characters')
    data['terminal_padded'] = data['window_51'].str.contains(r'\*')
    return data


def attach_provenance(data: pd.DataFrame, map_path: Path) -> pd.DataFrame:
    source = pd.read_csv(map_path, dtype=str).fillna('')
    source['position'] = source['position'].str.strip()
    source['pmid'] = source['pmid'].str.replace('411925556', '41192556', regex=False)
    grouped = source.groupby(['uniprot_id', 'position', 'species'])['pmid'].apply(
        lambda values: ','.join(sorted({p.strip() for value in values
                                        for p in value.split(',') if p.strip()}))
    ).rename('source_pmids').reset_index()
    result = data.merge(grouped, on=['uniprot_id', 'position', 'species'], how='left')
    result['source_pmids'] = result['source_pmids'].fillna('')
    result.loc[result['label'].eq(0), 'source_pmids'] = ''
    result.loc[result['label'].eq(0), 'source_record'] = 'candidate_unreported_lysine'
    return result


def stable_id(row: pd.Series) -> str:
    value = '|'.join(str(row[column]) for column in
                     ('species', 'uniprot_id', 'position', 'window_51', 'label'))
    return 'kla_' + hashlib.sha256(value.encode('utf-8')).hexdigest()[:20]


def remove_exact_ambiguity(data: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    labels = data.groupby('window_51')['label'].nunique()
    ambiguous = set(labels[labels.gt(1)].index)
    filtered = data[~data['window_51'].isin(ambiguous)].copy()
    before_dedup = len(filtered)
    filtered = filtered.sort_values(
        ['species', 'uniprot_id', 'position', 'source_record'], kind='stable'
    ).drop_duplicates(['window_51', 'label'], keep='first')
    return filtered, {
        'exact_conflicting_windows': len(ambiguous),
        'rows_removed_exact_conflict': int(data['window_51'].isin(ambiguous).sum()),
        'same_label_exact_duplicates_removed': before_dedup - len(filtered),
    }


def write_fasta(rows: pd.DataFrame, path: Path) -> None:
    with path.open('w', encoding='ascii', newline='\n') as handle:
        for row in rows.itertuples():
            handle.write(f'>{row.sample_id}\n{normalize_sequence(row.window_51)}\n')


def parse_clusters(path: Path) -> list[list[str]]:
    clusters = []
    current = []
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('>Cluster'):
            if current:
                clusters.append(current)
            current = []
            continue
        match = re.search(r'>([^.\s]+)\.\.\.', line)
        if match:
            current.append(match.group(1))
    if current:
        clusters.append(current)
    return clusters


def cluster_candidates(data: pd.DataFrame, identity: float, work: Path,
                       commands: list[dict]) -> tuple[pd.DataFrame, dict]:
    input_path = work / 'candidate.fa'
    output_path = work / 'candidate.clustered.fa'
    write_fasta(data, input_path)
    docker_cdhit('cd-hit', work, [
        '-i', input_path.name, '-o', output_path.name,
        '-c', f'{identity:.2f}', '-n', str(word_size(identity)),
        '-aS', '0.80', '-aL', '0.80', '-g', '1', '-d', '0', '-T', '0', '-M', '0',
    ], commands)
    by_id = data.set_index('sample_id')
    keep = []
    ambiguous_clusters = 0
    ambiguous_rows = 0
    for members in parse_clusters(Path(str(output_path) + '.clstr')):
        labels = set(by_id.loc[members, 'label'].tolist())
        if len(labels) > 1:
            ambiguous_clusters += 1
            ambiguous_rows += len(members)
            continue
        keep.append(members[0])
    result = by_id.loc[keep].reset_index()
    return result, {
        'homology_clusters': len(keep) + ambiguous_clusters,
        'mixed_label_clusters_removed': ambiguous_clusters,
        'mixed_label_rows_removed': ambiguous_rows,
        'same_label_cluster_redundancy_removed': len(data) - ambiguous_rows - len(result),
    }


def write_training_fasta(sequences: list[str], path: Path, prefix: str) -> None:
    with path.open('w', encoding='ascii', newline='\n') as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f'>{prefix}_{index}\n{sequence}\n')


def filter_against_training(data: pd.DataFrame, sequences: list[str], identity: float,
                            work: Path, name: str,
                            commands: list[dict]) -> pd.DataFrame:
    train_path = work / f'{name}.train.fa'
    query_path = work / f'{name}.candidate.fa'
    output_path = work / f'{name}.unmatched.fa'
    write_training_fasta(sequences, train_path, name)
    write_fasta(data, query_path)
    docker_cdhit('cd-hit-2d', work, [
        '-i', train_path.name, '-i2', query_path.name, '-o', output_path.name,
        '-c', f'{identity:.2f}', '-n', str(word_size(identity)),
        '-aS', '0.80', '-aL', '0.80', '-g', '1', '-d', '0', '-T', '0', '-M', '0',
    ], commands)
    unmatched = {line[1:].strip().split()[0]
                 for line in output_path.read_text(encoding='utf-8').splitlines()
                 if line.startswith('>')}
    return data[data['sample_id'].isin(unmatched)].copy()


def save_dataset(data: pd.DataFrame, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    ordered = data[OUTPUT_COLUMNS].sort_values(
        ['species', 'label', 'uniprot_id', 'position', 'sample_id'], kind='stable'
    )
    ordered.to_csv(directory / 'all_samples.tsv', sep='\t', index=False,
                   quoting=csv.QUOTE_MINIMAL)
    for species, group in ordered.groupby('species'):
        group.to_csv(directory / f'{slug(species)}.tsv', sep='\t', index=False)
    ordered[~ordered['terminal_padded']].to_csv(
        directory / 'all_samples_without_terminal_padding.tsv', sep='\t', index=False
    )


def count_row(stage: str, threshold: int, data: pd.DataFrame, scope='all') -> dict:
    counts = data['label'].value_counts().to_dict()
    return {
        'identity_threshold': threshold, 'scope': scope, 'stage': stage,
        'positive': int(counts.get(1, 0)), 'putative_negative': int(counts.get(0, 0)),
        'total': len(data), 'terminal_padded': int(data['terminal_padded'].sum()),
        'unique_proteins': data['uniprot_id'].nunique(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', type=Path, default=Path('inputs'))
    parser.add_argument('--training-dir', type=Path, default=Path('training'))
    parser.add_argument('--output-dir', type=Path, default=Path('outputs'))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    commands: list[dict] = []
    raw = attach_provenance(
        read_candidates(args.input_dir), args.input_dir / 'site_source_map.csv'
    )
    raw['sample_id'] = raw.apply(stable_id, axis=1)
    raw['label_name'] = raw['label'].map({1: 'experimentally_reported_positive',
                                          0: 'putative_unlabeled_negative'})
    exact, exact_audit = remove_exact_ambiguity(raw)

    training = {}
    training_manifest = []
    for model, filename in TRAINING_FILES.items():
        path = args.training_dir / filename
        training[model] = read_training(path, model)
        training_manifest.append({
            'model': model, 'source_file': filename, 'sha256': sha256(path),
            'unique_training_sequences': len(training[model]),
        })
    union_training = sorted({seq for values in training.values() for seq in values})

    audit = [count_row('raw_candidates', 0, raw), count_row('after_exact_rules', 0, exact)]
    details = {'exact_rules': exact_audit}
    with tempfile.TemporaryDirectory(prefix='kladb_strict_') as temp:
        temp_root = Path(temp)
        cluster_work = temp_root / 'fixed_global_cluster_70'
        cluster_work.mkdir()
        clustered, cluster_audit = cluster_candidates(
            exact, 0.70, cluster_work, commands
        )
        details['fixed_global_cluster_70'] = cluster_audit
        audit.append(count_row('after_global_label_aware_clustering', 70, clustered))
        previous_common = None
        previous_specific = {model: None for model in MODELS}
        nesting_adjustments = []
        for identity in THRESHOLDS:
            pct = int(identity * 100)
            work = temp_root / str(pct)
            work.mkdir()
            common = filter_against_training(
                clustered, union_training, identity, work, 'training_union', commands
            )
            if previous_common is not None:
                missing = previous_common[
                    ~previous_common['sample_id'].isin(set(common['sample_id']))
                ]
                if len(missing):
                    common = pd.concat([common, missing], ignore_index=True)
                nesting_adjustments.append({
                    'scope': 'common_strict', 'identity_threshold': pct,
                    'rows_restored_to_enforce_nested_sensitivity_set': len(missing),
                })
            previous_common = common.copy()
            save_dataset(common, args.output_dir / 'common_strict' / f'identity_{pct}')
            audit.append(count_row('common_strict', pct, common))
            for species, group in common.groupby('species'):
                audit.append(count_row('common_strict', pct, group, species))

            for model in MODELS:
                specific = filter_against_training(
                    clustered, training[model], identity, work, model, commands
                )
                if previous_specific[model] is not None:
                    missing = previous_specific[model][
                        ~previous_specific[model]['sample_id'].isin(set(specific['sample_id']))
                    ]
                    if len(missing):
                        specific = pd.concat([specific, missing], ignore_index=True)
                    nesting_adjustments.append({
                        'scope': model, 'identity_threshold': pct,
                        'rows_restored_to_enforce_nested_sensitivity_set': len(missing),
                    })
                previous_specific[model] = specific.copy()
                save_dataset(
                    specific,
                    args.output_dir / 'model_specific' / model / f'identity_{pct}',
                )
                audit.append(count_row('model_specific', pct, specific, model))
        details['sensitivity_nesting_adjustments'] = nesting_adjustments

    pd.DataFrame(audit).to_csv(args.output_dir / 'filtering_audit.csv', index=False)
    (args.output_dir / 'training_manifest.json').write_text(
        json.dumps(training_manifest, indent=2), encoding='utf-8'
    )
    (args.output_dir / 'build_manifest.json').write_text(json.dumps({
        'dataset_name': 'KlaDB external leakage-screened benchmark',
        'positive_definition': 'Experimentally reported Kla site',
        'negative_definition': 'Putative/unlabeled lysine not reported as Kla',
        'coverage_rule': 'Bidirectional alignment coverage >= 80%',
        'primary_identity_threshold': 0.70,
        'sensitivity_identity_thresholds': [0.80, 0.90],
        'sensitivity_nesting_rule': (
            'Each relaxed-threshold set explicitly includes the preceding stricter set; '
            'this removes rare non-monotonic omissions caused by CD-HIT heuristic prefilters.'
        ),
        'terminal_padding_symbol': '*',
        'candidate_source_sha256': {
            path.name: sha256(path) for path in sorted(args.input_dir.glob('*.xlsx'))
        },
        'training_union_unique_sequences': len(union_training),
        'details': details,
        'commands': commands,
    }, indent=2), encoding='utf-8')


if __name__ == '__main__':
    main()
