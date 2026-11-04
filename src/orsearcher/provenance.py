"""Guard reuse of assembly-backed structural/damage evidence."""
import json
from pathlib import Path


def verify_previous_genome(previous_sample, genome, allow_unverified=False):
    sample, genome = Path(previous_sample).resolve(), Path(genome).resolve()
    manifest=sample.parent.parent/'input_manifest.json'
    reason='previous_input_manifest_missing'
    if manifest.is_file():
        data=json.loads(manifest.read_text(encoding='utf-8'))
        match=next((item for item in data.get('inputs',[]) if Path(item['path']).resolve()==genome),None)
        reason='genome_not_recorded_in_previous_run'
        if match is not None:
            stat=genome.stat()
            if stat.st_size==match['size'] and stat.st_mtime_ns==match['mtime_ns']:
                return {'verified':True,'check':'same_resolved_path_size_mtime','assembly_sha256_verified':False}
            reason='genome_metadata_changed_since_structural_review'
    if not allow_unverified:
        raise ValueError(reason+': polish requires the unchanged original genome. '
                         'Use refine for a different assembly; --allow-unverified-previous explicitly accepts missing provenance, not proof of equivalence.')
    return {'verified':False,'reason':reason,'assembly_sha256_verified':False}
