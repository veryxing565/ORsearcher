"""Parameter overrides for the four single-run presets."""
from pathlib import Path
import re


DISCOVERY = dict(first_seed_hits=5, second_seed_hits=2, first_output_coverage=.20,
    second_output_coverage=.10, first_alignment_ratio=.25, second_alignment_ratio=.15,
    alignment_output_score=.15, candidate_min_aa=60, candidate_bits=30., candidate_query_coverage=.25)

RECONSTRUCTION = dict(adaptive_require_full_reference=False, adaptive_reference_coverage=.65,
    evidence_or_targets=100, evidence_gr_targets=20, identity_bits=45., identity_min_identity=20.,
    identity_query_coverage=.40, distant_bits=30., distant_query_coverage=.30, identity_use_tree=True,
    reconstruction_reference_policy='active-external', structure_max_loci=1000, structure_references=10,
    structure_rescue_budget=128, recall_max_regions=1000, recall_complementary_regions=200,
    translated_max_regions=7500, polish_repair_budget=512, repair_hit_bits=60., repair_hit_evalue=1e-8,
    repair_hit_identity=25., model_hit_evalue=1e-12, model_hit_bits=80., model_hit_identity=25.,
    model_hit_coverage=.65, partial_hit_identity=30.)

COMPLETENESS = dict(complete_evidence='external-homology', complete_hit_bits=100.,
    complete_query_coverage=.75, complete_reference_coverage=.80, complete_terminal_slack_aa=35.,
    complete_terminal_slack_fraction=.15)

PRESETS = {'G0': {}, 'G1': DISCOVERY, 'G2': DISCOVERY | RECONSTRUCTION,
           'G3': DISCOVERY | RECONSTRUCTION | COMPLETENESS}


def set_run_defaults(args):
    if args.sample is None:
        name = args.genome.name if args.genome else 'batch'
        name = re.sub(r'\.(?:gz|bgz)$', '', name, flags=re.I)
        name = re.sub(r'\.(?:fa|fna|fasta|fas)$', '', name, flags=re.I)
        args.sample = re.sub(r'[^A-Za-z0-9_.-]+', '_', name).strip('._-')[:80] or 'sample'
    if args.outdir is None:
        args.outdir = Path('orsearcher_' + args.sample)
    if args.jobs is None:
        args.jobs = min(4, args.threads)
    if args.workers is None:
        args.workers = min(18, args.threads)


def apply_preset(args, tokens, run_parser):
    options = run_parser._option_string_actions
    explicit = {options[flag].dest for token in tokens
                if (flag := token.split('=', 1)[0]) in options}
    for name, value in PRESETS[args.preset].items():
        if name not in explicit:
            setattr(args, name, value)
