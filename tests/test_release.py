"""Release CLI/reference contracts; do not require bioinformatics binaries."""
import json
import logging
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from orsearcher.cli import parser
from orsearcher.defaults import FILES, inventory, resolve_references
from orsearcher.logging_utils import ConciseFormatter
from orsearcher.engine import load_samples,discovery_reports
from orsearcher.provenance import verify_previous_genome


class ReferenceTests(unittest.TestCase):
    def bundle(self, root):
        for names in FILES.values():
            for name in names if isinstance(names, list) else [names]:
                (root / name).write_text('fixture\n')

    def args(self, root, **changes):
        args=parser().parse_args(['run','--genome','genome.fna','--outdir','new',
                                  '--reference-dir',str(root)])
        vars(args).update(changes)
        return args

    def test_defaults_are_portable(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.bundle(root);args=self.args(root)
            resolve_references(args)
            self.assertEqual(args.or_reference,root/'OR.faa')
            self.assertEqual(args.gr_hmm,[root/'GR_1.hmm',root/'GR_2.hmm'])

    def test_custom_or_does_not_inherit_full_reference_subset(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.bundle(root);custom=root/'new.faa';custom.write_text('new')
            args=self.args(root,or_reference=custom)
            resolve_references(args)
            self.assertIsNone(args.complete_reference)
            self.assertEqual(args.or_reference,custom)

    def test_gr_hmm_replacement_is_not_append(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.bundle(root);custom=root/'new.hmm';custom.write_text('new')
            args=self.args(root,gr_hmm=[custom])
            resolve_references(args)
            self.assertEqual(args.gr_hmm,[custom])

    def test_missing_defaults_fail_with_actionable_message(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'reference-dir'):
                resolve_references(self.args(Path(d)))

    def test_disable_background(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);self.bundle(root);(root/'background.faa').unlink()
            args=self.args(root,no_background=True);resolve_references(args)
            self.assertIsNone(args.background_reference)

    def test_inventory_does_not_require_full_bundle(self):
        with tempfile.TemporaryDirectory() as d:
            report=inventory(d)
            self.assertFalse(report['files']['OR.faa']['exists'])

    def test_polish_is_separate_subcommand(self):
        args=parser().parse_args(['polish','--genome','g.fna','--outdir','new',
                                 '--previous-sample','old','--polish-repair-budget','0'])
        self.assertEqual(args.command,'polish')
        self.assertEqual(args.polish_repair_budget,0)

    def test_manifest_paths_are_relative_to_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);manifest=root/'genomes.tsv'
            manifest.write_text('sample\tgenome\tannotation\na\tgenomes/a.fna\t\n')
            args=SimpleNamespace(manifest=manifest)
            rows=load_samples(args)
            self.assertEqual(rows[0]['genome'],str((root/'genomes/a.fna').resolve()))
            self.assertEqual(rows[0]['annotation'],'')

    def test_discovery_deferred_is_not_silently_completed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder=root/'discovery/recall';folder.mkdir(parents=True)
            (folder/'summary.json').write_text(json.dumps({'status':'completed_with_warnings','warnings':['genewise_budget_exceeded']}))
            reports,warnings=discovery_reports(root)
            self.assertEqual(warnings,['recall:genewise_budget_exceeded'])
            self.assertIn('recall',reports)


class LoggingTests(unittest.TestCase):
    def test_concise_log_suppresses_traceback_without_losing_debug_data(self):
        try:
            raise ValueError('fixture failure')
        except ValueError:
            import sys
            record=logging.LogRecord('test',logging.ERROR,'',1,'Failed: fixture',(),sys.exc_info())
        result=ConciseFormatter('%(message)s').format(record)
        self.assertEqual(result,'Failed: fixture')
        self.assertIsNotNone(record.exc_info)
        self.assertIn('ValueError',logging.Formatter().format(record))


class ProvenanceTests(unittest.TestCase):
    def test_missing_previous_manifest_is_not_silently_trusted(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            with self.assertRaisesRegex(ValueError,'unchanged original genome'):
                verify_previous_genome(root/'old/samples/a',root/'genome.fna')
            report=verify_previous_genome(root/'old/samples/a',root/'genome.fna',True)
            self.assertFalse(report['verified'])

    def test_changed_genome_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);genome=root/'genome.fna';genome.write_text('>a\nAAA\n')
            old=root/'old';old.mkdir()
            stat=genome.stat()
            manifest={'inputs':[{'path':str(genome),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns}]}
            (old/'input_manifest.json').write_text(json.dumps(manifest))
            self.assertTrue(verify_previous_genome(old/'samples/a',genome)['verified'])
            genome.write_text('>a\nAAAAAA\n')
            with self.assertRaisesRegex(ValueError,'metadata_changed'):
                verify_previous_genome(old/'samples/a',genome)


if __name__=='__main__':
    unittest.main()
