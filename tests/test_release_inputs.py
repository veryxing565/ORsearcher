import gzip
from pathlib import Path
import tempfile
import unittest

from orsearcher.cli import parser
from orsearcher.defaults import reference_directory, resolve_references
from orsearcher.fasta import read_fasta
from orsearcher.input_files import materialize_references, merge_profiles, user_proteins
from orsearcher.presets import apply_preset, set_run_defaults


class ReleaseInputsTests(unittest.TestCase):
    def args(self, *tokens):
        cli = parser()
        args = cli.parse_args(['run','--genome','example.fna',*tokens])
        apply_preset(args,tokens,cli._subparsers._group_actions[0].choices['run'])
        set_run_defaults(args)
        return args

    def test_genome_only_and_small_thread_budget(self):
        args = self.args()
        self.assertEqual((args.preset,args.threads,args.sample,args.outdir),('G0',36,'example',Path('orsearcher_example')))
        self.assertEqual((self.args('--threads','1').jobs,self.args('--threads','1').workers),(1,1))

    def test_preset_and_explicit_value_precedence(self):
        args = self.args('--candidate-min-aa=70','--preset','g3','--no-identity-use-tree')
        self.assertEqual((args.preset,args.candidate_min_aa,args.structure_max_loci,args.complete_evidence),
                         ('G3',70,1000,'external-homology'))
        self.assertFalse(args.identity_use_tree)

    def test_replacement_plus_addition_infers_templates_and_normalizes_ids(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            original = next(read_fasta(reference_directory()/'complete_OR.faa')).sequence
            replace = folder/'new'
            replace.write_text(f'>x first\n{original}\n>x second\n{original[20:160]}\n')
            extra = folder/'addition.txt'
            with gzip.open(extra,'wt') as out:
                out.write(f'>x repeated\n{original}\n>y third\n{original[30:160]}\n')
            args = self.args('--replace-or',str(replace),'--add-or',str(extra))
            resolve_references(args)
            self.assertIsNone(args.complete_reference)
            materialize_references(args,folder)
            records = list(read_fasta(args.or_reference))
            self.assertEqual(len(records),3)
            self.assertEqual(len({r.id for r in records}),3)
            self.assertEqual({r.id for r in records},args.user_or_ids)
            self.assertEqual(args.gr_reference,reference_directory()/'GR.faa')

    def test_family_hmm_options_keep_other_role(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            default = reference_directory()/'OR.hmm'
            extra = folder/'model'
            extra.write_text(default.read_text().replace('NAME  7tm_6','NAME  second_OR'))
            self.assertEqual(len(merge_profiles([default,default,extra])),2)
            args = self.args('--replace-hmm','OR',str(extra),'--add-hmm','OR',str(default))
            resolve_references(args)
            materialize_references(args,folder)
            self.assertEqual(len(merge_profiles([args.or_hmm])),2)
            self.assertEqual(len(args.gr_hmm),2)

    def test_rejects_nucleotide_reference_and_wrong_hmm_content(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'file'
            path.write_text('>dna\n'+'ATGC'*30+'\n')
            with self.assertRaisesRegex(ValueError,'looks like DNA'):
                list(user_proteins(path))
            with self.assertRaisesRegex(ValueError,'HMMER3'):
                merge_profiles([path])
