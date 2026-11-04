"""Build synthetic I/O fixtures from a reference; NOT an accuracy benchmark.

The planted OR exactly matches an input reference, so this deliberately tests
plumbing, coordinates, both strands and exports rather than discovery novelty.
"""
import argparse
from pathlib import Path
import random

from orsearcher.defaults import reference_directory
from orsearcher.fasta import CODON_TABLE_1, FastaRecord, read_fasta, reverse_complement, write_fasta


def make(reference, output):
    proteins=[r for r in read_fasta(Path(reference)) if r.sequence.startswith('M')
              and 300<=len(r.sequence)<=500 and set(r.sequence)<=set(CODON_TABLE_1.values())-{'*'}]
    if not proteins:
        raise ValueError('No suitable intact reference for synthetic fixture')
    record=proteins[0]
    codons={}
    for codon,aa in sorted(CODON_TABLE_1.items()):
        codons.setdefault(aa,codon)
    dna=''.join(codons[aa] for aa in record.sequence)+'TAA'
    middle=(len(record.sequence)//2)*3
    rng=random.Random(20261028)
    flank=lambda: ''.join(rng.choice('ACGT') for _ in range(500))
    spliced=dna[:middle]+'GT'+'A'*196+'AG'+dna[middle:]
    write_fasta(Path(output),[
        FastaRecord('synthetic_plus','synthetic reference-derived positive; NOT independent truth',flank()+spliced+flank()),
        FastaRecord('synthetic_minus','reverse-strand synthetic positive',flank()+reverse_complement(spliced)+flank()),
        FastaRecord('synthetic_background','random short negative fixture',flank()+flank())])
    print(f'Synthetic plumbing fixture written; planted reference={record.id}, protein_aa={len(record.sequence)}')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference',type=Path,default=reference_directory()/'complete_OR.faa')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    make(args.reference,args.out)
