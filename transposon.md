---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3.7.12 64-bit (conda)
    language: python
    name: python3712jvsc74a57bd0dc0fe05456373cce17991fe9e2e9264df4f4d1e972d77814a9700f21c9e7a8e2
---

# Globals

```python
import fastq
import pysam
import gzip
from collections import Counter
```

```python
import os
import pandas as pd
import numpy as np
```

```python
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from pyfaidx import Fasta
```

```python
from Bio import motifs
```

```python
import logomaker
```

```python
import matplotlib.pyplot as plt
```

```python
def read_vcf(df_path):
    vcf_cols = ['chrom1', 'pos1', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample']
    df = pd.read_table(df_path, comment='#', names=vcf_cols, dtype={'chrom':str, 'pos1':float})
    return df
```

```python
def str_to_fasta(seq):
    fastaseq = ''
    length = len(seq)
    for start in range(0, length, 80):
        end = min(start+80, length)
        subseq = seq[start:end]
        fastaseq += f'{subseq}\n'
    return fastaseq
```

```python
def get_secondaries(read):
    secondaries = [a[1].split(';') for a in read.tags if a[0] == 'SA']
    if len(secondaries) > 0:
        secondaries = secondaries[0]
        secondaries = list(filter(lambda x: x!='', secondaries))
        return secondaries
    return []
```

```python
def get_secondary_aligned_seq(secondary, ncseq):
    query_consumers = 'MIS=X' # https://samtools.github.io/hts-specs/SAMv1.pdf
    query_aligners = 'MI=X' # removing clip from query consumers
    clen_sum = 0
    secondary_aligned_seq = ''
    tname, tpos, tstrand, tcigar, _, _ = secondary.split(',')
    cigar_tuple = re.findall('(\d+)([A-Z])', tcigar)
    for ix, (clen, ctag) in enumerate(cigar_tuple):
        print(clen, ctag)
        clen = int(clen)
        if ctag in query_aligners:
            aligned_seq = ncseq[clen_sum:clen_sum+clen]
            secondary_aligned_seq += aligned_seq
        if ctag in query_consumers:
            clen_sum += clen
    return (tname, tpos, tstrand, tcigar, secondary_aligned_seq)
```

## get fqs seq

```python
def read_fastq(path):
    """
    path: path to fastq[.gz]
    returns: dict[rname] = seq
    """
    rname_to_seq = {}
    rname_to_qual = {}
    fos = fastq.read(path)
    for fo in fos:
        rname = fo.head.split(' ')[0][1:]
        seq = fo.body
        qual = fo.qstr
        rname_to_seq[rname] = seq
        rname_to_qual[rname] = qual
    return rname_to_seq, rname_to_qual
```

## get refgenome

```python
def read_chromosome_lengths(genome_fasta_index):
    fai = pd.read_csv(genome_fasta_index, sep='\t', header=None, names=['chrom', 'length', 'V3', 'V4', 'V5'])
    fai = fai.set_index('chrom')['length']
    # print(fai.to_dict())
    return fai.to_dict()

class RefGenomeInfo(object):
    def __init__(self, version):
        if version == 'hg19':
            self.chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']

            genome_fasta_index = pkg_resources.resource_filename('wgs_analysis', 'data/GRCh37-lite.fa.fai')

            self.chromosome_lengths = pd.Series(read_chromosome_lengths(genome_fasta_index)).reindex(self.chromosomes).astype(int)
            self.chromosome_lengths.index.name = 'chr'

            self.chromosome_end = np.cumsum(self.chromosome_lengths)
            self.chromosome_start = self.chromosome_end.shift(1)
            self.chromosome_start[0] = 0
            self.chromosome_start = self.chromosome_start.astype(int)
            self.chromosome_mid = (self.chromosome_start + self.chromosome_end) / 2.

            self.chromosome_info = pd.DataFrame({
                'chromosome_length': self.chromosome_lengths,
                'chromosome_end': self.chromosome_end,
                'chromosome_start': self.chromosome_start,
                'chromosome_mid': self.chromosome_mid,
            }).reset_index()

        elif version == 'hg38':
            self.chromosomes = [str(a) for a in range(1, 23)] + ['X', 'Y']
            self.chromosomes = ['chr'+c for c in self.chromosomes]

            # genome_fasta_index = pkg_resources.resource_filename('wgs_analysis', 'data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai')
            genome_fasta_index = '/juno/work/shah/users/chois7/packages/wgs_analysis/wgs_analysis/data/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai'

            self.chromosome_lengths = pd.Series(read_chromosome_lengths(genome_fasta_index)).reindex(self.chromosomes).astype(int)
            self.chromosome_lengths.index.name = 'chr'
            # self.chromosome_lengths.index = self.chromosome_lengths.index.str.replace('chr', '')

            self.chromosome_end = np.cumsum(self.chromosome_lengths)
            self.chromosome_start = self.chromosome_end.shift(1)
            self.chromosome_start[0] = 0
            self.chromosome_start = self.chromosome_start.astype(int)
            self.chromosome_mid = (self.chromosome_start + self.chromosome_end) / 2.

            self.chromosome_info = pd.DataFrame({
                'chromosome_length': self.chromosome_lengths,
                'chromosome_end': self.chromosome_end,
                'chromosome_start': self.chromosome_start,
                'chromosome_mid': self.chromosome_mid,
            }).reset_index()
```

## get faidx

```python
import pyfaidx
```

```python
genome = pyfaidx.Fasta('/juno/work/shah/users/chois7/tickets/sarcoma89/14474/resources/GRCh38_PBEF1NeoTransposon.fa')
```

# Test alignment

```python
df_path = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/vcf/14427_1_100.vcf.gz'
vcf_cols = ['chrom', 'pos', 'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample']
df = pd.read_table(df_path, comment='#', names=vcf_cols)
```

```python
ins = df[df['ID'].str.count('INS') > 0]
ins = ins[ins['alt'].str.len() > 3000]
```

```python
ins['alt'].str.len().hist(bins=30)
```

```python
ins['alt'].str.len().hist(bins=30)
```

```python
transposon = fa['PBEF1NeoTransposon'][:].seq.upper()
```

```python
alignment = pairwise2.align.globalxx(ins.loc[4810].alt, transposon)
```

```python
len(ins.loc[4810].alt)
```

```python
alignment = pairwise2.align.localms(ins.loc[4810].alt, transposon,
                                    1, -1, -1, -0.1)
```

```python
alignment[1].score / len(seq)
```

```python
for rix, row in ins.iterrows():
    alt = row['alt']
    seq, revcmp = str(Seq(row['alt'])), str(Seq(row['alt']).reverse_complement())
    if str(seq) in transposon:
        print(f'seq found in {rix}')
    if str(revcmp) in transposon:
        print(f'revcmp found in {rix}')
    # alignments = pairwise2.align.globalxx(seq, transposon)
    alignments = pairwise2.align.localms(seq, transposon,
                                         1, -1, -1, -0.1)
    if alignments[0].score / len(seq) > 0.7:
        print(f'seq found in {rix}')
    # alignments = pairwise2.align.globalxx(revcmp, transposon)
    # if alignments[0].score > 3000:
    #     print(f'revcmp found in {rix}')
```

```python
seq
```

```python

```

```python
fa_path = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/PBEF1NeoTransposon.fa'
fa = Fasta(fa_path)
```

```python
# Import pairwise2 module
from Bio import pairwise2

# Import format_alignment method
from Bio.pairwise2 import format_alignment

# Define two sequences to be aligned
X = "ACGGGT"
Y = "ACG"

# Get a list of the global alignments between the two sequences ACGGGT and ACG
# No parameters. Identical characters have score of 1, else 0.
# No gap penalties.
alignments = pairwise2.align.globalxx(X, Y)

# Use format_alignment method to format the alignments in the list
for a in alignments:
    print(format_alignment(*a))
```

```python
a.score
```

# Make report


## Load data


## vector genomes

```python
genomes = {
    '14472_1_100': pyfaidx.Fasta('/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/resources/GRCh38_PBEF1_GFPVector.fa'),
    '14472_2_200': pyfaidx.Fasta('/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/resources/GRCh38_PBEF1_PGBD5Vector.fa'),
    '14472_3_300': pyfaidx.Fasta('/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/resources/GRCh38_PBEF1_PiggyBacVector.fa'),
}
```

### load data per sample into dict *

```python
fqs = {}
for sample in genomes.keys():
    fq_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/fastq/{sample}.rnames.fastq.gz'
    fqs[sample] = read_fastq(fq_path)
```

```python
data_paths = {
    '14472_1_100': '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/vcf/14427_1_100.vcf.gz',
    '14472_2_200': '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/vcf/14427_2_200.vcf.gz',
    '14472_3_300': '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/vcf/14427_3_300.vcf.gz',
}
data = {k: read_vcf(data_paths[k]) for k in data_paths}
```

```python
dst_cols = ['chrom1', 'pos1', 'chrom2', 'pos2', 'type', 'length',
            'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO',
            'sample', 'rnames']
for key, df in data.items():
    df['chrom2'] = df['chrom1']
    df['type'] = 'NA'
    tra_brk_cnt = df['alt'].str.count('\]') + df['alt'].str.count('\[')
    tra = tra_brk_cnt > 0
    df['type'] = df['INFO'].str.extract('SVTYPE=([A-Z]+)')
    del_inv_dup = df['type'].isin({'DEL', 'INV', 'DUP'})
    ins = df['type'] == 'INS'
    df['length'] = abs(df['INFO'].str.extract('SVLEN=([\-\d]+)').astype(float))
    df['rnames'] = df['INFO'].str.extract('RNAMES=([a-z0-9\-,]+)')
    df.loc[tra, 'type'] = 'TRA'
    df.loc[tra, 'length'] = 0
    df.loc[tra, 'chrom2'] = df.loc[tra, 'alt'].str.extract("[\[\]]([_A-Za-z0-9]+):\d+[\[\]]").squeeze()
    df.loc[tra, 'pos2'] = df.loc[tra, 'alt'].str.extract("[\[\]][_A-Za-z0-9]+:(\d+)[\[\]]").squeeze()
    df.loc[~tra, 'chrom2'] = df.loc[~tra, 'chrom1']
    df.loc[del_inv_dup, 'pos2'] = df.loc[del_inv_dup, 'pos1'] + df.loc[del_inv_dup, 'length']
    df.loc[ins, 'pos2'] = df.loc[ins, 'pos1']
    df[['pos1', 'pos2']] = df[['pos1', 'pos2']].astype(int)
    data[key] = df[dst_cols]
```

```python
data['14472_1_100'].head(2)
```

### load fastq per sample into dict *


#### save in .fa format for blast

```python
out_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/fasta'
sv_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/svs'
for sample, (rname2seq, rname2qual) in fqs.items():
    svs = data[sample]
    svs = svs[(svs['chrom1']=='PBEF1NeoTransposon') | (svs['chrom2']=='PBEF1NeoTransposon')]
    svdf = svs[['chrom1', 'pos1', 'chrom2', 'pos2', 'type']]
    n_supports = []
    with open(f'{out_dir}/{sample}.rnames.fasta', 'w') as out:
        for riw, row in svs.iterrows(): # row: one SV
            chrom1, pos1, chrom2, pos2, svtype, svlength, _svid, _ref, _alt, _qual, _filter, _info, _sample, rnames = row.squeeze()
            sv_identity = [chrom1, pos1, chrom2, pos2, svtype]
            sv_identity_str = '|'.join([str(_) for _ in sv_identity])
            # print(sample, chrom1, pos1, chrom2, pos2, svtype)
            rnames = rnames.split(',')
            n_support = len(rnames)
            n_supports.append(n_support)
            
            for rname in rnames: # supporting & chimreric reads into fasta
                seq = rname2seq[rname]
                header = f'>{rname}|{sv_identity_str}'
                out.write(f'{header}\n')
                out.write(f'{str_to_fasta(seq)}\n')
                # break
    svdf['number_of_chimeric_reads'] = n_supports
    svdf.to_csv(f'{sv_dir}/{sample}.sv.csv', index=False)
```

### unfiltered vcf

```python
df = data['14472_1_100']
df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
```

```python
df = data['14472_2_200']
df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
```

```python
df = data['14472_3_300']
df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]#.loc[[19863, 2382, 5079, 9560]]
```

### filtered data

```python
df = data['14472_1_100']
df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
```

```python
df = data['14472_2_200']
df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')].loc[[1, 0]]
```

```python
df = data['14472_3_300']
df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')].loc[[19863, 2382, 5079, 9560]]
```

## PBEF1NeoTransposon analysis

```python
translocation_size = 6894
margin = 20
reg_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/region'
for key, df in data.items():
    reg_path = f'{reg_dir}/{key}.region.txt'
    with open(reg_path, 'w') as out:
        subdf = df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
        for rix, row in subdf.iterrows():
            chrom1, chrom2 = row['chrom1'], row['chrom2']
            pos1, pos2 = row['pos1'], row['pos2']
            start1, end1 = pos1-margin, pos1+margin
            start2, end2 = pos2-margin, pos2+margin
            svtype = row['type']
            svlength = row['length']
            print(f'{key} - {chrom1}:{pos1} - {chrom2}:{pos2} - {svtype} ({svlength:.0f} bp)')
            if chrom1 == 'PBEF1NeoTransposon':
                start1 = max(1, start1)
                end1 = min(translocation_size, end1)
            if chrom2 == 'PBEF1NeoTransposon':
                start2 = max(1, start2)
                end2 = min(translocation_size, end2)
            line = f"{chrom1}:{start1}-{end1} {chrom2}:{start2}-{end2} {row['type']}"
            out.write(line+'\n')
```

### Extract seq near brks

```python
import pyfaidx
```

```python
fa_path = '/juno/work/shah/users/grewald/tickets/SHAH-4349/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa'
fa = pyfaidx.Fasta(fa_path)
```

```python
fa['PBEF1NeoTransposon'][:100].seq.upper()
```

```python
translocation_size = 6894
margin = 200
seq_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/seq'
for key, df in data.items():
    seq_path = f'{seq_dir}/{key}.SV.fasta'
    with open(seq_path, 'w') as out:
        subdf = df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
        for i, (rix, row) in enumerate(subdf.iterrows()):
            chrom1, chrom2 = row['chrom1'], row['chrom2']
            pos1, pos2 = row['pos1'], row['pos2']
            start1, end1 = pos1-margin, pos1+margin # 0-based [,)
            start2, end2 = pos2-margin, pos2+margin # 0-based [,)
            if chrom1 == 'PBEF1NeoTransposon':
                start1 = max(0, start1)
                end1 = min(translocation_size, end1)
            if chrom2 == 'PBEF1NeoTransposon':
                start2 = max(0, start2)
                end2 = min(translocation_size, end2) 
            svtype = row['type']
            svlength = row['length']
            seq1 = fa[chrom1][start1:end1].seq.upper()
            seq2 = fa[chrom2][start2:end2].seq.upper()
            tag1 = f'{key}|SV{i+1}_brk1|{chrom1}|{pos1}|{svtype}|{svlength:.0f}bp'
            block1 = f'>{tag1}\n{seq1}'
            tag2 = f'{key}|SV{i+1}_brk2|{chrom2}|{pos2}|{svtype}|{svlength:.0f}bp'
            block2 = f'>{tag2}\n{seq2}'
            
            # print(f'{key}{chrom1}:{pos1} - {chrom2}:{pos2} - {svtype} ({svlength:.0f} bp) / seq: {len(seq1)}bp, {len(seq2)}bp')
            # line = f"{chrom1}:{start1}-{end1} {chrom2}:{start2}-{end2} {row['type']}"
            out.write(block1+'\n')
            out.write(block2+'\n')
```

### Save rnames for selected breakpoints

```python
for sample, df in data.items():
    tdf = df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
    
    rnames_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/rnames/{sample}.rnames.txt'
    rname_col = tdf['rnames']
    with open(rnames_path, 'w') as out:
        for rname in rname_col:
            line = rname.replace(',', '\n')
            out.write(line + '\n')
```

### Extract clipped sequences from supporting reads *

```python
[read.tags
```

#### save clipped reads

```python
clipped_cols = ['read_name', 'read_alignment_contig_name', 'read_alignment_start_pos', "read_5'-clipped_seq", "read_3'-clipped_seq", "qual_5'-clipped_seq", "qual_3'-clipped_seq"]
samples = list(data.keys())
wkdir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results'
clipped_tags = [4, 5]
sample2seqdf = {}
for sample in samples:
    seqdf = pd.DataFrame(columns=clipped_cols)
    sample_seqs, sample_quals = fqs[sample]
    six = 0
    df = data[sample]
    df = df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
    bam_path = f'{wkdir}/bam/{sample}.bam'
    bam = pysam.AlignmentFile(bam_path, 'rb')
    
    for rix, row in df.iterrows():
        chrom1, pos1, chrom2, pos2 = row['chrom1'], row['pos1'], row['chrom2'], row['pos2']
        rnames = row['rnames'].split(',') # SV-supporting reads for this adjacency
        regions = [(chrom1, pos1), (chrom2, pos2)]
        for chrom, pos in regions:
            if chrom != 'PBEF1NeoTransposon': continue # only take reads from transposon
            reads = bam.fetch(chrom, pos-1, pos) # all bam reads from this breakpoint
            for read in reads:
                if read.qname in rnames: # is this read a supporter?
                    ncseq = sample_seqs[read.qname] # non-clipped seq
                    ncqual = sample_quals[read.qname] # non-clipped seq
                    # need to revcomp fastq seq if read is reverse mapped
                    if read.is_reverse:
                        ncseq = str(Seq(ncseq).reverse_complement())
                        ncqual = ncqual[::-1]
                    # check 5'/3' of read
                    seq_5p, seq_3p, qual_5p, qual_3p = '', '', '', ''
                    cigar_5p = read.cigar[0]
                    cigar_3p = read.cigar[-1]
                    if cigar_5p[0] in clipped_tags:
                        # seq_5p = read.seq[:cigar_5p[1]]
                        seq_5p = ncseq[:cigar_5p[1]]
                        qual_5p = ncqual[:cigar_5p[1]]
                    if cigar_3p[0] in clipped_tags:
                        # seq_3p = read.seq[-cigar_3p[1]:]
                        seq_3p = ncseq[-cigar_3p[1]:] 
                        qual_3p = ncqual[-cigar_3p[1]:] 
                        
                    # if len(seq_5p) and len(seq_3p):
                    field = [read.qname, read.reference_name, read.pos+1, seq_5p, seq_3p, qual_5p, qual_3p]
                    seqdf.loc[six] = field
                    six += 1
    # clipped_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/clipped/{sample}.clipped.csv'
    # seqdf.to_csv(clipped_path, index=False)
    sample2seqdf[sample] = seqdf
```

#### save sv information + n_chimeric_reads

```python
secondary
```

```python
clipped_cols = ['read_name', 'read_alignment_contig_name', 'read_alignment_start_pos', "read_5'-clipped_seq", "read_3'-clipped_seq", "qual_5'-clipped_seq", "qual_3'-clipped_seq"]
samples = list(data.keys())
wkdir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results'
clipped_tags = [4, 5]
sample2seqdf = {}
for sample in samples:
    secondary_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/secondary/{sample}.2ndary.csv'
    with open(secondary_path, 'w') as out:
        header = ','.join(['sample', 'read_name', 'chrom1_breakpoint', 'pos1_breakpoint', 'chrom2_breakpoint', 'pos2_breakpoint', 
                           'type', 'chrom_secondary', 'pos_secondary', 'strand_secondary'])
        out.write(header + '\n')
        seqdf = pd.DataFrame(columns=clipped_cols)
        sample_seqs, sample_quals = fqs[sample]
        six = 0
        df = data[sample]
        df = df[(df['chrom1']=='PBEF1NeoTransposon') | (df['chrom2']=='PBEF1NeoTransposon')]
        svs = df[['chrom1', 'pos1', 'chrom2', 'pos2']]
        # n_supportings = []
        n_chimerics = []
        bam_path = f'{wkdir}/bam/{sample}.bam'
        bam = pysam.AlignmentFile(bam_path, 'rb')

        for rix, row in df.iterrows():
            chrom1, pos1, chrom2, pos2 = row['chrom1'], row['pos1'], row['chrom2'], row['pos2']
            rnames = row['rnames'].split(',') # SV-supporting reads for this adjacency
            # n_supporting = len(rnames)
            n_chimeric = 0
            regions = [(chrom1, pos1), (chrom2, pos2)]
            for chrom, pos in regions:
                if chrom != 'PBEF1NeoTransposon': continue # only take reads from transposon
                reads = bam.fetch(chrom, pos-1, pos) # all bam reads from this breakpoint
                for read in reads:
                    if read.qname in rnames: # is this read a supporter?
                        secondaries = get_secondaries(read)
                        if len(secondaries) > 0:
                            n_chimeric += 1
                            for secondary in secondaries:
                                chrom_2nd, pos_2nd, strand_2nd, cigar_2nd = secondary.split(',')[:4]
                                field = [sample, read.qname, chrom1, pos1, chrom2, pos2, row['type'], chrom_2nd, pos_2nd, strand_2nd]
                                field = [str(s) for s in field]
                                line = ','.join(field)
                                out.write(line + '\n')

                        ncseq = sample_seqs[read.qname] # non-clipped seq
                        ncqual = sample_quals[read.qname] # non-clipped seq
                        # need to revcomp fastq seq if read is reverse mapped
                        if read.is_reverse:
                            ncseq = str(Seq(ncseq).reverse_complement())
                            ncqual = ncqual[::-1]
                        # check 5'/3' of read
                        seq_5p, seq_3p, qual_5p, qual_3p = '', '', '', ''
                        cigar_5p = read.cigar[0]
                        cigar_3p = read.cigar[-1]
                        if cigar_5p[0] in clipped_tags:
                            # seq_5p = read.seq[:cigar_5p[1]]
                            seq_5p = ncseq[:cigar_5p[1]]
                            qual_5p = ncqual[:cigar_5p[1]]
                        if cigar_3p[0] in clipped_tags:
                            # seq_3p = read.seq[-cigar_3p[1]:]
                            seq_3p = ncseq[-cigar_3p[1]:] 
                            qual_3p = ncqual[-cigar_3p[1]:] 

                        # if len(seq_5p) and len(seq_3p):
                        field = [read.qname, read.reference_name, read.pos+1, seq_5p, seq_3p, qual_5p, qual_3p]
                        seqdf.loc[six] = field
                        six += 1
                    # fi
                # for chrom, pos
            n_chimerics.append(n_chimeric)
        # svs['n_supportings'] = n_supportings
        svs['n_chimeric'] = n_chimerics # chimeric reads, that support the sv
        # svs_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/svs/{sample}.svs.csv'
        # svs.to_csv(svs_path, index=False)
        # sample2seqdf[sample] = seqdf
```

```python
seqdf.head(2)
```

```python
out_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/clipped/fastq'
for sample, seqdf in sample2seqdf.items():
    with open(f'{out_dir}/{sample}.fastq', 'w') as out:
        for rix, row in seqdf.iterrows():
            rname, contig, start, seq5p, seq3p, qual5p, qual3p = row.squeeze()
            assert len(seq5p) == len(qual5p)
            assert len(seq3p) == len(qual3p)
            if len(seq5p) > 0:
                header = f'{rname}:{contig}:{start}:clipped5'
                out.write(f'{"@"+header}\n')
                out.write(f'{seq5p}\n')
                out.write(f'+\n')
                out.write(f'{qual5p}\n')
            if len(seq3p) > 0:
                header = f'{rname}:{contig}:{start}:clipped3'
                out.write(f'{"@"+header}\n')
                out.write(f'{seq3p}\n')
                out.write(f'+\n')
                out.write(f'{qual3p}\n')
```

```python
#                     secondaries = get_secondaries(read) # get 2ndary alignment infos
#                     cstats = read.get_cigar_stats()[0] # base counts
#                     ncseq = fqs[sample][read.qname] # non-clipped seq
#                     # [:1904] == read.seq[:1904]
                    
#                     tname, tpos, tstrand, secondary_aligned_seq = '', '', '', ''
#                     if len(secondaries):
#                         for secondary in secondaries:
#                             (tname, tpos, tstrand, tcigar, secondary_aligned_seq) = get_secondary_aligned_seq(secondary, ncseq)
#                             if len(secondary_aligned_seq) > 0:
#                                 print(tname, tpos, tstrand, tcigar, secondary_aligned_seq)
#                                 raise ValueError
#                     # field = [read.qname, 
```

## INS analysis

```python
chroms = ['chr'+str(c) for c in range(1, 22+1)] + ['chrX', 'chrY', 'chrM']
translocation_size = 6894
margin = 100
reg_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/region'
for key, df in data.items():
    ins = df[df['type'] == 'INS']
    ins = ins[ins['length'] > 6000]
    ins = ins[ins['length'] < 7000]
    ins = ins[ins['chrom1'].isin(chroms) | ins['chrom2'].isin(chroms)]
```

```python
ins[abs(ins['length'] - translocation_size) < 100]

seq = ins.iloc[1]['alt']
alignments = pairwise2.align.localms(seq, transposon,
                                         1, -1, -1, -0.1)
```

# Vector alignment analysis


## Load data


### load data per sample into dict *

```python
chromosomes = ['chr'+str(i) for i in range(1, 22+1)] + ['chrX', 'chrY']
vectors = ['PBEF1NeoTransposon', 'GFPVector', 'puro-GFP-PGBD5_seq', 'PiggyBacVector']
chromosomes += vectors
len(chromosomes)
```

```python
data_paths = {
    '14472_1_100': '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/GFPVector/14472_1_100/results/sniffles/14472_1_100.vcf',
    '14472_2_200': '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/PGBD5Vector/14472_2_200/results/sniffles/14472_2_200.vcf',
    '14472_3_300': '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/PiggyBacVector/14472_3_300/results/sniffles/14472_3_300.vcf',
}
vcf_data = {k: read_vcf(data_paths[k]) for k in data_paths}
```

```python
fqs = {}
for sample in data_paths.keys():
    fq_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results/fastq/{sample}.rnames.fastq.gz'
    fqs[sample] = read_fastq(fq_path)
```

```python
assert len(chromosomes) == 28, chromosomes
data = {}
dst_cols = ['chrom1', 'pos1', 'chrom2', 'pos2', 'type', 'length',
            'ID', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO',
            'sample', 'rnames']
for key, df in vcf_data.items():
    df['chrom2'] = df['chrom1']
    df['type'] = 'NA'
    tra_brk_cnt = df['alt'].str.count('\]') + df['alt'].str.count('\[')
    tra = tra_brk_cnt > 0
    df['type'] = df['INFO'].str.extract('SVTYPE=([A-Z]+)')
    del_inv_dup = df['type'].isin({'DEL', 'INV', 'DUP'})
    ins = df['type'] == 'INS'
    df['length'] = abs(df['INFO'].str.extract('SVLEN=([\-\d]+)').astype(float))
    df['rnames'] = df['INFO'].str.extract('RNAMES=([a-z0-9\-,]+)')
    df.loc[tra, 'type'] = 'TRA'
    df.loc[tra, 'length'] = 0
    df.loc[tra, 'chrom2'] = df.loc[tra, 'alt'].str.extract("[\[\]]([_A-Za-z0-9]+):\d+[\[\]]").squeeze()
    df.loc[tra, 'pos2'] = df.loc[tra, 'alt'].str.extract("[\[\]][_A-Za-z0-9]+:(\d+)[\[\]]").squeeze()
    df.loc[~tra, 'chrom2'] = df.loc[~tra, 'chrom1']
    df.loc[del_inv_dup, 'pos2'] = df.loc[del_inv_dup, 'pos1'] + df.loc[del_inv_dup, 'length']
    df.loc[ins, 'pos2'] = df.loc[ins, 'pos1']
    df[['pos1', 'pos2']] = df[['pos1', 'pos2']].astype(int)
    
    df = df[(df['chrom1'].isin(chromosomes)) & (df['chrom2'].isin(chromosomes))]
    df = df[(df['chrom1'].isin(vectors)) | (df['chrom2'].isin(vectors))]
    
    data[key] = df[dst_cols]
```

```python
data['14472_3_300']
```

### Save rnames for selected breakpoints

```python
chroms_proc = ['PBEF1NeoTransposon', 'GFPVector', 'puro-GFP-PGBD5_seq', 'PiggyBacVector']
for ix, sample in enumerate(samples):
    ref0 = chroms_proc[0]
    ref1 = chroms_proc[ix+1]
    bam_path = f'{wkdir}/../bams/{sample}.sorted.bam'
    bam = pysam.AlignmentFile(bam_path, 'rb')
    rnames = set()
    for ref in [ref0, ref1]:
        for read in bam.fetch(ref):
            rnames.add(read.qname)
    print(sample, 'done')
            
    with open(f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results/rnames/{sample}.rnames.txt', 'w') as out:
        content = '\n'.join(list(rnames)) + '\n'
        out.write(content)
```

## save clipped reads

```python
refgenome = RefGenomeInfo('hg38')
```

```python
clipped_cols = ['read_name', 'read_alignment_contig_name', 'read_alignment_start_pos', 'read_alignment_strand', "read_5'-clipped_seq", "read_3'-clipped_seq", "qual_5'-clipped_seq", "qual_3'-clipped_seq", "read_sequence"]
chroms_proc = ['PBEF1NeoTransposon', 'GFPVector', 'puro-GFP-PGBD5_seq', 'PiggyBacVector']
samples = list(data.keys())
wkdir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results'
clipped_tags = [4, 5]
sample2seqdf = {}
for sample in samples:
    secondary_path = f'{wkdir}/secondary/{sample}.2ndary.csv'
    with open(secondary_path, 'w') as out:
        # header = ','.join(['sample', 'read_name', 'chrom1_breakpoint', 'pos1_breakpoint', 'chrom2_breakpoint', 'pos2_breakpoint', 
        #                    'type', 'chrom_secondary', 'pos_secondary', 'strand_secondary'])
        # out.write(header + '\n')
        seqdf = pd.DataFrame(columns=clipped_cols)
        sample_seqs, sample_quals = fqs[sample]
        six = 0
        df = data[sample]
        svs = df[['chrom1', 'pos1', 'chrom2', 'pos2']]
        n_chimerics = []
        bam_path = f'{wkdir}/../bams/{sample}.sorted.bam'
        bam = pysam.AlignmentFile(bam_path, 'rb')

        for rix, row in df.iterrows():
            chrom1, pos1, chrom2, pos2 = row['chrom1'], row['pos1'], row['chrom2'], row['pos2']
            rnames = row['rnames'].split(',') # SV-supporting reads for this adjacency
            # n_supporting = len(rnames)
            n_chimeric = 0
            regions = [(chrom1, pos1), (chrom2, pos2)]
            for chrom, pos in regions:
                if chrom not in chroms_proc: continue
                reads = bam.fetch(chrom, pos-1, pos) # all bam reads from this breakpoint
                for read in reads:
                    if read.qname in rnames: # is this read a supporter?
                        secondaries = get_secondaries(read)
                        if len(secondaries) > 0:
                            n_chimeric += 1

                        ncseq = sample_seqs[read.qname] # non-clipped seq
                        ncqual = sample_quals[read.qname] # non-clipped seq
                        # need to revcomp fastq seq if read is reverse mapped
                        if read.is_reverse: # only change ncseq, ncqual from fastq
                            ncseq = str(Seq(ncseq).reverse_complement())
                            ncqual = ncqual[::-1]
                        # check 5'/3' of read
                        seq_5p, seq_3p, qual_5p, qual_3p = '', '', '', ''
                        cigar_5p_type, cigar_5p_len = read.cigar[0]
                        cigar_3p_type, cigar_3p_len = read.cigar[-1]
                        
                        if cigar_5p_type in clipped_tags:
                            # seq_5p = read.seq[:cigar_5p[1]]
                            seq_5p = ncseq[:cigar_5p_len]
                            qual_5p = ncqual[:cigar_5p_len]
                        if cigar_3p_type in clipped_tags:
                            # seq_3p = read.seq[-cigar_3p[1]:]
                            seq_3p = ncseq[-cigar_3p_len:] 
                            qual_3p = ncqual[-cigar_3p_len:] 

                        # if len(seq_5p) and len(seq_3p):
                        strand = '-' if read.is_reverse else '+'
                        field = [read.qname, read.reference_name, read.pos+1, strand, seq_5p[-50:], seq_3p[:50], qual_5p[-50:], qual_3p[:50], read.seq]
                        seqdf.loc[six] = field
                        six += 1
                        
                        # if read.qname == '8d49baf4-f7be-46e5-b495-3e6df968f6c1':
                        # if read.is_reverse:
                        # if cigar_5p_type == 4: 
                        if read.qname == 'b9315c0c-67bd-47f9-996b-5b12124a5824':
                            _sample_seqs = sample_seqs
                            _read = read
                            _field = field
                            raise ValueError
                    # fi
                # for chrom, pos
            n_chimerics.append(n_chimeric)
        # svs['n_supportings'] = n_supportings
        # svs['n_chimeric'] = n_chimerics # chimeric reads, that support the sv
        # svs_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/svs/{sample}.svs.csv'
        # svs.to_csv(svs_path, index=False)
        sample2seqdf[sample] = seqdf
```

```python
[t for t in read.tags if t[0]=='SA'][0][1].split(';')
```

```python
for sample, seqdf in sample2seqdf.items():
    out_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results/clipped/{sample}.clipped.csv'
    print(seqdf.shape, seqdf.drop_duplicates().shape)
    seqdf = seqdf.drop_duplicates()
    seqdf.to_csv(out_path, sep='\t', index=False)
```

#### 26655a81-a557-4ea9-9ef1-c48748672de5 (read.is_reverse == True)

```python
sample, read.qname, read.reference_name, read.pos, bam_path, read.is_reverse, read.cigarstring
```

```python
read.seq[958:968]
```

#### 26655a81-a557-4ea9-9ef1-c48748672de5 (read.is_reverse == True)

```python
sample, read.qname, read.reference_name, read.pos, bam_path, read.is_reverse, read.cigarstring
```

```python
ncseq[cigar_5p_len:cigar_5p_len+10]
```

```python
read.seq[:10]
```

```python
ncseq[-cigar_3p_len-10:-cigar_3p_len]
```

```python
read.seq[-10:]
```

#### 8d49baf4-f7be-46e5-b495-3e6df968f6c1 validation

```python
sample, read.qname, bam_path
```

```python
read.seq[:10]
```

```python
read.reference_name
```

## save breakpoint sequences **

```python
chrom_seq = genomes[sample][chrom]
```

```python
def add_info_to_df(df, rname, svtype, ref, pos, strand, end_type, seq_5p, seq_3p):
    ix = df.shape[0]
    field = [rname, svtype, ref, pos, strand, end_type, seq_5p, seq_3p]
    df.loc[ix] = field
```

```python
def get_cigar_tuples(cigarstring):
    cigar_types = 'MIDNSHP=X'
    converter = {
        cigar_types[i]:i for i in range(len(cigar_types))
    }
    cigar_tuples = []
    for match in re.finditer('(\d+)([A-Z\=])', cigarstring):
        clen, ctype = match.groups()
        cigar_tuple = (converter[ctype], int(clen))
        cigar_tuples.append(cigar_tuple)
    return cigar_tuples
```

```python
def append_breakpoint_seq(seqdf, cigar_tuples, qname, refname, strand, tpos, read_seq, ncseq, chrom_fasta):
    read_seq = None
    qpos = 0 # init qpos
    for cix, (cigar_type, cigar_len) in enumerate(cigar_tuples):
        seq_5p, seq_3p = '', ''
        pos = tpos+1
        if cigar_type == 5:
            svtype = 'HC'
            hc_brk, hc_end = qpos, qpos + cigar_len
            if cix == 0:
                seq_5p = ncseq[hc_brk:hc_end][-nbp:]
                seq_3p = ncseq[hc_end:hc_end+nbp]
                end_type = "5'"
            elif cix == len(cigar_tuples)-1:
                seq_5p = ncseq[hc_brk-nbp:hc_brk]
                seq_3p = ncseq[hc_brk:hc_end][:nbp]
                end_type = "3'"
            else:
                raise ValueError(qname, cigar_tuples, cix, cigar_type, cigar_len)
            add_info_to_df(seqdf, qname, svtype, refname, pos, strand, end_type, seq_5p, seq_3p)

        elif cigar_type == 4:
            svtype = 'SC'
            pos = tpos+1
            sc_brk, sc_end = qpos, qpos + cigar_len
            if cix == 0:
                seq_5p = ncseq[sc_brk:sc_end][-nbp:]
                seq_3p = ncseq[sc_end:sc_end+nbp]
                end_type = "5'"
            elif cix == len(cigar_tuples)-1:
                seq_5p = ncseq[sc_brk-nbp:sc_brk]
                seq_3p = ncseq[sc_brk:sc_end][:nbp]
                end_type = "3'"
            else:
                raise ValueError(qname, cigar_tuples, cix, cigar_type, cigar_len)
            add_info_to_df(seqdf, qname, svtype, refname, pos, strand, end_type, seq_5p, seq_3p)
            # if qname == 'ca4af55b-b014-4135-9996-571e5d7d2132':
            #     print(svtype, refname, pos, strand)
            #     print(f'seq_5p {seq_5p}')
            #     print(f'seq_3p {seq_3p}')
            #     print(f'qpos {qpos}, cigar_len:{cigar_len}')
            #     print(f'read_seq:\n{read_seq}')
            
        elif cigar_type == 1: # INS
            svtype = 'INS'
            # ins_seq = read_seq[qpos : qpos+cigar_len]
            ins_seq = ncseq[qpos : qpos+cigar_len]
            if cigar_len > 100:
                seq_5p_1 = ncseq[qpos-nbp:qpos] #chrom_fasta[tpos-nbp:tpos].seq.upper()
                seq_3p_1 = ins_seq[:nbp]
                seq_5p_2 = ins_seq[-nbp:]
                seq_3p_2 = ncseq[qpos+cigar_len:qpos+cigar_len+nbp] #chrom_fasta[tpos:tpos+nbp].seq.upper()

                add_info_to_df(seqdf, qname, svtype, refname, qpos+1, strand, "5'", seq_5p_1, seq_3p_1)
                add_info_to_df(seqdf, qname, svtype, refname, qpos+1, strand, "3'", seq_5p_2, seq_3p_2)
                
        elif cigar_type == 2: # DEL
            svtype = 'DEL'
            # del_seq = read_seq[qpos : qpos+cigar_len]
            del_seq = chrom_fasta[tpos : tpos+cigar_len].seq.upper()
            if cigar_len > 100:
                seq_5p_1 = ncseq[qpos-nbp:qpos] #chrom_fasta[tpos-nbp:tpos].seq.upper()
                seq_3p_1 = del_seq[:nbp]
                seq_5p_2 = del_seq[-nbp:]
                seq_3p_2 = ncseq[qpos:qpos+nbp] #chrom_fasta[tpos:tpos+nbp].seq.upper()

                add_info_to_df(seqdf, qname, svtype, refname, pos, strand, "5'", seq_5p_1, seq_3p_1)
                add_info_to_df(seqdf, qname, svtype, refname, pos+cigar_len, strand, "3'", seq_5p_2, seq_3p_2)
                
                if qname == '1d42bd2c-63a1-41c3-9e88-0bb900d5a517':
                    return qname, svtype, refname, qpos, strand, ncseq, del_seq
                
        if qname == '6f18bef6-1f57-4f9e-9914-8d633176d182':
            print('pos:', pos, 'cigar:', cigar_type, cigar_len, query_consumers)

        if cigar_type in reference_consumers:
            tpos += cigar_len
        if cigar_type in query_consumers:
            if qname == '6f18bef6-1f57-4f9e-9914-8d633176d182':
                print(f'adding in: cigar_type:{cigar_type}[{type(cigar_type)}], cigar_len:{cigar_len}, {query_consumers}')
            qpos += cigar_len
```

```python
# # clipped_cols = ['read_name', 'read_alignment_contig_name', 'read_alignment_start_pos', 'read_alignment_strand', "read_5'-clipped_seq", "read_3'-clipped_seq", "qual_5'-clipped_seq", "qual_3'-clipped_seq", "read_sequence"]
# brk_cols = ['read_name', 'breakpoint_type', 'contig_name', 'pos', 'strand', 'end_type', '5-prime', '3-prime']
# chroms_proc = ['PBEF1NeoTransposon', 'GFPVector', 'puro-GFP-PGBD5_seq', 'PiggyBacVector']
# query_consumers = {0, 1, 4, 7, 8} | {5} # SAM format + {5}
# reference_consumers = {0, 2, 3, 7, 8}
# nbp = 50 # 50 bp of near breakpoint seq
# samples = list(data.keys())
# wkdir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results'
# clipped_tags = [4, 5]
# sample2seqdf = {}

# for sample in samples:
#     seqdf = pd.DataFrame(columns=brk_cols)
#     sample_seqs, sample_quals = fqs[sample]
#     six = 0
#     df = data[sample]
#     svs = df[['chrom1', 'pos1', 'chrom2', 'pos2']]
#     n_chimerics = []
#     bam_path = f'{wkdir}/../bams/{sample}.sorted.bam'
#     bam = pysam.AlignmentFile(bam_path, 'rb')

#     for rix, row in df.iterrows():
#         chrom1, pos1, chrom2, pos2 = row['chrom1'], row['pos1'], row['chrom2'], row['pos2']
#         rnames = row['rnames'].split(',') # SV-supporting reads for this adjacency
#         # n_supporting = len(rnames)
#         n_chimeric = 0
#         regions = [(chrom1, pos1), (chrom2, pos2)]
#         for chrom, pos in regions:
#             chrom_fasta = genomes[sample][chrom]
#             if chrom not in chroms_proc: continue
#             reads = bam.fetch(chrom, pos-1, pos) # all bam reads from this breakpoint
#             for read in reads:
#                 secondaries = get_secondaries(read)
#                 if len(secondaries) < 1: # only consider reads with secondaries (chimeric reads)
#                     continue
#                 ncseq_pos = sample_seqs[read.qname].upper() # non-clipped seq
#                 ncseq_neg = str(Seq(ncseq_pos).reverse_complement())
#                 ncseq = ncseq_pos
#                 if read.is_reverse: # only change ncseq, ncqual from fastq
#                     ncseq = ncseq_neg

#                 strand = '-' if read.is_reverse else '+'

#                 tpos = read.pos
#                 cigarstring = read.cigarstring
#                 cigar_tuples = get_cigar_tuples(cigarstring)
#                 qname = read.qname
#                 refname = read.reference_name
#                 read_seq = None #read.seq
#                 obj = append_breakpoint_seq(seqdf, cigar_tuples, qname, refname, strand, tpos, read_seq, ncseq, chrom_fasta)
                
#                 # if qname == '1d42bd2c-63a1-41c3-9e88-0bb900d5a517':
#                 #     raise ValueError
                    
#                 for secondary in secondaries:
#                     refname_s, tpos_s, strand_s, cigarstring_s, _, _ = secondary.split(',')
#                     tpos_s = int(tpos_s) - 1
#                     cigar_tuples_s = get_cigar_tuples(cigarstring_s)
#                     # ['chr12,122054404,-,1086S517M41S,60,18']
#                     ncseq = ncseq_pos
#                     if strand_s == '-': # only change ncseq, ncqual from fastq
#                         ncseq = ncseq_neg
#                     obj = append_breakpoint_seq(seqdf, cigar_tuples_s, qname, refname_s, strand_s, tpos_s, read_seq, ncseq, genomes[sample][refname_s])
                    


#     brk_path = f'{wkdir}/breakpoints/{sample}.brk.csv'
#     seqdf.drop_duplicates().to_csv(brk_path, index=False)
```

## motif analysis


### modules

```python
def calc_bits(seq):
    counts = Counter(seq)
    keys = counts.keys()
    total = sum([counts[k] for k in keys])
    probs = {k: counts[k]/total for k in keys} 
    entropies = {k: -probs[k] * np.log2(probs[k]) for k in keys}
    total_entropy = sum(entropies.values())
    # (4⋅−0.25⋅log2(0.25))−(−0.7⋅log2(0.7)+−0.3⋅log2(0.3))=1.12
    expected_entropy = -4 * 1/4 * np.log2(1/4)
    # print(expected_entropy)
    total_bit = expected_entropy - total_entropy
    # print(total_bit)
    bits = {k: probs[k] * total_bit for k in keys}
    return entropies, bits
```

```python
def get_sequence_bit_data(cdf):
    bases = ['A', 'C', 'G', 'T']
    data = pd.DataFrame(columns=bases)
    
    up_min_len = cdf['5-prime'].str.len().min()
    down_min_len = cdf['3-prime'].str.len().min()
    up_seqs = [x[-up_min_len:][::-1] for x in cdf['5-prime'].tolist()]
    down_seqs = [x[:down_min_len] for x in cdf['3-prime'].tolist()]
    
    for i, seq in enumerate(zip(*up_seqs)):
        entropies, bits = calc_bits(seq)
        ix = -i - 1
        bit_list = [bits[k] if k in bits else 0 for k in bases]
        data.loc[ix] = bit_list
    for i, seq in enumerate(zip(*down_seqs)):
        entropies, bits = calc_bits(seq)
        ix = i + 1
        bit_list = [bits[k] if k in bits else 0 for k in bases]
        data.loc[ix] = bit_list
        
    return data
```

```python
def plot_logo(bits, title='', png_path=''):
    logo = logomaker.Logo(bits)
    logo.fig.suptitle
    logo.fig.set_figheight(2); logo.fig.set_figwidth(20);
    logo.ax.spines['top'].set_visible(False)
    logo.ax.spines['right'].set_visible(False)
    if title:
        logo.ax.set_title(title)
    if png_path:
        plt.tight_layout()
        plt.savefig(png_path)
    return logo
```

### process

```python
def get_motif_consensus(cdf):
    instances5 = [Seq(s[-5:]) for s in cdf["5-prime"]]
    m5 = motifs.create(instances5)
    m5r = m5.reverse_complement()

    instances3 = [Seq(s[:5]) for s in cdf["3-prime"]]
    m3 = motifs.create(instances3)
    m3r = m3.reverse_complement()
    
    return str(m5.consensus), str(m5r.consensus), str(m3.consensus), str(m3r.consensus)
```

```python
distances = pd.DataFrame(columns=['index_i', 'index_j', 'seq_i', 'seq_j', 'levenshtein'])
n_seqs = mseqs.shape[0]
for i in range(n_seqs):
    for j in range(i+1, n_seqs):
        indexi, seqi = mseqs.iloc[i].squeeze()
        indexj, seqj = mseqs.iloc[j].squeeze()
        distance = Levenshtein.distance(seqi, seqj)
        distances.loc[distances.shape[0]] = [indexi, indexj, seqi, seqj, distance]
```

```python
close_dists = distances[distances['levenshtein'] <= 1]
```

```python
def get_split_df(close_dists, index_id):
    split_cols = ['sample', 'breakpoint_type', 'contig_name', 'pos', 'end_type', 'motif_position']
    i_split = close_dists[f'index_{index_id}'].str.split('__', expand=True)
    i_split.columns = ['sample', 'breakpoint_type', 'chrpos', 'end_type', 'motif_position']
    i_chrpos = i_split['chrpos'].str.split(':', expand=True)
    i_chrpos.columns = ['contig_name', 'pos']
    i_split = i_split.join(i_chrpos)
    i_split = i_split[split_cols]
    i_split.columns = [f'{c}_{index_id}' for c in split_cols]
    return i_split

i_split = get_split_df(close_dists, 'i')
j_split = get_split_df(close_dists, 'j')
# j_split = close_dists['index_j'].str.split('__', expand=True)
```

```python
dist_cols = [
    'sample_i', 'breakpoint_type_i', 'contig_name_i', 'pos_i', 'end_type_i', 'motif_position_i', 
    'sample_j', 'breakpoint_type_j', 'contig_name_j', 'pos_j', 'end_type_j', 'motif_position_j',
    'seq_i', 'seq_j', 'levenshtein', 
]
close_dists = close_dists.join(i_split).join(j_split)
close_dists = close_dists[dist_cols]
```

```python
# close_dists.to_csv('~/chois7/tickets/sarcoma89/14474/vectors/results/motifs/motif_comparison.cutoff1.tsv', sep='\t', index=False)
```

```python
close_dists = pd.read_table('~/chois7/tickets/sarcoma89/14474/vectors/results/motifs/motif_comparison.cutoff1.tsv')
```

```python
motif_counts = pd.Series(close_dists['seq_i'].tolist() + close_dists['seq_j'].tolist()).value_counts()
motif_counts = pd.DataFrame(motif_counts).reset_index()
motif_counts.columns = ['motif', 'count']
```

```python
motif_counts.to_csv('~/chois7/tickets/sarcoma89/14474/vectors/results/motifs/motif_counts.cutoff1.tsv', sep='\t', index=False)
```

```python
samples = ['14472_1_100', '14472_2_200', '14472_3_300']
```

```python
for (si, sj), subdf in close_dists.groupby(['sample_i', 'sample_j']):
    motif_counts = pd.Series(subdf['seq_i'].tolist() + subdf['seq_j'].tolist()).value_counts()
    motif_counts = pd.DataFrame(motif_counts).reset_index()
    motif_counts.columns = ['motif', 'count']
    motif_counts.to_csv(f'~/chois7/tickets/sarcoma89/14474/vectors/results/motifs/counts/motif_counts.cutoff1.{si}__{sj}.tsv', sep='\t', index=False)
```

```python
pivot = distances[['index_i', 'index_j', 'levenshtein']].pivot(index='index_i', columns='index_j')
```

```python
import seaborn as sns
```

```python
fig, ax = plt.subplots(figsize=(20, 20))
sns.heatmap((5-pivot).fillna(0), annot=True, cbar=False)
```

```python
dist0s = distances[distances['levenshtein']==0]
```

```python
counts = dist0s['seq_i'].value_counts()
```

```python
dist0s[dist0s.index_i.str.count('2_200') > 0].to_csv('')
```

```python
dist0s[dist0s['seq_i'].isin((counts > 5).index)]
```

```python
wkdir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results'
samples = ['14472_1_100', '14472_2_200', '14472_3_300']

mseqs = pd.DataFrame(columns=['index', 'seq']) # motif seqs

cluster_count_cutoff = 1
for sample in samples:
    # if sample == '14472_1_100': continue ##@##
    brk_path = f'{wkdir}/breakpoints/{sample}.brk.csv'
    df = pd.read_csv(brk_path)
    
    for (svtype, contig, pos, end_type), cdf in df.groupby(['breakpoint_type', 'contig_name', 'pos', 'end_type']): # coordinate df
        if cdf.shape[0] < cluster_count_cutoff:
            continue
        # bits = get_sequence_bit_data(cdf)
        # end_tag = 'upstream' if end_type.startswith('5') else 'downstream'
        # png_path = f'./vectors/results/motifs/{sample}.{svtype}.{contig}__{pos}.{end_tag}.png'
        # logo = plot_logo(bits, title=f'{sample}: {svtype} at {contig}:{pos} {end_tag}', 
        #                  png_path = None)
                         # png_path=png_path)
        # raise ValueError
        m5_seq, m5r_seq, m3_seq, m3r_seq = get_motif_consensus(cdf)
        seq_tag = f'{sample}__{svtype}__{contig}:{pos}__{end_tag}'
        mseqs.loc[mseqs.shape[0]] = [f'{seq_tag}__5prime', m5_seq]
        mseqs.loc[mseqs.shape[0]] = [f'{seq_tag}__3prime', m3_seq]
```

```python
from Bio import motifs
```

```python
from collections import Counter
```

```python

```

```python
bits
```

```python
seq = 'AAAAA'
counts = Counter(seq)
keys = counts.keys()
total = sum([counts[k] for k in keys])
probs = {k:counts[k]/total for k in keys}
infos = {k: probs[k] * np.log2(probs[k] / 0.25) for k in keys} # ∑(P(i,j) * log2(P(i,j)/B(i)))
total_infos = sum(infos.values())
bits = {k: infos[k] / total_infos for k in keys}
```

```python
bits
```

## ---------

```python
qname, svtype, refname, qpos, strand, ncseq, del_seq = obj
```

```python
len(del_seq)
```

```python
ncseq[qpos : qpos+len(del_seq)]
```

```python
refname, qpos, qpos+len(del_seq)
```

```python
ncseq.index(obj[-1])
```

```python
ncseq[35376:35386]
```

```python
obj[-1]
```

```python
len(s) - 11
```

```python
# clipped_cols = ['read_name', 'read_alignment_contig_name', 'read_alignment_start_pos', "read_5'-clipped_seq", "read_3'-clipped_seq", "qual_5'-clipped_seq", "qual_3'-clipped_seq"]
# chroms_proc = ['PBEF1NeoTransposon', 'GFPVector', 'puro-GFP-PGBD5_seq', 'PiggyBacVector']
# samples = list(data.keys())
# wkdir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/vectors/results'
# clipped_tags = [4, 5]
# sample2seqdf = {}
# for sample in samples:
#     secondary_path = f'{wkdir}/secondary/{sample}.2ndary.csv'
#     with open(secondary_path, 'w') as out:
#         header = ','.join(['sample', 'read_name', 'chrom1_breakpoint', 'pos1_breakpoint', 'chrom2_breakpoint', 'pos2_breakpoint', 
#                            'type', 'chrom_secondary', 'pos_secondary', 'strand_secondary'])
#         out.write(header + '\n')
#         seqdf = pd.DataFrame(columns=clipped_cols)
#         sample_seqs, sample_quals = fqs[sample]
#         six = 0
#         df = data[sample]
#         df = df[(df['chrom1'].isin(chroms_proc)) | (df['chrom2'].isin(chroms_proc))]
#         svs = df[['chrom1', 'pos1', 'chrom2', 'pos2']]
#         # n_supportings = []
#         n_chimerics = []
#         bam_path = f'{wkdir}/../bams/{sample}.sorted.bam'
#         bam = pysam.AlignmentFile(bam_path, 'rb')

#         for rix, row in df.iterrows():
#             chrom1, pos1, chrom2, pos2 = row['chrom1'], row['pos1'], row['chrom2'], row['pos2']
#             rnames = row['rnames'].split(',') # SV-supporting reads for this adjacency
#             # n_supporting = len(rnames)
#             n_chimeric = 0
#             regions = [(chrom1, pos1), (chrom2, pos2)]
#             for chrom, pos in regions:
#                 if chrom not in chroms_proc: 
#                     continue # only take reads from transposon etc
#                 reads = bam.fetch(chrom, pos-1, pos) # all bam reads from this breakpoint
#                 for read in reads:
#                     if read.qname in rnames: # is this read a supporter?
#                         secondaries = get_secondaries(read)
#                         if len(secondaries) > 0:
#                             n_chimeric += 1
#                             for secondary in secondaries:
#                                 chrom_2nd, pos_2nd, strand_2nd, cigar_2nd = secondary.split(',')[:4]
#                                 field = [sample, read.qname, chrom1, pos1, chrom2, pos2, row['type'], chrom_2nd, pos_2nd, strand_2nd]
#                                 field = [str(s) for s in field]
#                                 line = ','.join(field)
#                                 out.write(line + '\n')

#                         ncseq = sample_seqs[read.qname] # non-clipped seq
#                         ncqual = sample_quals[read.qname] # non-clipped seq
#                         # need to revcomp fastq seq if read is reverse mapped
#                         if read.is_reverse:
#                             ncseq = str(Seq(ncseq).reverse_complement())
#                             ncqual = ncqual[::-1]
#                         # check 5'/3' of read
#                         seq_5p, seq_3p, qual_5p, qual_3p = '', '', '', ''
#                         cigar_5p = read.cigar[0]
#                         cigar_3p = read.cigar[-1]
#                         if cigar_5p[0] in clipped_tags:
#                             # seq_5p = read.seq[:cigar_5p[1]]
#                             seq_5p = ncseq[:cigar_5p[1]]
#                             qual_5p = ncqual[:cigar_5p[1]]
#                         if cigar_3p[0] in clipped_tags:
#                             # seq_3p = read.seq[-cigar_3p[1]:]
#                             seq_3p = ncseq[-cigar_3p[1]:] 
#                             qual_3p = ncqual[-cigar_3p[1]:] 

#                         # if len(seq_5p) and len(seq_3p):
#                         field = [read.qname, read.reference_name, read.pos+1, seq_5p, seq_3p, qual_5p, qual_3p]
#                         seqdf.loc[six] = field
#                         six += 1
                        
#                         if read.qname == '8d49baf4-f7be-46e5-b495-3e6df968f6c1':
#                             _sample_seqs = sample_seqs
#                             _read = read
#                             _field = field
#                     # fi
#                 # for chrom, pos
#             n_chimerics.append(n_chimeric)
#         # svs['n_supportings'] = n_supportings
#         svs['n_chimeric'] = n_chimerics # chimeric reads, that support the sv
#         # svs_path = f'/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/svs/{sample}.svs.csv'
#         # svs.to_csv(svs_path, index=False)
#         # sample2seqdf[sample] = seqdf
```

### debugging a seq mismatch

```python
_read.qname == '8d49baf4-f7be-46e5-b495-3e6df968f6c1'
```

```python
cigar_5p, cigar_3p = _read.cigar[0], _read.cigar[-1]
```

```python
cigar_5p[1]
```

```python
ncseq[:cigar_5p[1]][:20], ncseq[:cigar_5p[1]][-20:]
```

```python
_sample_seqs[_read.qname][:cigar_5p[1]][:20], _sample_seqs[_read.qname][:cigar_5p[1]][-20:]
```

```python
ncseq[:cigar_3p[1]][:20], ncseq[:cigar_3p[1]][-20:]
```

```python
_sample_seqs[_read.qname][:cigar_3p[1]][:20], _sample_seqs[_read.qname][:cigar_3p[1]][-20:]
```

```python
chrseq = genome['PBEF1NeoTransposon'][:].seq.upper()
```

```python
_read.pos
```

```python
_read.seq[:cigar_5p[1]][:20], _read.seq[:cigar_5p[1]][-20:]
```

```python
_read.seq[:20], _read.seq[-20:]
```

```python
chrseq[_read.pos+1-20:_read.pos+1], str(Seq(chrseq[_read.pos+1-20:_read.pos+1]).reverse_complement())
```

```python
_read.qlen
```

```python
len(ncseq)
```

```python
len(_read.seq)
```

```python
'ATTT'[:3:-1]
```

```python
_read.is_reverse
```

```python
ncseq == str(Seq(_read.seq).reverse_complement())
```

```python
read.seq[:cigar_5p[1]][:10]
```

```python
ncseq = _sample_seqs[_read.qname] # non-clipped seq
```

```python
_field
```

### load fastq per sample into dict *


#### save in .fa format for blast

```python
out_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/fasta'
sv_dir = '/juno/work/shah/users/chois7/tickets/sarcoma89/14474/results/svs'
for sample, (rname2seq, rname2qual) in fqs.items():
    svs = data[sample]
    svs = svs[(svs['chrom1']=='PBEF1NeoTransposon') | (svs['chrom2']=='PBEF1NeoTransposon')]
    svdf = svs[['chrom1', 'pos1', 'chrom2', 'pos2', 'type']]
    n_supports = []
    with open(f'{out_dir}/{sample}.rnames.fasta', 'w') as out:
        for riw, row in svs.iterrows(): # row: one SV
            chrom1, pos1, chrom2, pos2, svtype, svlength, _svid, _ref, _alt, _qual, _filter, _info, _sample, rnames = row.squeeze()
            sv_identity = [chrom1, pos1, chrom2, pos2, svtype]
            sv_identity_str = '|'.join([str(_) for _ in sv_identity])
            # print(sample, chrom1, pos1, chrom2, pos2, svtype)
            rnames = rnames.split(',')
            n_support = len(rnames)
            n_supports.append(n_support)
            
            for rname in rnames: # supporting & chimreric reads into fasta
                seq = rname2seq[rname]
                header = f'>{rname}|{sv_identity_str}'
                out.write(f'{header}\n')
                out.write(f'{str_to_fasta(seq)}\n')
                # break
    svdf['number_of_chimeric_reads'] = n_supports
    svdf.to_csv(f'{sv_dir}/{sample}.sv.csv', index=False)
```

# Motif discovery from clipped reads

```python
seqdf.head(1)
```

```python
from Bio import motifs
from Bio.Seq import Seq
```

```python
min_length = 20
seq5ps = seqdf["read_5'-clipped_seq"].tolist()
seq5ps = [s for s in seq5ps if len(s) >= min_length]
seq5ps = [s[::-1] for s in seq5ps]
seq5ps = [s[:20] for s in seq5ps]
seq3ps = seqdf["read_3'-clipped_seq"].tolist()
seq3ps = [s for s in seq3ps if len(s) >= min_length]
seq3ps = [s[:20] for s in seq3ps]
seqs = seq5ps + seq3ps
```

```python
sequences5 = [Seq(s) for s in seq5ps]
motif5 = motifs.create(sequences5)
consensus5 = motif5.consensus
# information_content = motif.pwm.mean_information(content=True)
```

```python
seq5ps
```

```python
consensus
```

```python
for i, s in enumerate(seq5ps):
    print(f'>seq{i+1}')
    print(s)
```

```python
for i, s in enumerate(seq3ps):
    print(f'>3\'-seq{i+1}')
    print(s)
```
