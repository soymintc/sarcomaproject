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

# Settings

```python
%matplotlib inline
import os
import re
import glob
import matplotlib
import numpy as np

from matplotlib.collections import LineCollection
from matplotlib import patches
from matplotlib.path import Path
import matplotlib.pyplot as plt
from collections import defaultdict

import pandas as pd
import pysam

import dvartk
from scgenome.cnplot import plot_cell_cn_profile
import wgs_analysis.plots.cnv as cnv
import wgs_analysis.plots.utils

matplotlib.rcParams.update({'font.size': 13})
```

```python
from matplotlib.patches import ConnectionPatch
```

## modules


```python
def get_gene_repr_exons(gtf, gene, lenient=False):
    transcript = get_repr_transcript_id(gtf, gene, lenient=lenient)
    exons = get_transcript_exons(gtf, transcript)
    return exons
```

```python
def get_transcript_exons(gtf, transcript_id):
    exons = gtf[
        (gtf['transcript_id']==transcript_id) &
        (gtf['Feature'].isin(['exon', 'CDS']))
    ]
    return exons
```

```python
def get_repr_transcript_id(gtf, gene_name, lenient=False):
    # gene_symbol = 'BCL2'
    if not lenient:
        transcripts = (
            gtf
            .query(f"gene_name == '{gene_name}'")
            .query("Feature == 'transcript'")
            .query("transcript_biotype == 'protein_coding'")
            .query("transcript_support_level == '1'")
        ).copy()
    else:
        transcripts = (
            gtf
            .query(f"gene_name == '{gene_name}'")
            .query("Feature == 'transcript'")
        ).copy()
    transcript_id = None
    if transcripts.shape[0] > 0:
        transcripts['length'] = transcripts['End'] - transcripts['Start']
        transcript = transcripts.sort_values(by=['length'], ascending=False).iloc[0]
        transcript_id = transcript['transcript_id']
    return transcript_id
```

```python
def get_transcript_exons(gtf, transcript_id):
    exons = gtf[
        (gtf['transcript_id']==transcript_id) &
        (gtf['Feature'].isin(['exon', 'CDS']))
    ]
    return exons
```

```python
def add_transcript_rectangles(ax, row, transcript_s, transcript_e, y_offset, strand, add_annot=False):
    box_color = '#090A76'
    strand_color = 'blue'
    x_direct, y_direct = -1, 0
    if strand == '+': 
        x_direct = 1
        strand_color = 'red'
    x_min, x_max = ax.get_xlim()
    y_offset_factor = 3
    y_offset *= y_offset_factor
    chrom = row['Chromosome']
    exon_s, exon_e = row['Start'], row['End']
    feature = row['Feature']
    rect_y, rect_height = 1+y_offset, 1
    transcript_line_y = rect_y + 0.5 * rect_height
    text_y = rect_y - rect_height + 3
    if feature == 'CDS':
        rect_y = 0.5 + y_offset
        rect_height = 2
    exon_len = exon_e - exon_s
    rect_x = exon_s #max(exon_s, start)
    
    # for each exon / CDS
    ax.add_patch(matplotlib.patches.Rectangle((rect_x, rect_y), exon_len, rect_height, 
                                              linewidth=1, edgecolor=strand_color, color=box_color, zorder=1))
    if add_annot: 
        quiver_interval = int((x_max - x_min) / 50)
        transcript_len = transcript_e - transcript_s
        # print(transcript_s, transcript_e, quiver_interval, type(transcript_s), type(transcript_e), type(quiver_interval))
        x_quiver_starts = np.arange(transcript_s, transcript_e, quiver_interval)[1:-1]
        n_quivers = len(x_quiver_starts)
        y_quiver_starts = [transcript_line_y] * n_quivers
        x_quiver_directs = [x_direct] * n_quivers
        y_quiver_directs = [y_direct] * n_quivers
        ax.quiver(x_quiver_starts, y_quiver_starts, x_quiver_directs, y_quiver_directs,
                  color=box_color,
                  width=0.001, 
                  headwidth=10, 
                  headaxislength=3,
                  headlength=6, 
                  pivot='tip',
                  scale=200,
                 )
        ax.plot([transcript_s+50, transcript_e-100], 
                [transcript_line_y, transcript_line_y], 
                linewidth=1, color=box_color, zorder=0)
        text_x = (transcript_s + transcript_e) / 2
        if text_x < x_min:
            text_x = x_min
        if text_x > x_max:
            text_x = x_max
        ax.text(x=text_x, y=text_y, s=row['gene_name'], color=box_color, horizontalalignment='center', fontdict={'size':11})
```

```python
def plot_gene_annotations(ax, gene_exons):
    gene2offset = {}
    transcript_start_ends = []
    for _gene, exons in gene_exons.items():
        strands = exons['Strand'].unique()
        assert len(strands) == 1, strands
        strand = strands[0]
        transcript_s, transcript_e = exons['Start'].min(), exons['End'].max()
        transcript_start_ends.append((transcript_s, transcript_e))
        transcript_offset = 0
        # if len(transcript_start_ends) >= 2:
        #     transcript_offset = calc_transcript_offset(transcript_start_ends)
        gene2offset[_gene] = transcript_offset
        for i, (rix, row) in enumerate(exons.iterrows()):
            add_annot = (i==0)
            # print(_gene, add_annot)
            add_transcript_rectangles(ax, row, transcript_s, transcript_e, transcript_offset, strand, add_annot=add_annot) # add rectangles and a line
    return gene2offset
```

#### savana modules

```python
class Breakpoint:
    def __init__(self, brk):
        self.chrom = brk['chrom']
        self.pos = brk['pos']
        self.self_id = brk['breakend'] # ID
        brk_ix = int(self.self_id[-1])
        assert brk_ix in {1, 2}, f'ERROR: brk_ix = {brk_ix}'
        self.ref = brk['ref']
        self.alt = brk['alt']
        self.info = brk['INFO']
        bp_notation = re.search('BP_NOTATION=([^;]+);', self.info).groups()[0]
        self.adjacency_id = brk['adjacency'] # ID
        self.mate_id = f'{self.adjacency_id}_{3-brk_ix}'
        
        _strand_combination = {'++', '--', '+-', '-+'}
        if bp_notation == '<INS>':
            self.strand = None
        elif bp_notation in _strand_combination:
            self.strand = bp_notation[brk_ix-1]
        else:
            raise ValueError(f'ERROR: bp_notation = {bp_notation}')

class Adjacency:
    def __init__(self, brks): # brks <- paired dataframe
        assert brks.shape[0] == 2, brks
        brk1, brk2 = brks.iloc[0], brks.iloc[1]
        self.brk1 = Breakpoint(brk1)
        self.brk2 = Breakpoint(brk2)
        self.type = 'n/a'
        
        self.type = self.get_svtype()
        self.length = abs(self.brk2.pos - self.brk1.pos)
        if self.type == 'translocation': 
            self.length = int(3e9)
    
    def get_svtype(self): # N-> <-N // <-N N-> // N<- N<- // ->N ->N
        if self.brk1.chrom != self.brk2.chrom: # - <TRA>
            return 'translocation'
        if (self.brk1.strand, self.brk2.strand) == ('+', '+'):
            return 'inversion'
        elif (self.brk1.strand, self.brk2.strand) == ('-', '-'):
            return 'inversion'
        elif (self.brk1.strand, self.brk2.strand) == ('+', '-'):
            return 'deletion'
        elif (self.brk1.strand, self.brk2.strand) == ('-', '+'):
            return 'duplication'
        else:
            raise ValueError(f'ERROR: (strand1, strand2) = ({self.brk1.strand}, {self.brk1.strand})')
```

```python
def resolve_breakpoints(df):
    df = df.copy()
    chroms = [str(c) for c in range(1, 22+1)] + ['X', 'Y']
    svsv_cols = ['chromosome_1', 'position_1', 'strand_1', 'chromosome_2', 'position_2', 'strand_2',
                 'type', 'length']
    svsv = pd.DataFrame(columns=svsv_cols) # savana sv

    df['adjacency'] = df['breakend'].str.rsplit('_', 1, expand=True)[0]
    
    svtype_cnt = defaultdict(int)
    for i, (ix, brks) in enumerate(df.groupby('adjacency')):
        brks_in_adj = brks.shape[0]
        if brks_in_adj == 2:
            brk1, brk2 = brks.iloc[0], brks.iloc[1]
            brk1.chrom = brk1.chrom.replace('chr', '')
            brk2.chrom = brk2.chrom.replace('chr', '')
            if brk1.chrom not in chroms: continue
            if brk2.chrom not in chroms: continue
            if brk1.chrom == brk2.chrom:
                assert brk1.pos < brk2.pos, (brk1.pos, brk2.pos)
            adj = Adjacency(brks)
            svtype_cnt[adj.type] += 1
            line = [adj.brk1.chrom, adj.brk1.pos, adj.brk1.strand, adj.brk2.chrom, adj.brk2.pos, adj.brk2.strand, adj.type, adj.length]
        elif brks_in_adj == 1:
            brk = brks.squeeze()
            brk1 = Breakpoint(brk)
            assert brk['alt'] == '<INS>', brk
            svtype = 'insertion'
            match = re.search('INSSEQ=([A-Z]+);', brk['INFO'])
            insseq = match.groups()[0]
            svlength = len(insseq)
            svtype_cnt[svtype] += 1
            line = [brk1.chrom, brk1.pos, brk1.strand, brk1.chrom, brk1.pos, brk1.strand, svtype, svlength]
        svsv.loc[i] = line

        # N-> <-N // <-N N-> // N-> N-> // <-N <-N
    # svsv['type'] = svsv['type'].replace('insertion', 'ins').replace('deletion', 'del').replace('inversion', 'inv').replace('duplication', 'dup')
    return svsv 
```

## plot

```python
def get_fusions_path(patient):
    paths_path = '/juno/work/shah/users/chois7/tickets/sarcoma89/sarcoma.ONT-NANOSEQ.tsv'
    paths = pd.read_table(paths_path)
    paths = paths[paths['sample_category'] == 'TUMOR']
    paths = paths[paths['technique_category'] == 'RNA']
    paths = paths[paths['isabl_patient_id'] == patient]
    paths = paths[paths['result_type'] == 'command_log']
    if paths.shape[0] == 0: 
        return None
    log_path = paths.iloc[0]['result_filepath']
    directory = os.path.split(log_path)[0]
    path = f'{directory}/results/jaffal/jaffa_results.csv'
    return path
```

```python
def get_svs_path(patient):
    fs = glob.glob(f'/juno/work/shah/users/chois7/tickets/savana/pipeline/results/{patient}/savana/*.sorted.somatic.sv_breakpoints.lenient.vcf')
    if len(fs) == 1:
        return fs[0]
    else:
        return None
```

```python
def get_fusions(fusions_path):
    fusions = pd.read_csv(fusions_path)
    fusions = fusions[
        (fusions['classification'] == 'HighConfidence') &
        (fusions['spanning reads'] >= 5)
    ]
    return fusions

def get_svs(svs_path):
    vcf_cols = ['chrom', 'pos', 'breakend', 'ref', 'alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample']
    proc_cols = ['chrom', 'pos', 'breakend', 'ref', 'alt']

    svs = pd.read_table(svs_path, comment='#', names=vcf_cols)
    svs = resolve_breakpoints(svs)
    _svs = svs.copy()
    svs = svs[
        svs['type'].isin({'inv', 'translocation'}) # filter in INV TRA only
    ]
    return svs
```

### get fusion data

```python
paths_path = '/juno/work/shah/users/chois7/tickets/sarcoma89/sarcoma.ONT-NANOSEQ.tsv'
paths = pd.read_table(paths_path)
patients = paths['isabl_patient_id'].unique()
```

```python
svs_patients = set()
for patient in patients:
    svs_path = get_svs_path(patient)
    if svs_path == None: continue
    if os.path.exists(svs_path):
        svs_patients.add(patient)
        
fusions_patients = set()
for patient in patients:
    fusions_path = get_fusions_path(patient)
    if fusions_path == None: 
        continue
    if os.path.exists(fusions_path):
        fusions_patients.add(patient)
```

```python
both_patients = svs_patients & fusions_patients
```

### get sv data


## plot fusion + matching sv


### gene plot

```python
def init_paired_axes(gene_exons, gene_chrom, gene1, gene2):
    gene1_start, gene1_end = gene_exons[gene1]['Start'].min(), gene_exons[gene1]['End'].max()
    gene2_start, gene2_end = gene_exons[gene2]['Start'].min(), gene_exons[gene2]['End'].max()
    chrom1, chrom2 = gene_chrom[gene1], gene_chrom[gene2]
    fig, axes = plt.subplots(2, 1, figsize=(10, 5))
    ax1, ax2 = axes
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.tick_top()
    ax1.set_xlabel(chrom1, fontsize=14)
    ax2.set_xlabel(chrom2, fontsize=14)

    offset1 = (gene1_end - gene1_start) / 30
    offset2 = (gene2_end - gene2_start) / 30
    ax1.set_xlim((gene1_start - offset1, gene1_end + offset1))
    ax2.set_xlim((gene2_start - offset2, gene2_end + offset2))
    ax1.set_ylim((-5, 8)); ax2.set_ylim((-5, 8));

    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.spines['top'].set_position(('outward', 5))
    ax2.spines['bottom'].set_position(('outward', 5))

    for ax in axes:
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False) #
        ax.set_yticks([])
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%d'))
        ax.tick_params(axis='both', which='major', labelsize=13)
        
    return fig, axes
```

```python
def get_fusion_matching_svs(svs, chrom1, x1_min, x1_max, chrom2, x2_min, x2_max):
    forward = svs[
        (
            ((svs['chromosome_1'] == chrom1) & (svs['chromosome_2'] == chrom2)) &
            ((svs['position_1'] < x1_max) & (svs['position_1'] >= x1_min) & (svs['position_2'] < x2_max) & (svs['position_2'] >= x2_min))
        )
    ].copy()
    backward = svs[
        (
            ((svs['chromosome_2'] == chrom1) & (svs['chromosome_1'] == chrom2)) &
            ((svs['position_2'] < x1_max) & (svs['position_2'] >= x1_min) & (svs['position_1'] < x2_max) & (svs['position_1'] >= x2_min))
        )
    ].copy()
    if backward.shape[0] > 0:
        backward = backward.rename(columns={
            'chromosome_1':'chromosome_2', 
            'chromosome_2':'chromosome_1', 
            'position_1':'position_2', 
            'position_2':'position_1', 
            'strand_1':'strand_2',
            'strand_2':'strand_1',
        })
    return pd.concat([forward, backward])
```

```python
def get_gene_exons_and_chroms(fusion):
    gene_exons = dict()
    gene_chrom = dict()
    for gene in [gene1, gene2]:
        exons = get_gene_repr_exons(gtf, gene)
        uniq_chroms = exons['Chromosome'].unique()
        if len(uniq_chroms) == 0:
            exons = get_gene_repr_exons(gtf, gene, lenient=True)
            uniq_chroms = exons['Chromosome'].unique()
            if len(uniq_chroms) == 0:
                return None, None
        assert len(uniq_chroms) == 1, (gene, uniq_chroms)
        chrom = uniq_chroms[0]
        gene_exons[gene] = exons
        gene_chrom[gene] = chrom
    return gene_exons, gene_chrom
```

```python
# gene1, gene2 = fusion['fusion genes'].split(':')
```

```python
# gene_exons = dict()
# gene_chrom = dict()
# for gene in [gene1, gene2]:
#     exons = get_gene_repr_exons(gtf, gene)
#     uniq_chroms = exons['Chromosome'].unique()
#     assert len(uniq_chroms) == 1, uniq_chroms
#     chrom = uniq_chroms[0]
#     gene_exons[gene] = exons
#     gene_chrom[gene] = chrom
```

## parse gtf

```python
gtf_path = '/juno/work/shah/reference/dlp-ref-GRCh38/human/Homo_sapiens.GRCh38.93.gtf'

```

```python
# gtf_path = '/home/chois7/chois7/datasets/hg38/GRCh38/gencode.v31.annotation.gtf'
# gtf_path = '/juno/work/shah/users/chois7/datasets/GRCh38/gencode.v43.annotation.gtf'
gtf_path = '/juno/work/shah/users/chois7/GDAN/DLBCL/breakpoint_qc/annotsv/resources/Homo_sapiens.GRCh38.93_chr.gtf'
gtf = pyranges.read_gtf(gtf_path, as_df=True)

```

```python
id2name = dict()
gene_ids_short = gtf.gene_id.str.extract('([^\.]+)').squeeze().values
gene_ids = gtf.gene_id.values
transcript_ids = gtf.transcript_id.values
gene_names = gtf.gene_name.values

id2name.update(dict(zip(gene_ids_short, gene_names)))
id2name.update(dict(zip(gene_ids, gene_names)))
id2name.update(dict(zip(transcript_ids, gene_names)))
```

```python
for patient in both_patients:
    # if patient != 'P-0068227': continue
    if patient != 'P-0021868': continue
    print(patient)
    svs_path = get_svs_path(patient)
    fusions_path = get_fusions_path(patient)
    svs = get_svs(svs_path)
    fusions = get_fusions(fusions_path)
    
    for fix, fusion in fusions.iterrows():
        gene1, gene2 = fusion['fusion genes'].split(':')
        gene_exons, gene_chrom = get_gene_exons_and_chroms(fusion)
        if gene_exons == None: continue

        chrom1, chrom2 = gene_chrom[gene1], gene_chrom[gene2]    
        gene1_start, gene1_end = gene_exons[gene1]['Start'].min(), gene_exons[gene1]['End'].max()
        gene2_start, gene2_end = gene_exons[gene2]['Start'].min(), gene_exons[gene2]['End'].max()
        offset1 = (gene1_end - gene1_start) / 30
        offset2 = (gene2_end - gene2_start) / 30
        x1_min, x1_max = gene1_start - offset1, gene1_end + offset1
        x2_min, x2_max = gene2_start - offset2, gene2_end + offset2

        tsvs = get_fusion_matching_svs(svs, chrom1, x1_min, x1_max, chrom2, x2_min, x2_max) # target SVs
        if tsvs.shape[0] > 0: 
        # if gene1.startswith('CMTM'):
            # print('tsv shape >0', gene1, gene2)
            fig, (ax1, ax2) = init_paired_axes(gene_exons, gene_chrom, gene1, gene2)
            # x1_min, x1_max = ax1.get_xlim()
            # x2_min, x2_max = ax2.get_xlim()

            # plot genes
            _ = plot_gene_annotations(ax1, {gene1: gene_exons[gene1]})
            _ = plot_gene_annotations(ax2, {gene2: gene_exons[gene2]})

            fusion_color = 'tab:green'
            sv_color = 'tab:red'

            # fusion
            fpos1, fpos2 = fusion['base1'], fusion['base2']
            fusion_y0, fusion_y1 = -1, 4
            ax1.plot([fpos1, fpos1], [fusion_y0, fusion_y1], linewidth=1, color=fusion_color)
            ax2.plot([fpos2, fpos2], [fusion_y0, fusion_y1], linewidth=1, color=fusion_color)
            fusion_conn = ConnectionPatch(xyA=[fpos1, fusion_y0], xyB=[fpos2, fusion_y1], coordsA="data", coordsB="data", 
                                          axesA=ax1, axesB=ax2, color=fusion_color, lw=1)
            ax2.add_artist(fusion_conn)

            fusion1_arrow_len = (x1_max - x1_min) / 40 * (-(fusion.strand1 == '-') + 0.5) * 2
            fusion2_arrow_len = (x2_max - x2_min) / 40 * (-('-' == '-') + 0.5) * 2
            ax1.annotate('', xytext=(fpos1-fusion1_arrow_len, fusion_y1), xy=(fpos1, fusion_y1),
                         arrowprops=dict(arrowstyle='->', color=fusion_color, lw=1, ls='-'))
            ax2.annotate('', xytext=(fpos2, fusion_y1), xy=(fpos2+fusion2_arrow_len, fusion_y1),
                         arrowprops=dict(arrowstyle='->', color=fusion_color, lw=1, ls='-'))

            # SV
            if tsvs.shape[0] > 0: 
                sv_y0, sv_y1 = -1.5, 4.5
                for rix, row in tsvs.iterrows():
                    sv_chrom1, sv_pos1, sv_strand1 = row['chromosome_1'], row['position_1'], row['strand_1']
                    sv_chrom2, sv_pos2, sv_strand2 = row['chromosome_2'], row['position_2'], row['strand_2']
                    sv_type, sv_length = row['type'], row['length']
                    ax1.plot([sv_pos1, sv_pos1], [sv_y0, sv_y1], linewidth=1, color=sv_color)
                    ax2.plot([sv_pos2, sv_pos2], [sv_y0, sv_y1], linewidth=1, color=sv_color)
                    sv_conn = ConnectionPatch(xyA=[sv_pos1, sv_y0], xyB=[sv_pos2, sv_y1], coordsA="data", coordsB="data", 
                                              axesA=ax1, axesB=ax2, color=sv_color, lw=1)
                    ax2.add_artist(sv_conn)
                    # print(row.strand_1, row.strand_2)
                    sv1_arrow_factor = -2 * (int(row.strand_1 == '+') - .5)
                    sv2_arrow_factor = -2 * (int(row.strand_2 == '+') - .5)
                    sv1_arrow_len = sv1_arrow_factor * (x1_max - x1_min) / 40
                    sv2_arrow_len = sv2_arrow_factor * (x2_max - x2_min) / 40 
                    ax1.annotate('', xytext=(sv_pos1, sv_y1), xy=(sv_pos1+sv1_arrow_len, sv_y1),
                                 arrowprops=dict(arrowstyle='<-', color=sv_color, lw=1, ls='-'))
                    ax2.annotate('', xytext=(sv_pos2, sv_y1), xy=(sv_pos2+sv2_arrow_len, sv_y1),
                                 arrowprops=dict(arrowstyle='<-', color=sv_color, lw=1, ls='-'))

                if gene1.startswith('CTCM'): raise ValueError()
```
