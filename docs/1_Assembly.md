Assembling complete pQBR plasmids using long and short reads
================
revised 2025

## Overview of approach

The plasmids were sequenced inside *P. putida* UWC1, a derivative of *P.
putida* KT2440. Initial Illumina sequencing was performed by MicrobesNG
and provided by Damian Rivett (Manchester Metropolitan University).
These data were analysed in a preliminary attempt to generate closed
plasmid sequences as follows:

1.  Map all reads against the reference Pseudomonas chromosomes using
    `bwa-mem` 0.7.17-r1188.
2.  Collect all reads that did not map to the chromosome and assemble
    *de novo* using SPAdes v3.15.5 with k-mer length of 21, 33, 55, 77.
3.  Long (\>5kb), deep (\>10x coverage) contigs were identified as
    putative pQBR plasmid contigs, and were used for initial comparative
    analysis.

Note that this analysis struggles with plasmids which have regions
matching the *P. putida* UWC1 chromosome (see [Hall et
al. 2015](https://pubmed.ncbi.nlm.nih.gov/25969927/)).

Alongside this analysis, putative plasmids were identified by the whole
genome *de novo* SPAdes assemblies that were provided by MicrobesNG by
visualising the assembly graphs in
[Bandage](https://rrwick.github.io/Bandage/) and identifying and
extracting closed circular contigs.

For plasmids that did not produced closed genomes by either of these
methods, we followed up with ONT sequencing of pools of samples,
grouping distinct plasmids based on the Illumina contigs, to provide
scaffolds that could then be polished with the Illumina reads.

Note: intermediate files are not included in the GitHub.

### Preliminary analysis: assemble unmapped reads

Map reads against reference Pseudomonas chromosomes, extract and
assemble unmapped reads.

``` bash
find ./reads/*_1_trimmed.fastq.gz > 1.tmp
find ./reads/*_2_trimmed.fastq.gz > 2.tmp
awk -v FS="_" '{print $2}' 1.tmp > n.tmp
paste n.tmp 1.tmp 2.tmp > ctrlsamp.txt
rm *.tmp

# make directories for outputs

mkdir ./unmapped_assemblies
mkdir ./unmapped_reads

REF=./ref/SBW25_KT2440_chromosomes.fasta
r1=./scratch/R1.fastq.gz
r2=./scratch/R2.fastq.gz

# use the iteration file to go through for each sample:
#  map using BWA, use samtools to take reads that don't map
#  assemble these reads using SPAdes

cat ctrlsamp.txt | while read SAMPLE R1 R2
do
echo "---- NEXT SAMPLE IS ${SAMPLE} ----"
bwa mem $REF $R1 $R2 \
  | samtools view -hb \
  | samtools sort \
  > mapped.bam

samtools index mapped.bam

samtools view -hf 4 mapped.bam \
  | awk '$3=="\*" {print $1}' \
  | uniq \
  > ./unmapped_reads/${SAMPLE}_unmapped.txt

# get the corresponding reads from the fastq files
seqtk subseq $R1 ./unmapped_reads/${SAMPLE}_unmapped.txt | gzip - > $r1
seqtk subseq $R2 ./unmapped_reads/${SAMPLE}_unmapped.txt | gzip - > $r2

# Assemble
spades.py -o ./unmapped_assemblies/${SAMPLE} \
  -k 21,33,55,77 --careful \
  -1 $r1 -2 $r2
done
```

The output files are in the scaffolds.fasta file in each output
directory.

First, see how many contigs (and what size) per plasmid. Make a summary
file that contains plasmid name, contig number, length, and coverage.

``` bash
grep "^>" ./unmapped_assemblies/UWC1*/scaffolds.fasta \
  | awk -v FS="[/_]" '{print $0, $4, $6, $8, $10}' > ./unmapped_assemblies/summary.txt
```

Plot data from this summary table.

``` r
summary <- read.table("./unmapped_assemblies/summary.txt", header=FALSE, col.names = c("seq","node","length","cov"))

ggplot(summary, aes(x=log10(length), y=log10(cov), colour=node)) + geom_point() + facet_wrap(~seq) +
  geom_hline(yintercept=log10(10), linetype="dotted") + geom_vline(xintercept=log10(5000), linetype="dotted")
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

This gives a good overview summary. It shows that for some samples,
there are just one or two contigs that are big, and have reasonable
coverage. These look most promising. Other samples have not assembled so
well. Note there were two non-pQBR plasmids included in the sequencing
run, UWC1UG451 and UWC1UG452, that are not analysed further.

How long are the assemblies for each sample?

``` r
summary %>% group_by(seq) %>% 
  summarise(n_contigs = n(),
    total_length = sum(length), 
            min_cov = min(cov), 
            max_cov = max(cov),
            mean_cov = mean(cov),
            median_cov = median(cov)) %>%
  arrange(total_length) %>%
  kable()
```

| seq         | n_contigs | total_length |   min_cov | max_cov |    mean_cov | median_cov |
|:------------|----------:|-------------:|----------:|--------:|------------:|-----------:|
| UWC1pQBR1   |        40 |        54554 |  1.050938 |    2707 |   70.527355 |   2.218342 |
| UWC1pQBR58  |         2 |        73781 | 13.864192 |     614 |  313.932096 | 313.932096 |
| UWC1pQBR132 |         4 |       139163 |  0.894602 |    3667 |  944.063399 |  54.179498 |
| UWC1pQBR55  |         4 |       140164 |  0.876033 |    1336 |  344.172217 |  19.906419 |
| UWC1pQBR28  |         3 |       140514 | 12.820513 |    2801 |  950.272285 |  36.996341 |
| UWC1pQBR53  |        39 |       154114 |  0.838384 |    1945 |   53.023001 |   1.118971 |
| UWC1pQBR105 |        72 |       164041 |  0.663399 |    6833 |  104.733654 |   1.044226 |
| UWC1pQBR57  |         7 |       308528 |  0.839721 |    6964 | 1014.584485 |  25.051282 |
| UWC1pQBR30  |         3 |       310667 |  0.980282 |    2951 | 1003.000029 |  57.019805 |
| UWC1pQBR150 |        12 |       310923 |  0.857143 |    2244 |  193.666878 |   1.278429 |
| UWC1pQBR56  |        13 |       311401 |  0.906915 |    5668 |  441.402890 |   1.048338 |
| UWC1pQBR102 |        11 |       315655 |  0.875445 |    2775 |  258.509623 |   1.038922 |
| UWC1pQBR47  |       122 |       337568 |  0.813620 |   10231 |   86.103346 |   1.447420 |
| UWC1pQBR24  |        19 |       386226 |  0.745575 |    2316 |  133.374671 |  15.022689 |
| UWC1pQBR23  |        16 |       386799 |  0.925532 |   10055 |  660.396746 |  40.833683 |
| UWC1pQBR43  |         4 |       393112 | 87.128439 |    4516 | 1221.934289 | 142.304359 |
| UWC1pQBR50  |         6 |       425289 |  0.860215 |    8759 | 1494.580102 |  54.817217 |
| UWC1pQBR5   |         5 |       425749 |  0.869674 |   12831 | 2583.335300 |  30.538462 |
| UWC1pQBR103 |         6 |       426251 |  1.167785 |   12950 | 2168.297131 |  12.821179 |
| UWC1pQBR4   |        14 |       428654 |  0.908078 |   11437 |  827.814143 |   1.166521 |
| UWC1pQBR11  |        21 |       438268 |  0.872396 |    7015 |  338.571414 |   0.953297 |
| UWC1pQBR106 |       105 |       441287 |  0.803150 |    4587 |   56.057636 |   1.857143 |
| UWC1pQBR49  |         8 |       449073 |  0.886957 |    8773 | 1151.128131 |  56.026832 |
| UWC1pQBR124 |       143 |       485914 |  0.848315 |    4741 |   37.140937 |   1.064417 |
| UWC1pQBR26  |        60 |       577309 |  0.946203 |    9410 |  222.874736 |  50.743658 |
| UWC1pQBR51  |       161 |       696326 |  0.711921 |    4885 |   32.958241 |   2.379043 |
| UWC1pQBR127 |       956 |      1115580 |  0.810277 |    6248 |    8.988482 |   1.365425 |
| UWC1pQBR8   |      1127 |      1380714 |  0.613636 |    3651 |    4.991678 |   1.464799 |

For each assembly, pull out the big (\>5 kb) high coverage (\>10x)
contigs.

Plot the size of the assemblies when just these contigs are included.

``` r
summary %>% filter(length > 5000 & cov > 10 & !(seq %in% c("UWC1UG451", "UWC1UG452"))) %>%
  mutate(label = gsub("UWC1(pQBR[0-9]+)", "\\1", seq)) %>%
  group_by(label) %>%
  summarise(total_length = sum(length)) %>%
  ggplot(aes(x=reorder(label, total_length), y=total_length)) +
  geom_bar(stat="identity") + scale_x_discrete(name = "") + theme_pub() +
  labs(y="total length of big high-coverage contigs") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Pull out these contigs.

``` bash
mkdir ./unmapped_assemblies/long_deep_unmapped

cat ctrlsamp.txt | while read SAMPLE R1 R2 
do
awk -v sam="$SAMPLE" '$1 == sam && $3 > 5000 && $4 > 10 {print $2}' \
  ./unmapped_assemblies/summary.txt \
  | sed 's/^/NODE_/' | sed 's/$/_/' > nodes.tmp
grep -f nodes.tmp ./unmapped_assemblies/${SAMPLE}/scaffolds.fasta \
  | sed 's/>//g' \
  > nodenames.tmp
seqtk subseq ./unmapped_assemblies/${SAMPLE}/scaffolds.fasta nodenames.tmp \
  | sed "s|>|>${SAMPLE}_|" \
  > ./unmapped_assemblies/long_deep_unmapped/${SAMPLE}_unmapped_assembly.fasta
done
```

Calculate mash distances between all sequences. For this, sequences
within each .fasta file need to be concatenated. Also output a file that
includes details on the number of nodes for each sequence.

``` bash
grep "UWC1pQBR" ctrlsamp.txt | while read SAMPLE R1 R2 
do
echo ">${SAMPLE}_unmapped_assembly" >> ./ref/ldu_pQBR.fasta
grep -v "^>" ./unmapped_assemblies/long_deep_unmapped/${SAMPLE}_unmapped_assembly.fasta \
| tr -d "\n" | fold >> ./unmapped_assemblies/ldu_pQBR.fasta
echo "\n" >> ./unmapped_assemblies/ldu_pQBR.fasta
done

grep "UWC1pQBR" ctrlsamp.txt | while read SAMPLE R1 R2 
do
grep "^>" ./unmapped_assemblies/long_deep_unmapped/${SAMPLE}_unmapped_assembly.fasta \
  | awk -v sam=$SAMPLE -v FS="_" '{print sam, $3, $5, $7}'
done > ./unmapped_assemblies/ldu_pQBR_details.txt
```

What are the details on the numbers of nodes per sequence?

``` r
pqbr_details <- read.table("./unmapped_assemblies/ldu_pQBR_details.txt",
                           col.names=c("plasmid","node_number","length","coverage"))

pqbr_details %>% group_by(plasmid) %>%
  summarise(number_of_nodes = n(),
            max_length = max(length),
            min_length = min(length),
            mean_coverage = mean(coverage)) %>%
  kable()
```

| plasmid     | number_of_nodes | max_length | min_length | mean_coverage |
|:------------|----------------:|-----------:|-----------:|--------------:|
| UWC1pQBR1   |               1 |       9878 |       9878 |      33.32211 |
| UWC1pQBR102 |               1 |     312323 |     312323 |      35.43129 |
| UWC1pQBR103 |               1 |     425168 |     425168 |      19.64236 |
| UWC1pQBR105 |               6 |      63130 |       7667 |      81.77337 |
| UWC1pQBR106 |               5 |      53387 |       5674 |      36.80394 |
| UWC1pQBR11  |               2 |     320086 |     110166 |      21.81896 |
| UWC1pQBR124 |               2 |     324264 |      99013 |      83.47026 |
| UWC1pQBR127 |               7 |     161015 |       5033 |      77.38717 |
| UWC1pQBR132 |               1 |     138541 |     138541 |      41.35900 |
| UWC1pQBR150 |               1 |     307116 |     307116 |      36.27496 |
| UWC1pQBR23  |               5 |     167665 |       7593 |      48.49636 |
| UWC1pQBR24  |               6 |     145724 |       7576 |      20.51830 |
| UWC1pQBR26  |              16 |     140175 |       5164 |      52.45696 |
| UWC1pQBR28  |               1 |     140281 |     140281 |      36.99634 |
| UWC1pQBR30  |               1 |     310157 |     310157 |      57.01980 |
| UWC1pQBR4   |               3 |     214260 |      38950 |      36.53968 |
| UWC1pQBR43  |               2 |     256230 |     136726 |      90.86858 |
| UWC1pQBR47  |               1 |     252601 |     252601 |      81.73491 |
| UWC1pQBR49  |               4 |     175972 |      42153 |      55.64029 |
| UWC1pQBR5   |               1 |     424612 |     424612 |      53.25670 |
| UWC1pQBR50  |               3 |     374297 |      11529 |      58.18543 |
| UWC1pQBR51  |               2 |     164905 |     142276 |      30.52411 |
| UWC1pQBR53  |               2 |      73086 |      66141 |      32.78773 |
| UWC1pQBR55  |               1 |     139190 |     139190 |      38.90255 |
| UWC1pQBR56  |               1 |     307407 |     307407 |      46.41393 |
| UWC1pQBR57  |               1 |     307276 |     307276 |      63.91373 |
| UWC1pQBR58  |               1 |      73703 |      73703 |      13.86419 |
| UWC1pQBR8   |               1 |       8488 |       8488 |      20.53846 |

Calculate mash distances. Use a high value for sketch size to detect
distant similarities between sequences.

- S = default seed function
- s = default sketch size
- k = default kmer size

``` bash
mash sketch ./unmapped_assemblies/ldu_pQBR.fasta -i -S 42 -s 100000 -k 21 -p 4 -o ./unmapped_assemblies/ldu_pQBR.msh
mash triangle ./unmapped_assemblies/ldu_pQBR.msh -i -k 21 -p 4 > ./unmapped_assemblies/ldu_pQBR_mash.dst
```

Plot as a heatmap in R. Make distance matrix square

``` r
pqbr_dist <- read.table("./unmapped_assemblies/ldu_pQBR_mash.dst", fill=TRUE, skip=1,
                        col.names=c("V0",paste("V", 1:28, sep="")))
pqbr_dist$V0 <- gsub("UWC1(pQBR[0-9]+)_unmapped_assembly", "\\1", pqbr_dist$V0)
pqbr_dist <- column_to_rownames(pqbr_dist, "V0")
pqbr_distmat <- as.dist(pqbr_dist, upper=TRUE, diag=TRUE)
pqbr_sqmat <- as.matrix(pqbr_distmat)

plasmid_order <- colnames(pqbr_sqmat)

pqbr_dist_df <- pqbr_sqmat %>% as.data.frame() %>%
  mutate(a = plasmid_order) %>%
  pivot_longer(-a, names_to = "b", values_to = "mash_distance") %>%
  filter(!is.na(mash_distance)) %>%
  mutate(a = factor(a, levels=plasmid_order), b=factor(b, levels=plasmid_order))

ggplot(data=pqbr_dist_df) +
  geom_tile(aes(x=a, y=b, fill=mash_distance)) +
  scale_fill_gradient(low = "dodgerblue", high = "black")
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Cluster according to similarities.

``` r
pqbr_dendro <- as.dendrogram(hclust(pqbr_distmat))

dendro_plot <- ggdendrogram(data = pqbr_dendro, rotate = TRUE)

plasmid_reorder <- order.dendrogram(pqbr_dendro)

(short_read_heatmap <- pqbr_dist_df %>%
  mutate(a = factor(a, levels=plasmid_order[plasmid_reorder]), 
         b = factor(b, levels=plasmid_order[plasmid_reorder])) %>%
  ggplot() +
  geom_tile(aes(x=a, y=b, fill=mash_distance)) +
  scale_fill_gradientn(colours = c("skyblue","dodgerblue","black"), values=c(0,0.25,1), name = "mash distance") +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.position="right"))
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
png(file = "./figs/short_read_heatmap.png", height=6, width=8, units="in", res=300)
short_read_heatmap
dev.off()
```

    ## quartz_off_screen 
    ##                 2

This gives a good overview of the plasmids.

#### Interpretation of preliminary analysis

There are five main clusters of sequences that show high internal
similarity: pQBR103-like (Group I), pQBR23-like (Group II), pQBR55-like
(Group III), pQBR57-like (Group IV), pQBR1-like (Group V), and pQBR105
(an outlier).

However, there are also some other observations, supported by manual
examination of the long, deep, unmapped assemblies, and the plots above.

**pQBR1**: Illumina sequencing indicated no plasmid present. Instead,
and in contrast to the ~70-kb plasmid referred to by [Lilley et
al. (1996)](https://doi.org/10.1111/j.1574-6941.1996.tb00320.x), an
integrated analysis of the assembly of all and unmapped reads shows a
9.9 kb sequence that has integrated into the chromosome. This suggests
that the stock strain has lost the plasmid but kept the transposon, as
described in various lab studies with other pQBR plasmids. We were not
able to detect conjugation in the lab, further supporting our inference
that the plasmid has been lost.

**pQBR8**: Illumina sequencing indicated low coverage for regions of
interest, except for a chromsomal-matching contig encoding an integrated
*mer* operon, suggesting that the plasmid was in the process of being
lost when the sample was sequenced. We were not able to detect
conjugation in the lab, further supporting our inference that the
plasmid is unstable/lost.

**pQBR58**: Illumina sequencing indicated a ~75 kb fragment, but this
was at lower coverage (~0.3x) compared with the chromosome. The *mer*
containing region (~8.5 kb) was connected to the chromosome and at ~1x,
suggesting that the plasmid was in the process of being lost when the
sample was sequenced. We were not able to detect conjugation in the lab,
further supporting our inference that the plasmid is unstable/lost.

**pQBR127**: MASH distance clustering indicated that this sample matched
both Group I and Group III sequences, which may indicate that more than
one plasmid is present in the sample.

### Resolving sequences

Sequences that were not fully resolved by either the assembly of all
reads, or the assembly of unmapped reads, were subject to ONT sequencing
by Plasmidsaurus.

The long-read contigs that matched each plasmid were identified by
BLASTing the ONT assemblies for each sample against the contigs
assembled by Illumina. The long-read contigs were extracted, named, and
placed in `./contigs` for polishing.

Some manual analysis was required for the following, possibly due to
sequence loss, chromosome integration, multiple plasmids, or some
combination of the above. This requires further analysis.

#### pQBR26 (Group II)

pQBR26 was previously identified as a Group II plasmid, like pQBR23 and
pQBR24. However, though transconjugants into *P. fluorescens* SBW25 were
easily generated, they did not give a positive PCR product for the
primers designed for Group II plasmids. Analysis of the Illumina reads
showed another contig that also contained a predicted *merA* gene, and
primers designed against this contig did give a positive product when
applied to transconjugant colonies. Long-read sequencing of UWC1(pQBR26)
showed two extrachromosomal replicons, one of which resembled pQBR23
(394 kb) but had a copy number of \<1, and the other with a copy number
of ~1 that matched the contig that was capable of conjugating (229 kb).
We infer that UWC1(pQBR26) contained two different *mer* plasmids, one
of which is less stable and/or less conjugative than the other, and the
229-kb conjugative plasmid is indicated as the canonical pQBR26 plasmid.

#### pQBR11 (Group I)

pQBR11 produced Group I-like contigs when the short reads were analysed.
However, the sample did not resolve in the long read sequencing. The
sequence therefore remains as a draft (pQBR11d) based on the Illumina
unmapped read assembly.

#### pQBR43 (Group I)

This sequence didn’t fully resolve owing to a duplicate sequence in the
chromosome. Two possibilities: both the chromosome and the plasmid
contain the same duplicated element (likely a Tn6290 transposon) but
remain separate replicons, or the plasmid has integrated into the
chromosome, flanked by copies of the duplication.

Resolving this is complex, as the duplicated element is 42 kb long and
there are insufficient reads overlapping the junction.

A path containing the plasmid contigs and Tn6290 transposon was
extracted for analysis.

#### pQBR50 (Group I)

This sequence didn’t fully resolve owing to duplicate sequences in the
chromosome and multiple copies of Tn6290. The mean depth of the
chromosome (from the Illumina sequencing) was ~21x, the mean depth of
Tn6290 was ~130x! Suggesting there are between 6 and 7 copies of Tn6290
(5-6 of which likely on or associated with the plasmid). Interestingly,
this plasmid seems to be at approximately 2x copy number of the
chromosome in the Illumina sequencing, but a bit lower in the ONT
sequencing.

As with pQBR43, this couldn’t be resolved due to the length of Tn6290,
even with the long-read sequencing.

The sequence therefore remains a draft (pQBR50d) based on conjoining the
largest contigs from the ONT sequencing.

``` bash
awk -v RS=">" '$1 ~ /^contig_3/ {print RS $0}' ./working/pQBR50_consensus_cov.fna \
  | revseq -filter > ./working/pQBR50_consensus_contigs.fasta
  
awk -v RS=">" '$1 ~ /^contig_6/ {print RS $0}' ./working/pQBR50_consensus_cov.fna \
  >> ./working/pQBR50_consensus_contigs.fasta

union ./working/pQBR50_consensus_contigs.fasta -filter \
  | sed 's/>.*/>pQBR50d/g' > ./working/pQBR50d_contig.fasta
```

These were polished as described below.

#### pQBR127 (“Group IV”)

Analysis of short-read sequencing suggested that this strain harboured
replicons that both resembled Group I and Group III plasmids. The
bandage plot shows matches to pQBR55 (Group III, blue) and pQBR103
(Group I, green).

<figure>
<img src="../figs/pQBR127_bandage.png" alt="pQBR127_bandage.png" />
<figcaption aria-hidden="true">pQBR127_bandage.png</figcaption>
</figure>

Performed long-read sequencing to resolve, which showed a distinct Group
III replicon and no Group I replicon. Infer that two plasmids were
present in the Illumina sequencing, of which only one remains in our
sample.

#### pQBR47 (Group I)

This sequence didn’t fully resolve owing to a duplicate Tn6290 sequence,
and the sequence is assembled as a single contig including plasmid and
chromosome. Perhaps the plasmid has also integrated into the chromosome,
but as above, resolving this with current sequencing data is not
straightforward.

The whole contig was used for downstream analysis, and the
plasmid-specific portion was extracted (see below).

#### pQBR106 (ND)

Similarly, this sequence didn’t fully resolve, in a similar manner to
pQBR47, with a single contig including the chromosome.

The whole contig was used for downstream analysis, and the
plasmid-specific portion was extracted (see below).

#### pQBR44 (Group I)

pQBR44 was previously sequenced in [Hall et
al. 2015](https://pubmed.ncbi.nlm.nih.gov/25969927/). It was not
resequenced here. The two sequences identified could not be resolved,
due to transposon insertions.

It was downloaded and the two contigs were concatenated (pQBR44d).

``` bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=CDLQ010000002.1&amp;rettype=fasta" \
  > ./working/pQBR44_1.fasta
  
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=CDLQ010000001.1&amp;rettype=fasta" \
  >> ./working/pQBR44_1.fasta
  
union -sequence ./working/pQBR44_1.fasta -filter \
  | sed 's/>.*/>pQBR44d/g' > ./polished/pQBR44d.fasta
```

### Polishing

Polished sequences using
[PolyPolish](https://github.com/rrwick/Polypolish).

Align reads to draft genomes with bwa and run PolyPolish.

``` bash
find ./contigs -name "pQBR*" -type f |
  sed 's/.*pQBR/pQBR/g' | sed 's/_contig.fasta//g' \
  > long_reads_draft.txt

find ./illumina/*_1_trimmed.fastq.gz > 1.tmp
find ./illumina/*_2_trimmed.fastq.gz > 2.tmp
awk -v FS="_" '{print $2}' 1.tmp > n.tmp
paste n.tmp 1.tmp 2.tmp > ctrlsamp.txt
rm *.tmp
  
cat long_reads_draft.txt | while read SAMPLE
do
grep "${SAMPLE}_" ./ctrlsamp.txt | while read STR R1 R2
do 
echo "---"
echo "Next plasmid is $SAMPLE"
echo "---"
bwa index ./contigs/${SAMPLE}_contig.fasta
bwa mem -t 16 -a ./contigs/${SAMPLE}_contig.fasta $R1 > alignments_1.sam
bwa mem -t 16 -a ./contigs/${SAMPLE}_contig.fasta $R2 > alignments_2.sam
/pub60/jamesh/bin/polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam \
  --out1 filtered_1.sam --out2 filtered_2.sam
polypolish ./contigs/${SAMPLE}_contig.fasta filtered_1.sam filtered_2.sam > ./polished/${SAMPLE}_polished.fasta
done
done
```

### Orienting/aligning

Checked with blastn and ACT to look for repeated sequences at the ends
of the ‘polished’ assemblies. Direct repeats, particularly from the
Illumina assemblies, likely indicates circularisation. Use
[ccfind](https://github.com/yosuken/ccfind) to identify and remove such
repeats.

``` bash
mkdir ./ccfind

find ./polished -name "*.fasta" -exec basename -- {} .fasta \; \
  | while read FASTA 
  do 
  ccfind ./polished/${FASTA}.fasta ./ccfind/${FASTA}
  done

mkdir ./archive

find ./ccfind/*/result -name "circ.detected.list" -exec wc -l {} \; \
  | awk '$1 == 1 {print $0}' | awk -v FS="/" '{print $3}' \
  | while read PLASMID
  do
  echo $PLASMID
  mv ./polished/${PLASMID}.fasta ./archive/${PLASMID}_withTR.fasta
  cp ./ccfind/${PLASMID}/result/circ.noTR.fasta ./polished/${PLASMID}.fasta
  done
```

Next, orient the sequences so they start at a similar/comparable
position. Use the first coding base of the putative replicase.

- For pQBR57-like (Group IV), use PQBR57_0001 in the reference sequence
  LN713926.
- For pQBR103-like (Group I) this is at 388361 in the reference sequence
  AM235768 (PQBR0445).
- For pQBR55-like (Group III), use the first of the dnaB genes referred
  to in Turner et al. (2002), which is PQBR55_0049.

For the others, do an analysis to check they are oriented correctly
relative to one another.

Create a blastn database of the three above CDS and blast all against it
to find the direction and first base of the corresponding gene in each
assembly.

``` bash
cat ./ref/PQBR55_0049.fasta ./ref/PQBR57_0001.fasta ./ref/pQBR0445.fasta > ./ref/rep.fasta

makeblastdb -in ./ref/rep.fasta -dbtype nucl

find ./polished -name "*.fasta" -exec blastn -query {} -db ./ref/rep.fasta -outfmt 6 \; \
  > ./working/blast_rep.blastn
```

Use this output to realign the sequences using EMBOSS.

``` bash
mkdir ./complete

cat ./working/blast_rep.blastn | while read SEQ MATCH PERC LEN GAP MM QSTART QFIN SSTART SFIN EVAL BITSC
do
let endposf=$QSTART-1
let endposr=$QFIN+1
filename=`echo $SEQ | sed 's/_polypolish/_polished/g'`
seqname=`echo $filename | sed 's/_polished//g'`
if [ $SFIN -lt $SSTART ]
then
echo "$SEQ is in reverse direction"
seqret -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send ${QFIN} \
  -srev \
  -sequence ./polished/${filename}.fasta \
  -filter > tmp1.fasta
seqret -sformat1 fasta \
  -osformat2 fasta \
  -srev \
  -sbegin1 $endposr \
  -sequence ./polished/${filename}.fasta \
  -filter > tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 ${seqname} \
  | sed 's/_polypolish//g' \
  > ./complete/${seqname}.fasta  
else
echo "$SEQ is in forward direction"
seqret -sequence ./polished/${filename}.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin ${QSTART} \
  -outseq tmp1.fasta
seqret -sequence ./polished/${filename}.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 $endposf \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 ${seqname} \
  | sed 's/_polypolish//g' \
  > ./complete/${seqname}.fasta
fi
done
```

This gives complete, oriented assemblies for most sequences. Some are
missing due to a lack of blastn matches: the Group II plasmids pQBR23
and pQBR24, pQBR26, pQBR105 (no similar plasmids amongst the groups).

For pQBR23 and pQBR24, these should be oriented to the same point. Run a
quick Prokka annotation and identify a putative replicase.

``` bash
prokka ./polished/pQBR23_polished.fasta -out ./working/pQBR23_annot
```

Finds a putative pair of plasmid replication genes, the first of which
is at position 147660..149189. Extract and blast against pQBR24.

``` bash
awk -v seq="PROKKA_00168" -v RS='>' '$1 == seq {print RS $0}' ./working/pQBR23_annot/PROKKA_12212023.ffn \
  > ./working/pQBR23_repA.fasta
  
makeblastdb -in ./working/pQBR23_repA.fasta -dbtype nucl

blastn -query ./polished/pQBR24_polished.fasta -db ./working/pQBR23_repA.fasta -outfmt 6 \
  >  ./working/pQBR24_rep.blastn
```

Identifies the sequence in pQBR24 at 384040..385569, in reverse. So
385569 is the start of the gene, and 384040 is the end.

``` bash
seqret -sequence ./polished/pQBR23_polished.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 147660 \
  -outseq tmp1.fasta
seqret -sequence ./polished/pQBR23_polished.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 147659 \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 pQBR23 \
  | sed 's/_polypolish//g' \
  > ./complete/pQBR23.fasta
  
seqret -sequence ./polished/pQBR24_polished.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 385569 \
  -srev \
  -outseq tmp1.fasta
seqret -sequence ./polished/pQBR24_polished.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 385570 \
  -srev \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 pQBR24 \
  | sed 's/_polypolish//g' \
  > ./complete/pQBR24.fasta
```

Do similar for pQBR26.

``` bash
prokka ./polished/pQBR26.fasta -out ./working/pQBR26_annot
```

Finds a repE-like replicase gene starting as position 112373.

``` bash
seqret -sequence ./polished/pQBR26.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 112373 \
  -outseq tmp1.fasta
seqret -sequence ./polished/pQBR26.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 112373 \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 pQBR26 \
  | sed 's/_polypolish//g' \
  > ./complete/pQBR26.fasta
```

…and pQBR105.

``` bash
prokka ./polished/pQBR105_polished.fasta -out ./working/pQBR105_annot
```

This couldn’t find a putative replicase, so the sequence was just copied
across.

``` bash
sed 's/_polypolish//g' ./polished/pQBR105_polished.fasta > ./complete/pQBR105.fasta
```

Add ‘R’ to pQBR103, pQBR55, and pQBR57 to indicate these variants have
been resequenced (this was done manually).

Download and re-centre the originally-sequenced pQBR103, pQBR55, and
pQBR57 to analyse alongside.

pQBR103:

``` bash
mkdir ./originals

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=AM235768&amp;rettype=fasta" \
  > ./originals/pQBR103.fasta

blastn -query ./originals/pQBR103.fasta -db ./ref/rep.fasta -outfmt 6 \
  > pQBR103_rep.blastn
  
# check and align
# 388361    389491

seqret -sequence ./originals/pQBR103.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 388361 \
  -outseq tmp1.fasta
seqret -sequence ./originals/pQBR103.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 388360 \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 pQBR103 \
  > ./originals/pQBR103.fasta
  
gzip ./originals/pQBR103.fasta
```

pQBR57:

``` bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=LN713926&amp;rettype=fasta" \
  > ./originals/pQBR57.fasta

blastn -query ./originals/pQBR57.fasta -db ./ref/rep.fasta -outfmt 6 \
  > pQBR57_rep.blastn
  
# 542

seqret -sequence ./originals/pQBR57.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 542 \
  -outseq tmp1.fasta
seqret -sequence ./originals/pQBR57.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 541 \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 pQBR57 \
  > ./originals/pQBR57.fasta
  
gzip ./originals/pQBR57.fasta
```

pQBR55:

``` bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=LN713927&amp;rettype=fasta" \
  > ./originals/pQBR55.fasta

blastn -query ./originals/pQBR55.fasta -db ./ref/rep.fasta -outfmt 6 \
  > pQBR55_rep.blastn
  
# 34006

seqret -sequence ./originals/pQBR55.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 34006 \
  -outseq tmp1.fasta
seqret -sequence ./originals/pQBR55.fasta \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send1 34005 \
  -outseq tmp2.fasta
cat tmp1.fasta tmp2.fasta \
  | union -filter -osname2 pQBR55 \
  > ./originals/pQBR55.fasta
  
gzip ./originals/pQBR55.fasta
```

Annotate with Bakta version 1.8.2. Database is version 5.0 full.

Use metagenome mode to ensure that the same CDS prediction model is used
across all of the plasmids, to aid identification of homologues on
related sequences.

``` bash
find ./complete -name "pQBR*.fasta" -exec gzip {} \;

cp ./originals/*.fasta.gz ./complete

find ./complete -name "pQBR*.fasta.gz" -exec basename -- {} .fasta.gz \; \
  > ./list_of_pQBR_plasmids.txt
  
cat list_of_pQBR_plasmids.txt |
while read PLASMID
do
bakta --db /pub60/jamesh/db --prefix ${PLASMID} \
  --complete \
  --locus ${PLASMID}_contig \
  --verbose \
  --output ./bakta_a/${PLASMID} \
  --plasmid ${PLASMID} \
  --threads 128 \
  --locus-tag ${PLASMID^^} \
  --meta \
  ./complete/${PLASMID}.fasta.gz
done
```

These annotations (bakta_annotated) used for subsequent analyses.

Make a table of information for each sequence:

``` bash
grep "Length: " ./bakta_annotated/*/*.txt \
  | sed 's:./bakta_annotated/::g' | sed 's;/pQBR[0-9Rpd]*.txt:Length: ;\t;g' | sort -n -k2
```

    ## pQBR106p 98802
    ## pQBR132  139938
    ## pQBR127  140415
    ## pQBR55R  140432
    ## pQBR28   141505
    ## pQBR44d  143197
    ## pQBR53   157450
    ## pQBR55   157450
    ## pQBR105  161330
    ## pQBR26   228742
    ## pQBR47p  293223
    ## pQBR56   307330
    ## pQBR57   307330
    ## pQBR30   310323
    ## pQBR57R  324348
    ## pQBR102  334884
    ## pQBR150  366385
    ## pQBR51   366385
    ## pQBR23   393597
    ## pQBR24   393604
    ## pQBR124  424568
    ## pQBR103  425094
    ## pQBR103R 425094
    ## pQBR11d  430820
    ## pQBR5    466604
    ## pQBR43   510643
    ## pQBR4    525650
    ## pQBR50d  563271
    ## pQBR49   587656
    ## pQBR106  6320937
    ## pQBR47   6515328

Shows that pQBR106 and pQBR47 are abnormally large. The specific
plasmid-encoding portions were identified using ACT with comparison to
the KT2440 chromosome.

Note that the concatenation between pQBR44 contigs is at position
92628-92629 in pQBR44d.

#### pQBR47 plasmid-specific portion

The plasmid-specific portion, with a single integrated copy of Tn6290,
goes from position `6501702..6515328` and `1..279596`, as identified by
ACT.

``` bash
seqret -sequence ./bakta_annotated/pQBR47/pQBR47.fna \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send 279596 \
  -outseq ./working/tmp1.fasta
seqret -sequence ./bakta_annotated/pQBR47/pQBR47.fna \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 6501702 \
  -send1 6515328 \
  -outseq ./working/tmp2.fasta

cat ./working/tmp1.fasta ./working/tmp2.fasta \
 > ./working/pQBR47p.fasta
```

Sequence was manually renamed and run through bakta.

``` bash
PLASMID=pQBR47p

bakta --db /pub60/jamesh/db --prefix ${PLASMID} \
  --locus ${PLASMID}_contig \
  --verbose \
  --output ./bakta_a/${PLASMID} \
  --plasmid ${PLASMID} \
  --threads 128 \
  --locus-tag ${PLASMID^^} \
  --meta \
  pQBR47p.fasta
```

#### pQBR106 plasmid-specific portion

The plasmid-specific portion goes from position `1..14487` (running into
a copy of Tn6290) and `6236623` (running from a copy of Tn6290) to the
end of the contig at `6320937`. The remaining plasmid-specific part of
pQBR106 is much smaller than it is for pQBR47.

``` bash
seqret -sequence ./bakta_annotated/pQBR106/pQBR106.fna \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 1 \
  -send 14487 \
  -outseq ./working/tmp1.fasta
seqret -sequence ./bakta_annotated/pQBR106/pQBR106.fna \
  -sformat1 fasta \
  -osformat2 fasta \
  -sbegin 6236623 \
  -send1 6320937 \
  -outseq ./working/tmp2.fasta

cat ./working/tmp1.fasta ./working/tmp2.fasta \
 > ./working/pQBR106p.fasta
```

Sequence was manually renamed and run through bakta.

``` bash
PLASMID=pQBR106p

bakta --db /pub60/jamesh/db --prefix ${PLASMID} \
  --locus ${PLASMID}_contig \
  --verbose \
  --output ./bakta_a/${PLASMID} \
  --plasmid ${PLASMID} \
  --threads 128 \
  --locus-tag ${PLASMID^^} \
  --meta \
  pQBR106p.fasta
```

Annotated sequences for downstream analysis are in the subdirectory
`./bakta_annotated`.

Note that the suffix ‘R’ indicates ‘resequenced’, ‘d’ indicates ‘draft’
(i.e. manual editing/concatenation of contigs necessary to generate the
final sequence, as described above), and ‘p’ indicates the
plasmid-specific portion of a (potentially) chromosomally-integrated
assembly.

### Sequence distance/similarity plot for assembled plasmids

Generate a heatmap like that above, but just for the assembled/annotated
plasmids used in downstream analyses.

``` bash
find ./bakta_annotated/*/*.fna | awk -v FS="/" '{print $3}' > ./1_sketches/annotated_plasmids.list

cat ./1_sketches/annotated_plasmids.list | while read PLASMID
do
mash sketch ./bakta_annotated/${PLASMID}/${PLASMID}.fna -S 42 -s 100000 -k 21 -p 4 -o ./1_sketches/${PLASMID}.msh
done

mash triangle ./1_sketches/pQBR*.msh -i -k 21 -p 4 > ./1_sketches/bakta_annotated_mash.dst
```

Plot as a heatmap in R. Make distance matrix square

``` r
pqbra_dist <- read.table("./1_sketches/bakta_annotated_mash.dst", fill=TRUE, skip=1,
                        col.names=c("V0",paste("V", 1:31, sep=""))) 
pqbra_dist$V0 <- gsub(".*/(pQBR[0-9Rpd]+)/.*.fna", "\\1", pqbra_dist$V0)
pqbra_dist <- column_to_rownames(pqbra_dist, "V0")
pqbra_distmat <- as.dist(pqbra_dist, upper=TRUE, diag=TRUE)
pqbra_sqmat <- as.matrix(pqbra_distmat)

plasmida_order <- colnames(pqbra_sqmat)

pqbra_dist_df <- pqbra_sqmat %>% as.data.frame() %>%
  mutate(a = plasmida_order) %>%
  pivot_longer(-a, names_to = "b", values_to = "mash_distance") %>%
  filter(!is.na(mash_distance)) %>%
  mutate(a = factor(a, levels=plasmida_order), b=factor(b, levels=plasmida_order))

ggplot(data=pqbra_dist_df) +
  geom_tile(aes(x=a, y=b, fill=mash_distance)) +
  scale_fill_gradient(low = "dodgerblue", high = "black")
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

Cluster according to similarities.

``` r
pqbra_dendro <- as.dendrogram(hclust(pqbra_distmat))

plasmida_reorder <- order.dendrogram(pqbra_dendro)

(annotated_heatmap <- pqbra_dist_df %>%
  mutate(a = factor(a, levels=plasmida_order[plasmida_reorder]), 
         b = factor(b, levels=plasmida_order[plasmida_reorder])) %>%
  ggplot() +
  geom_tile(aes(x=a, y=b, fill=mash_distance)) +
  scale_fill_gradientn(colours = c("skyblue","dodgerblue","black"), values=c(0,0.25,1), name = "mash distance") +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(), axis.title.y = element_blank(), 
        legend.position="right"))
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

Remove pQBR47 and pQBR106 from the heatmap, as these erroneously include
the whole chromosome as well.

``` r
(annotated_heatmap_edit <- pqbra_dist_df %>%
   filter(!(a %in% c("pQBR47", "pQBR106")) & !(b %in% c("pQBR47", "pQBR106"))) %>%
   mutate(a = factor(a, levels=plasmida_order[plasmida_reorder]), 
          b = factor(b, levels=plasmida_order[plasmida_reorder])) %>%
   ggplot() +
   geom_tile(aes(x=a, y=b, fill=mash_distance)) +
   scale_fill_gradientn(colours = c("skyblue","dodgerblue","black"), values=c(0,0.25,1), name = "mash distance") +
   theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(), axis.title.y = element_blank(), 
         legend.position="right"))
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
png(file = "./figs/annotated_heatmap_edit.png", height=8, width=11, units="in", res=300)
annotated_heatmap_edit
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Comparison of all-against-all using average nucleotide identity (ANI) rather than mash

Average nucleotide identity (ANI) is a more intuitive measure than mash
distance.

Perform all-vs-all comparison using fastANI. Remove the whole-genome
plasmids from the comparison.

``` bash
find ./bakta_annotated/*/ -name "*.fna" \
  | grep -v "pQBR106.fna" | grep -v "pQBR47.fna" > ./working/genome_list.txt

fastani --ql ./working/genome_list.txt \
        --rl ./working/genome_list.txt \
        --output ./1_sketches/ani.tsv --fragLen 500 \
        --minFraction 0.2 --threads 1 \
        --matrix
```

## Summary of plasmid features

``` bash
find ./bakta_annotated/ -name "pQBR*.txt" > ./3_summaries/summary_filenames.txt

cat ./3_summaries/summary_filenames.txt | while read FILE
do
PLASMID=`echo $FILE | sed -E 's/.*(pQBR[0-9Rdp]+).txt/\1/'`
LENGTH=`grep "Length" $FILE | sed 's/Length: //g'`
GC=`grep "GC:" $FILE | sed 's/GC: //g'`
DENSITY=`grep "coding density:" $FILE | sed 's/coding density: //g'`
CDS=`grep "CDSs:" $FILE | sed 's/CDSs: //g'`
echo -e "$PLASMID\t$LENGTH\t$GC\t$DENSITY\t$CDS"
done > ./3_summaries/summaries.txt
```

Input to R and format.

``` r
summaries <- read.table("./3_summaries/summaries.txt", sep="\t", 
                        col.names=c("plasmid","length","GC","coding_density","n_cds")) %>%
  mutate(last_char = substr(plasmid, nchar(plasmid), nchar(plasmid)),
         assembly_details = case_when(last_char == "R" ~ "resequenced",
                                      last_char == "d" ~ "draft",
                                      last_char == "p" ~ "nonchromosomal_only",
                                      .default = "assembled"),
         Plasmid = gsub("([0-9])[a-zA-Z]", "\\1", plasmid)) %>%
  select(-last_char)

lilley <- read.table("./ref/Lilley.csv", sep=",", header=TRUE)

(full_summaries <- left_join(summaries, lilley, by="Plasmid") %>% 
  arrange(Group, Plasmid) %>% 
    select(Group, Plasmid, assembly_details, length, GC, coding_density, n_cds)) %>% kable()
```

| Group | Plasmid | assembly_details    |  length |   GC | coding_density | n_cds |
|:------|:--------|:--------------------|--------:|-----:|---------------:|------:|
| I     | pQBR103 | resequenced         |  425094 | 53.2 |           86.0 |   550 |
| I     | pQBR103 | assembled           |  425094 | 53.2 |           86.0 |   550 |
| I     | pQBR11  | draft               |  430820 | 53.2 |           86.3 |   558 |
| I     | pQBR124 | assembled           |  424568 | 53.2 |           85.9 |   547 |
| I     | pQBR4   | assembled           |  525650 | 53.7 |           85.8 |   658 |
| I     | pQBR43  | assembled           |  510643 | 53.8 |           86.1 |   614 |
| I     | pQBR44  | draft               |  143197 | 53.7 |           85.3 |   191 |
| I     | pQBR47  | assembled           | 6515328 | 61.1 |           89.6 |  5985 |
| I     | pQBR47  | nonchromosomal_only |  293223 | 53.3 |           86.7 |   338 |
| I     | pQBR49  | assembled           |  587656 | 53.7 |           85.5 |   712 |
| I     | pQBR5   | assembled           |  466604 | 53.4 |           85.8 |   590 |
| I     | pQBR50  | draft               |  563271 | 54.0 |           85.9 |   684 |
| I     | pQBR51  | assembled           |  366385 | 54.2 |           84.9 |   462 |
| II    | pQBR23  | assembled           |  393597 | 57.5 |           81.5 |   499 |
| II    | pQBR24  | assembled           |  393604 | 57.5 |           81.6 |   499 |
| II    | pQBR26  | assembled           |  228742 | 57.1 |           87.3 |   273 |
| III   | pQBR28  | assembled           |  141505 | 51.5 |           84.9 |   176 |
| III   | pQBR53  | assembled           |  157450 | 52.2 |           85.3 |   199 |
| III   | pQBR55  | assembled           |  157450 | 52.2 |           85.3 |   198 |
| III   | pQBR55  | resequenced         |  140432 | 51.6 |           84.9 |   174 |
| IV    | pQBR102 | assembled           |  334884 | 53.6 |           85.1 |   421 |
| IV    | pQBR127 | assembled           |  140415 | 51.6 |           85.3 |   177 |
| IV    | pQBR30  | assembled           |  310323 | 53.7 |           84.6 |   390 |
| IV    | pQBR56  | assembled           |  307330 | 53.8 |           84.6 |   395 |
| IV    | pQBR57  | resequenced         |  324348 | 53.9 |           84.8 |   419 |
| IV    | pQBR57  | assembled           |  307330 | 53.8 |           84.6 |   395 |
| ND    | pQBR105 | assembled           |  161330 | 56.4 |           88.3 |   172 |
| ND    | pQBR106 | nonchromosomal_only |   98802 | 53.7 |           82.4 |   144 |
| ND    | pQBR106 | assembled           | 6320937 | 61.4 |           89.7 |  5803 |
| ND    | pQBR132 | assembled           |  139938 | 51.6 |           85.8 |   176 |
| ND    | pQBR150 | assembled           |  366385 | 54.2 |           84.9 |   462 |

``` r
problem_assemblies <- c("pQBR47","pQBR106")

(cropped_summaries <- full_summaries %>%
  filter(!(Plasmid %in% problem_assemblies) & assembly_details == "assembled" | assembly_details == "draft") %>%
  mutate(revised_group = case_when(Plasmid == "pQBR106" ~ "I",
                                   Plasmid == "pQBR127" ~ "III",
                                   Plasmid == "pQBR132" ~ "III",
                                   Plasmid == "pQBR150" ~ "IV",
                                   Plasmid == "pQBR51" ~ "IV",
                                   .default = Group)) %>%
    select(Group, revised_group, Plasmid, assembly_details, length, GC, coding_density, n_cds) %>%
    arrange(revised_group, Plasmid)) %>% kable()
```

| Group | revised_group | Plasmid | assembly_details | length |   GC | coding_density | n_cds |
|:------|:--------------|:--------|:-----------------|-------:|-----:|---------------:|------:|
| I     | I             | pQBR103 | assembled        | 425094 | 53.2 |           86.0 |   550 |
| I     | I             | pQBR11  | draft            | 430820 | 53.2 |           86.3 |   558 |
| I     | I             | pQBR124 | assembled        | 424568 | 53.2 |           85.9 |   547 |
| I     | I             | pQBR4   | assembled        | 525650 | 53.7 |           85.8 |   658 |
| I     | I             | pQBR43  | assembled        | 510643 | 53.8 |           86.1 |   614 |
| I     | I             | pQBR44  | draft            | 143197 | 53.7 |           85.3 |   191 |
| I     | I             | pQBR49  | assembled        | 587656 | 53.7 |           85.5 |   712 |
| I     | I             | pQBR5   | assembled        | 466604 | 53.4 |           85.8 |   590 |
| I     | I             | pQBR50  | draft            | 563271 | 54.0 |           85.9 |   684 |
| II    | II            | pQBR23  | assembled        | 393597 | 57.5 |           81.5 |   499 |
| II    | II            | pQBR24  | assembled        | 393604 | 57.5 |           81.6 |   499 |
| II    | II            | pQBR26  | assembled        | 228742 | 57.1 |           87.3 |   273 |
| IV    | III           | pQBR127 | assembled        | 140415 | 51.6 |           85.3 |   177 |
| ND    | III           | pQBR132 | assembled        | 139938 | 51.6 |           85.8 |   176 |
| III   | III           | pQBR28  | assembled        | 141505 | 51.5 |           84.9 |   176 |
| III   | III           | pQBR53  | assembled        | 157450 | 52.2 |           85.3 |   199 |
| III   | III           | pQBR55  | assembled        | 157450 | 52.2 |           85.3 |   198 |
| IV    | IV            | pQBR102 | assembled        | 334884 | 53.6 |           85.1 |   421 |
| ND    | IV            | pQBR150 | assembled        | 366385 | 54.2 |           84.9 |   462 |
| IV    | IV            | pQBR30  | assembled        | 310323 | 53.7 |           84.6 |   390 |
| I     | IV            | pQBR51  | assembled        | 366385 | 54.2 |           84.9 |   462 |
| IV    | IV            | pQBR56  | assembled        | 307330 | 53.8 |           84.6 |   395 |
| IV    | IV            | pQBR57  | assembled        | 307330 | 53.8 |           84.6 |   395 |
| ND    | ND            | pQBR105 | assembled        | 161330 | 56.4 |           88.3 |   172 |

``` r
cropped_summaries %>% filter(assembly_details == "assembled" & (length == min(length) | length == max(length)))
```

    ##   Group revised_group Plasmid assembly_details length   GC coding_density n_cds
    ## 1     I             I  pQBR49        assembled 587656 53.7           85.5   712
    ## 2    ND           III pQBR132        assembled 139938 51.6           85.8   176

``` r
(cropped_summaries %>% group_by(revised_group) %>% 
  summarise(n = n(),
            mean_length_kb = mean(length/1000),
            max_length = max(length),
            min_length = min(length),
            mean_gc = mean(GC),
            max_gc = max(GC),
            min_gc = min(GC),
            mean_cd = mean(coding_density),
            max_cd = max(coding_density),
            min_cd = min(coding_density))) %>% kable()
```

| revised_group | n | mean_length_kb | max_length | min_length | mean_gc | max_gc | min_gc | mean_cd | max_cd | min_cd |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| I | 9 | 453.0559 | 587656 | 143197 | 53.54444 | 54.0 | 53.2 | 85.84444 | 86.3 | 85.3 |
| II | 3 | 338.6477 | 393604 | 228742 | 57.36667 | 57.5 | 57.1 | 83.46667 | 87.3 | 81.5 |
| III | 5 | 147.3516 | 157450 | 139938 | 51.82000 | 52.2 | 51.5 | 85.32000 | 85.8 | 84.9 |
| IV | 6 | 332.1062 | 366385 | 307330 | 53.88333 | 54.2 | 53.6 | 84.78333 | 85.1 | 84.6 |
| ND | 1 | 161.3300 | 161330 | 161330 | 56.40000 | 56.4 | 56.4 | 88.30000 | 88.3 | 88.3 |

### Compare with ‘complete’ Pseudomonasdb genomes.

Genome length and GC content of PseudomonasDB (DB version 22.1
(2023-10-06))

``` bash
echo -e "file\tlength\tgc_content" > fna_complete_summary.tsv

INPUT_DIR="/Volumes/bottlenose/DATABASES/pseudomonasdb/fna_complete"

for file in "$INPUT_DIR"/*.fna; do
    if [[ -f "$file" ]]; then
        # Extract filename
        filename=$(basename "$file")
        
        echo $file
        
        # Remove FASTA headers and concatenate all sequence lines
        sequence=$(grep -v '^>' "$file" | tr -d '\n' | tr '[:lower:]' '[:upper:]')
        
        # Calculate total length
        total_length=${#sequence}
        
        # Count G and C bases
        gc_count=$(echo "$sequence" | grep -o '[GC]' | wc -l)

        # Calculate GC content percentage
        if [[ $total_length -gt 0 ]]; then
            gc_content=$(awk "BEGIN { printf \"%.2f\", ($gc_count / $total_length) * 100 }")
        else
            gc_content="0.00"
        fi
        
        echo -e "$filename\t$total_length\t$gc_content" >> fna_complete_summary.tsv
    fi
done 
```

Get the summary from this file. Add a filter to remove unlikely
sequences.

``` r
pseudomonasdb_complete <- read.table("/Volumes/bottlenose/DATABASES/pseudomonasdb/fna_complete_summary.tsv",
                                     header=TRUE) %>% filter(length > 2e6)

(psdb_complete_summ <- pseudomonasdb_complete %>% 
  summarise(n = n(),
            mean_length = mean(length),
            median_length = median(length),
            max_length = max(length),
            min_length = min(length),
            sd_length = sd(length),
            mean_gc = mean(gc_content),
            median_gc = median(gc_content),
            max_gc = max(gc_content),
            min_gc = min(gc_content),
            sd_gc = sd(gc_content))) %>%
  kable()
```

| n | mean_length | median_length | max_length | min_length | sd_length | mean_gc | median_gc | max_gc | min_gc | sd_gc |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1321 | 6430819 | 6504659 | 7840830 | 3217366 | 621661.5 | 63.61528 | 65.16 | 67.43 | 38.05 | 2.940273 |

``` r
(pqbr132_ratio <- (filter(cropped_summaries, Plasmid == "pQBR132") %>% select(length) %>% pull()) /
  psdb_complete_summ$median_length) # 0.0215135
```

    ## [1] 0.0215135

``` r
(pqbr49_ratio <- (filter(cropped_summaries, Plasmid == "pQBR49") %>% select(length) %>% pull()) /
  psdb_complete_summ$median_length) # 0.09034386
```

    ## [1] 0.09034386

Calculate whether the GC content of each plasmid was significantly lower
than the chromosome.

``` r
ggplot(data = pseudomonasdb_complete, aes(x=gc_content)) + geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](1_Assembly_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

``` r
## Not normally distributed; use a nonparametric test

testGCvsChr <- function(x){
  test <- wilcox.test(pseudomonasdb_complete$gc_content, mu = x)
  return(test$p.value)
}

cropped_summaries %>% 
  rowwise() %>%
  mutate(wilcox_gc_vs_chromosome = testGCvsChr(GC)) %>% kable()
```

| Group | revised_group | Plasmid | assembly_details | length | GC | coding_density | n_cds | wilcox_gc_vs_chromosome |
|:---|:---|:---|:---|---:|---:|---:|---:|---:|
| I | I | pQBR103 | assembled | 425094 | 53.2 | 86.0 | 550 | 0 |
| I | I | pQBR11 | draft | 430820 | 53.2 | 86.3 | 558 | 0 |
| I | I | pQBR124 | assembled | 424568 | 53.2 | 85.9 | 547 | 0 |
| I | I | pQBR4 | assembled | 525650 | 53.7 | 85.8 | 658 | 0 |
| I | I | pQBR43 | assembled | 510643 | 53.8 | 86.1 | 614 | 0 |
| I | I | pQBR44 | draft | 143197 | 53.7 | 85.3 | 191 | 0 |
| I | I | pQBR49 | assembled | 587656 | 53.7 | 85.5 | 712 | 0 |
| I | I | pQBR5 | assembled | 466604 | 53.4 | 85.8 | 590 | 0 |
| I | I | pQBR50 | draft | 563271 | 54.0 | 85.9 | 684 | 0 |
| II | II | pQBR23 | assembled | 393597 | 57.5 | 81.5 | 499 | 0 |
| II | II | pQBR24 | assembled | 393604 | 57.5 | 81.6 | 499 | 0 |
| II | II | pQBR26 | assembled | 228742 | 57.1 | 87.3 | 273 | 0 |
| IV | III | pQBR127 | assembled | 140415 | 51.6 | 85.3 | 177 | 0 |
| ND | III | pQBR132 | assembled | 139938 | 51.6 | 85.8 | 176 | 0 |
| III | III | pQBR28 | assembled | 141505 | 51.5 | 84.9 | 176 | 0 |
| III | III | pQBR53 | assembled | 157450 | 52.2 | 85.3 | 199 | 0 |
| III | III | pQBR55 | assembled | 157450 | 52.2 | 85.3 | 198 | 0 |
| IV | IV | pQBR102 | assembled | 334884 | 53.6 | 85.1 | 421 | 0 |
| ND | IV | pQBR150 | assembled | 366385 | 54.2 | 84.9 | 462 | 0 |
| IV | IV | pQBR30 | assembled | 310323 | 53.7 | 84.6 | 390 | 0 |
| I | IV | pQBR51 | assembled | 366385 | 54.2 | 84.9 | 462 | 0 |
| IV | IV | pQBR56 | assembled | 307330 | 53.8 | 84.6 | 395 | 0 |
| IV | IV | pQBR57 | assembled | 307330 | 53.8 | 84.6 | 395 | 0 |
| ND | ND | pQBR105 | assembled | 161330 | 56.4 | 88.3 | 172 | 0 |

All plasmids have a GC content significantly lower than the chromosome.

Examine ANI data for the different Groups.

``` r
ani <- read.table("./1_sketches/ani.tsv", header=FALSE,
                 col.names = c("q", "r", "ani", "frags", "total_frags")) %>%
  mutate(query = sub(".*/(pQBR[0-9]+).*.fna", "\\1", q),
         ref = sub(".*/(pQBR[0-9]+).*.fna", "\\1", r),
         q = basename(q),
         r = basename(r))

match_groups <- cropped_summaries %>% 
  mutate(ref = Plasmid, ref_group = revised_group,
         query = Plasmid, query_group = revised_group) %>%
  select(ref, query, ref_group, query_group)

ani_annot <- ani %>% 
  left_join(select(match_groups, query, query_group), by="query") %>%
  left_join(select(match_groups, ref, ref_group), by="ref") %>%
  mutate(frag_rat = frags/total_frags) %>%
  select(query, query_group, ref, ref_group, ani, frag_rat)
  
(within_groups_summ <- ani_annot %>% filter(ref_group == query_group) %>%
  group_by(ref_group) %>% summarise(n = sqrt(n()),
                                    min_ani = min(ani),
                                    mean_ani = mean(ani),
                                    max_ani = max(ani),
                                    minfrags = min(frag_rat),
                                    meanfrags = mean(frag_rat),
                                    maxfrags = max(frag_rat))) %>% kable()
```

| ref_group |   n |  min_ani |  mean_ani | max_ani |  minfrags | meanfrags |  maxfrags |
|:----------|----:|---------:|----------:|--------:|----------:|----------:|----------:|
| I         |  10 |  98.5342 |  99.58990 |     100 | 0.2178723 | 0.8389526 | 1.0000000 |
| II        |   3 |  93.5298 |  97.26161 |     100 | 0.1410419 | 0.6406145 | 0.9974587 |
| III       |   6 |  95.5494 |  98.62566 |     100 | 0.7961783 | 0.9390469 | 1.0000000 |
| IV        |   7 |  99.1892 |  99.69788 |     100 | 0.8237705 | 0.9480383 | 1.0000000 |
| ND        |   1 | 100.0000 | 100.00000 |     100 | 0.8757764 | 0.8757764 | 0.8757764 |

``` r
(between_groups <- ani_annot %>% filter(ref_group != query_group))
```

    ##     query query_group     ref ref_group     ani   frag_rat
    ## 1 pQBR150          IV   pQBR4         I 95.1008 0.22540984
    ## 2   pQBR4           I pQBR150        IV 95.4101 0.15699334
    ## 3   pQBR4           I  pQBR51        IV 95.4101 0.15699334
    ## 4   pQBR4           I  pQBR53       III 94.1893 0.06089439
    ## 5   pQBR4           I  pQBR55       III 94.1893 0.06089439
    ## 6  pQBR51          IV   pQBR4         I 95.1008 0.22540984
    ## 7  pQBR55         III   pQBR4         I 94.4375 0.20063694

``` r
ani_annot %>% filter(ref_group == "I" & query_group == "I") %>% arrange(-frag_rat)
```

    ##       query query_group     ref ref_group      ani  frag_rat
    ## 1   pQBR103           I pQBR103         I 100.0000 1.0000000
    ## 2   pQBR103           I pQBR103         I 100.0000 1.0000000
    ## 3   pQBR103           I pQBR103         I 100.0000 1.0000000
    ## 4   pQBR103           I pQBR103         I 100.0000 1.0000000
    ## 5    pQBR11           I  pQBR11         I 100.0000 1.0000000
    ## 6   pQBR124           I pQBR124         I 100.0000 1.0000000
    ## 7   pQBR124           I   pQBR5         I  99.9984 1.0000000
    ## 8   pQBR124           I   pQBR4         I  99.9887 1.0000000
    ## 9   pQBR124           I  pQBR49         I  99.9812 1.0000000
    ## 10   pQBR44           I  pQBR44         I 100.0000 1.0000000
    ## 11    pQBR5           I   pQBR5         I 100.0000 1.0000000
    ## 12    pQBR5           I   pQBR4         I  99.9864 1.0000000
    ## 13    pQBR5           I  pQBR49         I  99.9804 1.0000000
    ## 14   pQBR44           I   pQBR4         I  99.9445 0.9965035
    ## 15  pQBR103           I pQBR124         I  99.5550 0.9964706
    ## 16  pQBR103           I   pQBR5         I  99.5533 0.9964706
    ## 17  pQBR103           I pQBR124         I  99.5550 0.9964706
    ## 18  pQBR103           I   pQBR5         I  99.5533 0.9964706
    ## 19  pQBR103           I   pQBR4         I  99.5408 0.9952941
    ## 20  pQBR103           I   pQBR4         I  99.5408 0.9952941
    ## 21  pQBR124           I pQBR103         I  99.5689 0.9952886
    ## 22  pQBR124           I pQBR103         I  99.5689 0.9952886
    ## 23  pQBR103           I  pQBR49         I  99.5565 0.9941176
    ## 24  pQBR103           I  pQBR49         I  99.5565 0.9941176
    ## 25  pQBR103           I  pQBR11         I  98.9089 0.9800000
    ## 26  pQBR103           I  pQBR11         I  98.9089 0.9800000
    ## 27  pQBR124           I  pQBR11         I  98.6480 0.9776207
    ## 28    pQBR5           I  pQBR50         I  99.9894 0.9753483
    ## 29  pQBR124           I  pQBR50         I  99.9921 0.9729093
    ## 30   pQBR11           I pQBR103         I  98.7993 0.9721254
    ## 31   pQBR11           I pQBR103         I  98.7993 0.9721254
    ## 32  pQBR103           I  pQBR50         I  99.5325 0.9705882
    ## 33  pQBR103           I  pQBR50         I  99.5325 0.9705882
    ## 34   pQBR11           I pQBR124         I  98.6698 0.9651568
    ## 35   pQBR11           I  pQBR49         I  98.6737 0.9639954
    ## 36   pQBR11           I   pQBR4         I  98.6643 0.9639954
    ## 37   pQBR11           I   pQBR5         I  98.6725 0.9628339
    ## 38   pQBR11           I  pQBR50         I  98.6370 0.9372822
    ## 39    pQBR5           I  pQBR43         I  99.9964 0.9303323
    ## 40  pQBR124           I  pQBR43         I  99.9907 0.9257951
    ## 41    pQBR4           I   pQBR4         I 100.0000 0.9257850
    ## 42  pQBR103           I  pQBR43         I  99.5192 0.9200000
    ## 43  pQBR103           I  pQBR43         I  99.5192 0.9200000
    ## 44    pQBR5           I pQBR124         I  99.9978 0.9099678
    ## 45    pQBR5           I pQBR103         I  99.5245 0.9078242
    ## 46    pQBR5           I pQBR103         I  99.5245 0.9078242
    ## 47   pQBR11           I  pQBR43         I  98.5543 0.8931475
    ## 48    pQBR5           I  pQBR11         I  98.6367 0.8906752
    ## 49    pQBR4           I   pQBR5         I  99.9833 0.8905804
    ## 50    pQBR4           I  pQBR49         I  99.9839 0.8896289
    ## 51   pQBR44           I  pQBR49         I  99.8283 0.8846154
    ## 52   pQBR44           I pQBR124         I  99.8746 0.8811189
    ## 53   pQBR44           I  pQBR50         I  99.8579 0.8811189
    ## 54   pQBR44           I   pQBR5         I  99.8859 0.8776224
    ## 55   pQBR44           I  pQBR43         I  99.8530 0.8776224
    ## 56   pQBR44           I  pQBR11         I  98.7026 0.8741259
    ## 57    pQBR4           I  pQBR50         I  99.9935 0.8696480
    ## 58   pQBR44           I pQBR103         I  99.2941 0.8636364
    ## 59   pQBR44           I pQBR103         I  99.2941 0.8636364
    ## 60   pQBR49           I  pQBR49         I 100.0000 0.8621277
    ## 61   pQBR43           I  pQBR50         I  99.9934 0.8589618
    ## 62   pQBR43           I  pQBR43         I 100.0000 0.8579824
    ## 63   pQBR43           I   pQBR5         I  99.9904 0.8550441
    ## 64   pQBR43           I   pQBR4         I  99.9891 0.8550441
    ## 65   pQBR43           I  pQBR49         I  99.9824 0.8550441
    ## 66    pQBR4           I  pQBR43         I  99.9950 0.8296860
    ## 67   pQBR50           I  pQBR50         I 100.0000 0.8214920
    ## 68   pQBR50           I   pQBR4         I  99.9915 0.8108348
    ## 69   pQBR50           I  pQBR49         I  99.9806 0.8099467
    ## 70   pQBR50           I   pQBR5         I  99.9939 0.8090586
    ## 71    pQBR4           I pQBR124         I  99.9876 0.8078021
    ## 72    pQBR4           I pQBR103         I  99.5469 0.8058991
    ## 73    pQBR4           I pQBR103         I  99.5469 0.8058991
    ## 74   pQBR49           I   pQBR4         I  99.9888 0.7957447
    ## 75   pQBR49           I   pQBR5         I  99.9861 0.7957447
    ## 76    pQBR4           I  pQBR11         I  98.6053 0.7916270
    ## 77   pQBR49           I  pQBR50         I  99.9845 0.7812766
    ## 78   pQBR50           I  pQBR43         I  99.9959 0.7753108
    ## 79   pQBR43           I pQBR124         I  99.9393 0.7717924
    ## 80   pQBR43           I pQBR103         I  99.5070 0.7688541
    ## 81   pQBR43           I pQBR103         I  99.5070 0.7688541
    ## 82   pQBR43           I  pQBR11         I  98.5342 0.7522037
    ## 83   pQBR49           I  pQBR43         I  99.9891 0.7421277
    ## 84   pQBR50           I pQBR124         I  99.9903 0.7335702
    ## 85   pQBR50           I pQBR103         I  99.5607 0.7317940
    ## 86   pQBR50           I pQBR103         I  99.5607 0.7317940
    ## 87   pQBR49           I pQBR124         I  99.9506 0.7234043
    ## 88   pQBR49           I pQBR103         I  99.5227 0.7217021
    ## 89   pQBR49           I pQBR103         I  99.5227 0.7217021
    ## 90   pQBR50           I  pQBR11         I  98.6054 0.7193606
    ## 91   pQBR49           I  pQBR11         I  98.6329 0.7072340
    ## 92  pQBR124           I  pQBR44         I  99.7527 0.2968198
    ## 93  pQBR103           I  pQBR44         I  99.1167 0.2952941
    ## 94  pQBR103           I  pQBR44         I  99.1167 0.2952941
    ## 95   pQBR11           I  pQBR44         I  98.6863 0.2903600
    ## 96    pQBR4           I  pQBR44         I  99.7393 0.2740247
    ## 97    pQBR5           I  pQBR44         I  99.7500 0.2711683
    ## 98   pQBR43           I  pQBR44         I  99.5474 0.2497551
    ## 99   pQBR50           I  pQBR44         I  99.7461 0.2246892
    ## 100  pQBR49           I  pQBR44         I  99.4830 0.2178723

Output a nicely formatted table.

``` r
cropped_summaries %>% 
  select(Plasmid, assembly_details, length, GC, n_cds, coding_density, revised_group) %>%
  arrange(revised_group, Plasmid) %>% write.table("./tabs/tab1.tsv", sep="\t",
                                                  quote=FALSE, row.names=FALSE)
```

### Extract transposon sequences for the plasmids that were lost

pQBR1, pQBR8 and pQBR58 did not give complete plasmid assemblies, and
instead were indicative of plasmid loss (but transposon maintenance).

``` r
filter(summary, seq %in% c("UWC1pQBR1","UWC1pQBR8","UWC1pQBR58")) %>%
  ggplot(aes(x=log10(length), y=log10(cov), colour=node)) + geom_point() + facet_wrap(~seq) +
  geom_hline(yintercept=log10(10), linetype="dotted") + geom_vline(xintercept=log10(5000), linetype="dotted")
```

![](1_Assembly_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

Extract the transposon-encoding sequences from each of these. Blastn
analysis indicated that the assembled sequence from pQBR1 exactly
matched the higher-depth region from pQBR58. There’s a region at the end
of the pQBR1 (~850 bp) sequence that matches two other regions from
pQBR58. pQBR8 exactly matched that from pQBR1 and pQBR58 (from position
15-8463 in pQBR8), without this additional region.

Annotate with Bakta.

``` bash
for TRANSPOSON in pQBR1 pQBR8 pQBR58; do
bakta --db /pub60/jamesh/db --prefix ${TRANSPOSON}tn \
  --complete \
  --locus ${TRANSPOSON}tn_contig \
  --verbose \
  --output ${TRANSPOSON}tn \
  --threads 128 \
  --locus-tag ${TRANSPOSON}tn \
  --meta \
  ./${TRANSPOSON}_transposon.fasta
done
```

The transposon is 8,449 bp, with detactable transposase and a merRTPFADE
*mer* operon. Blast on TnCentral shows high identity (\>97%) and
coverage (\>98%) to transposon Tn512, a known *Pseudomonas* mercury
resistance transposon.

------------------------------------------------------------------------

**[Back to index.](../README.md)**
