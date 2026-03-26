Attempted classification of pQBR plasmids
================

## PInc classification

Attempt to classify the plasmids according to the IncP scheme outlined
by [Nishimura et al. 2024](https://doi.org/10.1101/2024.09.03.610885)

Download sequences described in the paper.

``` r
incp_groups <- read.csv("./ref/Nishimura2024.csv", header=TRUE, sep=",")

incp_groups %>% kable()
```

| Inc_group  | Alt_group  | Representative | Accession | MOB | MPF |
|:-----------|:-----------|:---------------|:----------|:----|:----|
| IncP-1     | IncP       | RP4            | BN000925  | P   | T   |
| IncP-2     |            | Rms139         | LC653116  | P   | I   |
| IncP-3A/C  | IncA/C     | pRA1           | FJ705807  | H   | F   |
| IncP-3A/C  | IncA/C     | pR55           | JQ010984  | H   | F   |
| IncP-4     | IncQ       | RSF1010        | M28829    | Q   | NA  |
| IncP-5     |            | Rms163         | LC685027  | NA  | NA  |
| IncP-6     | IncG       | Rms149         | AJ877225  | P   | NA  |
| IncP-7     |            | pCAR1          | AB088420  | H   | F   |
| IncP-9     |            | pWW0           | AJ344068  | F   | F   |
| IncP-9     |            | NAH7           | AB237655  | F   | F   |
| IncP-9     |            | Rsu2           | LC685593  | F   | F   |
| IncP-10    |            | R91-5          | X54695    | NA  | NA  |
| IncP-10    |            | pPAB546        | MN433456  | NA  | NA  |
| IncP-11    |            | RP1-1          | LC700336  | NA  | NA  |
| IncP-12    |            | R716           | LC685026  | NA  | NA  |
| IncP-13    |            | pMG26d         | LC685025  | F   | G   |
| IncP-15    | PromA      | pSN1104-11     | AP018707  | P   | T   |
| IncP-16    | pSN1216-29 | pSN1216-29     | AP018710  | P   | T   |
| unassigned | pQBR103    | pQBR103        | AM235768  | NA  | NA  |
| IncPSTY    | pSTY       | pSTY           | NC_022739 | F   | F   |

Download the sequences.

``` bash
mkdir -p 3_incp_groups/seqs

tail -n+2  ./ref/Nishimura2024.csv \
  | awk -v FS="," -v OFS="\t" '{print $4}' | while read ACCESSION
  do
  curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=${ACCESSION}&amp;rettype=fasta" > ./3_incp_groups/seqs/${ACCESSION}.fasta
  done
  
cat ./3_incp_groups/seqs/*.fasta > ./3_incp_groups/seqs.fasta
```

Do a MASH comparison.

``` bash
mkdir ./3_incp_groups/mash_results

mash sketch -i -s 100000 ./3_incp_groups/seqs.fasta -o ./3_incp_groups/seqs.msh

cat ./1_sketches/annotated_plasmids.list | while read PLASMID
do
mash dist -p 64 \
  ./1_sketches/${PLASMID}.msh \
  ./3_incp_groups/seqs.msh \
  | sort -gk3 \
  > ./3_incp_groups/mash_results/${PLASMID}_incpref_mash_results.txt
head -n1 ./3_incp_groups/mash_results/${PLASMID}_incpref_mash_results.txt >> ./3_incp_groups/mash_tophits.tsv
done
```

Examine in R.

``` r
mash_tophits <- read.table("./3_incp_groups/mash_tophits.tsv", header=FALSE,
           col.names=c("query","subject","dist","e_val","num_hits")) %>%
  mutate(plasmid = gsub(".*(pQBR[0-9Rpd]+).fna", "\\1", query),
         Accession = gsub("\\.[0-9]", "", subject),
         N_hits = gsub("/100000", "", num_hits))

mash_tophits %>% left_join(incp_groups, by = "Accession") %>% 
  select(plasmid, dist, N_hits, Inc_group, Alt_group, Representative, Accession, MOB, MPF) %>% 
  arrange(-dist) %>% kable()
```

| plasmid  |      dist | N_hits | Inc_group  | Alt_group | Representative | Accession | MOB | MPF |
|:---------|----------:|:-------|:-----------|:----------|:---------------|:----------|:----|:----|
| pQBR102  | 0.1904890 | 924    | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR30   | 0.1888330 | 957    | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR56   | 0.1886860 | 960    | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR57   | 0.1886860 | 960    | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR28   | 0.1762150 | 1251   | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR55R  | 0.1760650 | 1255   | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR127  | 0.1759900 | 1257   | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR132  | 0.1757290 | 1264   | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR106  | 0.1707920 | 1404   | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR105  | 0.1616020 | 1708   | IncP-12    |           | R716           | LC685026  | NA  | NA  |
| pQBR150  | 0.1265550 | 3633   | IncP-9     |           | pWW0           | AJ344068  | F   | F   |
| pQBR51   | 0.1265550 | 3633   | IncP-9     |           | pWW0           | AJ344068  | F   | F   |
| pQBR47   | 0.1262770 | 3655   | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR57R  | 0.1224510 | 3973   | IncP-9     |           | pWW0           | AJ344068  | F   | F   |
| pQBR53   | 0.0997823 | 6554   | IncP-9     |           | pWW0           | AJ344068  | F   | F   |
| pQBR55   | 0.0997755 | 6555   | IncP-9     |           | pWW0           | AJ344068  | F   | F   |
| pQBR26   | 0.0725135 | 12240  | IncP-7     |           | pCAR1          | AB088420  | H   | F   |
| pQBR106p | 0.0508376 | 20761  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR44d  | 0.0426369 | 25664  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR24   | 0.0296269 | 36685  | IncPSTY    | pSTY      | pSTY           | NC_022739 | F   | F   |
| pQBR23   | 0.0296221 | 36690  | IncPSTY    | pSTY      | pSTY           | NC_022739 | F   | F   |
| pQBR47p  | 0.0201117 | 48755  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR43   | 0.0073717 | 74914  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR49   | 0.0072453 | 75263  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR11d  | 0.0068294 | 76428  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR4    | 0.0061848 | 78285  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR50d  | 0.0060468 | 78691  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR5    | 0.0052690 | 81037  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR124  | 0.0029878 | 88534  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR103  | 0.0000024 | 99990  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |
| pQBR103R | 0.0000024 | 99990  | unassigned | pQBR103   | pQBR103        | AM235768  | NA  | NA  |

This table indicates:

- The Group I plasmids match the pQBR103 IncP group defined by Nishimura
  et al. (pQBR103), as expected.
- The Group II plasmids (roughly) match the IncPSTY group defined by
  Nishimura et al. (pSTY)
- The Group IIb plasmid pQBR26 roughly matches the IncP-7 group (pCAR1)
- The non-determined plasmid pQBR105 matches IncP-12 group (R716)
- The plasmids matching the IncP-9 plasmid (pWW0) are likely matching
  through the Tn4652 transposon, since they include representatives of
  Group III and Group IV that have Tn4652, and pWW0 has the complete
  version of this *tol* transposon.
- The Group III plasmids and Group IV plasmids do not match one of the
  defined IncP groups.

Run a BLAST comparison to visualise.

``` bash
makeblastdb -dbtype nucl -in ./3_incp_groups/seqs/NC_022739.fasta
makeblastdb -dbtype nucl -in ./3_incp_groups/seqs/AB088420.fasta
makeblastdb -dbtype nucl -in ./3_incp_groups/seqs/LC685026.fasta


mkdir ./3_incp_groups/blast

tblastx -query ./bakta_annotated/pQBR23/pQBR23.fna -db ./3_incp_groups/seqs/NC_022739.fasta \
  -outfmt 6 > ./3_incp_groups/blast/pQBR23_pSTY.tblastx
  
tblastx -query ./bakta_annotated/pQBR26/pQBR26.fna -db ./3_incp_groups/seqs/AB088420.fasta \
  -outfmt 6 > ./3_incp_groups/blast/pQBR26_pCAR1.tblastx
  
tblastx -query ./bakta_annotated/pQBR105/pQBR105.fna -db ./3_incp_groups/seqs/LC685026.fasta \
  -outfmt 6 > ./3_incp_groups/blast/pQBR105_R716.tblastx
```

Examined using ACT. Rough plots here with ggplot.

``` r
blast_columns <- c("queryId","subjectId","percIdentity", "alnLength",
                   "mismatchCount", "gapOpenCount", "queryStart", "queryEnd", 
                   "subjectStart", "subjectEnd", "eVal", "bitScore")

pQBR23_pSTY_tblastx <- read.table("./3_incp_groups/blast/pQBR23_pSTY.tblastx",
                                  col.names=blast_columns)

makeComparison <- function(x){
  vals <- x %>% mutate(id = row_number()) %>%
    filter(eVal < 1e-100 & alnLength > 100) %>%
    select(id, percIdentity, queryStart, queryEnd, subjectEnd, subjectStart) %>%
    mutate(direction = case_when(queryEnd > queryStart & subjectEnd > subjectStart ~ "F",
                                 queryEnd < queryStart & subjectEnd > subjectStart ~ "R",
                                 queryEnd > queryStart & subjectEnd < subjectStart ~ "R",
                                 queryEnd < queryStart & subjectEnd < subjectStart ~ "F")) %>%
    pivot_longer(cols = c("queryStart","queryEnd","subjectEnd","subjectStart"),
                 names_to = "type",
                 values_to = "x_coord")
  
  vals <- vals %>%
    mutate(y_coord = rep(c(1,1,2,2), length.out = nrow(vals)))
  
  ggplot(data = vals) +
    geom_polygon(aes(x = x_coord, y=(0-y_coord), group=id, fill = direction, alpha=percIdentity))
}

(pQBR23_pSTY_comparisons <- makeComparison(pQBR23_pSTY_tblastx))
```

![](3_Classification_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
pQBR26_pCAR1_tblastx <- read.table("./3_incp_groups/blast/pQBR26_pCAR1.tblastx",
                                  col.names=blast_columns)
(pQBR26_pCAR1_comparisons <- makeComparison(pQBR26_pCAR1_tblastx))
```

![](3_Classification_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
pQBR105_R716_tblastx <- read.table("./3_incp_groups/blast/pQBR105_R716.tblastx",
                                  col.names=blast_columns)
(pQBR105_R716_comparisons <- makeComparison(pQBR105_R716_tblastx))
```

![](3_Classification_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

pQBR23 and pQBR26 match their corresponding groups well, but pQBR105
matches its corresponding plasmid R716 only through a transposon.

``` bash
blastn -query ./bakta_annotated/pQBR23/pQBR23.fna -db ./3_incp_groups/seqs/NC_022739.fasta \
  -perc_identity 80 \
  -outfmt "6 std qcovus" # 56% coverage at >80% identity
  
blastn -query ./bakta_annotated/pQBR24/pQBR24.fna -db ./3_incp_groups/seqs/NC_022739.fasta \
  -perc_identity 80 \
  -outfmt "6 std qcovus" # 56% coverage at >80% identity

blastn -query ./bakta_annotated/pQBR26/pQBR26.fna -db ./3_incp_groups/seqs/AB088420.fasta \
  -perc_identity 80 \
  -outfmt "6 std qcovus" # 42% coverage at >80% identity

blastn -query ./bakta_annotated/pQBR105/pQBR105.fna -db ./3_incp_groups/seqs/LC685026.fasta  \
  -perc_identity 80 \
  -outfmt "6 std qcovus" # only 4% coverage at >80% identity
```

------------------------------------------------------------------------

**[Back to index.](../README.md)**
