Creating tree for Group I
================

Data analysis for manuscript “The pQBR mercury resistance plasmids: a
model set of sympatric environmental mobile genetic elements”. Code
written by Victoria T Orr (<Victoria.Orr@liverpool.ac.uk>), supervised
by [James P. J. Hall](mailto:j.p.j.hall@liverpool.ac.uk).

Code for creating trees and heatmaps for Group I plasmids.

``` r
library(ggplot2)
library(tidyr)
library(readr)
library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)
library(tidytree)
library(patchwork)
```

Prerequisite file: - PIRATE.gene_families.ordered.tsv

``` r
group_i_rel_pirate <- read.table("PIRATE.gene_families.ordered.tsv",
                                header=TRUE, sep="\t")

giprel_matrix <- data.frame(cbind(group_i_rel_pirate[2], ifelse(group_i_rel_pirate[23:length(group_i_rel_pirate)]=="", 0, 1))) %>%
  column_to_rownames("gene_family")

giprel_genefam_dendro <- hclust(dist(giprel_matrix))  
giprel_plasmid_dendro <- hclust(dist(t(giprel_matrix)))

giprel_genefam_order <- rownames(giprel_matrix)[giprel_genefam_dendro$order]
giprel_plasmid_order <- colnames(giprel_matrix)[giprel_plasmid_dendro$order]

giprel_long <- group_i_rel_pirate %>% select(gene_family, threshold, starts_with("pQ")) %>%
  pivot_longer(cols = starts_with("pQ"), names_to = "plasmid", values_to = "gene") %>% 
  filter(gene != "") %>%
  mutate(gene_family = factor(gene_family, levels=rev(giprel_genefam_order)),
         plasmid = factor(plasmid, levels=giprel_plasmid_order))

ggplot(giprel_long, aes(x=gene_family, y=plasmid)) + geom_tile()
```

![](GroupI_tree_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

#### Create a core genome tree for Group I and relatives

Check alignment from PIRATE. Trim gappy/misaligned columns. Align using
IQtree3. Remove pQGrel01 as it is known to be identical to pQBR103.

    -m MFP               Extended model selection followed by tree inference
    -B 1000              1000 bootstrap replicates
    --alrt 1000          1000 replicates for SH approximate likelihood ratio test
    -T AUTO              Automatic thread assignment

Need Seqkit, TrimAL and IQ-TREE

``` bash
seqkit grep -v -p "pQG1rel01" core_alignment.fasta \
> core_alignment_crop.fasta

trimal -in core_alignment_crop.fasta \
-out core_rel_aln.trimmed.fasta \
-gt 0.9

iqtree3 -s core_rel_aln.trimmed.fasta  \
-m MFP \
-B 1000 --alrt 1000 \
-T AUTO \
-pre group_i_tree

mv group_i_tree* .trees/
```

Note: several sequences were identical in the tree.

Plot tree.

``` r
if (!exists("is.waive")) {
  is.waive <- \(x) inherits(x, "waiver") # nolint
}

gi_tr <- read.tree("./trees/group_i_tree.treefile")

tr_mid <- midpoint(gi_tr)

tr_mid$tip.label <- recode(tr_mid$tip.label,
                           "pQG1rel02" = "pPma4326F(1)",
                           "pQG1rel03" = "pPma4326F(2)",
                           "pQG1rel04" = "unnamed1",
                           "pQG1rel05" = "pPlaYM7902A",
                           "pQG1rel06" = "unnamed2",
                           "pQG1rel07" = "pPlaYM7902A*",
                           "pQG1rel08" = "unnamed3",
                           "pQG1rel09" = "unnamed4")

tb <- as_tibble(tr_mid)

supports <- tb %>%
  transmute(node,
            sh_alrt = as.numeric(str_extract(label, "^[^/]+")),
            ufboot  = as.numeric(str_extract(label, "(?<=/).+$")))
```

    ## Warning: There was 1 warning in `transmute()`.
    ## ℹ In argument: `sh_alrt = as.numeric(str_extract(label, "^[^/]+"))`.
    ## Caused by warning:
    ## ! NAs introduced by coercion

``` r
p_sup <- ggtree(tr_mid) +
  theme_tree2() +
  geom_text2(
    aes(subset = !isTip, label =
          ifelse(!is.na(supports$ufboot),
                paste0("SH=", round(supports$sh_alrt),
                       " / UF=", round(supports$ufboot)),
                 round(as.numeric(label)))),
    hjust = -0.1, size = 2.8, color = "grey20") +
  ggtitle("Tree with SH-aLRT / UFBoot support")
```

    ## Warning: `aes_()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`
    ## ℹ The deprecated feature was likely used in the ggtree package.
    ##   Please report the issue at <https://github.com/YuLab-SMU/ggtree/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

``` r
p_sup
```

    ## Warning in ifelse(!is.na(supports$ufboot), paste0("SH=",
    ## round(supports$sh_alrt), : NAs introduced by coercion

![](GroupI_tree_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
tree <- ggtree(tr_mid) + geom_tiplab(size=6, align=TRUE) +xlim(0,0.2)
```

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the ggtree package.
    ##   Please report the issue at <https://github.com/YuLab-SMU/ggtree/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
tree_order <- get_taxa_name(p_sup)


giprel_long_filtered <- giprel_long[giprel_long$plasmid != "pQG1rel01", ]


giprel_long_filtered$plasmid <- recode(giprel_long_filtered$plasmid,
                           "pQG1rel02" = "pPma4326F(1)",
                           "pQG1rel03" = "pPma4326F(2)",
                           "pQG1rel04" = "unnamed1",
                           "pQG1rel05" = "pPlaYM7902A",
                           "pQG1rel06" = "unnamed2",
                           "pQG1rel07" = "pPlaYM7902A*",
                           "pQG1rel08" = "unnamed3",
                           "pQG1rel09" = "unnamed4")

heatmap <- ggplot(giprel_long_filtered,
                  aes(x=gene_family,y=factor(plasmid,
                                            levels=rev(tree_order)))) +
  geom_tile() +
  ylab(NULL) +
  xlab("Gene Family") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.y = element_blank())

tree + heatmap + plot_layout(widths = c(0.9,1))
```

![](GroupI_tree_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
ggsave("GroupI_heattree.png", dpi = 500, height=10, width=11)
```

------------------------------------------------------------------------

**[Back to index.](../4_Analysis.md)**
