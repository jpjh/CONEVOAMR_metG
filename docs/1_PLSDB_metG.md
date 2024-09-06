Plasmid-encoded *metG* homologues part 1: bacteria
================
JPJ Hall and L Santos
September 2024

## Identifying plasmid-borne metG homologues

Genomic epidemiology indicated a plasmid-borne *metG* homologue as a key
genetic factor associated with an epidemic strain of Shigella.

Here, we attempt to address two main questions:

1.  How widespread are plasmid-encoded *metG* homologues in the PLSDB
    database ([Schwartz et
    al. 2022](https://academic.oup.com/nar/article/50/D1/D273/6439675))
    database?
2.  What taxa harbour such plasmids?

### Annotation of PLSDB

The [PLSDB database version 2023_11_23 was
downloaded](https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/),
and annotated using Bakta ([Schwengers et
al. 2021](https://doi.org/10.1099/mgen.0.000685)) version 1.8.2 with the
‘full’ database.

Annotation was run in batches. Example command:

``` bash
bakta -v --threads 64 --db /pub60/jamesh/db \
  --complete \
  --skip-crispr \
  --locus-tag plsdb_${GROUP} \
  --skip-plot \
  --meta \
  --output ./PLSDB_${GROUP}_bakta \
  --keep-contig-headers plsdb_${GROUP}.fna.gz
```

### Extraction of metG homologues from the database

Homologues of metG were detected by searching the `.gff3` annotation
file for genes matching “\[Mm\]etG”.

``` bash
grep "[Mm]etG" ./PLSDB_*_bakta/plsdb_*.fna.gff3 \
  | sed 's/.*.fna.gff3://g' \
  > ./1_data/PLSDB_bakta_MetG_matches.gff
```

The corresponding gene names were extracted for examination

``` bash
awk -v FS=";" '{print $2}' ./1_data/PLSDB_bakta_MetG_matches.gff \
  | sed 's/Name=//g' | sort | uniq
```

- Methionine–tRNA ligase
- Methionine–tRNA ligase 1
- Methionyl-tRNA synthetase
- methionine–tRNA ligase
- tRNA-binding protein

All makes sense.

Extract `locus_tag` IDs for subsequent analysis, as well as the amino
acid and nucleic acid sequences.

``` bash
awk -v FS=";" '{print $1}' ./1_data/PLSDB_bakta_MetG_matches.gff \
  | sed 's/.*ID=//g' > ./1_data/PLSDB_bakta_MetG_matches.txt
  
cat ./PLSDB_*_bakta/plsdb_*.fna.faa \
  | seqtk subseq - \
  ./1_data/PLSDB_bakta_MetG_matches.txt \
  | awk '{print $1}' \
  | tr '\n' '\t' | tr '>' '\n' \
  > ./1_data/PLSDB_bakta_MetG_matches.tsv
  
cat ./PLSDB_*_bakta/plsdb_*.fna.ffn \
  | seqtk subseq - \
  ./1_data/PLSDB_bakta_MetG_matches.txt \
  | awk '{print $1}' \
  | tr '\n' '\t' | tr '>' '\n' \
  > ./1_data/PLSDB_bakta_MetG_matches_na.tsv
```

### Analysis of PLSDB plasmids

The corresponding PLSDB metadata was downloaded from the database
website.

Load up PLSDB metadata.

``` r
plsdb <- read.table("./1_data/plsdb.tsv", 
                    header = TRUE, 
                    sep = "\t", 
                    fill = TRUE, 
                    quote = "", 
                    na.strings = "") %>%
  filter(!is.na(NUCCORE_UID)) %>%
  rename(Accession = NUCCORE_ACC)
```

Note that preliminary analyses showed there was at least one one
sequence in PLSDB that was annotated with genus but not class. Correct
this.

``` r
plsdb <- plsdb %>%
  mutate(TAXONOMY_class = ifelse(NUCCORE_UID == 1563790493, "Mollicutes", TAXONOMY_class))
```

Load up functions for extracting information from `.gff` filesm and load
up the gff file, which contains information on metG-encoding CDS.

``` r
source("./func/gff_functions.R")

metg_gff_file <- gffRead("./1_data/PLSDB_bakta_MetG_matches.gff")
```

    ## Reading ./1_data/PLSDB_bakta_MetG_matches.gff: found 316 rows with classes: character, character, character, integer, integer, character, character, character, character

``` r
metg_gff <- data.frame(Accession = metg_gff_file$seqname,
                       start = metg_gff_file$start,
                       end = metg_gff_file$end,
                       strand = metg_gff_file$strand,
                       ID = getAttributeField(metg_gff_file$attributes, "ID"))
```

Merge with the PLSDB information.

``` r
plasmid_metg <- metg_gff %>%
  left_join(plsdb, by = "Accession")

nrow(plasmid_metg)
```

    ## [1] 316

``` r
length(unique(plasmid_metg$Accession))
```

    ## [1] 309

`plasmid_metg` contains details of plasmid-encoded metG genes.

316 matches for plasmid-encoded metG, and 309 unique plasmids with metG
homologues.

#### Proportion plot

Make a plot of counts for different taxa.

``` r
class_colours <- c("Actinobacteria" = "deeppink",
                   "Alphaproteobacteria" = "blueviolet",   
                   "Bacilli" = "darkseagreen",
                   "Clostridia" = "burlywood",
                   "Gammaproteobacteria" = "cornflowerblue",   
                   "Planctomycetia" = "darkorange2",
                   "Thermomicrobia" = "yellow",
                   "Actinomycetes" = "deeppink",
                   "Bacteroidia" = "salmon",
                   "Betaproteobacteria" = "brown",
                   "Coriobacteriia" = "darkslategray2",
                   "Cyanophyceae" = "lavender",
                   "Mollicutes" = "darkcyan",
                   "NA" = "grey",
                   "Other" = "grey")   

ggplot(data=plasmid_metg, aes(x=TAXONOMY_genus, fill=TAXONOMY_class)) + 
  geom_bar() + 
  scale_fill_manual(values = class_colours, name="") +
  ylim(0, 200) +
  labs(title = "Number of plasmid-borne metG homologues in PLSDB",
    x = "Genus",
    y = "Plasmid Count") +
  theme(axis.text.x = element_text(face="italic", angle=45, hjust=1),
        legend.position = "bottom",
        legend.direction="horizontal")
```

![](1_PLSDB_metG_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Create proportions (proportion of plasmids with metG by Genus).

``` r
plsdb_total_genus <- plsdb %>% group_by(TAXONOMY_genus) %>% count()

plasmid_metg_proportions <- plasmid_metg %>%
  group_by(TAXONOMY_class, TAXONOMY_genus, Accession) %>%
  summarise(n = n()) %>%  # first 'peels off' the counts per Accession  
  summarise(n = n()) %>%  # next peels of the counts of #Accession per genus
  left_join(plsdb_total_genus, by = "TAXONOMY_genus") %>%
  mutate(proportion = n.x/n.y)
```

    ## `summarise()` has grouped output by 'TAXONOMY_class', 'TAXONOMY_genus'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'TAXONOMY_class'. You can override using
    ## the `.groups` argument.

Note: there is one sequence with NA class and NA genus.

Generate proportion plot.

``` r
plsdb_metg_proportions_text <- plasmid_metg_proportions %>%
  group_by(TAXONOMY_class, TAXONOMY_genus) %>%
  summarise(positive = sum(n.x), total = sum(n.y)) 
```

    ## `summarise()` has grouped output by 'TAXONOMY_class'. You can override using
    ## the `.groups` argument.

``` r
ggplot(data=plasmid_metg_proportions, aes(x=TAXONOMY_genus, y=proportion, fill=TAXONOMY_class)) +  
  geom_bar(stat="identity") + 
  scale_fill_manual(values = class_colours, name="") +
  labs(title = "The Proportion of Plasmids per Genus Containing ≥1 MetG Homologues in Bakta-annotated PLSDB", 
       x = "Genus", 
       y = "Proportion of Plasmids With ≥1 MetG Homologue", 
       fill = "Class") + 
  geom_text(data=plsdb_metg_proportions_text, 
            aes(label=paste(positive, " / ", total, sep=""), y=-0.05), 
            hjust=0.5, size=2.8, angle=60) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0))) + 
  theme(legend.position="bottom", 
        axis.text.x=element_text(face="italic", angle=45, hjust=1))
```

![](1_PLSDB_metG_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Joy Plot

Make a ‘joy plot’ with size histogram density for different families,
annotated with the size of plasmids harbouring metG homologues.

Put ‘encodes metG’ information into the big PLSDB table, to facilitate
filtering and analysis.

``` r
plsdb <- plsdb %>%
  mutate(encodes_metg = ifelse(Accession %in% plasmid_metg$Accession, "Y", "N"))
```

Group by family and calculate counts…

``` r
metg_families <- plasmid_metg$TAXONOMY_family %>% unique()

plsdb_subset_fam <- filter(plsdb, TAXONOMY_family %in% metg_families)

(counts_family <- plsdb_subset_fam %>%
    group_by(TAXONOMY_family, encodes_metg) %>%
    summarise(total = n()) %>%
    pivot_wider(names_from = encodes_metg,
                values_from = total,
                values_fill = 0) %>%            
    mutate(totals = N + Y)) %>%                  
  kable()
```

    ## `summarise()` has grouped output by 'TAXONOMY_family'. You can override using
    ## the `.groups` argument.

| TAXONOMY_family      |     N |   Y | totals |
|:---------------------|------:|----:|-------:|
| Acaryochloridaceae   |    27 |   2 |     29 |
| Actinomycetaceae     |     6 |   1 |      7 |
| Aeromonadaceae       |   464 |   1 |    465 |
| Bacillaceae          |  2448 |   3 |   2451 |
| Bacteroidaceae       |   436 |   2 |    438 |
| Burkholderiaceae     |   558 |  99 |    657 |
| Coriobacteriaceae    |    15 |   1 |     16 |
| Enterobacteriaceae   | 29078 | 176 |  29254 |
| Lactobacillaceae     |  2186 |   1 |   2187 |
| Legionellaceae       |    93 |   1 |     94 |
| Metamycoplasmataceae |    71 |   1 |     72 |
| Moraxellaceae        |  2342 |   1 |   2343 |
| Morganellaceae       |   512 |   2 |    514 |
| Neisseriaceae        |   209 |   2 |    211 |
| Nocardiaceae         |   329 |   1 |    330 |
| Oscillospiraceae     |    75 |   3 |     78 |
| Rhizobiaceae         |  1263 |   5 |   1268 |
| Streptosporangiaceae |     8 |   1 |      9 |
| Tannerellaceae       |   113 |   1 |    114 |
| Tsukamurellaceae     |    26 |   1 |     27 |
| Vibrionaceae         |   557 |   1 |    558 |
| Xanthomonadaceae     |   525 |   2 |    527 |
| NA                   |  1035 |   1 |   1036 |

Plot only Families with \>20 plasmids in PLSDB, and \>1 metG-encoding
plasmid, to avoid misannotations.

``` r
family_count_20 <- filter(counts_family, totals > 20 & Y > 1)

plsdb_subset_fam_20 <- filter(plsdb_subset_fam, TAXONOMY_family %in% family_count_20$TAXONOMY_family)
```

Create the Joy Plot.

``` r
metg_joy_plot_fam <- ggplot(plsdb_subset_fam_20, aes(x = log10(NUCCORE_Length))) +
  geom_density(linewidth = 0, aes(fill = TAXONOMY_class), alpha = 0.8) +
   geom_text(data = family_count_20, 
             aes(x = 2.5, y = 0, 
                 label = paste("n = ", Y, " / ", totals, sep = "")), 
             hjust = 0, vjust = 0, size = 3) +
  scale_fill_manual(values = class_colours, name = "") +
  geom_point(data = subset(plsdb_subset_fam_20, encodes_metg == "Y"), 
             aes(y = 0.1), shape = 17, alpha = 0.6, colour = "black") + 
  facet_grid(TAXONOMY_family ~ ., scales = "free_y") +
  labs(x = expression(paste("log"[10], " plasmid size in bp")), y = "") +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),         
    axis.title.x = element_text(size = 14),         
    axis.title.y = element_text(size = 14),         
    axis.text.x = element_text(size = 12),          
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

metg_joy_plot_fam
```

![](1_PLSDB_metG_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

In many families, the *metG* plasmids are extremely large, and possibly
represent secondary chromosomes and/or misannotated chromosomes.

Examine these. First the non-Burkholderiaceae.

``` r
plsdb %>% filter(TAXONOMY_family != "Burkholderiaceae" & 
                                 NUCCORE_Length > 1e6 & 
                                 encodes_metg == "Y") %>%
  select(Accession, NUCCORE_Description, NUCCORE_Length) %>% kable()
```

| Accession     | NUCCORE_Description                                                                | NUCCORE_Length |
|:--------------|:-----------------------------------------------------------------------------------|---------------:|
| NZ_LT559118.1 | Nonomuraea gerenzanensis isolate nono1 plasmid I, complete sequence                |       11850240 |
| NZ_CP119424.1 | Photobacterium sp. DA100 plasmid pDA100, complete sequence                         |        2139191 |
| NZ_CP093325.1 | Bacillus sp. N5-665 strain N5/665 plasmid p1, complete sequence                    |        5519593 |
| NZ_CP119142.1 | Neisseria gonorrhoeae strain 2022NG-0032 plasmid p2022NG-0032_1, complete sequence |        2167441 |
| NZ_CP119140.1 | Neisseria gonorrhoeae strain 2022NG-0076 plasmid p2022NG-0076_1, complete sequence |        2168363 |
| CP104135.1    | Sinorhizobium sp. C101 plasmid pC101a, complete sequence                           |        3521547 |
| CP104131.1    | Sinorhizobium sp. K101 plasmid pK101a, complete sequence                           |        3521393 |
| CP104127.1    | Sinorhizobium sp. M103 plasmid pM103a, complete sequence                           |        3521513 |
| NZ_OL790812.1 | Proteus mirabilis strain H110 plasmid unnamed, complete sequence                   |        4088254 |
| NC_008269.1   | Rhodococcus jostii RHA1 plasmid pRHL1, complete sequence                           |        1123075 |
| NZ_AP026896.1 | Aeromonas caviae strain KAM479 plasmid pKAM479_3, complete sequence                |        4944632 |
| CP061244.1    | Escherichia coli strain STEC411 plasmid pSTEC411_1, complete sequence              |        1412998 |

- [The NCBI BioProject page for
  NZ_LT559118.1](https://www.ncbi.nlm.nih.gov/bioproject/357576) implies
  that this 12 Mb ‘plasmid’ is the chromosomal sequence.
- The Photobacterium plasmid NZ_CP119424.1 is likely a chromid or
  secondary chromosome. *Photobacterium* is a Vibrionaceae, which are
  known for secondary chromosomes ([Urbanczyk et
  al. 2011](https://academic.oup.com/femsre/article/35/2/324/660773))
- [The NCBI assembly page for
  NZ_CP093325.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022637375.2/)
  shows that the huge 5 Mb Bacillus ‘plasmid’ is a misannotation (the
  largest replicon in the isolate).
- [The NCBI assembly page for
  NZ_CP119142.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029378455.1/)
  and [NZ_CP1119140.1](Z_CP119140.1) shows that these 2.2 Mb Neisseria
  plasmids are misannotations (the largest replicon in the isolate).
- [The NCBI assembly pages for
  CP10435.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029201665.1/),
  [CP104131.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029201645.1/),
  [CP104127.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029201625.1/)
  show that these 3.5 Mb Sinorhizobium plasmids are misannotations (the
  largest replicon in the isolate).
- There was no Assembly page for NZ_OL790812.1 but it resembles a
  chromosome with rRNA gene sequences.
- [The NCBI assembly page for
  NC_008269.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014565.1/)
  and the published genome sequence indicates that this is a (linear)
  megaplasmid ([McLeod et
  al. 2006](https://www.pnas.org/doi/full/10.1073/pnas.0607048103)).
- [The NCBI assembly page for
  NZ_AP026896.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026000055.1/)
  shows that this is a misannotated chromosome.
- CP061244 could be a possible megaplasmid (1.4 Mb), smaller than the
  annotated chromosome at 4.2 Mb.

Overall, these sequences are predominantly misannotations, and thus
should be removed from further analysis.

``` r
plsdb %>% 
  filter(TAXONOMY_family == "Burkholderiaceae" & encodes_metg == "Y") %>%
  group_by(TAXONOMY_species) %>%
  summarise(n = n(),
            min_size = min(NUCCORE_Length)/1000,
            max_size = max(NUCCORE_Length)/1000,
            mean_size = mean(NUCCORE_Length)/1000) %>%
  kable()
```

| TAXONOMY_species             |   n | min_size | max_size | mean_size |
|:-----------------------------|----:|---------:|---------:|----------:|
| Ralstonia_nicotianae         |   1 | 1994.515 | 1994.515 |  1994.515 |
| Ralstonia_pseudosolanacearum |   7 | 2008.534 | 2094.509 |  2070.806 |
| Ralstonia_solanacearum       |  86 | 1781.111 | 2326.525 |  2044.370 |
| Ralstonia_sp.\_RS642         |   1 | 1996.921 | 1996.921 |  1996.921 |
| Ralstonia_sp.\_RS647         |   1 | 2069.415 | 2069.415 |  2069.415 |
| Ralstonia_sp.\_RS650         |   1 | 2073.962 | 2073.962 |  2073.962 |
| Ralstonia_syzygii            |   2 | 1578.131 | 1614.911 |  1596.521 |

The Burkholderiaceae plasmids are mostly *Ralstonia solanacearum* and
related *Ralstonia*s. These organisms are known to have a multipartite
genome consisting of a chromosome and a chromid, as described by
[DiCenzo et
al. 2019](https://academic.oup.com/mbe/article/36/3/562/5273485), and
several of the plasmids in the dataset here were characterised as
chromids by DiCenzo et al.

Output a list of all the IDs of the metG encoding replicons in PLSDB.

``` r
filter(plsdb, encodes_metg == "Y") %>%
  select(NUCCORE_UID) %>% pull() %>%
  write.table(file = "./1_data/plsdb_metG_uid.txt",
            row.names=FALSE, sep=",", quote=FALSE, col.names=FALSE)

filter(plsdb, encodes_metg == "Y") %>%
  select(Accession) %>% pull() %>%
  write.table(file = "./1_data/plsdb_metG.txt",
            row.names=FALSE, sep=",", quote=FALSE, col.names=FALSE)
```

#### Filtered proportion plot

To focus on plasmids that represent likely horizontal gene transfers
(rather than core genome chromids), edit the plot to remove plasmid-rare
genera (genera with fewer than 10 plasmids in PLSDB), and remove
‘plasmids’ which are equal or above 1Mb.

``` r
metg_genera <- filter(plsdb, encodes_metg == "Y" & NUCCORE_Length < 1e6) %>%
  select(TAXONOMY_genus) %>% pull() %>% unique()

plsdb_subset_genera <- filter(plsdb, 
                              TAXONOMY_genus %in% metg_genera &
                                NUCCORE_Length < 1e6)

(metg_counts_genera <- plsdb_subset_genera %>% 
    group_by(TAXONOMY_class, TAXONOMY_genus, encodes_metg) %>%
    summarise(total = n()) %>% 
    pivot_wider(names_from=encodes_metg, values_from=total, values_fill=0) %>%
    mutate(totals = N + Y, proportion = Y/totals) %>%
    filter(Y > 0 & totals > 9 & !is.na(TAXONOMY_genus))) %>% kable()
```

    ## `summarise()` has grouped output by 'TAXONOMY_class', 'TAXONOMY_genus'. You can
    ## override using the `.groups` argument.

| TAXONOMY_class      | TAXONOMY_genus   |     N |   Y | totals | proportion |
|:--------------------|:-----------------|------:|----:|-------:|-----------:|
| Actinomycetes       | Tsukamurella     |    26 |   1 |     27 |  0.0370370 |
| Alphaproteobacteria | Rhizobium        |   571 |   2 |    573 |  0.0034904 |
| Bacilli             | Lactobacillus    |   148 |   1 |    149 |  0.0067114 |
| Bacilli             | Lysinibacillus   |    21 |   1 |     22 |  0.0454545 |
| Bacteroidia         | Bacteroides      |   310 |   1 |    311 |  0.0032154 |
| Bacteroidia         | Parabacteroides  |   113 |   1 |    114 |  0.0087719 |
| Bacteroidia         | Phocaeicola      |   126 |   1 |    127 |  0.0078740 |
| Clostridia          | Monoglobus       |    19 |   1 |     20 |  0.0500000 |
| Coriobacteriia      | Collinsella      |    15 |   1 |     16 |  0.0625000 |
| Cyanophyceae        | Acaryochloris    |    27 |   2 |     29 |  0.0689655 |
| Gammaproteobacteria | Acinetobacter    |  2074 |   1 |   2075 |  0.0004819 |
| Gammaproteobacteria | Escherichia      | 12110 | 163 |  12273 |  0.0132812 |
| Gammaproteobacteria | Klebsiella       |  9851 |   4 |   9855 |  0.0004059 |
| Gammaproteobacteria | Legionella       |    92 |   1 |     93 |  0.0107527 |
| Gammaproteobacteria | Providencia      |   158 |   1 |    159 |  0.0062893 |
| Gammaproteobacteria | Salmonella       |  3034 |   2 |   3036 |  0.0006588 |
| Gammaproteobacteria | Shigella         |   651 |   6 |    657 |  0.0091324 |
| Gammaproteobacteria | Stenotrophomonas |    17 |   2 |     19 |  0.1052632 |

``` r
(plasmid_proportion_filtered <- ggplot(data=metg_counts_genera,
                                       aes(x=TAXONOMY_genus, y=proportion, fill=TAXONOMY_class)) +  
    geom_bar(stat="identity") + 
    scale_fill_manual(values = class_colours, name="") + 
    labs(x = "Genus", 
         y = "Proportion of Filtered Plasmids With ≥1 MetG Homologue", 
         fill = "Class") + 
    ggtitle("The Proportion of Plasmids per Genus Containing ≥1 \nmetG Homologues in BAKTA-annotated PLSDB (filtered)") +
    geom_text(aes(label=paste(Y, " / ", totals, sep=""), y=-0.014), 
              hjust=0.5, size=3.3, angle=60) +
    scale_y_continuous(limits=c(-0.019, 0.12),
                       breaks = seq(0, 1, by = 0.04)) +
    theme(axis.text.x = element_text(face = "italic", angle=60, hjust=1),
          legend.position="bottom"))
```

![](1_PLSDB_metG_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
png("./figs/proportion_plot.png", height=7.2, width=6.2, units="in", res=300)
plasmid_proportion_filtered + ggtitle("")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

How many plasmids in this plot?

``` r
metg_counts_genera %>% ungroup() %>% summarise(metg_encoding = sum(Y))
```

    ## # A tibble: 1 × 1
    ##   metg_encoding
    ##           <int>
    ## 1           192

192 plasmids.

Output a list of the plasmids that fulfil the filtering criteria for
analysis in subsequent scripts.

``` r
(metg_plasmids_uid <- filter(plsdb_subset_genera, encodes_metg == "Y") %>%
  select(NUCCORE_UID) %>% pull()) %>% length()
```

    ## [1] 198

``` r
write.table(filter(plsdb_subset_genera, encodes_metg == "Y"),
            file = "./1_data/plsdb_subset.tsv",
            row.names=FALSE, sep="\t", quote=FALSE, col.names=TRUE)

write.table(metg_plasmids_uid, file = "./1_data/plsdb_subset_uid.txt",
            row.names=FALSE, sep=",", quote=FALSE, col.names=FALSE)
```

Note that this list includes plasmids from plasmid-rare genera that were
not included in the proportions plot above, specifically:

``` r
filter(plsdb_subset_genera, encodes_metg == "Y" &
         !(TAXONOMY_genus %in% metg_counts_genera$TAXONOMY_genus)) %>%
  select(NUCCORE_UID, Accession, TAXONOMY_class, TAXONOMY_genus, TAXONOMY_species, NUCCORE_Length)
```

    ##   NUCCORE_UID     Accession TAXONOMY_class TAXONOMY_genus
    ## 1   319788691   NC_014825.1     Clostridia   Ruminococcus
    ## 2  2431558231 NZ_CP116225.1  Actinomycetes    Actinomyces
    ## 3  2161066756 NZ_CP088021.1        Bacilli  Radiobacillus
    ## 4  2039922910 NZ_AP023416.1     Clostridia     Vescimonas
    ## 5  1563790493    LR214964.1     Mollicutes Mycoplasmopsis
    ## 6  1314784850    MF554642.1           <NA>           <NA>
    ##            TAXONOMY_species NUCCORE_Length
    ## 1        Ruminococcus_albus         352646
    ## 2          Actinomyces_oris         415503
    ## 3  Radiobacillus_kanasensis          66608
    ## 4     Vescimonas_fastidiosa         860917
    ## 5 Mycoplasmopsis_fermentans          79596
    ## 6      uncultured_bacterium         113109

------------------------------------------------------------------------

**[Back to index.](../README.md)**
