UKBB deCODE Protein Analysis
================
Lily Andrews
2023-04-20

This analysis is to see if UK Biobank (UKBB) proteins published on
pre-print
(<https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1.abstract>)
are replicated in deCODE. UKBB published pQTLs which reached a threshold
of 5x10-8/(number of proteins n=1463)=3.3x10-11. Similarily I calculated
the same threshold for deCODE 5x10-8/(number of proteins
n=4907)=1.0x10-11. This generated 1021 SNPs which have been replicated
in the same protein in UKBB and deCODE.

Import data

``` r
library(ggplot2)
ukbb_decode <- read.table(file="merged_ukbb_decode_cleaned.txt", header=T, fill=T)
```

Double check alleles are flipped correctly.

``` r
allele <- data.frame(do.call("rbind", strsplit(as.character(ukbb_decode$MarkerName_ukbb), ":", fixed = TRUE)))
names(allele)[3] <- "other_allele_ukbb"
names(allele)[4] <- "effect_allele_ukbb"
ukbb_decode['other_allele_ukbb'] <- allele$other_allele_ukbb
ukbb_decode['effect_allele_ukbb'] <- allele$effect_allele_ukbb
subset(ukbb_decode, effect_allele_decode!=effect_allele_ukbb)
```

    ##            SNP chr_decode position_decode other_allele_decode
    ## 4   rs10993994      chr10        46046326                   A
    ## 48    rs367070      chr19        54296648                   G
    ## 141 rs11668882      chr19        54171328                   C
    ## 201  rs2511241      chr11        73234296                   C
    ## 228  rs6420185       chr8       143220901                   A
    ## 314 rs11668882      chr19        54171328                   C
    ## 487   rs367070      chr19        54296648                   G
    ## 703  rs9728526       chr1       145718303                   C
    ## 812  rs3752404       chr7       142762725                   G
    ## 898 rs10993994      chr10        46046326                   A
    ## 906 rs11868696      chr17        79732608                   G
    ## 996   rs712048      chr17        35999179                   C
    ## 998  rs7255591      chr19        54904812                   C
    ##     effect_allele_decode beta_decode se_decode pval_decode
    ## 4                      G     -0.0521  0.007648    9.64e-12
    ## 48                     A     -0.0691  0.009839    2.17e-12
    ## 141                    T     -0.0632  0.008427    6.39e-14
    ## 201                 TRUE      0.1097  0.013906    3.05e-15
    ## 228                    G      0.0687  0.008505    6.61e-16
    ## 314                    T     -0.0737  0.008393    1.62e-18
    ## 487                    A     -0.1092  0.010148    5.27e-27
    ## 703                    T     -0.1308  0.007926    3.54e-61
    ## 812                    A      0.1954  0.008533   4.76e-116
    ## 898                    G      1.0228  0.004614   1.00e-200
    ## 906                    A     -0.9136  0.005380   1.00e-200
    ## 996                    A     -0.7401  0.011037   1.00e-200
    ## 998                    G      0.4191  0.012759   1.00e-200
    ##     effect_allele_freq_decode samplesize.exposure protein
    ## 4                     0.38221               35328  CRISP2
    ## 48                    0.21919               35367   KITLG
    ## 141                   0.41554               35371     CPM
    ## 201                   0.08753               35366     OSM
    ## 228                   0.42520               35264   SFRP1
    ## 314                   0.41554               35359   SMPD1
    ## 487                   0.21919               35367  LILRB1
    ## 703                   0.46076               35361   CD160
    ## 812                   0.40391               35372   PRSS2
    ## 898                   0.38221               35316    MSMB
    ## 906                   0.39653               35371   ENPP7
    ## 996                   0.09333               35368   CCL23
    ## 998                   0.10814               35371    NCR1
    ##                      exposure_file        MarkerName_ukbb chr_ukbb
    ## 4      9282_12_CRISP2_CRIS2.txt.gz 10:51549496:T:C:imp:v1       10
    ## 48        9377_25_KITLG_SCF.txt.gz 19:54800500:A:G:imp:v1       19
    ## 141        9416_77_CPM_CBPM.txt.gz 19:54675097:T:C:imp:v1       19
    ## 201        14063_17_OSM_OSM.txt.gz 11:72945341:C:T:imp:v1       11
    ## 228    3221_54_SFRP1_SARP_2.txt.gz 8:144303071:G:A:imp:v1        8
    ## 314      10818_36_SMPD1_ASM.txt.gz 19:54675097:T:C:imp:v1       19
    ## 487    5090_49_LILRB1_ILT_2.txt.gz 19:54800500:A:G:imp:v1       19
    ## 703     17677_47_CD160_BY55.txt.gz 1:145716763:G:A:imp:v1        1
    ## 812 5034_79_PRSS2_Trypsin_2.txt.gz 7:142470574:A:G:imp:v1        7
    ## 898    10620_21_MSMB_PSP_94.txt.gz 10:51549496:T:C:imp:v1       10
    ## 906     4435_66_ENPP7_ENPP7.txt.gz 17:77706406:A:G:imp:v1       17
    ## 996  3028_36_CCL23_Ck_b_8_1.txt.gz 17:34326215:A:C:imp:v1       17
    ## 998     8360_169_NCR1_NKp46.txt.gz 19:55416170:G:C:imp:v1       19
    ##     position_ukbb region_id_ukbb region_start_ukbb region_end_ukbb mhc_ukbb
    ## 4        46046326             13          44463408        47046326        0
    ## 48       54296648             15          53296648        55302450        0
    ## 141      54171328           1123          53042277        55173495        0
    ## 201      73234296            521          72234296        74234296        0
    ## 228     143220901            552         142220901       144220901        0
    ## 314      54171328           1123          53042277        55173495        0
    ## 487      54296648             15          53296648        55302450        0
    ## 703     145718303            189         144718303       146718303        0
    ## 812     142762725            796         141752462       143954267        0
    ## 898      46046326             13          44463408        47046326        0
    ## 906      79732608             23          78732608        80732608        0
    ## 996      35999179            612          34766707        37001607        0
    ## 998      54904812            390          53904812        55904812        0
    ##               protein_id_ukbb assay_target_ukbb uniprot_target_ukbb
    ## 4   CRISP2:P16562:OID21456:v1            CRISP2              P16562
    ## 48   KITLG:P21583:OID20196:v1             KITLG              P21583
    ## 141    CPM:P14384:OID21054:v1               CPM              P14384
    ## 201    OSM:P13725:OID20574:v1               OSM              P13725
    ## 228  SFRP1:Q8N474:OID20984:v1             SFRP1              Q8N474
    ## 314  SMPD1:P17405:OID21028:v1             SMPD1              P17405
    ## 487 LILRB1:Q8NHL6:OID20323:v1            LILRB1              Q8NHL6
    ## 703  CD160:O95971:OID20647:v1             CD160              O95971
    ## 812  PRSS2:P07478:OID20364:v1             PRSS2              P07478
    ## 898   MSMB:P08118:OID20275:v1              MSMB              P08118
    ## 906  ENPP7:Q6UWV6:OID20689:v1             ENPP7              Q6UWV6
    ## 996  CCL23:P55773:OID20693:v1             CCL23              P55773
    ## 998   NCR1:O76036:OID20566:v1              NCR1              O76036
    ##     A1_freq_discovery_ukbb beta_discovery_ukbb se_discovery_ukbb
    ## 4                0.6083000          -0.0400485        0.00539774
    ## 48               0.2304170           0.1064750        0.00860863
    ## 141              0.4408190           0.0536416        0.00717886
    ## 201              0.9295160           0.3633080        0.01439080
    ## 228              0.3872050          -0.0641971        0.00737050
    ## 314              0.4409760           0.0579802        0.00682831
    ## 487              0.2304450          -0.6664440        0.00851100
    ## 703              0.4393840          -0.3924850        0.00742716
    ## 812              0.4162090          -0.1373940        0.00743728
    ## 898              0.6086310           0.9287290        0.00715392
    ## 906              0.3664510           0.9425100        0.00786969
    ## 996              0.8692290           0.5619500        0.01093110
    ## 998              0.0888507          -0.4279300        0.01279910
    ##     log10_pval_discovery_ukbb A1_freq_replication_ukbb beta_replication_ukbb
    ## 4                     12.9297                 0.568076            -0.0414217
    ## 48                    34.4119                 0.211336             0.1150300
    ## 141                   13.1030                 0.433821             0.0499719
    ## 201                  139.9000                 0.924484             0.2802670
    ## 228                   17.5173                 0.395985            -0.0684623
    ## 314                   16.6891                 0.434288             0.0544676
    ## 487                 1333.4300                 0.211547            -0.5952830
    ## 703                  608.2140                 0.469836            -0.3574110
    ## 812                   75.4728                 0.420460            -0.1270010
    ## 898                 3661.9100                 0.568265             0.8848000
    ## 906                 3116.8400                 0.405327             0.8956640
    ## 996                  575.6920                 0.793098             0.4643720
    ## 998                  244.3640                 0.105457            -0.3530790
    ##     se_replication_ukbb log10_pval_replication_ukbb cistrans_ukbb cis_gene_ukbb
    ## 4            0.00780281                     6.95667         trans             -
    ## 48           0.01251400                    19.41410         trans             -
    ## 141          0.00953707                     6.79373         trans             -
    ## 201          0.01916360                    47.71050         trans             -
    ## 228          0.01014340                    10.82840         trans             -
    ## 314          0.00986372                     7.47477         trans             -
    ## 487          0.01218220                   520.29000           cis        LILRB1
    ## 703          0.01015730                   270.51000           cis         CD160
    ## 812          0.01050550                    32.91820           cis         PRSS2
    ## 898          0.00983181                  1760.70000           cis          MSMB
    ## 906          0.01098050                  1446.80000           cis         ENPP7
    ## 996          0.01336660                   263.72600           cis         CCL23
    ## 998          0.01609570                   105.93100           cis          NCR1
    ##     bioinformatic_annotated_gene_ukbb ensembl_gene_id_ukbb
    ## 4                                MSMB      ENSG00000263639
    ## 48                         AC245884.7      ENSG00000251431
    ## 141                              TMC4      ENSG00000167608
    ## 201                             P2RY2      ENSG00000175591
    ## 228                           GPIHBP1      ENSG00000277494
    ## 314                              TMC4      ENSG00000167608
    ## 487                        AC245884.7      ENSG00000251431
    ## 703                             CD160      ENSG00000117281
    ## 812                           PRSS3P1      ENSG00000250591
    ## 898                              MSMB      ENSG00000263639
    ## 906                             ENPP7      ENSG00000182156
    ## 996                       CCL15-CCL14      ENSG00000275688
    ## 998                              NCR1      ENSG00000189430
    ##           annotated_gene_consequence_ukbb            biotype_ukbb
    ## 4                     5_prime_UTR_variant          protein_coding
    ## 48                  upstream_gene_variant  unprocessed_pseudogene
    ## 141                        intron_variant          protein_coding
    ## 201                      missense_variant          protein_coding
    ## 228               downstream_gene_variant          protein_coding
    ## 314                        intron_variant          protein_coding
    ## 487                 upstream_gene_variant  unprocessed_pseudogene
    ## 703                 upstream_gene_variant          protein_coding
    ## 812    non_coding_transcript_exon_variant  unprocessed_pseudogene
    ## 898                   5_prime_UTR_variant          protein_coding
    ## 906                        intron_variant          protein_coding
    ## 996 intron_variant&NMD_transcript_variant nonsense_mediated_decay
    ## 998                 upstream_gene_variant          protein_coding
    ##     cadd_phred_ukbb   sift_ukbb polyphen_ukbb phast_phylop_ukbb
    ## 4             2.259        <NA>          <NA>              <NA>
    ## 48            0.579        <NA>          <NA>              <NA>
    ## 141           2.609        <NA>          <NA>              <NA>
    ## 201          17.560 Deleterious          <NA>               375
    ## 228           2.368        <NA>          <NA>              <NA>
    ## 314           2.609        <NA>          <NA>              <NA>
    ## 487           0.579        <NA>          <NA>              <NA>
    ## 703           1.952        <NA>          <NA>              <NA>
    ## 812          10.310        <NA>          <NA>               375
    ## 898           2.259        <NA>          <NA>              <NA>
    ## 906           5.519        <NA>          <NA>              <NA>
    ## 996           0.743        <NA>          <NA>              <NA>
    ## 998           0.159        <NA>          <NA>              <NA>
    ##     other_allele_ukbb effect_allele_ukbb
    ## 4                   T                  C
    ## 48                  A                  G
    ## 141                 T                  C
    ## 201                 C                  T
    ## 228                 G                  A
    ## 314                 T                  C
    ## 487                 A                  G
    ## 703                 G                  A
    ## 812                 A                  G
    ## 898                 T                  C
    ## 906                 A                  G
    ## 996                 A                  C
    ## 998                 G                  C

Sort out alleles between UKBB and deCODE so they match

``` r
ukbb_decode <- ukbb_decode[-4,] #remove line 4 as alleles weren't matching and couldn't be flipped
ukbb_decode[47,4] <- "A"
ukbb_decode[47,5] <- "G"
ukbb_decode[47,6] <- ukbb_decode[47,6]*-1
ukbb_decode[47,9] <- 1-ukbb_decode[47,9]
ukbb_decode <- ukbb_decode[-140,] #remove line 141 now 140 as alleles weren't matching and couldn't be flipped
ukbb_decode[226,4] <- "G"
ukbb_decode[226,5] <- "A"
ukbb_decode[226,6] <- ukbb_decode[226,6]*-1
ukbb_decode[226,9] <- 1-ukbb_decode[226,9]
ukbb_decode[312,4] <- "T"
ukbb_decode[312,5] <- "C"
ukbb_decode[312,6] <- ukbb_decode[312,6]*-1
ukbb_decode[312,9] <- 1-ukbb_decode[312,9]
ukbb_decode[485,4] <- "A"
ukbb_decode[485,5] <- "G"
ukbb_decode[485,6] <- ukbb_decode[485,6]*-1
ukbb_decode[485,9] <- 1-ukbb_decode[485,9]
ukbb_decode <- ukbb_decode[-701,] #remove line 703 now 701 as alleles weren't matching and couldn't be flipped
ukbb_decode[809,4] <- "A"
ukbb_decode[809,5] <- "G"
ukbb_decode[809,6] <- ukbb_decode[809,6]*-1
ukbb_decode[809,9] <- 1-ukbb_decode[809,9]
ukbb_decode <- ukbb_decode[-895,] #remove line 898 now 895 as alleles weren't matching and couldn't be flipped
ukbb_decode[902,4] <- "A"
ukbb_decode[902,5] <- "G"
ukbb_decode[902,6] <- ukbb_decode[902,6]*-1
ukbb_decode[902,9] <- 1-ukbb_decode[902,9]
ukbb_decode[992,4] <- "A"
ukbb_decode[992,5] <- "C"
ukbb_decode[992,6] <- ukbb_decode[992,6]*-1
ukbb_decode[992,9] <- 1-ukbb_decode[992,9]
ukbb_decode[994,4] <- "G"
ukbb_decode[994,5] <- "C"
ukbb_decode[994,6] <- ukbb_decode[994,6]*-1
ukbb_decode[994,9] <- 1-ukbb_decode[994,9]
ukbb_decode$other_allele_decode[ukbb_decode$other_allele_decode=="TRUE"] <- "T"
ukbb_decode$effect_allele_decode[ukbb_decode$effect_allele_decode=="TRUE"] <- "T"
allele <- data.frame(do.call("rbind", strsplit(as.character(ukbb_decode$MarkerName_ukbb), ":", fixed = TRUE)))
names(allele)[3] <- "other_allele_ukbb"
names(allele)[4] <- "effect_allele_ukbb"
ukbb_decode['other_allele_ukbb'] <- allele$other_allele_ukbb
ukbb_decode['effect_allele_ukbb'] <- allele$effect_allele_ukbb
subset(ukbb_decode, effect_allele_decode!=effect_allele_ukbb)
```

    ##  [1] SNP                               chr_decode                       
    ##  [3] position_decode                   other_allele_decode              
    ##  [5] effect_allele_decode              beta_decode                      
    ##  [7] se_decode                         pval_decode                      
    ##  [9] effect_allele_freq_decode         samplesize.exposure              
    ## [11] protein                           exposure_file                    
    ## [13] MarkerName_ukbb                   chr_ukbb                         
    ## [15] position_ukbb                     region_id_ukbb                   
    ## [17] region_start_ukbb                 region_end_ukbb                  
    ## [19] mhc_ukbb                          protein_id_ukbb                  
    ## [21] assay_target_ukbb                 uniprot_target_ukbb              
    ## [23] A1_freq_discovery_ukbb            beta_discovery_ukbb              
    ## [25] se_discovery_ukbb                 log10_pval_discovery_ukbb        
    ## [27] A1_freq_replication_ukbb          beta_replication_ukbb            
    ## [29] se_replication_ukbb               log10_pval_replication_ukbb      
    ## [31] cistrans_ukbb                     cis_gene_ukbb                    
    ## [33] bioinformatic_annotated_gene_ukbb ensembl_gene_id_ukbb             
    ## [35] annotated_gene_consequence_ukbb   biotype_ukbb                     
    ## [37] cadd_phred_ukbb                   sift_ukbb                        
    ## [39] polyphen_ukbb                     phast_phylop_ukbb                
    ## [41] other_allele_ukbb                 effect_allele_ukbb               
    ## <0 rows> (or 0-length row.names)

Flipped alleles, beta and allele frequency and deleted rows with
mismatching alleles e.g. CG in deCODE and AT in UKBB for same SNP.

Number of pQTLs which meet p-value threshold fot both UKBB and deCODE

``` r
num_pqtl <- nrow(ukbb_decode)
print(num_pqtl)
```

    ## [1] 1017

Compare betas of UKBB SNPs and deCODE SNPs

``` r
ggplot(aes(x=beta_decode, y=beta_discovery_ukbb), data=ukbb_decode)+geom_point() #create plot to compare betas
```

![](ukbb_decode_protein_analysis_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Log transform deCODE p-values to plot against UKBB.

``` r
ukbb_decode['log10_pval_decode'] <- -log10(ukbb_decode$pval_decode)
ggplot(aes(x=log10_pval_decode, y=log10_pval_discovery_ukbb), data=ukbb_decode)+geom_point() #create plot to compare p-values
```

![](ukbb_decode_protein_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Deeper look into the protein data to see which proteins are associated
with pQTLs found in bothe UKBB and deCODE.

``` r
prot_ukbb_decode <- unique(ukbb_decode$protein)
print(prot_ukbb_decode)
```

    ##   [1] "LGALS9"    "EGLN1"     "CR2"       "SCG3"      "ADAMTS13"  "CTSZ"     
    ##   [7] "MDGA1"     "KIT"       "MFGE8"     "USO1"      "PLAU"      "SMPD1"    
    ##  [13] "L1CAM"     "CPM"       "SELP"      "LPO"       "FCN2"      "GSTA1"    
    ##  [19] "FAM3B"     "CD55"      "FLT3LG"    "LILRA5"    "NOTCH3"    "PPY"      
    ##  [25] "PTGDS"     "IL5RA"     "SERPINE1"  "SFRP1"     "RTN4R"     "HYAL1"    
    ##  [31] "STX8"      "TINAGL1"   "SHMT1"     "CCL23"     "ICAM4"     "IL3RA"    
    ##  [37] "ANG"       "WIF1"      "CD58"      "CD163"     "ST3GAL1"   "DSG2"     
    ##  [43] "AMIGO2"    "KITLG"     "PLXNB2"    "WFDC2"     "IL1RL1"    "IL1R1"    
    ##  [49] "CST6"      "CTSO"      "PRSS2"     "SMAD1"     "ART3"      "NCAM1"    
    ##  [55] "TNFSF13B"  "GDF2"      "TNFRSF11B" "CLPS"      "CLEC4C"    "ARHGAP25" 
    ##  [61] "ITIH3"     "CDH5"      "MANF"      "TFRC"      "HMOX1"     "IL1R2"    
    ##  [67] "SPP1"      "CD14"      "CDKN2D"    "ROBO2"     "PDGFB"     "KLK13"    
    ##  [73] "CD59"      "FCRLB"     "IGFBP2"    "IGFBP6"    "GNLY"      "CLEC11A"  
    ##  [79] "HBEGF"     "C2"        "SOST"      "FLRT2"     "NELL1"     "IGF2R"    
    ##  [85] "LRP11"     "VEGFA"     "PRSS8"     "CREG1"     "IDUA"      "GRN"      
    ##  [91] "DKK3"      "CA6"       "CHI3L1"    "MATN2"     "CSF2RA"    "TIMD4"    
    ##  [97] "ENTPD5"    "CRISP2"    "LEPR"      "FKBP4"     "CCL21"     "OLR1"     
    ## [103] "FGFBP1"    "APEX1"     "CD177"     "LY75"      "VNN2"      "HAVCR2"   
    ## [109] "CD160"     "SPON1"     "WFIKKN1"   "PRTN3"     "MMP9"      "MAD1L1"   
    ## [115] "PLTP"      "TNFRSF1B"  "TCN2"      "LILRB1"    "CDH17"     "CLEC1B"   
    ## [121] "TNFRSF19"  "TBCB"      "QPCT"      "MSMB"      "CXCL9"     "HPGDS"    
    ## [127] "IL18R1"    "NCAM2"     "MMP12"     "DUSP3"     "NTRK2"     "SFTPD"    
    ## [133] "ANGPT1"    "SERPINA11" "PLA2G7"    "ANGPTL1"   "ROR1"      "VASN"     
    ## [139] "TPSAB1"    "MANSC1"    "DOK2"      "PEAR1"     "DNER"      "AXL"      
    ## [145] "PRCP"      "TNXB"      "INHBC"     "OSM"       "THPO"      "LCN2"     
    ## [151] "S100A4"    "CTSS"      "KRT5"      "HS6ST1"    "ADH4"      "FGF19"    
    ## [157] "MIF"       "DPY30"     "NPPB"      "ICAM1"     "LYPD3"     "TNFSF12"  
    ## [163] "NADK"      "CXCL10"    "CD244"     "EIF4B"     "LAG3"      "PDGFA"    
    ## [169] "TGFBR3"    "PLXDC1"    "LEP"       "RHOC"      "GGH"       "CRTAM"    
    ## [175] "TPP1"      "FKBP1B"    "NPTX1"     "USP8"      "DBNL"      "APOM"     
    ## [181] "HDGF"      "CCL18"     "IL1RAP"    "TNR"       "WARS"      "TFF2"     
    ## [187] "TNC"       "UXS1"      "SELE"      "NTF3"      "MIA"       "SIGLEC5"  
    ## [193] "ADAM22"    "TREML2"    "IL12B"     "BMP4"      "FBP1"      "IGFBP1"   
    ## [199] "PSME2"     "CD46"      "CD93"      "LGALS3"    "FHIT"      "CST5"     
    ## [205] "TFPI2"     "BST1"      "NCR1"      "ADGRE2"    "CBLN4"     "SCARA5"   
    ## [211] "CHRDL1"    "B4GALT1"   "EPHB4"     "EPHB6"     "AGR3"      "ANGPTL3"  
    ## [217] "LILRB5"    "IGF1R"     "SEMA4C"    "FLT1"      "CDON"      "REG1B"    
    ## [223] "CHL1"      "PLAUR"     "SEMA7A"    "ACP6"      "ENO2"      "FCGR3B"   
    ## [229] "DAPP1"     "CD109"     "THBS2"     "FUT8"      "COMP"      "EDIL3"    
    ## [235] "RRM2B"     "ARHGAP1"   "SMOC2"     "BOC"       "KYNU"      "GPNMB"    
    ## [241] "STX4"      "GHRL"      "ASGR1"     "CNTNAP2"   "ITGA5"     "HAVCR1"   
    ## [247] "HSPB1"     "PLA2G1B"   "CTSF"      "IGSF3"     "CPXM1"     "PDGFRA"   
    ## [253] "F9"        "CCDC80"    "TFF1"      "CNTN4"     "ADAM23"    "NUDC"     
    ## [259] "OMG"       "CHGB"      "CCL8"      "ALDH1A1"   "CTSB"      "SNAP29"   
    ## [265] "SCARF2"    "FCRL3"     "ARG1"      "STAMBP"    "PDCD6"     "BID"      
    ## [271] "SMOC1"     "COL1A1"    "REG1A"     "DCTPP1"    "TGFBI"     "SLIT2"    
    ## [277] "ICAM5"     "FETUB"     "FRZB"      "TXNDC15"   "ICAM2"     "SPARCL1"  
    ## [283] "GLB1"      "COL6A3"    "SPINK5"    "XRCC4"     "ANXA5"     "REG4"     
    ## [289] "FABP1"     "NFASC"     "ERP44"     "CLEC4G"    "DECR1"     "TYRO3"    
    ## [295] "SMPDL3A"   "SELPLG"    "PGLYRP1"   "LILRA2"    "ISLR2"     "HPCAL1"   
    ## [301] "KLB"       "SIAE"      "CCS"       "APLP1"     "SPINK1"    "HGF"      
    ## [307] "DLK1"      "FLT4"      "ENTPD6"    "FCRL6"     "GPC1"      "CAPG"     
    ## [313] "CDCP1"     "ANGPT2"    "PAPPA"     "TXNDC5"    "IDI2"      "BTC"      
    ## [319] "FAS"       "ALDH3A1"   "AGER"      "PDCD1LG2"  "MUC16"     "AIF1"     
    ## [325] "MMP10"     "IGSF8"     "CRNN"      "CRHBP"     "EPHA2"     "CTSD"     
    ## [331] "EFEMP1"    "FGR"       "ARSB"      "TFF3"      "DPEP1"     "PILRA"    
    ## [337] "DPEP2"     "SIGLEC1"   "NID1"      "DSC2"      "BMP6"      "PCOLCE"   
    ## [343] "SERPINA9"  "FIS1"      "IL1RL2"    "RGMA"      "SULT2A1"   "IDS"      
    ## [349] "GAS6"      "TNFAIP8"   "TFPI"      "MAX"       "PMVK"      "SERPINB8" 
    ## [355] "CTSH"      "IFNLR1"    "AOC1"      "BCAN"      "BLVRB"     "CCL5"     
    ## [361] "CLMP"      "MMP8"      "RSPO3"     "PLAT"      "GDF15"     "ESAM"     
    ## [367] "CRTAC1"    "SORD"      "DBI"       "CD274"     "NUDT2"     "NT5E"     
    ## [373] "SPON2"     "PDCD5"     "KIR2DL3"   "ENPP7"     "CXCL6"     "FABP2"    
    ## [379] "SERPINA12" "MYOC"      "PI3"       "APOH"      "TIMP4"     "UMOD"     
    ## [385] "HEBP1"     "IL15RA"    "TLR3"      "LHB"       "SIRPA"     "RARRES1"  
    ## [391] "OGN"       "SIRPB1"    "CHIT1"     "OXT"       "ENPP5"

Next, looking at the number of SNPs found in both UKBB and deCODE for
each protein.

``` r
df <- data.frame(matrix(ncol = 2, nrow = 0))
for (protein in prot_ukbb_decode){
  total <- sum(ukbb_decode$protein == protein)
  print(paste("Protein", protein, total, "SNP(s)"))
  df_temp = c(protein, total)
  df = rbind(df, df_temp)
  df
}
```

    ## [1] "Protein LGALS9 2 SNP(s)"
    ## [1] "Protein EGLN1 3 SNP(s)"
    ## [1] "Protein CR2 9 SNP(s)"
    ## [1] "Protein SCG3 2 SNP(s)"
    ## [1] "Protein ADAMTS13 4 SNP(s)"
    ## [1] "Protein CTSZ 5 SNP(s)"
    ## [1] "Protein MDGA1 3 SNP(s)"
    ## [1] "Protein KIT 8 SNP(s)"
    ## [1] "Protein MFGE8 5 SNP(s)"
    ## [1] "Protein USO1 2 SNP(s)"
    ## [1] "Protein PLAU 5 SNP(s)"
    ## [1] "Protein SMPD1 7 SNP(s)"
    ## [1] "Protein L1CAM 2 SNP(s)"
    ## [1] "Protein CPM 4 SNP(s)"
    ## [1] "Protein SELP 2 SNP(s)"
    ## [1] "Protein LPO 11 SNP(s)"
    ## [1] "Protein FCN2 5 SNP(s)"
    ## [1] "Protein GSTA1 2 SNP(s)"
    ## [1] "Protein FAM3B 7 SNP(s)"
    ## [1] "Protein CD55 3 SNP(s)"
    ## [1] "Protein FLT3LG 1 SNP(s)"
    ## [1] "Protein LILRA5 9 SNP(s)"
    ## [1] "Protein NOTCH3 3 SNP(s)"
    ## [1] "Protein PPY 6 SNP(s)"
    ## [1] "Protein PTGDS 1 SNP(s)"
    ## [1] "Protein IL5RA 9 SNP(s)"
    ## [1] "Protein SERPINE1 3 SNP(s)"
    ## [1] "Protein SFRP1 3 SNP(s)"
    ## [1] "Protein RTN4R 2 SNP(s)"
    ## [1] "Protein HYAL1 3 SNP(s)"
    ## [1] "Protein STX8 1 SNP(s)"
    ## [1] "Protein TINAGL1 3 SNP(s)"
    ## [1] "Protein SHMT1 2 SNP(s)"
    ## [1] "Protein CCL23 3 SNP(s)"
    ## [1] "Protein ICAM4 6 SNP(s)"
    ## [1] "Protein IL3RA 3 SNP(s)"
    ## [1] "Protein ANG 2 SNP(s)"
    ## [1] "Protein WIF1 3 SNP(s)"
    ## [1] "Protein CD58 1 SNP(s)"
    ## [1] "Protein CD163 11 SNP(s)"
    ## [1] "Protein ST3GAL1 4 SNP(s)"
    ## [1] "Protein DSG2 7 SNP(s)"
    ## [1] "Protein AMIGO2 4 SNP(s)"
    ## [1] "Protein KITLG 4 SNP(s)"
    ## [1] "Protein PLXNB2 4 SNP(s)"
    ## [1] "Protein WFDC2 2 SNP(s)"
    ## [1] "Protein IL1RL1 4 SNP(s)"
    ## [1] "Protein IL1R1 8 SNP(s)"
    ## [1] "Protein CST6 5 SNP(s)"
    ## [1] "Protein CTSO 4 SNP(s)"
    ## [1] "Protein PRSS2 6 SNP(s)"
    ## [1] "Protein SMAD1 2 SNP(s)"
    ## [1] "Protein ART3 2 SNP(s)"
    ## [1] "Protein NCAM1 3 SNP(s)"
    ## [1] "Protein TNFSF13B 1 SNP(s)"
    ## [1] "Protein GDF2 4 SNP(s)"
    ## [1] "Protein TNFRSF11B 3 SNP(s)"
    ## [1] "Protein CLPS 3 SNP(s)"
    ## [1] "Protein CLEC4C 9 SNP(s)"
    ## [1] "Protein ARHGAP25 2 SNP(s)"
    ## [1] "Protein ITIH3 2 SNP(s)"
    ## [1] "Protein CDH5 4 SNP(s)"
    ## [1] "Protein MANF 1 SNP(s)"
    ## [1] "Protein TFRC 5 SNP(s)"
    ## [1] "Protein HMOX1 6 SNP(s)"
    ## [1] "Protein IL1R2 6 SNP(s)"
    ## [1] "Protein SPP1 1 SNP(s)"
    ## [1] "Protein CD14 6 SNP(s)"
    ## [1] "Protein CDKN2D 2 SNP(s)"
    ## [1] "Protein ROBO2 2 SNP(s)"
    ## [1] "Protein PDGFB 3 SNP(s)"
    ## [1] "Protein KLK13 6 SNP(s)"
    ## [1] "Protein CD59 3 SNP(s)"
    ## [1] "Protein FCRLB 2 SNP(s)"
    ## [1] "Protein IGFBP2 1 SNP(s)"
    ## [1] "Protein IGFBP6 2 SNP(s)"
    ## [1] "Protein GNLY 6 SNP(s)"
    ## [1] "Protein CLEC11A 5 SNP(s)"
    ## [1] "Protein HBEGF 6 SNP(s)"
    ## [1] "Protein C2 1 SNP(s)"
    ## [1] "Protein SOST 3 SNP(s)"
    ## [1] "Protein FLRT2 5 SNP(s)"
    ## [1] "Protein NELL1 7 SNP(s)"
    ## [1] "Protein IGF2R 5 SNP(s)"
    ## [1] "Protein LRP11 5 SNP(s)"
    ## [1] "Protein VEGFA 5 SNP(s)"
    ## [1] "Protein PRSS8 1 SNP(s)"
    ## [1] "Protein CREG1 3 SNP(s)"
    ## [1] "Protein IDUA 4 SNP(s)"
    ## [1] "Protein GRN 7 SNP(s)"
    ## [1] "Protein DKK3 5 SNP(s)"
    ## [1] "Protein CA6 6 SNP(s)"
    ## [1] "Protein CHI3L1 2 SNP(s)"
    ## [1] "Protein MATN2 1 SNP(s)"
    ## [1] "Protein CSF2RA 1 SNP(s)"
    ## [1] "Protein TIMD4 10 SNP(s)"
    ## [1] "Protein ENTPD5 5 SNP(s)"
    ## [1] "Protein CRISP2 6 SNP(s)"
    ## [1] "Protein LEPR 4 SNP(s)"
    ## [1] "Protein FKBP4 2 SNP(s)"
    ## [1] "Protein CCL21 4 SNP(s)"
    ## [1] "Protein OLR1 1 SNP(s)"
    ## [1] "Protein FGFBP1 1 SNP(s)"
    ## [1] "Protein APEX1 1 SNP(s)"
    ## [1] "Protein CD177 2 SNP(s)"
    ## [1] "Protein LY75 4 SNP(s)"
    ## [1] "Protein VNN2 6 SNP(s)"
    ## [1] "Protein HAVCR2 3 SNP(s)"
    ## [1] "Protein CD160 1 SNP(s)"
    ## [1] "Protein SPON1 4 SNP(s)"
    ## [1] "Protein WFIKKN1 2 SNP(s)"
    ## [1] "Protein PRTN3 5 SNP(s)"
    ## [1] "Protein MMP9 3 SNP(s)"
    ## [1] "Protein MAD1L1 1 SNP(s)"
    ## [1] "Protein PLTP 4 SNP(s)"
    ## [1] "Protein TNFRSF1B 7 SNP(s)"
    ## [1] "Protein TCN2 5 SNP(s)"
    ## [1] "Protein LILRB1 5 SNP(s)"
    ## [1] "Protein CDH17 6 SNP(s)"
    ## [1] "Protein CLEC1B 2 SNP(s)"
    ## [1] "Protein TNFRSF19 2 SNP(s)"
    ## [1] "Protein TBCB 1 SNP(s)"
    ## [1] "Protein QPCT 7 SNP(s)"
    ## [1] "Protein MSMB 1 SNP(s)"
    ## [1] "Protein CXCL9 2 SNP(s)"
    ## [1] "Protein HPGDS 2 SNP(s)"
    ## [1] "Protein IL18R1 3 SNP(s)"
    ## [1] "Protein NCAM2 3 SNP(s)"
    ## [1] "Protein MMP12 3 SNP(s)"
    ## [1] "Protein DUSP3 1 SNP(s)"
    ## [1] "Protein NTRK2 5 SNP(s)"
    ## [1] "Protein SFTPD 5 SNP(s)"
    ## [1] "Protein ANGPT1 2 SNP(s)"
    ## [1] "Protein SERPINA11 3 SNP(s)"
    ## [1] "Protein PLA2G7 5 SNP(s)"
    ## [1] "Protein ANGPTL1 2 SNP(s)"
    ## [1] "Protein ROR1 5 SNP(s)"
    ## [1] "Protein VASN 3 SNP(s)"
    ## [1] "Protein TPSAB1 7 SNP(s)"
    ## [1] "Protein MANSC1 4 SNP(s)"
    ## [1] "Protein DOK2 1 SNP(s)"
    ## [1] "Protein PEAR1 5 SNP(s)"
    ## [1] "Protein DNER 2 SNP(s)"
    ## [1] "Protein AXL 1 SNP(s)"
    ## [1] "Protein PRCP 5 SNP(s)"
    ## [1] "Protein TNXB 3 SNP(s)"
    ## [1] "Protein INHBC 1 SNP(s)"
    ## [1] "Protein OSM 1 SNP(s)"
    ## [1] "Protein THPO 1 SNP(s)"
    ## [1] "Protein LCN2 5 SNP(s)"
    ## [1] "Protein S100A4 1 SNP(s)"
    ## [1] "Protein CTSS 3 SNP(s)"
    ## [1] "Protein KRT5 3 SNP(s)"
    ## [1] "Protein HS6ST1 5 SNP(s)"
    ## [1] "Protein ADH4 2 SNP(s)"
    ## [1] "Protein FGF19 4 SNP(s)"
    ## [1] "Protein MIF 1 SNP(s)"
    ## [1] "Protein DPY30 1 SNP(s)"
    ## [1] "Protein NPPB 3 SNP(s)"
    ## [1] "Protein ICAM1 2 SNP(s)"
    ## [1] "Protein LYPD3 3 SNP(s)"
    ## [1] "Protein TNFSF12 2 SNP(s)"
    ## [1] "Protein NADK 2 SNP(s)"
    ## [1] "Protein CXCL10 3 SNP(s)"
    ## [1] "Protein CD244 2 SNP(s)"
    ## [1] "Protein EIF4B 1 SNP(s)"
    ## [1] "Protein LAG3 7 SNP(s)"
    ## [1] "Protein PDGFA 4 SNP(s)"
    ## [1] "Protein TGFBR3 3 SNP(s)"
    ## [1] "Protein PLXDC1 1 SNP(s)"
    ## [1] "Protein LEP 1 SNP(s)"
    ## [1] "Protein RHOC 1 SNP(s)"
    ## [1] "Protein GGH 4 SNP(s)"
    ## [1] "Protein CRTAM 5 SNP(s)"
    ## [1] "Protein TPP1 3 SNP(s)"
    ## [1] "Protein FKBP1B 2 SNP(s)"
    ## [1] "Protein NPTX1 5 SNP(s)"
    ## [1] "Protein USP8 3 SNP(s)"
    ## [1] "Protein DBNL 1 SNP(s)"
    ## [1] "Protein APOM 3 SNP(s)"
    ## [1] "Protein HDGF 2 SNP(s)"
    ## [1] "Protein CCL18 3 SNP(s)"
    ## [1] "Protein IL1RAP 2 SNP(s)"
    ## [1] "Protein TNR 2 SNP(s)"
    ## [1] "Protein WARS 2 SNP(s)"
    ## [1] "Protein TFF2 4 SNP(s)"
    ## [1] "Protein TNC 3 SNP(s)"
    ## [1] "Protein UXS1 2 SNP(s)"
    ## [1] "Protein SELE 6 SNP(s)"
    ## [1] "Protein NTF3 2 SNP(s)"
    ## [1] "Protein MIA 2 SNP(s)"
    ## [1] "Protein SIGLEC5 2 SNP(s)"
    ## [1] "Protein ADAM22 2 SNP(s)"
    ## [1] "Protein TREML2 1 SNP(s)"
    ## [1] "Protein IL12B 5 SNP(s)"
    ## [1] "Protein BMP4 2 SNP(s)"
    ## [1] "Protein FBP1 2 SNP(s)"
    ## [1] "Protein IGFBP1 2 SNP(s)"
    ## [1] "Protein PSME2 1 SNP(s)"
    ## [1] "Protein CD46 3 SNP(s)"
    ## [1] "Protein CD93 1 SNP(s)"
    ## [1] "Protein LGALS3 3 SNP(s)"
    ## [1] "Protein FHIT 1 SNP(s)"
    ## [1] "Protein CST5 3 SNP(s)"
    ## [1] "Protein TFPI2 1 SNP(s)"
    ## [1] "Protein BST1 4 SNP(s)"
    ## [1] "Protein NCR1 3 SNP(s)"
    ## [1] "Protein ADGRE2 3 SNP(s)"
    ## [1] "Protein CBLN4 4 SNP(s)"
    ## [1] "Protein SCARA5 1 SNP(s)"
    ## [1] "Protein CHRDL1 1 SNP(s)"
    ## [1] "Protein B4GALT1 3 SNP(s)"
    ## [1] "Protein EPHB4 2 SNP(s)"
    ## [1] "Protein EPHB6 1 SNP(s)"
    ## [1] "Protein AGR3 2 SNP(s)"
    ## [1] "Protein ANGPTL3 1 SNP(s)"
    ## [1] "Protein LILRB5 2 SNP(s)"
    ## [1] "Protein IGF1R 2 SNP(s)"
    ## [1] "Protein SEMA4C 2 SNP(s)"
    ## [1] "Protein FLT1 2 SNP(s)"
    ## [1] "Protein CDON 2 SNP(s)"
    ## [1] "Protein REG1B 4 SNP(s)"
    ## [1] "Protein CHL1 4 SNP(s)"
    ## [1] "Protein PLAUR 5 SNP(s)"
    ## [1] "Protein SEMA7A 5 SNP(s)"
    ## [1] "Protein ACP6 2 SNP(s)"
    ## [1] "Protein ENO2 2 SNP(s)"
    ## [1] "Protein FCGR3B 5 SNP(s)"
    ## [1] "Protein DAPP1 1 SNP(s)"
    ## [1] "Protein CD109 2 SNP(s)"
    ## [1] "Protein THBS2 1 SNP(s)"
    ## [1] "Protein FUT8 2 SNP(s)"
    ## [1] "Protein COMP 4 SNP(s)"
    ## [1] "Protein EDIL3 2 SNP(s)"
    ## [1] "Protein RRM2B 1 SNP(s)"
    ## [1] "Protein ARHGAP1 2 SNP(s)"
    ## [1] "Protein SMOC2 2 SNP(s)"
    ## [1] "Protein BOC 3 SNP(s)"
    ## [1] "Protein KYNU 1 SNP(s)"
    ## [1] "Protein GPNMB 2 SNP(s)"
    ## [1] "Protein STX4 1 SNP(s)"
    ## [1] "Protein GHRL 1 SNP(s)"
    ## [1] "Protein ASGR1 3 SNP(s)"
    ## [1] "Protein CNTNAP2 2 SNP(s)"
    ## [1] "Protein ITGA5 1 SNP(s)"
    ## [1] "Protein HAVCR1 3 SNP(s)"
    ## [1] "Protein HSPB1 1 SNP(s)"
    ## [1] "Protein PLA2G1B 1 SNP(s)"
    ## [1] "Protein CTSF 4 SNP(s)"
    ## [1] "Protein IGSF3 2 SNP(s)"
    ## [1] "Protein CPXM1 2 SNP(s)"
    ## [1] "Protein PDGFRA 1 SNP(s)"
    ## [1] "Protein F9 1 SNP(s)"
    ## [1] "Protein CCDC80 3 SNP(s)"
    ## [1] "Protein TFF1 3 SNP(s)"
    ## [1] "Protein CNTN4 3 SNP(s)"
    ## [1] "Protein ADAM23 3 SNP(s)"
    ## [1] "Protein NUDC 1 SNP(s)"
    ## [1] "Protein OMG 2 SNP(s)"
    ## [1] "Protein CHGB 1 SNP(s)"
    ## [1] "Protein CCL8 2 SNP(s)"
    ## [1] "Protein ALDH1A1 1 SNP(s)"
    ## [1] "Protein CTSB 2 SNP(s)"
    ## [1] "Protein SNAP29 1 SNP(s)"
    ## [1] "Protein SCARF2 2 SNP(s)"
    ## [1] "Protein FCRL3 2 SNP(s)"
    ## [1] "Protein ARG1 1 SNP(s)"
    ## [1] "Protein STAMBP 1 SNP(s)"
    ## [1] "Protein PDCD6 1 SNP(s)"
    ## [1] "Protein BID 1 SNP(s)"
    ## [1] "Protein SMOC1 2 SNP(s)"
    ## [1] "Protein COL1A1 3 SNP(s)"
    ## [1] "Protein REG1A 4 SNP(s)"
    ## [1] "Protein DCTPP1 1 SNP(s)"
    ## [1] "Protein TGFBI 2 SNP(s)"
    ## [1] "Protein SLIT2 1 SNP(s)"
    ## [1] "Protein ICAM5 6 SNP(s)"
    ## [1] "Protein FETUB 3 SNP(s)"
    ## [1] "Protein FRZB 1 SNP(s)"
    ## [1] "Protein TXNDC15 3 SNP(s)"
    ## [1] "Protein ICAM2 1 SNP(s)"
    ## [1] "Protein SPARCL1 3 SNP(s)"
    ## [1] "Protein GLB1 1 SNP(s)"
    ## [1] "Protein COL6A3 1 SNP(s)"
    ## [1] "Protein SPINK5 1 SNP(s)"
    ## [1] "Protein XRCC4 1 SNP(s)"
    ## [1] "Protein ANXA5 1 SNP(s)"
    ## [1] "Protein REG4 1 SNP(s)"
    ## [1] "Protein FABP1 1 SNP(s)"
    ## [1] "Protein NFASC 2 SNP(s)"
    ## [1] "Protein ERP44 1 SNP(s)"
    ## [1] "Protein CLEC4G 1 SNP(s)"
    ## [1] "Protein DECR1 1 SNP(s)"
    ## [1] "Protein TYRO3 1 SNP(s)"
    ## [1] "Protein SMPDL3A 2 SNP(s)"
    ## [1] "Protein SELPLG 2 SNP(s)"
    ## [1] "Protein PGLYRP1 2 SNP(s)"
    ## [1] "Protein LILRA2 1 SNP(s)"
    ## [1] "Protein ISLR2 3 SNP(s)"
    ## [1] "Protein HPCAL1 1 SNP(s)"
    ## [1] "Protein KLB 2 SNP(s)"
    ## [1] "Protein SIAE 1 SNP(s)"
    ## [1] "Protein CCS 1 SNP(s)"
    ## [1] "Protein APLP1 1 SNP(s)"
    ## [1] "Protein SPINK1 1 SNP(s)"
    ## [1] "Protein HGF 1 SNP(s)"
    ## [1] "Protein DLK1 2 SNP(s)"
    ## [1] "Protein FLT4 3 SNP(s)"
    ## [1] "Protein ENTPD6 1 SNP(s)"
    ## [1] "Protein FCRL6 1 SNP(s)"
    ## [1] "Protein GPC1 2 SNP(s)"
    ## [1] "Protein CAPG 1 SNP(s)"
    ## [1] "Protein CDCP1 3 SNP(s)"
    ## [1] "Protein ANGPT2 2 SNP(s)"
    ## [1] "Protein PAPPA 2 SNP(s)"
    ## [1] "Protein TXNDC5 1 SNP(s)"
    ## [1] "Protein IDI2 1 SNP(s)"
    ## [1] "Protein BTC 1 SNP(s)"
    ## [1] "Protein FAS 2 SNP(s)"
    ## [1] "Protein ALDH3A1 1 SNP(s)"
    ## [1] "Protein AGER 3 SNP(s)"
    ## [1] "Protein PDCD1LG2 2 SNP(s)"
    ## [1] "Protein MUC16 2 SNP(s)"
    ## [1] "Protein AIF1 1 SNP(s)"
    ## [1] "Protein MMP10 2 SNP(s)"
    ## [1] "Protein IGSF8 1 SNP(s)"
    ## [1] "Protein CRNN 1 SNP(s)"
    ## [1] "Protein CRHBP 2 SNP(s)"
    ## [1] "Protein EPHA2 1 SNP(s)"
    ## [1] "Protein CTSD 2 SNP(s)"
    ## [1] "Protein EFEMP1 1 SNP(s)"
    ## [1] "Protein FGR 1 SNP(s)"
    ## [1] "Protein ARSB 1 SNP(s)"
    ## [1] "Protein TFF3 2 SNP(s)"
    ## [1] "Protein DPEP1 2 SNP(s)"
    ## [1] "Protein PILRA 2 SNP(s)"
    ## [1] "Protein DPEP2 2 SNP(s)"
    ## [1] "Protein SIGLEC1 1 SNP(s)"
    ## [1] "Protein NID1 2 SNP(s)"
    ## [1] "Protein DSC2 1 SNP(s)"
    ## [1] "Protein BMP6 1 SNP(s)"
    ## [1] "Protein PCOLCE 1 SNP(s)"
    ## [1] "Protein SERPINA9 1 SNP(s)"
    ## [1] "Protein FIS1 1 SNP(s)"
    ## [1] "Protein IL1RL2 1 SNP(s)"
    ## [1] "Protein RGMA 2 SNP(s)"
    ## [1] "Protein SULT2A1 1 SNP(s)"
    ## [1] "Protein IDS 1 SNP(s)"
    ## [1] "Protein GAS6 3 SNP(s)"
    ## [1] "Protein TNFAIP8 1 SNP(s)"
    ## [1] "Protein TFPI 1 SNP(s)"
    ## [1] "Protein MAX 1 SNP(s)"
    ## [1] "Protein PMVK 1 SNP(s)"
    ## [1] "Protein SERPINB8 1 SNP(s)"
    ## [1] "Protein CTSH 1 SNP(s)"
    ## [1] "Protein IFNLR1 1 SNP(s)"
    ## [1] "Protein AOC1 1 SNP(s)"
    ## [1] "Protein BCAN 1 SNP(s)"
    ## [1] "Protein BLVRB 1 SNP(s)"
    ## [1] "Protein CCL5 1 SNP(s)"
    ## [1] "Protein CLMP 1 SNP(s)"
    ## [1] "Protein MMP8 1 SNP(s)"
    ## [1] "Protein RSPO3 2 SNP(s)"
    ## [1] "Protein PLAT 1 SNP(s)"
    ## [1] "Protein GDF15 1 SNP(s)"
    ## [1] "Protein ESAM 1 SNP(s)"
    ## [1] "Protein CRTAC1 2 SNP(s)"
    ## [1] "Protein SORD 1 SNP(s)"
    ## [1] "Protein DBI 1 SNP(s)"
    ## [1] "Protein CD274 1 SNP(s)"
    ## [1] "Protein NUDT2 1 SNP(s)"
    ## [1] "Protein NT5E 1 SNP(s)"
    ## [1] "Protein SPON2 1 SNP(s)"
    ## [1] "Protein PDCD5 1 SNP(s)"
    ## [1] "Protein KIR2DL3 2 SNP(s)"
    ## [1] "Protein ENPP7 1 SNP(s)"
    ## [1] "Protein CXCL6 1 SNP(s)"
    ## [1] "Protein FABP2 1 SNP(s)"
    ## [1] "Protein SERPINA12 1 SNP(s)"
    ## [1] "Protein MYOC 1 SNP(s)"
    ## [1] "Protein PI3 1 SNP(s)"
    ## [1] "Protein APOH 1 SNP(s)"
    ## [1] "Protein TIMP4 1 SNP(s)"
    ## [1] "Protein UMOD 1 SNP(s)"
    ## [1] "Protein HEBP1 1 SNP(s)"
    ## [1] "Protein IL15RA 1 SNP(s)"
    ## [1] "Protein TLR3 1 SNP(s)"
    ## [1] "Protein LHB 1 SNP(s)"
    ## [1] "Protein SIRPA 1 SNP(s)"
    ## [1] "Protein RARRES1 1 SNP(s)"
    ## [1] "Protein OGN 1 SNP(s)"
    ## [1] "Protein SIRPB1 1 SNP(s)"
    ## [1] "Protein CHIT1 1 SNP(s)"
    ## [1] "Protein OXT 1 SNP(s)"
    ## [1] "Protein ENPP5 1 SNP(s)"

``` r
column_title <- c("protein_name", "n_snps")
colnames(df) <- column_title
```

Order the number of SNPs from largest to smallest.

``` r
prot_order <- df[order(df$n_snps, decreasing=T), ]
print(prot_order)
```

    ##     protein_name n_snps
    ## 3            CR2      9
    ## 22        LILRA5      9
    ## 26         IL5RA      9
    ## 59        CLEC4C      9
    ## 8            KIT      8
    ## 48         IL1R1      8
    ## 12         SMPD1      7
    ## 19         FAM3B      7
    ## 42          DSG2      7
    ## 83         NELL1      7
    ## 90           GRN      7
    ## 116     TNFRSF1B      7
    ## 123         QPCT      7
    ## 139       TPSAB1      7
    ## 167         LAG3      7
    ## 24           PPY      6
    ## 35         ICAM4      6
    ## 51         PRSS2      6
    ## 65         HMOX1      6
    ## 66         IL1R2      6
    ## 68          CD14      6
    ## 72         KLK13      6
    ## 77          GNLY      6
    ## 79         HBEGF      6
    ## 92           CA6      6
    ## 98        CRISP2      6
    ## 107         VNN2      6
    ## 119        CDH17      6
    ## 189         SELE      6
    ## 277        ICAM5      6
    ## 6           CTSZ      5
    ## 9          MFGE8      5
    ## 11          PLAU      5
    ## 17          FCN2      5
    ## 49          CST6      5
    ## 64          TFRC      5
    ## 78       CLEC11A      5
    ## 82         FLRT2      5
    ## 84         IGF2R      5
    ## 85         LRP11      5
    ## 86         VEGFA      5
    ## 91          DKK3      5
    ## 97        ENTPD5      5
    ## 112        PRTN3      5
    ## 117         TCN2      5
    ## 118       LILRB1      5
    ## 131        NTRK2      5
    ## 132        SFTPD      5
    ## 135       PLA2G7      5
    ## 137         ROR1      5
    ## 142        PEAR1      5
    ## 145         PRCP      5
    ## 150         LCN2      5
    ## 154       HS6ST1      5
    ## 174        CRTAM      5
    ## 177        NPTX1      5
    ## 195        IL12B      5
    ## 224        PLAUR      5
    ## 225       SEMA7A      5
    ## 228       FCGR3B      5
    ## 5       ADAMTS13      4
    ## 14           CPM      4
    ## 41       ST3GAL1      4
    ## 43        AMIGO2      4
    ## 44         KITLG      4
    ## 45        PLXNB2      4
    ## 47        IL1RL1      4
    ## 50          CTSO      4
    ## 56          GDF2      4
    ## 62          CDH5      4
    ## 89          IDUA      4
    ## 99          LEPR      4
    ## 101        CCL21      4
    ## 106         LY75      4
    ## 110        SPON1      4
    ## 115         PLTP      4
    ## 140       MANSC1      4
    ## 156        FGF19      4
    ## 168        PDGFA      4
    ## 173          GGH      4
    ## 186         TFF2      4
    ## 206         BST1      4
    ## 209        CBLN4      4
    ## 222        REG1B      4
    ## 223         CHL1      4
    ## 233         COMP      4
    ## 249         CTSF      4
    ## 273        REG1A      4
    ## 2          EGLN1      3
    ## 7          MDGA1      3
    ## 20          CD55      3
    ## 23        NOTCH3      3
    ## 27      SERPINE1      3
    ## 28         SFRP1      3
    ## 30         HYAL1      3
    ## 32       TINAGL1      3
    ## 34         CCL23      3
    ## 36         IL3RA      3
    ## 38          WIF1      3
    ## 54         NCAM1      3
    ## 57     TNFRSF11B      3
    ## 58          CLPS      3
    ## 71         PDGFB      3
    ## 73          CD59      3
    ## 81          SOST      3
    ## 88         CREG1      3
    ## 108       HAVCR2      3
    ## 113         MMP9      3
    ## 127       IL18R1      3
    ## 128        NCAM2      3
    ## 129        MMP12      3
    ## 134    SERPINA11      3
    ## 138         VASN      3
    ## 146         TNXB      3
    ## 152         CTSS      3
    ## 153         KRT5      3
    ## 159         NPPB      3
    ## 161        LYPD3      3
    ## 164       CXCL10      3
    ## 169       TGFBR3      3
    ## 175         TPP1      3
    ## 178         USP8      3
    ## 180         APOM      3
    ## 182        CCL18      3
    ## 187          TNC      3
    ## 200         CD46      3
    ## 202       LGALS3      3
    ## 204         CST5      3
    ## 207         NCR1      3
    ## 208       ADGRE2      3
    ## 212      B4GALT1      3
    ## 238          BOC      3
    ## 243        ASGR1      3
    ## 246       HAVCR1      3
    ## 254       CCDC80      3
    ## 255         TFF1      3
    ## 256        CNTN4      3
    ## 257       ADAM23      3
    ## 272       COL1A1      3
    ## 278        FETUB      3
    ## 280      TXNDC15      3
    ## 282      SPARCL1      3
    ## 299        ISLR2      3
    ## 308         FLT4      3
    ## 313        CDCP1      3
    ## 321         AGER      3
    ## 349         GAS6      3
    ## 1         LGALS9      2
    ## 4           SCG3      2
    ## 10          USO1      2
    ## 13         L1CAM      2
    ## 15          SELP      2
    ## 18         GSTA1      2
    ## 29         RTN4R      2
    ## 33         SHMT1      2
    ## 37           ANG      2
    ## 46         WFDC2      2
    ## 52         SMAD1      2
    ## 53          ART3      2
    ## 60      ARHGAP25      2
    ## 61         ITIH3      2
    ## 69        CDKN2D      2
    ## 70         ROBO2      2
    ## 74         FCRLB      2
    ## 76        IGFBP6      2
    ## 93        CHI3L1      2
    ## 100        FKBP4      2
    ## 105        CD177      2
    ## 111      WFIKKN1      2
    ## 120       CLEC1B      2
    ## 121     TNFRSF19      2
    ## 125        CXCL9      2
    ## 126        HPGDS      2
    ## 133       ANGPT1      2
    ## 136      ANGPTL1      2
    ## 143         DNER      2
    ## 155         ADH4      2
    ## 160        ICAM1      2
    ## 162      TNFSF12      2
    ## 163         NADK      2
    ## 165        CD244      2
    ## 176       FKBP1B      2
    ## 181         HDGF      2
    ## 183       IL1RAP      2
    ## 184          TNR      2
    ## 185         WARS      2
    ## 188         UXS1      2
    ## 190         NTF3      2
    ## 191          MIA      2
    ## 192      SIGLEC5      2
    ## 193       ADAM22      2
    ## 196         BMP4      2
    ## 197         FBP1      2
    ## 198       IGFBP1      2
    ## 213        EPHB4      2
    ## 215         AGR3      2
    ## 217       LILRB5      2
    ## 218        IGF1R      2
    ## 219       SEMA4C      2
    ## 220         FLT1      2
    ## 221         CDON      2
    ## 226         ACP6      2
    ## 227         ENO2      2
    ## 230        CD109      2
    ## 232         FUT8      2
    ## 234        EDIL3      2
    ## 236      ARHGAP1      2
    ## 237        SMOC2      2
    ## 240        GPNMB      2
    ## 244      CNTNAP2      2
    ## 250        IGSF3      2
    ## 251        CPXM1      2
    ## 259          OMG      2
    ## 261         CCL8      2
    ## 263         CTSB      2
    ## 265       SCARF2      2
    ## 266        FCRL3      2
    ## 271        SMOC1      2
    ## 275        TGFBI      2
    ## 290        NFASC      2
    ## 295      SMPDL3A      2
    ## 296       SELPLG      2
    ## 297      PGLYRP1      2
    ## 301          KLB      2
    ## 307         DLK1      2
    ## 311         GPC1      2
    ## 314       ANGPT2      2
    ## 315        PAPPA      2
    ## 319          FAS      2
    ## 322     PDCD1LG2      2
    ## 323        MUC16      2
    ## 325        MMP10      2
    ## 328        CRHBP      2
    ## 330         CTSD      2
    ## 334         TFF3      2
    ## 335        DPEP1      2
    ## 336        PILRA      2
    ## 337        DPEP2      2
    ## 339         NID1      2
    ## 346         RGMA      2
    ## 363        RSPO3      2
    ## 367       CRTAC1      2
    ## 375      KIR2DL3      2
    ## 16           LPO     11
    ## 40         CD163     11
    ## 96         TIMD4     10
    ## 21        FLT3LG      1
    ## 25         PTGDS      1
    ## 31          STX8      1
    ## 39          CD58      1
    ## 55      TNFSF13B      1
    ## 63          MANF      1
    ## 67          SPP1      1
    ## 75        IGFBP2      1
    ## 80            C2      1
    ## 87         PRSS8      1
    ## 94         MATN2      1
    ## 95        CSF2RA      1
    ## 102         OLR1      1
    ## 103       FGFBP1      1
    ## 104        APEX1      1
    ## 109        CD160      1
    ## 114       MAD1L1      1
    ## 122         TBCB      1
    ## 124         MSMB      1
    ## 130        DUSP3      1
    ## 141         DOK2      1
    ## 144          AXL      1
    ## 147        INHBC      1
    ## 148          OSM      1
    ## 149         THPO      1
    ## 151       S100A4      1
    ## 157          MIF      1
    ## 158        DPY30      1
    ## 166        EIF4B      1
    ## 170       PLXDC1      1
    ## 171          LEP      1
    ## 172         RHOC      1
    ## 179         DBNL      1
    ## 194       TREML2      1
    ## 199        PSME2      1
    ## 201         CD93      1
    ## 203         FHIT      1
    ## 205        TFPI2      1
    ## 210       SCARA5      1
    ## 211       CHRDL1      1
    ## 214        EPHB6      1
    ## 216      ANGPTL3      1
    ## 229        DAPP1      1
    ## 231        THBS2      1
    ## 235        RRM2B      1
    ## 239         KYNU      1
    ## 241         STX4      1
    ## 242         GHRL      1
    ## 245        ITGA5      1
    ## 247        HSPB1      1
    ## 248      PLA2G1B      1
    ## 252       PDGFRA      1
    ## 253           F9      1
    ## 258         NUDC      1
    ## 260         CHGB      1
    ## 262      ALDH1A1      1
    ## 264       SNAP29      1
    ## 267         ARG1      1
    ## 268       STAMBP      1
    ## 269        PDCD6      1
    ## 270          BID      1
    ## 274       DCTPP1      1
    ## 276        SLIT2      1
    ## 279         FRZB      1
    ## 281        ICAM2      1
    ## 283         GLB1      1
    ## 284       COL6A3      1
    ## 285       SPINK5      1
    ## 286        XRCC4      1
    ## 287        ANXA5      1
    ## 288         REG4      1
    ## 289        FABP1      1
    ## 291        ERP44      1
    ## 292       CLEC4G      1
    ## 293        DECR1      1
    ## 294        TYRO3      1
    ## 298       LILRA2      1
    ## 300       HPCAL1      1
    ## 302         SIAE      1
    ## 303          CCS      1
    ## 304        APLP1      1
    ## 305       SPINK1      1
    ## 306          HGF      1
    ## 309       ENTPD6      1
    ## 310        FCRL6      1
    ## 312         CAPG      1
    ## 316       TXNDC5      1
    ## 317         IDI2      1
    ## 318          BTC      1
    ## 320      ALDH3A1      1
    ## 324         AIF1      1
    ## 326        IGSF8      1
    ## 327         CRNN      1
    ## 329        EPHA2      1
    ## 331       EFEMP1      1
    ## 332          FGR      1
    ## 333         ARSB      1
    ## 338      SIGLEC1      1
    ## 340         DSC2      1
    ## 341         BMP6      1
    ## 342       PCOLCE      1
    ## 343     SERPINA9      1
    ## 344         FIS1      1
    ## 345       IL1RL2      1
    ## 347      SULT2A1      1
    ## 348          IDS      1
    ## 350      TNFAIP8      1
    ## 351         TFPI      1
    ## 352          MAX      1
    ## 353         PMVK      1
    ## 354     SERPINB8      1
    ## 355         CTSH      1
    ## 356       IFNLR1      1
    ## 357         AOC1      1
    ## 358         BCAN      1
    ## 359        BLVRB      1
    ## 360         CCL5      1
    ## 361         CLMP      1
    ## 362         MMP8      1
    ## 364         PLAT      1
    ## 365        GDF15      1
    ## 366         ESAM      1
    ## 368         SORD      1
    ## 369          DBI      1
    ## 370        CD274      1
    ## 371        NUDT2      1
    ## 372         NT5E      1
    ## 373        SPON2      1
    ## 374        PDCD5      1
    ## 376        ENPP7      1
    ## 377        CXCL6      1
    ## 378        FABP2      1
    ## 379    SERPINA12      1
    ## 380         MYOC      1
    ## 381          PI3      1
    ## 382         APOH      1
    ## 383        TIMP4      1
    ## 384         UMOD      1
    ## 385        HEBP1      1
    ## 386       IL15RA      1
    ## 387         TLR3      1
    ## 388          LHB      1
    ## 389        SIRPA      1
    ## 390      RARRES1      1
    ## 391          OGN      1
    ## 392       SIRPB1      1
    ## 393        CHIT1      1
    ## 394          OXT      1
    ## 395        ENPP5      1

A total of 1017 SNPs and 395 proteins were found in both UKBB and
deCODE.
