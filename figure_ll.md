Feature Plots and Differential Expression
================
2025-04-13

``` r
p1 <- DimPlot(HuLu, label = TRUE) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  coord_fixed(ratio = 1)

p2 <- FeaturePlot(HuLu, features = "CD69", order = TRUE) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  coord_fixed(ratio = 1)

p3 <- FeaturePlot(HuLu, features = "ITGAE", order = TRUE) +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  coord_fixed(ratio = 1)

p1 | p2 | p3
```

![](/Users/aj397/Documents/Data/HULU%2010X%202023/R-analysis/github2/figure_ll_files/figure-gfm/FeaturePlots-1.png)<!-- -->

``` r
HuLu$Group <- paste(Idents(HuLu), HuLu$orig.ident, sep = "_")
Idents(HuLu) <- "Group"
levels(HuLu)
```

    ##  [1] "CD8_T_cells_Un_Inf"        "Endothelial_cells_Un_Inf"  "Epithelial_Un_Inf"        
    ##  [4] "CD4_T_cells_Un_Inf"        "NK_cells_Un_Inf"           "Mast_cells_Un_Inf"        
    ##  [7] "T_cell_mito_high_Un_Inf"   "Macrophages_Un_Inf"        "Multiciliated_Un_Inf"     
    ## [10] "Monocytes_Un_Inf"          "B_cells_Un_Inf"            "NK_cells_IAV_Inf"         
    ## [13] "CD8_T_cells_IAV_Inf"       "Macrophages_IAV_Inf"       "CD4_T_cells_IAV_Inf"      
    ## [16] "Epithelial_IAV_Inf"        "B_cells_IAV_Inf"           "T_cell_mito_high_IAV_Inf" 
    ## [19] "Monocytes_IAV_Inf"         "Endothelial_cells_IAV_Inf" "Mast_cells_IAV_Inf"       
    ## [22] "Multiciliated_IAV_Inf"

``` r
my_levels <- c("CD8_T_cells_Un_Inf"   ,     "Endothelial_cells_Un_Inf",  "Epithelial_Un_Inf"   ,      "CD4_T_cells_Un_Inf"    ,   
               "NK_cells_Un_Inf"      ,     "Mast_cells_Un_Inf"    ,     "T_cell_mito_high_Un_Inf" ,  "Macrophages_Un_Inf" ,      
              "Multiciliated_Un_Inf"  ,    "Monocytes_Un_Inf"    ,      "B_cells_Un_Inf",
              "CD8_T_cells_IAV_Inf"   ,     "Endothelial_cells_IAV_Inf",  "Epithelial_IAV_Inf"   ,      "CD4_T_cells_IAV_Inf"    ,   
               "NK_cells_IAV_Inf"      ,     "Mast_cells_IAV_Inf"    ,     "T_cell_mito_high_IAV_Inf" ,  "Macrophages_IAV_Inf" ,      
               "Multiciliated_IAV_Inf"  ,    "Monocytes_IAV_Inf"    ,      "B_cells_IAV_Inf")


# Re-level object@ident
HuLu@active.ident <- factor(x = HuLu@active.ident, levels = my_levels)


### Viral transcripts
DotPlot(HuLu, features = c("fluHA","fluNA","fluNS","fluPB1","fluPA",
                     "fluPB2", "fluNP"),) + theme(axis.text.x = element_text(angle = 45, hjust=1))
```

![](/Users/aj397/Documents/Data/HULU%2010X%202023/R-analysis/github2/figure_ll_files/figure-gfm/DE_setup-1.png)<!-- -->

``` r
CD8_T_cells <-  FindMarkers(HuLu, ident.1 = "CD8_T_cells_IAV_Inf", 
                            ident.2 = "CD8_T_cells_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
CD4_T_cells <-  FindMarkers(HuLu, ident.1 = "CD4_T_cells_IAV_Inf", 
                            ident.2 = "CD4_T_cells_Un_Inf",   min.pct = 0.01,   min.diff.pct = -Inf)
Endothelial_cells <-  FindMarkers(HuLu, ident.1 = "Endothelial_cells_IAV_Inf", 
                                  ident.2 = "Endothelial_cells_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
Epithelial <-  FindMarkers(HuLu, ident.1 = "Epithelial_IAV_Inf", 
                           ident.2 = "Epithelial_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
NK_cells <-  FindMarkers(HuLu, ident.1 = "NK_cells_IAV_Inf", 
                         ident.2 = "NK_cells_Un_Inf", min.pct = 0.01,   min.diff.pct = -Inf)
Mast_cells <-  FindMarkers(HuLu, ident.1 = "Mast_cells_IAV_Inf", 
                           ident.2 = "Mast_cells_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
Macrophages <-  FindMarkers(HuLu, ident.1 = "Macrophages_IAV_Inf", 
                            ident.2 = "Macrophages_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
Multiciliated <-  FindMarkers(HuLu, ident.1 = "Multiciliated_IAV_Inf", 
                              ident.2 = "Multiciliated_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
Monocytes <-  FindMarkers(HuLu, ident.1 = "Monocytes_IAV_Inf", 
                          ident.2 = "Monocytes_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
B_cells <-  FindMarkers(HuLu, ident.1 = "B_cells_IAV_Inf", 
                        ident.2 = "B_cells_Un_Inf",  min.pct = 0.01,   min.diff.pct = -Inf)
```

``` r
vp1 <- EnhancedVolcano(CD8_T_cells, lab = rownames(CD8_T_cells),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "CD8_T_cells", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-15)))

vp2 <- EnhancedVolcano(CD4_T_cells, lab = rownames(CD4_T_cells),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "CD4_T_cells", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-26)))

vp3 <- EnhancedVolcano(Endothelial_cells, lab = rownames(Endothelial_cells),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "Endothelial_cells", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-3)))

vp4 <- EnhancedVolcano(Epithelial, lab = rownames(Epithelial),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "Epithelial", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-3)))

vp5 <- EnhancedVolcano(NK_cells, lab = rownames(NK_cells),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "NK_cells", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-10)))

vp6 <- EnhancedVolcano(Mast_cells, lab = rownames(Mast_cells),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "Mast_cells", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-10)))

vp7 <- EnhancedVolcano(Macrophages, lab = rownames(Macrophages),
                       x = "avg_log2FC", y = 'p_val_adj',
                       title = "Macrophages", legendPosition = "none",
                       pCutoff = 0.1, ylim = c(0, -log10(10e-9)))


grid.arrange(grobs = list(vp1, vp2, vp3, vp5, vp6, vp7), ncol = 3)
```

![](/Users/aj397/Documents/Data/HULU%2010X%202023/R-analysis/github2/figure_ll_files/figure-gfm/VolcanoPlots-1.png)<!-- -->
