# vui-202012-1

University of Exeter outcomes analysis for the UK "variant of concern" VOC-202012/1, also referred to a lineage B.1.1.7

This code is the aimed to be reproducible given a set of data files. These data files are not generally shareable due to the potential for re-identification but can
be obtained upon request and signing a data sharing agreement from public health england. The particular files used in our analysis are:

* 20210129 COVID19 Deaths.xlsx
* Anonymised Combined Line List 20210129.csv
* SGTF_linelist_20210129.csv

This code depends on the following libraries:
* tidyverse
* data.table
* survival
* knitr / kable
* DiagrammeR
* DiagrammeRsvg

There is also an optional dependency on a library for producing formatted outputs: [standardPrintOutput](https://github.com/terminological/standard-print-output). 
This has a lot of its own dependencies, and may not be particularly useful.

The main analysis is in the new-variant-mortality.Rmd markdown file which depends on the functions defined in new-variant-mortality.R file, the code is not 
heavily optimised and relies on a substantial mount of memory being available - 16Gb minimum.

The environment the analysis was run on was:

## sessionInfo():

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Linux Mint 18.2

Matrix products: default
BLAS:   /usr/lib/openblas-base/libblas.so.3
LAPACK: /usr/lib/libopenblasp-r0.2.18.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8      
 [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] extrafont_0.17                 standardPrintOutput_0.0.0.9000 testthat_2.3.2                 survminer_0.4.8                ggpubr_0.2.5                   magrittr_1.5                  
 [7] survival_3.1-11                patchwork_1.0.0.9000           forcats_0.5.0                  stringr_1.4.0                  dplyr_1.0.2                    purrr_0.3.4                   
[13] readr_1.3.1                    tidyr_1.1.2                    tibble_3.0.4                   ggplot2_3.3.2                  tidyverse_1.3.0               

loaded via a namespace (and not attached):
  [1] clipr_0.7.0           tidyselect_1.1.0      lme4_1.1-23           htmlwidgets_1.5.1     grid_3.6.3            devtools_2.2.2        EpiEstim_2.2-2        polyCub_0.7.1        
  [9] munsell_0.5.0         units_0.6-5           captioner_2.2.3       statmod_1.4.34        miniUI_0.1.1.1        withr_2.3.0           colorspace_1.4-1      config_0.3           
 [17] knitr_1.30            uuid_0.1-4            rstudioapi_0.13       stats4_3.6.3          incidence_1.7.3       tensor_1.5            ggsignif_0.6.0        officer_0.3.15       
 [25] Rttf2pt1_1.3.8        labeling_0.3          bbmle_1.0.23.1        polyclip_1.10-0       KMsurv_0.1-5          MCMCpack_1.4-6        farver_2.0.3          rprojroot_1.3-2      
 [33] coda_0.19-3           vctrs_0.3.4           generics_0.1.0        coarseDataTools_0.6-5 metafor_2.4-0         xfun_0.18             R6_2.4.1              rsvg_2.1             
 [41] spatstat.utils_1.17-0 assertthat_0.2.1      promises_1.1.0        scales_1.1.1          DiagrammeRsvg_0.1     ggExtra_0.9           gtable_0.3.0          npsurv_0.4-0         
 [49] goftest_1.2-2         processx_3.4.4        mcmc_0.9-6.1          rlang_0.4.8           MatrixModels_0.4-1    systemfonts_0.3.2     splines_3.6.3         extrafontdb_1.0      
 [57] checkmate_2.0.0       broom_0.7.2           abind_1.4-5           yaml_2.2.1            reshape2_1.4.3        modelr_0.1.6          crosstalk_1.0.0       backports_1.1.10     
 [65] httpuv_1.5.2          DiagrammeR_1.0.6.1    surveillance_1.18.0   tools_3.6.3           tcltk_3.6.3           usethis_1.5.1         ellipsis_0.3.1        RColorBrewer_1.1-2   
 [73] sessioninfo_1.1.1     Rcpp_1.0.4            plyr_1.8.6            base64enc_0.1-3       visNetwork_2.0.9      classInt_0.4-2        ps_1.4.0              prettyunits_1.1.1    
 [81] deldir_0.1-28         rpart_4.1-15          zoo_1.8-7             haven_2.2.0           fs_1.3.2              data.table_1.12.8     magick_2.3            flextable_0.5.11     
 [89] SparseM_1.78          reprex_0.3.0          ukcovidtools_0.0.1    mvtnorm_1.1-0         fitdistrplus_1.0-14   pkgload_1.1.0         hms_0.5.3             lsei_1.2-0           
 [97] mime_0.9              evaluate_0.14         xtable_1.8-4          leaflet_2.0.3         readxl_1.3.1          gridExtra_2.3         compiler_3.6.3        bdsmatrix_1.3-4      
[105] KernSmooth_2.23-18    V8_3.2.0              crayon_1.3.4          minqa_1.2.4           htmltools_0.4.0       mgcv_1.8-33           staplr_2.9.0          later_1.0.0          
[113] expm_0.999-4          lubridate_1.7.4       DBI_1.1.0             dbplyr_1.4.2          MASS_7.3-53           sf_0.9-3              boot_1.3-24           Matrix_1.2-18        
[121] cli_2.1.0             huxtable_5.1.0        pkgconfig_2.0.3       km.ci_0.5-2           metaplus_0.7-11       sp_1.4-1              numDeriv_2016.8-1.1   binom_1.1-1          
[129] xml2_1.3.2            webshot_0.5.2         gghighlight_0.3.0     rvest_0.3.5           callr_3.5.1           digest_0.6.26         fastGHQuad_1.0        spatstat.data_1.4-3  
[137] msm_1.6.8             rmarkdown_2.5         cellranger_1.1.0      survMisc_0.5.5        gdtools_0.2.1         curl_4.3              shiny_1.4.0.2         quantreg_5.54        
[145] nloptr_1.2.2.2        lifecycle_0.2.0       nlme_3.1-145          jsonlite_1.6.1        desc_1.2.0            fansi_0.4.1           pillar_1.4.6          lattice_0.20-41      
[153] fastmap_1.0.1         httr_1.4.2            pkgbuild_1.1.0        glue_1.4.2            remotes_2.1.1         spatstat_1.64-1       zip_2.1.1             png_0.1-7            
[161] class_7.3-15          stringi_1.4.6         metR_0.7.0            memoise_1.1.0         e1071_1.7-3          
```