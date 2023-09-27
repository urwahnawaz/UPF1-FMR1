libs = c("dplyr", "ggplot2", "reshape2", "tools", "magrittr", "tibble", "readxl",
         "data.table", "scales", "tidyr", "reshape2", "stringr", "tidyverse", "readxl",
         "corrplot", "purrr", "DESeq2", "edgeR", "biomaRt", "EnsDb.Mmusculus.v79",
         "org.Mm.eg.db", "goseq", "tximport", "tximeta", "viridis", "venn", "ggpubr",
         "MetBrewer", "tidyr", "SummarizedExperiment", "ggvenn", "pheatmap", "GSEABase", "msigdb",
        "DT", "ggiraph", "openxlsx", "readxl","AnnotationHub","tximeta", "ngsReports" )

#"fgsea", "UpSetR", "VennDiagram", "cowplot", "msigdbr", "pander",
   #      "gridExtra", "ggrepel", "IsoformSwitchAnalyzeR")
libsLoaded = lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})
