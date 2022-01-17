library(here)
library(tidyverse)
library(magrittr)
library(zeallot)
library(Seurat)
library(SeuratDisk)

src_dir  <- here("code")
data_dir <- here("data")
source(here(src_dir, "genes.R"))
source(here(src_dir, "functions.R"))

rar2020_ages_all     <- c("E15", "E17", "P00", "P02", "P10", "P23")
rar2020_ages_postnat <-               c("P00", "P02", "P10", "P23")

samples_df <- read_tsv(here("data/samples.tsv"))
rar2020_srt_pub <- read_rds(here(data_dir, "oldCCA_nae_srt.rds"))
colours <- read_lines(here(data_dir, "colours_wtree.tsv"))
clrlev <- levels(rar2020_srt_pub$wtree)
names(colours) <- clrlev
rm(rar2020_srt_pub)
gc()

# neuro  <- LoadH5Seurat(here(data_dir, "rar2020.srt.neuro.raw.h5seurat"))
onecut <- LoadH5Seurat(here(data_dir, "rar2020.srt.cont.oc2or3.raw.h5seurat"))
onecut3 <- subset(onecut, subset = Onecut3 > 0)

c(cge15, coe15) %<-% sort_heatmap_dimension(onecut,
        genes = c("Onecut2", "Onecut3", glutr),
        the_gene = "Onecut3",
        key = "age",
        value = "E15",
        use_assay = "RNA")

p15 <- DoHeatmap(onecut,
        cells = coe15,
        features = cge15,
        size = 3,
        assay = "RNA",
        slot = "scale.data",
        group.by = "age")

# pd15 <- DotPlot(subset(neuro, subset = age == "E15"), features = sorted_genes, group.by = "wtree")

c(cge17, coe17) %<-% sort_heatmap_dimension(onecut3,
        genes = c("Onecut2", "Onecut3", glutr),
        the_gene = "Onecut3",
        key = "age",
        value = "E17",
        use_assay = "RNA")

p17 <- DoHeatmap(onecut, cells = coe17, features = cge17, size = 3, assay = "RNA", slot = "scale.data", group.by = "age") 
# pd17 <- DotPlot(subset(neuro, subset = age == "E17"), features = cge17, group.by = "wtree")


cgp00 <- GetAssayData(subset(onecut, 
                                          subset = age == "P00"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      corrr::correlate() %>% 
                      arrange(Onecut3) %>% 
                      .$term

cop00 <- GetAssayData(subset(onecut, 
                                          subset = age == "P00"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      rownames_to_column(var = "cname") %>% 
                      select(cname, Onecut3) %>% 
                      arrange(-Onecut3) %>% 
                      .$cname

p00 <- DoHeatmap(onecut, cells = cop00, features = cgp00, size = 3, assay = "RNA", slot = "scale.data", group.by = "age")
# pd00 <- DotPlot(subset(neuro, subset = age == "P00"), features = cgp00, group.by = "wtree")

cgp02 <- GetAssayData(subset(onecut, 
                                          subset = age == "P02"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      corrr::correlate() %>% 
                      arrange(Onecut3) %>% 
                      .$term

cop02 <- GetAssayData(subset(onecut, 
                                          subset = age == "P02"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      rownames_to_column(var = "cname") %>% 
                      select(cname, Onecut3) %>% 
                      arrange(-Onecut3) %>% 
                      .$cname

p02 <- DoHeatmap(onecut, cells = cop02, features = cgp02, size = 3, assay = "RNA", slot = "scale.data", group.by = "age")
# pd02 <- DotPlot(subset(neuro, subset = age == "P02"), features = cgp02, group.by = "wtree")

cgp10 <- GetAssayData(subset(onecut, 
                                          subset = age == "P10"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      corrr::correlate() %>% 
                      arrange(Onecut3) %>% 
                      .$term

cop10 <- GetAssayData(subset(onecut, 
                                          subset = age == "P10"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      rownames_to_column(var = "cname") %>% 
                      select(cname, Onecut3) %>% 
                      arrange(-Onecut3) %>% 
                      .$cname

p10 <- DoHeatmap(onecut, cells = cop10, features = cgp10, size = 3, assay = "RNA", slot = "scale.data", group.by = "age")
# pd10 <- DotPlot(subset(neuro, subset = age == "P10"), features = cgp10, group.by = "wtree")

cgp23 <- GetAssayData(subset(onecut, 
                                          subset = age == "P23"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      corrr::correlate() %>% 
                      arrange(Onecut3) %>% 
                      .$term

cop23 <- GetAssayData(subset(onecut, 
                                          subset = age == "P23"), 
                                   "data", 
                                   "RNA") %>% 
                      as.data.frame() %>% 
                      .[unique(c("Onecut2", "Onecut3", glutr)) %>% 
                              .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
                      as.data.frame() %>% 
                      t() %>% 
                      as.data.frame() %>% 
                      rownames_to_column(var = "cname") %>% 
                      select(cname, Onecut3) %>% 
                      arrange(-Onecut3) %>% 
                      .$cname

p23 <- DoHeatmap(onecut, cells = cop23, features = cgp23, size = 3, assay = "RNA", slot = "scale.data", group.by = "age")
# pd23 <- DotPlot(subset(neuro, subset = age == "P23"), features = cgp23, group.by = "wtree")

cgp <- GetAssayData(onecut3, "data", "RNA") %>%
        as.data.frame() %>%
        .[unique(c("Onecut2", "Onecut3", glutr)) %>%
                .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>%
        as.data.frame() %>%
        t() %>%
        corrr::correlate() %>%
        arrange(Onecut3) %>%
        .$term

#TODO Plot only Glutamatergic cells and remove progenitors
pd <- DotPlot(onecut, features = cgp, group.by = "age")
pd3 <- DotPlot(onecut3, features = cgp, group.by = "age")
ph <- DoHeatmap(onecut3, features = cgp, size = 3, assay = "RNA", slot = "scale.data", group.by = "age")

#TODO 1) split for Glut and GABA separetly, 2) remove progenitors (place only cells with Slc17a6 or Gad1/2 or Slc32a1). Make plot with Th and Trh
pdnt3 <- DotPlot(onecut3, features = c("Onecut2", "Onecut3", neurotrans), group.by = "age")
phnt <- DoHeatmap(onecut3, features = c("Onecut2", "Onecut3", neurotrans), size = 3, assay = "RNA", slot = "scale.data", group.by = "age")
phntr <- DoHeatmap(onecut3, features = c("Onecut2", "Onecut3", neurotrans), size = 3, assay = "RNA", slot = "data", group.by = "age")


corn <- GetAssayData(onecut, "data", "RNA") %>% 
            as.data.frame() %>% 
            .[unique(c("Onecut2", "Onecut3", glutr, neurotrans, "Th", "Trh")) %>% 
                    .[. %in% row.names(onecut@assays$RNA@data)[onecut@assays$RNA@data %>% rowSums() > 1]], ] %>% 
            as.data.frame() %>% 
            t() %>% 
            corrr::correlate()
write_csv(corn, file = here("output/cor_oc2or3.csv"))

corn3 <- GetAssayData(onecut3, "data", "RNA") %>% 
            as.data.frame() %>% 
            .[unique(c("Onecut2", "Onecut3", glutr, neurotrans, "Th", "Trh")) %>% 
                    .[. %in% row.names(onecut3@assays$RNA@data)[onecut3@assays$RNA@data %>% rowSums() > 1]], ] %>%
            as.data.frame() %>% 
            t() %>% 
            corrr::correlate()
write_csv(corn3, file = here("output/cor_oc3.csv"))

.vsc.attach()

p15
pd15 + RotatedAxis()
p17
pd17 + RotatedAxis()
p00
pd00 + RotatedAxis()
p02
pd02 + RotatedAxis()
p10
pd10 + RotatedAxis()
p23
pd23 + RotatedAxis()

pd + RotatedAxis()
pd3 + RotatedAxis()
ph

pdnt3 + RotatedAxis()
phnt
phntr


