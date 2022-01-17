plot_sorted_heatmap <- function(srt, feature, sub) {
        # Specific for Onecut3 and makes subset by age
        srt_filtered <-
                GetAssayData(subset(srt,
                                    subset = age == sub),
                             "data",
                             "RNA") %>%
                as.data.frame() %>%
                .[unique(c(feature)) %>%
                        .[. %in% row.names(srt@assays$RNA@data)[srt@assays$RNA@data %>% rowSums() > 1]], ] %>%
                as.data.frame() %>%
                t()
        custom_genes <- c(
                corrr::correlate(srt_filtered) %>%
                arrange(Onecut3) %>%
                select(term, Onecut3) %>%
                .$term
                )
        custom_order <-
                srt_filtered %>%
                as.data.frame() %>%
                rownames_to_column(var = "cname") %>%
                select(cname, Onecut3) %>%
                arrange(-Onecut3 # Ascending # Descending
                        ) %>%
                .$cname
        (p <- DoHeatmap(srt,
                       cells = custom_order,
                       features = custom_genes,
                       size = 3,
                       assay = "RNA",
                       slot = "scale.data"))
        return(p)
}
