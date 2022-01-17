sort_heatmap_dimension <-
    function(srt, genes, the_gene, key = NULL, value = NULL, use_assay ="RNA") {

        if (!is.null(key) && !is.null(value)) {
            temp <- FetchData(srt, key) == value
            srt <- subset(srt, cells = attr(temp, 'dimnames')[[1]][temp])
        }

        expressed_genes <- srt %>%
            GetAssayData("data", use_assay) %>%
            rowSums() %>%
            .[. > 1] %>%
            names()

        srt <- Seurat::DietSeurat(srt, features = expressed_genes)

        use_genes <- unique(genes) %>%
                .[. %in% expressed_genes]
        holder <- srt %>%
            GetAssayData("data", use_assay) %>%
            as.data.frame() %>%
            .[use_genes, ] %>%
            as.data.frame() %>%
            t()

        sorted_genes <- holder %>%
            corrr::correlate() %>%
            arrange(get(the_gene)) %>%
            .$term

        sorted_cells <- holder %>%
            as.data.frame() %>%
            rownames_to_column(var = "cname") %>%
            select(cname, rev(sorted_genes)) %>%
            arrange(!!!str_c("-", rev(sorted_genes))) %>%
            .$cname

        return(list(sorted_genes, sorted_cells))
}
