`%>%` <- magrittr::`%>%`

#' Available RAM in kB
check_ram <- function() {
  as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))
}

#' Available cores
available_cores <- function(prop2use = .9) {
    max(1, floor(parallel::detectCores(logical = FALSE) * prop2use))
}

#' Get number of cores to fit RAM needs
cores4ram <- function(need) {
  max(1, min(available_cores(), floor(check_ram() / need)))
}

sort_heatmap_dimension <-
  function(srt, genes, the_gene, key = NULL, value = NULL, use_assay = "RNA") {
    stopifnot(require(Seurat))
    if (!is.null(key) && !is.null(value)) {
      temp <- Seurat::FetchData(srt, key) == value
      srt <- subset(srt, cells = attr(temp, "dimnames")[[1]][temp])
    }

    expressed_genes <- srt %>%
      Seurat::GetAssayData("data", use_assay) %>%
      rowSums() %>%
      .[. > 1] %>%
      names()

    srt <- Seurat::DietSeurat(srt, features = expressed_genes)

    use_genes <- unique(genes) %>%
      .[. %in% expressed_genes]
    holder <- srt %>%
      Seurat::GetAssayData("data", use_assay) %>%
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

#' Save plot
save_my_plot <- function(name,
                         plt,
                         type,
                         h = 10,
                         asp = 1.618,
                         path = plots_dir,
                         format = ".pdf") {
  cowplot::save_plot(
    filename = here::here(
      path,
      stringr::str_glue(type,
        as.character(name),
        format,
        .sep = "_"
      )
    ),
    plot = plt,
    base_height = h,
    base_asp = asp,
    limitsize = FALSE
  )
}

occurence_ident <- function(srt) {
  stopifnot(require(Seurat))
  hold_idents <- rep_len("negative", ncol(srt))
  names(hold_idents) <- colnames(srt)
  hold_idents[Seurat::WhichCells(srt, expression = Th > 0)] <- "positive"
  return(hold_idents)
}

occurence_markers <- function(srt) {
  stopifnot(require(Seurat))
  if (!Seurat::Idents(srt) == srt$status) {
    srt$status <- occurence_ident(srt)
    Seurat::Idents(srt) <- srt$status
  }
  markers <- Seurat::FindAllMarkers(srt)
  return(markers)
}

matrix_corr <- function(mydata, glist) {
  stopifnot(require(magrittr))
  mydata <- tibble::as_tibble(mydata) %>% dplyr::select(!!glist)
  plot <- (ggstatsplot::ggcorrmat(
    data = mydata,
    ggcorrplot.args = list(outline.color = "black", hc.order = TRUE)
  ))
}
