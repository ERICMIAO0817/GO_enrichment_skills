#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

if (requireNamespace("BiocParallel", quietly = TRUE)) {
  try({
    BiocParallel::register(BiocParallel::SerialParam())
  }, silent = TRUE)
}

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    token <- args[[i]]
    if (startsWith(token, "--")) {
      key <- substring(token, 3)
      if (i == length(args) || startsWith(args[[i + 1]], "--")) {
        out[[key]] <- TRUE
        i <- i + 1
      } else {
        out[[key]] <- args[[i + 1]]
        i <- i + 2
      }
    } else {
      i <- i + 1
    }
  }
  out
}

write_empty <- function(path, fields) {
  df <- as.data.frame(matrix(nrow = 0, ncol = length(fields)))
  names(df) <- fields
  write.csv(df, file = path, row.names = FALSE)
}

write_no_data_plot <- function(path, title_txt, msg_txt) {
  if (is.null(path) || path == "") return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(filename = path, width = 1440, height = 900, res = 180)
  par(mar = c(2, 2, 2, 2))
  plot.new()
  title(main = title_txt)
  text(0.5, 0.5, msg_txt, cex = 1)
  dev.off()
}

render_ora_plots <- function(ora_df, dotplot_out = NULL, barplot_out = NULL) {
  has_dot <- !(is.null(dotplot_out) || dotplot_out == "")
  has_bar <- !(is.null(barplot_out) || barplot_out == "")
  if (!(has_dot || has_bar)) return(invisible(NULL))

  top_df <- ora_df
  if (!is.null(top_df) && nrow(top_df) > 0) {
    suppressWarnings({
      top_df$padj_num <- as.numeric(top_df$padj)
      top_df$hits_num <- as.numeric(top_df$gene_hits)
    })
    top_df$padj_num[is.na(top_df$padj_num)] <- 1
    top_df$hits_num[is.na(top_df$hits_num) | top_df$hits_num <= 0] <- 1
    top_df <- top_df[order(top_df$padj_num, decreasing = FALSE), , drop = FALSE]
    top_df <- head(top_df, 15)
  }

  if (is.null(top_df) || nrow(top_df) == 0) {
    write_no_data_plot(dotplot_out, "GO Enrichment Dotplot", "No ORA terms available")
    write_no_data_plot(barplot_out, "GO Enrichment Top Terms", "No ORA terms available")
    return(invisible(NULL))
  }

  x <- -log10(pmax(top_df$padj_num, 1e-300))
  hits <- top_df$hits_num
  labels <- substr(as.character(top_df$term_name), 1, 48)
  ontology <- toupper(trimws(as.character(top_df$ontology)))
  ontology <- sub("^GO:", "", ontology)
  ontology[!(ontology %in% c("BP", "CC", "MF"))] <- "UNKNOWN"
  ontology_palette <- c(BP = "#4C78A8", CC = "#54A24B", MF = "#F58518", UNKNOWN = "#9E9E9E")
  ontology_labels <- c(BP = "BP", CC = "CC", MF = "MF", UNKNOWN = "Other")
  point_fill <- unname(ontology_palette[ontology])
  legend_order <- c("BP", "CC", "MF", "UNKNOWN")
  present_ontology <- legend_order[legend_order %in% unique(ontology)]

  if (has_dot) {
    dir.create(dirname(dotplot_out), recursive = TRUE, showWarnings = FALSE)
    png(filename = dotplot_out, width = 1440, height = 900, res = 180)
    par(mar = c(5, 5, 4, 8), xpd = NA)
    cex_vals <- pmax(0.8, pmin(3.0, sqrt(hits)))
    plot(
      x,
      hits,
      pch = 21,
      bg = point_fill,
      col = "#1F2D3D",
      cex = cex_vals,
      xlab = "-log10(FDR)",
      ylab = "Gene hits",
      main = "GO Enrichment Dotplot"
    )
    if (length(present_ontology) > 0) {
      legend(
        "topright",
        inset = c(-0.24, 0),
        legend = ontology_labels[present_ontology],
        pt.bg = ontology_palette[present_ontology],
        pch = 21,
        title = "Ontology",
        bty = "n",
        pt.cex = 1.2
      )
    }
    dev.off()
  }

  if (has_bar) {
    dir.create(dirname(barplot_out), recursive = TRUE, showWarnings = FALSE)
    png(filename = barplot_out, width = 1620, height = 900, res = 180)
    par(mar = c(5, 16, 4, 8), xpd = NA)
    idx <- rev(seq_len(nrow(top_df)))
    barplot(
      x[idx],
      names.arg = labels[idx],
      horiz = TRUE,
      las = 1,
      col = point_fill[idx],
      border = NA,
      xlab = "-log10(FDR)",
      main = "GO Enrichment Top Terms",
      cex.names = 0.7
    )
    if (length(present_ontology) > 0) {
      legend(
        "topright",
        inset = c(-0.24, 0),
        legend = ontology_labels[present_ontology],
        fill = ontology_palette[present_ontology],
        title = "Ontology",
        bty = "n"
      )
    }
    dev.off()
  }

  invisible(NULL)
}

as_bool <- function(x) {
  if (is.null(x)) return(FALSE)
  tolower(as.character(x)) %in% c("1", "true", "yes")
}

non_empty_count <- function(x) {
  if (is.null(x)) return(0L)
  y <- trimws(as.character(x))
  sum(!is.na(y) & y != "")
}

extract_non_empty <- function(x) {
  if (is.null(x)) return(character(0))
  y <- trimws(as.character(x))
  unique(y[!is.na(y) & y != ""])
}

pick_gene_keytype <- function(df) {
  entrez_vals <- if ("entrez" %in% names(df)) df[["entrez"]] else NULL
  symbol_vals <- if ("symbol" %in% names(df)) df[["symbol"]] else NULL
  ensembl_vals <- if ("ensembl" %in% names(df)) df[["ensembl"]] else NULL
  input_vals <- if ("input_id" %in% names(df)) df[["input_id"]] else NULL

  entrez_n <- non_empty_count(entrez_vals)
  symbol_n <- non_empty_count(symbol_vals)
  ensembl_n <- non_empty_count(ensembl_vals)

  if (entrez_n >= symbol_n && entrez_n >= ensembl_n && entrez_n > 0) {
    return(list(keyType = "ENTREZID", genes = extract_non_empty(entrez_vals)))
  }
  if (symbol_n >= ensembl_n && symbol_n > 0) {
    return(list(keyType = "SYMBOL", genes = extract_non_empty(symbol_vals)))
  }
  if (ensembl_n > 0) {
    return(list(keyType = "ENSEMBL", genes = extract_non_empty(ensembl_vals)))
  }
  list(keyType = "SYMBOL", genes = extract_non_empty(input_vals))
}

to_ont_list <- function(x) {
  if (is.null(x) || toupper(x) == "ALL") return(c("BP", "MF", "CC"))
  items <- unlist(strsplit(toupper(x), "[,[:space:]]+"))
  items <- items[items %in% c("BP", "MF", "CC")]
  if (length(items) == 0) return(c("BP", "MF", "CC"))
  unique(items)
}

parse_scores <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  raw <- read.table(path, sep = "", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(raw) < 2) return(NULL)
  genes <- as.character(raw[[1]])
  vals <- suppressWarnings(as.numeric(raw[[2]]))
  keep <- !is.na(vals) & genes != ""
  genes <- genes[keep]
  vals <- vals[keep]
  if (length(genes) == 0) return(NULL)
  dedup <- !duplicated(genes)
  genes <- genes[dedup]
  vals <- vals[dedup]
  out <- vals
  names(out) <- genes
  out <- sort(out, decreasing = TRUE)
  out
}

guess_rank_keytype <- function(ids) {
  x <- trimws(as.character(ids))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return("SYMBOL")

  frac_entrez <- mean(grepl("^[0-9]+$", x))
  frac_ensembl <- mean(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", toupper(x)))

  if (frac_entrez >= 0.8) return("ENTREZID")
  if (frac_ensembl >= 0.8) return("ENSEMBL")
  "SYMBOL"
}

read_normalized_input <- function(path) {
  tab_df <- tryCatch(
    read.table(path, sep = "\t", header = TRUE, quote = "", stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) NULL
  )
  need <- c("input_id", "symbol", "entrez", "ensembl")
  if (is.null(tab_df) || !all(need %in% names(tab_df))) {
    csv_df <- tryCatch(
      read.csv(path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(csv_df)) {
      tab_df <- csv_df
    }
  }
  if (is.null(tab_df)) {
    tab_df <- data.frame(stringsAsFactors = FALSE)
  }
  for (nm in need) {
    if (!(nm %in% names(tab_df))) {
      tab_df[[nm]] <- character(nrow(tab_df))
    }
    tab_df[[nm]] <- as.character(tab_df[[nm]])
  }
  tab_df
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

ora_fields <- c("term_id", "term_name", "ontology", "source", "pvalue", "padj", "gene_hits", "gene_ratio", "engine")
gsea_fields <- c("term_id", "term_name", "ontology", "source", "nes", "pvalue", "padj", "leading_edge", "direction", "engine")

if (is.null(args[["ora-out"]])) stop("--ora-out is required")
if (is.null(args[["normalized-input"]])) stop("--normalized-input is required")

ora_out <- normalizePath(args[["ora-out"]], mustWork = FALSE)
gsea_out <- if (!is.null(args[["gsea-out"]])) normalizePath(args[["gsea-out"]], mustWork = FALSE) else NULL
dotplot_out <- if (!is.null(args[["dotplot-out"]])) normalizePath(args[["dotplot-out"]], mustWork = FALSE) else NULL
barplot_out <- if (!is.null(args[["barplot-out"]])) normalizePath(args[["barplot-out"]], mustWork = FALSE) else NULL
mode <- if (!is.null(args[["mode"]])) tolower(args[["mode"]]) else "auto"
species <- if (!is.null(args[["species"]])) tolower(args[["species"]]) else "human"
threshold <- if (!is.null(args[["sig-threshold"]])) as.numeric(args[["sig-threshold"]]) else 0.05
onts <- to_ont_list(args[["ontology"]])

has_clusterprofiler <- requireNamespace("clusterProfiler", quietly = TRUE)
has_annotation <- requireNamespace("AnnotationDbi", quietly = TRUE)
orgdb <- NULL
if (species == "human" && requireNamespace("org.Hs.eg.db", quietly = TRUE)) orgdb <- get("org.Hs.eg.db", asNamespace("org.Hs.eg.db"))
if (species == "mouse" && requireNamespace("org.Mm.eg.db", quietly = TRUE)) orgdb <- get("org.Mm.eg.db", asNamespace("org.Mm.eg.db"))

if (!(has_clusterprofiler && has_annotation && !is.null(orgdb))) {
  write_empty(ora_out, ora_fields)
  if (!is.null(gsea_out)) write_empty(gsea_out, gsea_fields)
  write_no_data_plot(dotplot_out, "GO Enrichment Dotplot", "R dependencies missing")
  write_no_data_plot(barplot_out, "GO Enrichment Top Terms", "R dependencies missing")
  message("Missing R dependencies: clusterProfiler/AnnotationDbi/org.*.eg.db")
  quit(status = 2)
}

norm <- read_normalized_input(args[["normalized-input"]])
picked <- pick_gene_keytype(norm)
gene_ids <- picked$genes
key_type <- picked$keyType

if (length(gene_ids) == 0) {
  write_empty(ora_out, ora_fields)
  if (!is.null(gsea_out)) write_empty(gsea_out, gsea_fields)
  write_no_data_plot(dotplot_out, "GO Enrichment Dotplot", "No normalized genes available")
  write_no_data_plot(barplot_out, "GO Enrichment Top Terms", "No normalized genes available")
  message("No normalized genes available")
  quit(status = 3)
}

run_ora <- mode %in% c("auto", "ora", "gsea")
run_gsea <- mode %in% c("auto", "gsea") && !is.null(args[["scores"]])

ora_rows <- list()
if (run_ora) {
  for (ont in onts) {
    res <- tryCatch(
      {
        clusterProfiler::enrichGO(
          gene = gene_ids,
          OrgDb = orgdb,
          keyType = key_type,
          ont = ont,
          pAdjustMethod = "BH",
          pvalueCutoff = threshold,
          qvalueCutoff = threshold,
          readable = FALSE
        )
      },
      error = function(e) NULL
    )

    if (!is.null(res)) {
      df <- as.data.frame(res)
      if (nrow(df) > 0) {
        for (i in seq_len(nrow(df))) {
          ora_rows[[length(ora_rows) + 1]] <- data.frame(
            term_id = as.character(df$ID[[i]]),
            term_name = as.character(df$Description[[i]]),
            ontology = ont,
            source = "clusterProfiler",
            pvalue = as.numeric(df$pvalue[[i]]),
            padj = as.numeric(df$p.adjust[[i]]),
            gene_hits = as.integer(df$Count[[i]]),
            gene_ratio = as.character(df$GeneRatio[[i]]),
            engine = "r",
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if (length(ora_rows) > 0) {
  ora_df <- do.call(rbind, ora_rows)
  write.csv(ora_df[, ora_fields], file = ora_out, row.names = FALSE)
} else {
  ora_df <- as.data.frame(matrix(nrow = 0, ncol = length(ora_fields)))
  names(ora_df) <- ora_fields
  write_empty(ora_out, ora_fields)
}

render_ora_plots(ora_df, dotplot_out = dotplot_out, barplot_out = barplot_out)

if (!is.null(gsea_out)) {
  if (run_gsea) {
    ranked <- parse_scores(args[["scores"]])
    if (is.null(ranked)) {
      write_empty(gsea_out, gsea_fields)
    } else {
      gsea_rows <- list()
      for (ont in onts) {
        gsea_key_type <- guess_rank_keytype(names(ranked))
        max_gs <- max(10, length(ranked) - 1)
        gres <- tryCatch(
          {
            clusterProfiler::gseGO(
              geneList = ranked,
              OrgDb = orgdb,
              keyType = gsea_key_type,
              ont = ont,
              pvalueCutoff = threshold,
              pAdjustMethod = "BH",
              minGSSize = 5,
              maxGSSize = max_gs,
              verbose = FALSE
            )
          },
          error = function(e) NULL
        )

        if (!is.null(gres)) {
          gdf <- as.data.frame(gres)
          if (nrow(gdf) > 0) {
            for (i in seq_len(nrow(gdf))) {
              nes <- as.numeric(gdf$NES[[i]])
              gsea_rows[[length(gsea_rows) + 1]] <- data.frame(
                term_id = as.character(gdf$ID[[i]]),
                term_name = as.character(gdf$Description[[i]]),
                ontology = ont,
                source = "clusterProfiler",
                nes = nes,
                pvalue = as.numeric(gdf$pvalue[[i]]),
                padj = as.numeric(gdf$p.adjust[[i]]),
                leading_edge = as.character(gdf$core_enrichment[[i]]),
                direction = ifelse(is.na(nes), "unknown", ifelse(nes >= 0, "up", "down")),
                engine = "r",
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }

      if (length(gsea_rows) > 0) {
        gsea_df <- do.call(rbind, gsea_rows)
        write.csv(gsea_df[, gsea_fields], file = gsea_out, row.names = FALSE)
      } else {
        write_empty(gsea_out, gsea_fields)
      }
    }
  } else {
    write_empty(gsea_out, gsea_fields)
  }
}

message(sprintf("R enrichment completed (keyType=%s, genes=%d)", key_type, length(gene_ids)))
quit(status = 0)
