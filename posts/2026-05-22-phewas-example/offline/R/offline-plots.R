# Adapted from MRCIEU/gbc-mr 3bbac89e254acfbfa3385aae9b4520e4e8b5932c.
# Base graphics exports. Display metadata never changes the scientific inputs.
format_probability <- function(x) ifelse(is.finite(x), formatC(x, digits = 2, format = "g"), "NA")
wrap_label <- function(x, width = 88) paste(strwrap(x, width), collapse = "\n")

.plot_value <- function(x, name, fallback) {
  if (!name %in% names(x) || !nrow(x)) return(fallback)
  value <- as.character(x[[name]][1])
  if (is.na(value) || !nzchar(trimws(value))) fallback else value
}
.plot_units <- function(x, column, kind) {
  value <- .plot_value(x, column, paste("Saved", kind, "GWAS units; scale unavailable"))
  if (grepl("unavailable|unknown", value, ignore.case = TRUE))
    paste("Saved", kind, "GWAS units (scale unknown)") else value
}
.plot_outcome <- function(x) .plot_value(x, "outcome_name", "Outcome")
.plot_limits <- function(lo, hi) {
  limits <- range(c(0, lo[is.finite(lo)], hi[is.finite(hi)]))
  if (diff(limits) == 0) limits <- c(-1, 1)
  limits + c(-1, 1) * diff(limits) * 0.08
}
.plot_family <- function(p, adjusted, name) {
  ok <- is.finite(p) & p >= 0 & p <= 1
  sprintf("%s: tested %d/%d; unavailable %d; P<0.05: %d; BH<0.05: %d",
    name, sum(ok), length(p), sum(!ok), sum(p[ok] < 0.05), sum(adjusted[ok] < 0.05, na.rm = TRUE))
}
.plot_exclusions <- function(x) {
  leads <- .plot_value(x, "n_exclusion_loci", "all known")
  window <- if ("window_bp" %in% names(x)) format(x$window_bp[1] / 1e6, trim = TRUE) else "1"
  paste0("Shared comparator excludes ", leads, " outcome leads and inclusive +/-", window, " Mb regions")
}

draw_estimates <- function(b, lo, hi, labels, colours, title, subtitle, xlab,
                           left_margin = 7, label_cex = 0.8, top_margin = 6,
                           bottom_margin = 4, title_cex = 1, subtitle_cex = 0.85) {
  available <- is.finite(b) & is.finite(lo) & is.finite(hi)
  limits <- .plot_limits(lo, hi)
  y <- rev(seq_along(b))
  par(mar = c(bottom_margin, left_margin, top_margin, 1.5),
      mgp = c(2.0, 0.5, 0), tcl = -0.2)
  plot(NA, xlim = limits, ylim = c(0.5, length(b) + 0.5), yaxt = "n",
    ylab = "", xlab = xlab, cex.lab = 0.9, cex.axis = 0.9, bty = "n")
  abline(v = 0, lty = 3, col = "grey65")
  axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = label_cex)
  segments(lo[available], y[available], hi[available], y[available], col = colours[available], lwd = 1.8)
  points(b[available], y[available], col = colours[available], pch = 19, cex = 0.95)
  if (any(!available)) text(mean(limits), y[!available], "Unavailable (reason in CSV)",
    cex = 0.85, col = "grey40")
  mtext(wrap_label(title, 78), side = 3, line = top_margin - 1.8, adj = 0,
    cex = title_cex, font = 2)
  mtext(subtitle, side = 3, line = 0.5, adj = 0, cex = subtitle_cex)
}

# One multipage PDF; PNG names are stem.png then stem-page-02.png, etc.
# PDF page size is uniform; optional PNG heights adapt to page content.
export_pages <- function(path, n, per_page, draw_page, width = 10, height = 11,
                         png_heights = NULL) {
  if (!n) stop("Cannot export an empty figure inventory")
  pages <- split(seq_len(n), ceiling(seq_len(n) / per_page))
  if (is.null(png_heights)) png_heights <- rep(height, length(pages))
  if (length(png_heights) != length(pages)) stop("One PNG height required per page")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  tryCatch(for (page in seq_along(pages)) draw_page(pages[[page]], page, length(pages)),
    finally = grDevices::dev.off())
  pngs <- character(length(pages))
  for (page in seq_along(pages)) {
    suffix <- if (page == 1) ".png" else sprintf("-page-%02d.png", page)
    pngs[page] <- sub("\\.pdf$", suffix, path)
    grDevices::png(pngs[page], width = width, height = png_heights[page], units = "in", res = 150, type = "cairo")
    tryCatch(draw_page(pages[[page]], page, length(pages)), finally = grDevices::dev.off())
  }
  pngs
}

export_comparison_plot <- function(x, path) {
  mediator_ids <- unique(x$mediator_id)
  counts <- vapply(mediator_ids, function(id) sum(x$mediator_id == id), integer(1))
  heights <- pmax(6.5, 3.8 + 0.64 * (counts + 1))
  q_family <- .plot_family(x$q_pval, x$q_padj, "Q family (all unique pairs)")
  export_pages(path, length(mediator_ids), 1, function(rows, page, pages) {
    par(mfrow = c(1, 1), oma = c(1.5, 0, 0, 0))
    z <- x[x$mediator_id == mediator_ids[rows], , drop = FALSE]
    for (field in c("comparator_b", "comparator_se", "comparator_nsnp"))
      if (length(unique(z[[field]])) != 1L) stop("Comparison panel must use one shared comparator per mediator")
    exact_labels <- paste0(z$target_snp, " | EA ", ifelse(is.na(z$effect_allele), "?", z$effect_allele),
      "\nQ=", format_probability(z$q), "; P=", format_probability(z$q_pval), "; BH=", format_probability(z$q_padj))
    draw_estimates(c(z$target_b, z$comparator_b[1]), c(z$target_lo, z$comparator_lo[1]),
      c(z$target_hi, z$comparator_hi[1]),
      c(exact_labels, paste0("Other SNPs | n=", z$comparator_nsnp[1], "\nOne shared mediator comparator")),
      c(rep("#2166ac", nrow(z)), "#b35806"),
      paste0(z$mediator[1], " | ", z$mediator_id[1], " | page ", page, "/", pages),
      paste(.plot_exclusions(z), q_family, sep = "\n"),
      paste0(.plot_units(z, "outcome_units", "outcome"), "\nper ", .plot_units(z, "mediator_scale", "mediator")),
      left_margin = 21, label_cex = 0.88, top_margin = 6.5, bottom_margin = 4.5)
    mtext(paste0(.plot_outcome(z), ": exact hits (blue) versus comparator (orange); 95% CIs. Heterogeneity is model incompatibility."),
      outer = TRUE, side = 1, line = 0.25, cex = 0.8)
  }, width = 14, height = max(heights), png_heights = heights)
}

export_decomposition_plot <- function(x, path) {
  colours <- c("#2166ac", "#1b7837", "#b35806")
  # Every page has the same outcome-unit axis, so pages can be compared visually.
  limits <- .plot_limits(c(x$t_lo, x$i_lo, x$d_lo), c(x$t_hi, x$i_hi, x$d_hi))
  families <- paste(.plot_family(x$i_pval, x$i_padj, "Indirect family"),
                    .plot_family(x$d_pval, x$d_padj, "Residual family"), sep = "\n")
  export_pages(path, nrow(x), 6, function(rows, page, pages) {
    z <- x[rows, , drop = FALSE]
    par(mfrow = c(1, 1), mar = c(7, 23, 7.3, 2), mgp = c(2.2, 0.6, 0))
    y <- rev(seq_len(nrow(z)))
    plot(NA, xlim = limits, ylim = c(0.4, nrow(z) + 0.6), yaxt = "n", ylab = "",
      xlab = paste0(.plot_units(z, "outcome_units", "outcome"), " per recorded target effect allele"),
      cex.lab = 0.9, cex.axis = 0.9, bty = "n")
    abline(v = 0, lty = 3, col = "grey60")
    labels <- paste0(z$mediator, " | ", z$mediator_id,
      "\n", z$target_snp, " | EA ", ifelse(is.na(z$effect_allele), "?", z$effect_allele), " | other n=", z$comparator_nsnp,
      "\nIndirect BH=", format_probability(z$i_padj), "; residual BH=", format_probability(z$d_padj))
    axis(2, at = y, labels = labels, las = 1, tick = FALSE, cex.axis = 0.85)
    for (component in seq_along(colours)) {
      field <- c("t", "i", "d")[component]
      available <- is.finite(z[[field]]) & is.finite(z[[paste0(field, "_lo")]]) & is.finite(z[[paste0(field, "_hi")]])
      yp <- y + c(0.23, 0, -0.23)[component]
      segments(z[[paste0(field, "_lo")]][available], yp[available],
        z[[paste0(field, "_hi")]][available], yp[available], col = colours[component], lwd = 1.7)
      points(z[[field]][available], yp[available], pch = 19, col = colours[component], cex = 0.85)
      if (any(!available)) text(limits[2], yp[!available], paste(c("Total", "Indirect", "Residual")[component], "unavailable"),
        adj = 1, cex = 0.72, col = colours[component])
    }
    mtext(sprintf("%s: model-based components | page %d/%d", .plot_outcome(z), page, pages),
      side = 3, line = 5.7, adj = 0, font = 2, cex = 1.05)
    mtext(families, side = 3, line = 2.8, adj = 0, cex = 0.8)
    legend("top", inset = c(0, -0.08), xpd = TRUE, horiz = TRUE, bty = "n", cex = 0.9,
      legend = c("Total", "Indirect", "Model-based residual"), col = colours, pch = 19)
    mtext("95% first-order delta-method CIs; zero covariance working assumption. Each mediator is analysed separately.\nResiduals are model-based components; indirect effects are not summed across mediators.",
      side = 1, line = 5.3, cex = 0.8)
  }, width = 15, height = 11)
}

export_pooled_plot <- function(x, path) {
  family <- .plot_family(x$pval, x$pval_bh, "Pooled mediator family")
  export_pages(path, nrow(x), 5, function(rows, page, pages) {
    par(mfrow = c(length(rows), 1), oma = c(5, 0, 2.5, 0), cex = 1)
    for (row in rows) {
      z <- x[row, , drop = FALSE]
      subtitle <- paste0(z$method, "; P=", format_probability(z$pval), "; BH=", format_probability(z$pval_bh),
        "; excluded lead/region n=", z$n_excluded_target + z$n_excluded_region)
      draw_estimates(z$b, z$lo, z$hi, paste0("Shared comparator\nn=", z$nsnp), "#1b7837",
        paste(z$mediator, "|", z$mediator_id), subtitle,
        paste0(.plot_units(z, "outcome_units", "outcome"), " / ", .plot_units(z, "mediator_scale", "mediator")),
        left_margin = 10, top_margin = 3.5, bottom_margin = 3, label_cex = 0.85,
        title_cex = 0.95, subtitle_cex = 0.8)
    }
    mtext(sprintf("%s: pooled mediator MR | page %d/%d", .plot_outcome(x), page, pages), outer = TRUE,
      side = 3, line = 0.7, font = 2, cex = 1.1)
    mtext(paste0(family, "\n95% CIs; one shared comparator per mediator; all outcome lead loci excluded before harmonisation. Panel scales differ."),
      outer = TRUE, side = 1, line = 2.3, cex = 0.82)
  }, width = 14, height = 15)
}
