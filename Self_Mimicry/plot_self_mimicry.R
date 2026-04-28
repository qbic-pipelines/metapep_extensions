# -------------------------------
# Setup
# -------------------------------
my_lib <- "~/R/library"
dir.create(my_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(my_lib)

cran_mirror <- "https://cloud.r-project.org/"

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE, lib.loc = my_lib)) {
    install.packages(pkg, lib = my_lib, repos = cran_mirror)
    library(pkg, character.only = TRUE, lib.loc = my_lib)
  }
}

install_if_missing("ggplot2")
install_if_missing("dplyr")
install_if_missing("ggExtra")



# -------------------------------
# Load mimicry data (shared for all)
# -------------------------------
mimicry_files <- sort(list.files(pattern = "^entities_self_mimicry_ratios\\.allele_[0-9]+\\.tsv$"))

# -------------------------------
# Load samplesheet (for names)
# -------------------------------
samplesheet <- read.table("samplesheet.csv", header = TRUE, sep = ",")
alleles_str <- samplesheet$alleles[1]
alleles <- unlist(strsplit(alleles_str, " +"))

cat("Alleles found in samplesheet:\n")
print(alleles)

# -------------------------------
# Loop over all allele files
# -------------------------------
allele_files <- list.files(pattern = "^entity_binding_ratios\\.allele_[0-9]+\\.tsv$")

cat("Found allele files:\n")
print(allele_files)

if (length(allele_files) != length(alleles)) {
  warning("Number of allele files does not match number of alleles in samplesheet!")
}

for (i in seq_len(min(length(allele_files), length(alleles)))) {
  file <- allele_files[i]
  allele <- alleles[i]
  mimicry_file <- mimicry_files[i]
  mimicry_data <- read.table(mimicry_file, header = TRUE, sep = "\t")
  cat("  Using mimicry file:", mimicry_file, "\n")

  cat("\nProcessing:", file, "→", allele, "\n")

  # Load binding data
  binding_data <- read.table(file, header = TRUE, sep = "\t")
  binding_data$entity_id <- seq_len(nrow(binding_data))

  # Combine with mimicry data
  if (nrow(binding_data) != nrow(mimicry_data)) {
    warning(paste("Row count mismatch in", file, "- skipping."))
    next
  }
  merged_data <- merge(binding_data, mimicry_data, by = "entity_id")

  # Calculate correlation
  cor_values <- merged_data %>%
    group_by(condition_name) %>%
    summarise(
    r = cor(self_mimicry_ratio, binding_rate, method = "pearson"),
    p_value = cor.test(self_mimicry_ratio, binding_rate)$p.value,
    r_label = sprintf("r = %.2f", r),
    p_label = sprintf("p = %.4f", p_value),
    .groups = "drop"
  )



  # -------------------------------
  # Plot
  # -------------------------------
  p <- ggplot(merged_data, aes(x = self_mimicry_ratio, y = binding_rate, color = condition_name)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    ylab("Entity-wise binding rate") +
    xlab("Self-mimicry ratio") +
    ggtitle(paste("Binding vs. Self-mimicry (", allele, ")", sep = "")) +
    scale_color_brewer(palette = "Dark2") +
    theme_classic() +
    theme(
      legend.position = "top",
      plot.title = element_text(
        hjust = 0.5,
        size = 16,
        face = "bold",
        margin = margin(t = 10, b = 10)
      ),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.border = element_blank()
    )

  # Add correlation description
  p <- p +
    geom_text(
      data = cor_values,
      aes(x = Inf, y = Inf, label = paste0(r_label, ", ", p_label), color = condition_name),
      hjust = 1.1,
      vjust = seq(2, 0.5 + 1.5 * (nrow(cor_values)), length.out = nrow(cor_values)) + 1,
      size = 4,
      show.legend = FALSE
    ) +
    annotate("text", x = Inf, y = Inf, label = "Correlation:", hjust = 1.1, vjust = 1.2, size = 4.2, fontface = "bold")


  # Attach Boxplots to axes
  p_box <- ggMarginal(p, type = "boxplot", groupColour = TRUE, groupFill = TRUE)
  p_box$layout$t[1] <- 1
  p_box$layout$r[1] <- max(p_box$layout$r)
  
  # Store results
  clean_name <- gsub("\\*", "", allele)
  clean_name <- gsub(":", "", clean_name)
  out_file <- paste0(clean_name, "_entity_binding_ratio_self_mimicry.png")

  ggsave(out_file, p_box, width = 8, height = 6, dpi = 300)
  cat("Saved plot to", out_file, "\n")


}

cat("\nAll allele plots finished!\n")
