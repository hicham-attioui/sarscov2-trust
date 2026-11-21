# Does the pipeline's QC catch bad samples and keep good ones?
#
# Applies the QC rules of the benchmark paper (Xiaoli et al. 2022, PeerJ, Table 2) to both
# datasets, processed by main.nf:
#   VOI/VOC   16 samples that passed the authors' QC   -> every sample should PASS
#   Failed QC 24 samples that failed the authors' QC   -> every sample should FAIL
#
# Rules (Illumina), as in the paper, with two adaptations stated below:
#   Mean base quality (Phred, before filtering)  > 25
#   Reads mapped to the reference                > 65%
#   Mean depth                                   > 25   (paper: mean depth > 10 AND mean coverage > 25)
#   Positions with depth < 10                    < 3,000
#   Consensus length                             > 29,400
#   Ns in the consensus                          < 10%
#   Nextclade frameshift check = good            (replaces the paper's VADR check, not in this pipeline)
# A missing value (for example no reads at all) counts as a failure.
#
# Input : results_nf/ and results_failedqc/ (main.nf outputs), data/samples.tsv,
#         data/samples_failedqc.tsv
# Output: results/qc_metrics.tsv, figures/fig2_qc_rules.png/.pdf
#
# Run from the repository root, in the R environment:
#   conda activate sarscov2-r
#   Rscript scripts/04_qc_check.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

DATASETS <- tribble(
  ~dataset,                          ~dir,               ~sheet,
  "VOI/VOC (expected: pass)",        "results_nf",       "data/samples.tsv",
  "Failed QC (expected: fail)",      "results_failedqc", "data/samples_failedqc.tsv"
)

# ---------------------------------------------------------------------------
# Metrics per sample
# ---------------------------------------------------------------------------
read_sheet <- function(path) {
  x <- read_tsv(path, col_names = FALSE, show_col_types = FALSE)
  if (ncol(x) >= 3) tibble(sample = x$X1, description = x$X3)
  else tibble(sample = x$X1, description = paste("lineage", x$X2))
}

sample_metrics <- function(dir, s) {
  d <- file.path(dir, s)

  fp <- fromJSON(file.path(d, "fastp.json"))
  qualities <- c(fp$read1_before_filtering$quality_curves$mean,
                 fp$read2_before_filtering$quality_curves$mean)

  mapped_line <- str_subset(read_lines(file.path(d, "flagstat.txt")), "primary mapped")
  pct_mapped <- suppressWarnings(as.numeric(str_match(mapped_line, "\\(([0-9.]+)%")[, 2]))

  depth <- read_tsv(file.path(d, paste0(s, ".depth.txt")),
                    col_names = c("chrom", "pos", "depth"), col_types = "cii")

  consensus <- paste(read_lines(file.path(d, paste0(s, ".consensus.fa")))[-1], collapse = "")
  n_len <- nchar(consensus)

  tibble(sample              = s,
         reads               = fp$summary$before_filtering$total_reads,
         mean_phred          = if (length(qualities)) mean(qualities) else NA_real_,
         pct_mapped          = if (length(pct_mapped)) pct_mapped[1] else NA_real_,
         mean_depth          = if (nrow(depth)) mean(depth$depth) else 0,
         low_depth_positions = if (nrow(depth)) sum(depth$depth < 10) else 29903,
         consensus_length    = n_len,
         pct_N               = if (n_len) 100 * str_count(consensus, "N") / n_len else 100)
}

read_nextclade <- function(dir) {
  read_tsv(file.path(dir, "nextclade.tsv"), show_col_types = FALSE) |>
    transmute(sample = seqName, pango = Nextclade_pango,
              nextclade_qc = `qc.overallStatus`, frameshift_qc = `qc.frameShifts.status`)
}

metrics <- pmap_dfr(DATASETS, function(dataset, dir, sheet) {
  sheet_df <- read_sheet(sheet)
  map_dfr(sheet_df$sample, \(s) sample_metrics(dir, s)) |>
    left_join(sheet_df, by = "sample") |>
    left_join(read_nextclade(dir), by = "sample") |>
    mutate(dataset = dataset)
})

# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------
RULES <- c("Mean base quality > 25", "Reads mapped > 65%", "Mean depth > 25",
           "Positions with depth < 10: under 3,000", "Consensus length > 29,400",
           "Ns < 10%", "Nextclade: no frameshift")

checks <- metrics |>
  transmute(dataset, sample, description,
            `Mean base quality > 25`                  = mean_phred > 25,
            `Reads mapped > 65%`                      = pct_mapped > 65,
            `Mean depth > 25`                         = mean_depth > 25,
            `Positions with depth < 10: under 3,000`  = low_depth_positions < 3000,
            `Consensus length > 29,400`               = consensus_length > 29400,
            `Ns < 10%`                                = pct_N < 10,
            `Nextclade: no frameshift`                = frameshift_qc == "good") |>
  mutate(across(all_of(RULES), \(x) coalesce(x, FALSE)))     # missing value = fail

checks$pass_all <- apply(checks[RULES], 1, all)
write_tsv(left_join(metrics, checks |> select(sample, pass_all), by = "sample"), "results/qc_metrics.tsv")

# ---------------------------------------------------------------------------
# Figure: one row per sample, one column per rule
# ---------------------------------------------------------------------------
long <- checks |>
  pivot_longer(all_of(RULES), names_to = "rule", values_to = "passed") |>
  mutate(rule = factor(rule, levels = RULES),
         label = paste0(sample, "  ", str_trunc(description, 55)),
         outcome = if_else(passed, "Pass", "Fail"))

fig <- ggplot(long, aes(rule, label, fill = outcome)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c(Pass = "#1baf7a", Fail = "#c0392b"), name = NULL) +
  facet_grid(dataset ~ ., scales = "free_y", space = "free_y") +
  labs(title = "QC rules of the benchmark paper applied to the pipeline outputs",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"), legend.position = "top")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig2_qc_rules.png", fig, width = 11, height = 10, dpi = 300, bg = "white")
ggsave("figures/fig2_qc_rules.pdf", fig, width = 11, height = 10)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
options(width = 160)
cat("\n=== Overall result per dataset ===\n")
print(checks |> group_by(dataset) |>
        summarise(samples = n(), passed_all = sum(pass_all), failed = sum(!pass_all)))

cat("\n=== Failed QC samples: which rules caught them ===\n")
caught <- checks |> filter(str_detect(dataset, "Failed")) |>
  pivot_longer(all_of(RULES), names_to = "rule", values_to = "passed") |>
  filter(!passed) |>
  group_by(sample, description) |>
  summarise(failed_rules = paste(rule, collapse = "; "), .groups = "drop")
missed <- checks |> filter(str_detect(dataset, "Failed"), pass_all) |> select(sample, description)
print(caught, n = Inf)
cat("\nFailed-QC samples NOT caught by any rule:\n")
print(missed, n = Inf)

cat("\n=== VOI/VOC samples wrongly rejected ===\n")
print(checks |> filter(str_detect(dataset, "VOI"), !pass_all) |> select(sample, description), n = Inf)

cat("\n=== Key metrics ===\n")
print(metrics |> select(dataset, sample, reads, mean_phred, pct_mapped, mean_depth,
                        low_depth_positions, pct_N, frameshift_qc) |>
        mutate(across(where(is.numeric), \(x) round(x, 1))), n = Inf)
cat("\nSaved results/qc_metrics.tsv and figures/fig2_qc_rules.png\n")
