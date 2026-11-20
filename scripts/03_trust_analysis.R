suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

MIN_DEPTH  <- 10                         
DROPOUT    <- 0.10                       
TOLERANCE  <- 3                          
THRESHOLDS <- c(1, 5, 10, 20, 50, 100)

WHO <- c("B.1.1.7" = "Alpha", "B.1.351" = "Beta", "P.1" = "Gamma", "B.1.617.2" = "Delta",
         "B.1.427" = "Epsilon", "B.1.429" = "Epsilon", "B.1.525" = "Eta", "B.1.526" = "Iota",
         "B.1.617.1" = "Kappa", "P.2" = "Zeta")

dir.create("figures", showWarnings = FALSE)

# Samples
samples <- read_tsv("data/samples.tsv", col_names = c("sample", "lineage"), show_col_types = FALSE) |>
  mutate(who = WHO[lineage],
         label = paste0(sample, "  ", who, " (", lineage, ")")) |>
  arrange(who, lineage, sample) |>
  mutate(label = factor(label, levels = rev(label)))      # grouped by lineage in figures

# A. Reads per amplicon (iVar trim logs)
read_ivar <- function(s) {
  lines <- read_lines(file.path("results", s, "ivar_trim.log"))
  m <- str_match(lines, "^nCoV-2019_(\\d+)_(LEFT|RIGHT)(?:_alt\\d+)?\\s+(\\d+)\\s*$")
  tibble(sample = s, amplicon = as.integer(m[, 2]), reads = as.numeric(m[, 4])) |>
    filter(!is.na(amplicon))
}

amplicons <- map_dfr(samples$sample, read_ivar) |>
  group_by(sample, amplicon) |>
  summarise(reads = sum(reads), .groups = "drop") |>      
  group_by(sample) |>
  mutate(relative = reads / median(reads),
         log2_relative = log2(pmax(relative, 1 / 64)),     
         dropout = relative < DROPOUT) |>
  ungroup() |>
  left_join(samples, by = "sample")

dropouts <- amplicons |> filter(dropout) |>
  select(sample, who, lineage, amplicon, reads, relative) |> arrange(who, sample, amplicon)

# B. Positions below the depth threshold, by cause
bed <- read_tsv("ref/nCoV-2019.primer.bed",
                col_names = c("chrom", "start", "end", "name", "pool", "strand"),
                col_types = "ciicic")
scheme_start <- max(bed$end[str_detect(bed$name, "_1_LEFT")]) + 1
scheme_end   <- min(bed$start[str_detect(bed$name, "_98_RIGHT")])

read_depth <- function(s) {
  read_tsv(file.path("results", s, paste0(s, ".depth.txt")),
           col_names = c("chrom", "pos", "depth"), col_types = "cii") |>
    mutate(sample = s)
}
depth <- map_dfr(samples$sample, read_depth)

deletions <- read_tsv("results/nextclade.tsv", show_col_types = FALSE) |>
  select(sample = seqName, deletions) |>
  filter(!is.na(deletions)) |>
  separate_longer_delim(deletions, ",") |>                 
  separate_wider_delim(deletions, "-", names = c("from", "to"), too_few = "align_start") |>
  mutate(from = as.integer(from),
         to = as.integer(coalesce(to, as.character(from))),   
         lo = from - TOLERANCE, hi = to + TOLERANCE)

in_deletion <- depth |>
  select(sample, pos) |>
  inner_join(deletions, by = join_by(sample, between(pos, lo, hi))) |>
  distinct(sample, pos) |>
  mutate(is_deletion = TRUE)

depth <- depth |>
  left_join(in_deletion, by = c("sample", "pos")) |>
  mutate(is_deletion = coalesce(is_deletion, FALSE),
         outside = pos < scheme_start | pos > scheme_end)

missing_by_cause <- depth |>
  filter(depth < MIN_DEPTH) |>
  mutate(cause = case_when(outside ~ "Outside amplicon scheme",
                           is_deletion ~ "Viral deletion",
                           TRUE ~ "Missing data")) |>
  count(sample, cause, name = "positions") |>
  complete(sample = samples$sample, cause, fill = list(positions = 0)) |>
  left_join(samples, by = "sample")

# C. Threshold sensitivity
thresholds <- map_dfr(THRESHOLDS, function(t) {
  depth |>
    filter(!is_deletion) |>
    group_by(sample) |>
    summarise(positions_n = sum(depth < t), .groups = "drop") |>
    mutate(min_depth = t)
}) |>
  left_join(samples, by = "sample")


check <- thresholds |> filter(min_depth == MIN_DEPTH) |>
  left_join(read_tsv("results/summary.tsv", show_col_types = FALSE) |> select(sample, n_bases),
            by = "sample") |>
  select(sample, who, predicted_N = positions_n, consensus_N = n_bases)

# Save tables
write_tsv(amplicons |> select(sample, who, lineage, amplicon, reads, relative, dropout),
          "results/trust_amplicons.tsv")
write_tsv(dropouts, "results/trust_dropouts.tsv")
write_tsv(missing_by_cause |> select(sample, who, lineage, cause, positions),
          "results/trust_missing_by_cause.tsv")
write_tsv(thresholds |> select(sample, who, lineage, min_depth, positions_n),
          "results/trust_thresholds.tsv")

# Figure
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank()))

pA <- ggplot(amplicons, aes(amplicon, label, fill = log2_relative)) +
  geom_tile() +
  geom_point(data = filter(amplicons, dropout), shape = 4, size = 1.2, colour = "black") +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#2a78d6", midpoint = 0,
                       limits = c(-6, 2), oob = scales::squish,
                       name = "Reads vs sample\nmedian (log2)") +
  scale_x_continuous(breaks = c(1, seq(10, 90, 10), 98), expand = c(0, 0)) +
  labs(title = "A. Reads per amplicon (x = dropout, below 10% of the sample median)",
       x = "ARTIC V3 amplicon", y = NULL)

pB <- ggplot(missing_by_cause, aes(positions, label, fill = cause)) +
  geom_col() +
  scale_fill_manual(values = c("Missing data" = "#c0392b", "Viral deletion" = "#2a78d6",
                               "Outside amplicon scheme" = "#b4b2a9"), name = NULL) +
  labs(title = paste0("B. Positions with depth below ", MIN_DEPTH, ", by cause"),
       x = "Genome positions", y = NULL) +
  theme(legend.position = "bottom")

pC <- ggplot(thresholds, aes(min_depth, positions_n, group = sample, colour = who)) +
  geom_line(alpha = 0.8) + geom_point(size = 1) +
  scale_x_log10(breaks = THRESHOLDS) + scale_y_log10() +
  labs(title = "C. Untrusted positions vs minimum depth",
       x = "Minimum depth to call a base", y = "Positions that would be N", colour = NULL) +
  theme(legend.position = "bottom")

fig <- pA / (pB | pC) + plot_layout(heights = c(1.1, 1))
ggsave("figures/fig1_trust.png", fig, width = 13, height = 11, dpi = 300, bg = "white")
ggsave("figures/fig1_trust.pdf", fig, width = 13, height = 11)

# Print the key results
options(width = 150)
cat("\n=== Amplicon dropouts (below 10% of the sample median) ===\n")
print(dropouts |> group_by(who, sample) |>
        summarise(amplicons = paste(amplicon, collapse = ", "), .groups = "drop"), n = Inf)
cat("\n=== Positions below depth", MIN_DEPTH, "by cause ===\n")
print(missing_by_cause |> select(sample, who, cause, positions) |>
        pivot_wider(names_from = cause, values_from = positions), n = Inf)
cat("\n=== Check: predicted vs real Ns at depth", MIN_DEPTH, "===\n")
print(check, n = Inf)
cat("\nSaved figures/fig1_trust.png and results/trust_*.tsv\n")
