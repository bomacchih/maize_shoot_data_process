# --- Packages ---
library(dplyr)
library(broom)
library(emmeans)
library(car)
library(rstatix)

# =========================
# Data: use your in-memory object
# =========================
# Expect a data.frame named SV_phenotyping_20251010 with columns:
#   part (trait name), line (A10/GRF6/GRF8/GRF10), measurement (numeric)
df <- SV_phenotyping_20251010 %>%
  mutate(
    part = factor(part),
    line = factor(line)
  ) %>%
  filter(line %in% c("A10","GRF6","GRF8","GRF10")) %>%
  droplevels()

# Traits to test
focus_parts <- c("SAM_ratio","Height","leaf_length","leaf_width","SAM_height","SAM_width")
alpha <- 0.05

# --- Helper to build compact-letter display from pairwise p-values ---
make_letters <- function(levels_vec, pairs_df, from = c("gh","tukey"), alpha = 0.05) {
  from <- match.arg(from)
  if (!requireNamespace("multcompView", quietly = TRUE)) {
    return(tibble(line = levels_vec, .group = NA_character_))
  }
  if (from == "gh") {
    # rstatix::games_howell_test output
    pairs_df <- pairs_df %>%
      transmute(g1 = as.character(group1), g2 = as.character(group2), p = p.adj)
  } else {
    # emmeans Tukey summary: contrast "A - B", p.value
    pairs_df <- pairs_df %>%
      tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
      transmute(g1, g2, p = p.value)
  }
  pairs_df <- pairs_df %>%
    mutate(a = pmin(g1,g2), b = pmax(g1,g2), key = paste(a,b,sep="-")) %>%
    distinct(key, .keep_all = TRUE)
  pvec <- pairs_df$p; names(pvec) <- pairs_df$key
  letters <- multcompView::multcompLetters(pvec, threshold = alpha)$Letters
  tibble(line = names(letters), .group = unname(letters))
}

# --- Main analyzer for one trait ---
analyze_part <- function(dat, part_name, alpha = 0.05) {
  sub <- dat %>% filter(part == part_name, is.finite(measurement)) %>% droplevels()
  if (n_distinct(sub$line) < 2L) {
    return(list(
      anova = tibble(part = part_name, method = NA, stat = NA, p.value = NA_real_,
                     df1 = NA, df2 = NA, note = "Fewer than 2 groups"),
      posthoc_pairwise = tibble(), posthoc_vsA10 = tibble(),
      grouping = tibble(part = part_name, line = levels(sub$line), .group = NA_character_)
    ))
  }

  # Set A10 as reference (for Dunnett)
  if ("A10" %in% levels(sub$line)) sub$line <- stats::relevel(sub$line, ref = "A10")

  fit_aov <- aov(measurement ~ line, data = sub)

  # Assumptions
  sh  <- tryCatch(shapiro.test(residuals(fit_aov)), error = function(e) NULL)
  lev <- tryCatch(car::leveneTest(measurement ~ line, data = sub), error = function(e) NULL)
  lev_p <- if (!is.null(lev)) lev[1,"Pr(>F)"] else NA_real_
  hetero <- is.finite(lev_p) && lev_p < alpha

  post_pair <- tibble(); post_vsA10 <- tibble()
  grouping_tbl <- tibble(part = part_name, line = levels(sub$line), .group = NA_character_)

  if (hetero) {
    # Welch + Games–Howell (robust)
    welch <- oneway.test(measurement ~ line, data = sub, var.equal = FALSE)
    anova_row <- tibble(
      part   = part_name,
      method = "Welch ANOVA (heteroscedastic)",
      stat   = unname(welch$statistic),
      p.value= unname(welch$p.value),
      df1    = unname(welch$parameter[1]),
      df2    = unname(welch$parameter[2]),
      note   = paste0("Levene p=", signif(lev_p,3),
                      if (!is.null(sh)) paste0("; Shapiro p=", signif(sh$p.value,3)) else "")
    )
    # Compute GH for letters regardless; show pairwise table only if omnibus sig
    gh_all <- rstatix::games_howell_test(sub, measurement ~ line)
    if (anova_row$p.value < alpha) {
      post_pair <- gh_all %>%
        transmute(
          part = part_name, method = "Games-Howell",
          contrast = paste(group1, "vs", group2),
          estimate = estimate, conf.low = conf.low, conf.high = conf.high,
          p.value  = p.adj, p.signif = p.adj.signif
        )
      post_vsA10 <- post_pair %>% filter(grepl("^A10 vs | vs A10$", contrast))
    }
    grouping_tbl <- make_letters(levels(sub$line), gh_all, from = "gh", alpha = alpha) %>%
      mutate(part = part_name, .before = 1)

  } else {
    # Classical ANOVA + Tukey; Dunnett vs A10
    aov_tidy <- broom::tidy(fit_aov)
    aov_line <- aov_tidy %>% filter(term == "line")
    anova_row <- tibble(
      part   = part_name,
      method = "Classical ANOVA",
      stat   = aov_line$statistic,
      p.value= aov_line$p.value,
      df1    = aov_line$df,
      df2    = df.residual(fit_aov),
      note   = paste0("Levene p=", signif(lev_p,3),
                      if (!is.null(sh)) paste0("; Shapiro p=", signif(sh$p.value,3)) else "")
    )

    em <- emmeans(fit_aov, ~ line)

    # Tukey pairs for grouping (compute always), report only if omnibus sig
    tuk_s <- summary(pairs(em, adjust = "tukey")) %>% as.data.frame()
    grouping_tbl <- make_letters(levels(sub$line), tuk_s, from = "tukey", alpha = alpha) %>%
      mutate(part = part_name, .before = 1)

    if (anova_row$p.value < alpha) {
      tuk_c <- confint(pairs(em, adjust = "tukey")) %>% as.data.frame()
      post_pair <- dplyr::left_join(
        tuk_s, tuk_c %>% dplyr::select(contrast, lower.CL, upper.CL), by = "contrast"
      ) %>%
        transmute(
          part = part_name, method = "Tukey HSD",
          contrast, estimate, conf.low = lower.CL, conf.high = upper.CL,
          p.value,
          p.signif = case_when(
            is.na(p.value) ~ "",
            p.value < 0.001 ~ "***",
            p.value < 0.01  ~ "**",
            p.value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
        )

      # Dunnett vs A10
      if ("A10" %in% levels(sub$line)) {
        dun <- contrast(em, method = "trt.vs.ctrl", ref = "A10", adjust = "dunnett")
        dun_s <- summary(dun) %>% as.data.frame()
        dun_c <- confint(dun)   %>% as.data.frame()
        post_vsA10 <- dplyr::left_join(
          dun_s, dun_c %>% dplyr::select(contrast, lower.CL, upper.CL), by = "contrast"
        ) %>%
          mutate(treatment = sub(" - A10", "", contrast)) %>%
          transmute(
            part = part_name, method = "Dunnett vs A10",
            contrast = paste0(treatment, " vs A10"),
            estimate, conf.low = lower.CL, conf.high = upper.CL,
            p.value,
            p.signif = case_when(
              is.na(p.value) ~ "",
              p.value < 0.001 ~ "***",
              p.value < 0.01  ~ "**",
              p.value < 0.05  ~ "*",
              TRUE ~ "ns"
            )
          )
      }
    }
  }

  list(anova = anova_row, posthoc_pairwise = post_pair, posthoc_vsA10 = post_vsA10,
       grouping = grouping_tbl)
}

# --- Run all tests and print outputs ---
res_list <- lapply(focus_parts, function(pn) analyze_part(df, pn, alpha = alpha))

anova_table      <- bind_rows(lapply(res_list, `[[`, "anova"))
posthoc_pairwise <- bind_rows(lapply(res_list, `[[`, "posthoc_pairwise"))
posthoc_vsA10    <- bind_rows(lapply(res_list, `[[`, "posthoc_vsA10"))
grouping_table   <- bind_rows(lapply(res_list, `[[`, "grouping"))

cat("\n=== Omnibus test (ANOVA or Welch) ===\n")
print(anova_table)

cat("\n=== Pairwise post hoc (reported only if omnibus significant) ===\n")
print(posthoc_pairwise %>% arrange(part, p.value))

cat("\n=== Post hoc vs A10 (Dunnett or GH filtered vs A10) ===\n")
print(posthoc_vsA10 %>% arrange(part, p.value))

cat("\n=== Grouping letters (compact letter display) ===\n")
print(grouping_table %>% arrange(part, .group, line))

grouping_table %>%
  arrange(part, .group, line) %>%
  print(n = Inf, width = Inf)




++++++++++++++++++++++++++++++++++++++++++
“ab” is from the compact letter display (CLD). It means that group’s mean is not significantly different from any group that has “a” or “b”. In other words, it overlaps both letter sets.
Quick rules of thumb:


Groups that share at least one letter are not significantly different at your α (e.g., 0.05).


Groups that share no letters are significantly different.


Letters aren’t ordered and don’t indicate effect size or direction—just overlap of non-significance sets.


Example interpretation
If you see:
LineLetterA10aGRF6abGRF8bGRF10b


A10 vs GRF6: share “a” → not significant


GRF6 vs GRF8/GRF10: share “b” → not significant


A10 vs GRF8/GRF10: no shared letter (a vs b) → significant


So “ab” often acts as a bridge group between “a” and “b”.
