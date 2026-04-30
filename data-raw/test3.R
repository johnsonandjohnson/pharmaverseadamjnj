
library(envsetup)
library(tern)
library(dplyr)
library(rtables)
library(junco)

################################################################################
# Define script level parameters:
################################################################################

tblid <- "TSFVIT04"
popfl <- "SAFFL"
trtvar <- "TRT01A"
ctrl_grp <- "Placebo"

# flag to indclude timepoint rows
timepoint_rows <- TRUE

# specify in the order you want to print on table
selparamcd <- c("ORTHYP", "SYSBPO", "DIABPO" )

selparamcdN <- tibble(PARAMCD = selparamcd, PARAMCDN = seq_along(selparamcd))

### Per email June 12: DAS/SDS confirmed to NOT restrict to on-treatment values

################################################################################
# Process Data:
################################################################################

adsl <- adsl1 %>%
  filter(.data[[popfl]] == "Y") %>%
  mutate(
    !!rlang::sym(trtvar) := factor(
      .data[[trtvar]],
      levels = c(
        "Xanomeline Low Dose",
        "Xanomeline High Dose",
        "Placebo"
      )
    )
  ) %>%
  select(USUBJID, all_of(c(popfl, trtvar)), SEX, AGEGR1, RACE, ETHNIC)


adsl$colspan_trt <- factor(
  ifelse(adsl[[trtvar]] == ctrl_grp, " ", "Active Study Agent"),
  levels = c("Active Study Agent", " ")
)


colspan_trt_map <- create_colspan_map(
  adsl,
  non_active_grp = ctrl_grp,
  non_active_grp_span_lbl = " ",
  active_grp_span_lbl = "Active Study Agent",
  colspan_var = "colspan_trt",
  trt_var = trtvar
)
ref_path <- c("colspan_trt", " ", trtvar, ctrl_grp)

### N is the number of subjects with postbaseline orthostatic measurements and without orthostatic hypotension at baseline
## do not use TRTEMFL as filter, as this will only select AVALC = Y records per definition of TRTEMFL
## instead : start from post-baseline records and retain one record per subject
## for those subjects with both Y and N records, keep the Y record
## for those subjects with only Y records or only N records, keep one Y, N record respectively
filtered_advs <- advs1 %>%
  filter(PARAMCD %in% selparamcd & APOBLFL == "Y") %>%
  filter(PARAMCD %in% selparamcd) %>%
  mutate(
    AVALC = case_when(
      (PARAMCD == "SYSBPO" | PARAMCD == "DIABPO") & CRIT1FL == "Y" ~ "Y",
      (PARAMCD == "SYSBPO" | PARAMCD == "DIABPO") &
        (CRIT1FL == "N" |
          CRIT1FL == NA) ~
        "N",
      PARAMCD == "ORTHYP" ~ AVALC
    ),
    PARAM = case_when(
      PARAMCD == "SYSBPO" ~ "SBP (STD-SUP) <-20",
      PARAMCD == "DIABPO" ~ "DBP (STD-SUP) <-10",
      PARAMCD == "ORTHYP" ~ PARAM
    )
  ) %>%
  select(
    STUDYID,
    USUBJID,
    PARAMCD,
    PARAM,
    AVALC,
    AVISIT,
    APOBLFL,
    TRTEMFL,
    ONTRTFL
  ) %>%
  inner_join(adsl) %>%
  ### ensure to keep only 1 result per subject, keep N only in case no Y was observed
  arrange(USUBJID, PARAMCD, AVALC) %>%
  group_by(USUBJID, PARAMCD) %>%
  mutate(navalc = n_distinct(AVALC)) %>%
  filter(!(navalc > 1 & AVALC == "N")) %>%
  ## only keep one record
  slice_head(n = 1) %>%
  ungroup()

#### remove subjects abnormal for "ORTHYP" at baseline
bl_abn_orthyp <- advs1 %>%
  filter(PARAMCD == "ORTHYP" & ABLFL == "Y" & AVALC == "Y")

### actually remove the subjects with AVALC = Y for ORTHYP
### N is the number of subjects with postbaseline orthostatic measurements and without orthostatic hypotension at baseline
filtered_advs <- filtered_advs %>%
  filter(!(USUBJID %in% unique(bl_abn_orthyp$USUBJID)))

### get sorting as per order in selparamcdN
selparamcdN <- selparamcdN %>%
  left_join(unique(filtered_advs %>% select(PARAMCD, PARAM))) %>%
  arrange(PARAMCDN)

param_levels <- unique(as.character(selparamcdN$PARAM))

filtered_advs$PARAM <- factor(
  as.character(filtered_advs$PARAM),
  levels = param_levels
)


### Mapping for AVALC
### alternative approach to retrieve from metadata iso dataset
xlabel_map <- unique(filtered_advs %>% select(PARAM, PARAMCD)) %>%
  mutate(
    value = "Y",
    label = as.character(PARAM)
  ) %>%
  mutate(
    label = case_when(
      label == "Orthostatic Hypotension" ~
        "Total number of subjects with orthostatic hypotension",
      TRUE ~ label
    )
  )

map <- unique(filtered_advs %>% select(PARAM, PARAMCD)) %>%
  mutate(
         PARAMCD=as.character(PARAMCD),
         PARAM = as.character(PARAM)
         ) %>% cross_join(test <- filtered_advs %>% select(AVISIT) %>% distinct() %>% mutate(AVISIT=as.character(AVISIT)))


################################################################################
# Define layout and build table:
################################################################################

extra_args_rr <- list(
  method = "wald",
  denom = "n_df",
  ref_path = ref_path,
  .stats = c("count_unique_fraction")
)

extra_args_rr1 <- list(
  method = "wald",
  denom = "n_df",
  denom_by = var,
  ref_path = ref_path,
  .stats = c("denom", "count_unique_fraction")
)

lyt0 <- basic_table(
  show_colcounts = TRUE,
  colcount_format = "N=xx"
) %>%
  split_cols_by(
    "colspan_trt",
    split_fun = trim_levels_to_map(map = colspan_trt_map)
  ) %>%
  split_cols_by(trtvar) %>%
  analyze(
    "AVALC",
    var_labels = "On-treatment",
    a_freq_j,
    show_labels = "visible",
    table_names = "AVALC_N",
    extra_args = list(denom = "n_df", .stats = c("n_df"))
  ) %>%
  split_rows_by(
    "PARAM",
    label_pos = "topleft",
    split_fun = drop_split_levels,
    child_labels = "hidden",
    split_label = "Orthostatic hypotension, n (%)"
  ) %>%
  # as in shell, do not show denom in count/denom (%)
  ### indent will be fixed to 1, will be updated later in post-processing
  analyze(
    "AVALC",
    a_freq_j,
    extra_args = append(
      extra_args_rr,
      list(
        val = c("Y"),
        label_map = xlabel_map
      )
    ),
    indent_mod = 2L,
    show_labels = "hidden"
  )

if(timepoint_rows){
  lyt <- lyt0 %>%
    split_rows_by(
    "AVISIT",
    label_pos = "topleft",
    split_fun = trim_levels_in_group("AVISIT"),
    child_labels = "visible",
    split_label = " ",
    section_div = " "
  ) %>%
    summarize_row_groups(
      "AVISIT",
      cfun = a_freq_j,
      extra_args = list(
        .stats = "n_df",
        label = "N",
        riskdiff = FALSE
      )
    ) %>%
    split_rows_by(
    "PARAM",
    label_pos = "topleft",
    split_fun = drop_split_levels,
    child_labels = "hidden",
    split_label = " "
  ) %>%
    analyze(
    "AVALC",
    a_freq_j,
    table_names = "AVALC_V1",
    extra_args = append(
      extra_args_rr,
      list(
        val = c("Y"),
        label_map = xlabel_map
      )
    ),
    indent_mod = 0L,
    show_labels = "hidden"
  )

} else{
  lyt <- lyt0
}

result <- build_table(lyt, filtered_advs, alt_counts_df = adsl)

