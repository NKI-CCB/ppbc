# Functions for outcome analyses in PPBC

read_outcome <- function(fn) {
  metadata <- purrr::map(rlang::set_names(c('sample_data', 'patient_data')), function (sheet) {
    readxl::read_excel(fn, sheet = sheet, col_types = 'text')
  })
  # Only select the data fields we know about to prevent hidden bugs and catch errors here
  metadata$sample_data <- metadata$sample_data %>%
    dplyr::filter(sample_type == 'slide', !is.na(sample_ID), !is.na(batch_HALO)) %>%
    dplyr::transmute(
      sample_ID,
      patient_ID,
      # Use parse_factor to check that the levels are what we expect and catch errors early
      batch_HALO = readr::parse_factor(batch_HALO, levels = paste0('batch ', 1:7)) %>%
        stringr::str_replace(' ', ''),
      experimental_platform = readr::parse_factor(experimental_platform,
                                                  levels = c("MPIF26", "MPIF27")))
  metadata$patient_data <- metadata$patient_data %>%
    dplyr::semi_join(metadata$sample_data, by = 'patient_ID') %>%
    dplyr::transmute(
      patient_ID,
      time_OS_months = readr::parse_double(time_OS_months),
      death = readr::parse_logical(death),
      time_DRS_months = readr::parse_double(time_DRS_months),
      distant_recurrence = readr::parse_logical(distant_recurrence),
      study_group = readr::parse_factor(study_group,
                                        levels = c("non_prbc", "prbc", "ppbc_inv", "ppbc_lac")),
      stage = readr::parse_factor(stage,
                                  levels = c("stage I", "stage II", "stage III", "stage IV")),
      grade = readr::parse_factor(grade, levels = c("grade I", "grade II", "grade III")),
    )
  full_join(metadata$sample_data, metadata$patient_data, by = 'patient_ID')
}
