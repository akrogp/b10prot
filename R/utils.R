# R/utils.R

#' Check for required columns in a data frame
#'
#' @description
#' This helper function verifies that all required columns exist in a data frame
#' or tibble. If any are missing, it throws an error.
#'
#' @param data A data frame or tibble.
#' @param required A character vector of required column names.
#'
#' @return The same tibble if all required columns exist.
#' @keywords internal
check_required_cols <- function(data, required) {
    missing <- setdiff(required, names(data))
    if (length(missing) > 0) {
        rlang::abort(
            message = paste0(
                "Missing required columns: ",
                paste(missing, collapse = ", ")
            ),
            class = "tda_missing_cols"
        )
    }
    data
}

