# lib_logging.R — canonical log_msg() shared across production/diagnostic
# entry-point scripts.
#
# Consolidates five previously copy-pasted local `log_msg()` definitions
# that had drifted into three distinct behaviors:
#   - run_microclimate_site.R, run_colonization.R: message() to console AND
#     append the same stamped line to a per-run log file
#   - height_resolution_experiment.R, resolution_diagnostics.R: message()
#     to console only, no file
#   - tests/test.R: variadic args, full date-time stamp (not just HH:MM:SS),
#     writes via cat() to stdout instead of message(), always flushes
#
# make_log_msg() returns a closure configured for one of those behaviors, so
# each call site below sources this file once and binds `log_msg` to a
# closure carrying its own log_file/stamp-format/flush choice — the actual
# call sites (log_msg("...")) are unchanged, only the five local function
# *definitions* are gone.
make_log_msg <- function(log_file = NULL, timestamp_fmt = c("short", "long"),
                          flush_output = FALSE) {
  timestamp_fmt <- match.arg(timestamp_fmt)

  function(...) {
    # Default sep = " " on purpose, matching tests/test.R's original
    # paste(..., collapse = ""); the "short" call sites all pass a single
    # already-formatted string, so sep is a no-op for them.
    msg <- paste(..., collapse = "")
    stamped <- if (timestamp_fmt == "long") {
      sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
    } else {
      paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
    }

    if (timestamp_fmt == "long") {
      cat(stamped, "\n", sep = "", file = stdout())
    } else {
      message(stamped)
    }

    if (!is.null(log_file)) cat(stamped, "\n", file = log_file, append = TRUE)
    if (flush_output) flush(stdout())
    invisible(stamped)
  }
}
