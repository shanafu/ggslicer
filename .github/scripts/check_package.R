# Run R CMD check on ggslicer inside the uvr environment.

if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("TAR", unset = ""))) {
  Sys.setenv(TAR = "internal")
}

check_dir <- Sys.getenv("GGSLICER_CHECK_DIR", unset = "")
if (!nzchar(check_dir)) {
  tmp_base <- Sys.getenv("TMPDIR", unset = "")
  if (!nzchar(tmp_base) && .Platform$OS.type == "unix" && dir.exists("/var/tmp")) {
    tmp_base <- "/var/tmp"
    Sys.setenv(TMPDIR = tmp_base)
  }
  if (!nzchar(tmp_base)) {
    tmp_base <- tempdir()
  }
  check_dir <- tempfile(pattern = "ggslicer-check-", tmpdir = tmp_base)
  on.exit(unlink(check_dir, recursive = TRUE), add = TRUE)
}
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

res <- rcmdcheck::rcmdcheck(
  path = ".",
  args = c("--no-manual"),
  build_args = c("--no-manual"),
  error_on = "warning",
  check_dir = check_dir
)

print(res)
