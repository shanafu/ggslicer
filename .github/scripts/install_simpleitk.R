# Install pre-built SimpleITK binary into the active R library (.uvr/library).

patch_macos_r_libs <- function(pkg_dir) {
  libs_dir <- file.path(pkg_dir, "libs")
  if (!dir.exists(libs_dir)) {
    return(invisible(NULL))
  }
  r_lib_dir <- R.home("lib")
  so_files <- list.files(
    libs_dir,
    pattern = "\\.(so|dylib)$",
    recursive = TRUE,
    full.names = TRUE
  )
  for (so in so_files) {
    Sys.chmod(so, mode = "0755")
    lines <- suppressWarnings(system2("otool", c("-L", shQuote(so)), stdout = TRUE))
    if (length(lines) <= 1) next
    changed <- FALSE
    for (line in lines[-1]) {
      old_dep <- strsplit(trimws(line), "\\s+")[[1]][1]
      if (!nzchar(old_dep)) next
      if (!grepl("R\\.framework|libR\\.dylib", old_dep)) next
      new_dep <- file.path(r_lib_dir, basename(old_dep))
      if (old_dep == new_dep || !file.exists(new_dep)) next
      status <- system2(
        "install_name_tool",
        c("-change", shQuote(old_dep), shQuote(new_dep), shQuote(so))
      )
      if (identical(status, 0L)) {
        changed <- TRUE
      }
    }
    if (changed) {
      system2("codesign", c("--force", "--sign", "-", shQuote(so)))
    }
  }
  invisible(NULL)
}

lib <- .libPaths()[1]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("SimpleITK", lib.loc = lib, quietly = TRUE)) {
  ios <- SimpleITK::ImageFileReader()$GetRegisteredImageIOs()
  if ("MINCImageIO" %in% ios) {
    message("SimpleITK is already installed in ", lib)
    quit(save = "no", status = 0)
  }
}

sitk_version <- Sys.getenv("SIMPLEITK_VERSION", unset = "2.5.7")
r_minor <- strsplit(R.version$minor, "\\.")[[1]][1]
r_ver <- paste0(R.version$major, ".", r_minor)

sysname <- Sys.info()[["sysname"]]
machine <- Sys.info()[["machine"]]
is_x86_64 <- grepl("x86[_-]64|amd64", machine, ignore.case = TRUE)
is_arm64 <- grepl("arm64|aarch64", machine, ignore.case = TRUE)

if (sysname == "Linux" && is_x86_64) {
  platform <- "linux-x86_64"
  ext <- "tar.gz"
} else if (sysname == "Darwin" && is_arm64) {
  platform <- "macos-arm64"
  ext <- "tgz"
} else if (sysname == "Darwin" && is_x86_64) {
  platform <- "macos-x86_64"
  ext <- "tgz"
} else if (sysname == "Windows" && is_x86_64) {
  platform <- "windows-x86_64"
  ext <- "zip"
} else {
  stop(
    "Unsupported platform for pre-built SimpleITK binaries: ",
    sysname, "/", machine,
    call. = FALSE
  )
}

asset <- sprintf("SimpleITK_%s_R%s_%s.%s", sitk_version, r_ver, platform, ext)
url <- sprintf(
  "https://github.com/SimpleITK/SimpleITKRInstaller/releases/download/v%s/%s",
  sitk_version, asset
)

tmp_base <- Sys.getenv("TMPDIR", unset = "")
if (!nzchar(tmp_base) && .Platform$OS.type == "unix" && dir.exists("/var/tmp")) {
  tmp_base <- "/var/tmp"
  Sys.setenv(TMPDIR = tmp_base)
}
if (!nzchar(tmp_base)) {
  tmp_base <- tempdir()
}

work_dir <- tempfile(pattern = "sitk-install-", tmpdir = tmp_base)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

old_timeout <- getOption("timeout")
on.exit(options(timeout = old_timeout), add = TRUE)
options(timeout = max(600, old_timeout))

archive_path <- file.path(work_dir, sprintf("SimpleITK_%s.%s", sitk_version, ext))

message("Downloading ", url)
utils::download.file(
  url,
  archive_path,
  mode = "wb",
  headers = c("User-Agent" = "ggslicer-ci/1.0 R")
)

if (sysname == "Windows") {
  utils::install.packages(archive_path, repos = NULL, type = "win.binary", lib = lib)
} else if (sysname == "Darwin") {
  pkg_dest <- file.path(lib, "SimpleITK")
  if (dir.exists(pkg_dest)) {
    unlink(pkg_dest, recursive = TRUE)
  }
  rc <- utils::untar(archive_path, exdir = lib)
  if (!identical(rc, 0L)) {
    stop("Failed to extract SimpleITK archive (exit ", rc, ").", call. = FALSE)
  }
  patch_macos_r_libs(pkg_dest)
} else {
  utils::install.packages(archive_path, repos = NULL, type = "source", lib = lib)
}

library(SimpleITK, lib.loc = lib)
ios <- SimpleITK::ImageFileReader()$GetRegisteredImageIOs()
if (!("MINCImageIO" %in% ios)) {
  stop("Installed SimpleITK binary does not include MINCImageIO support.", call. = FALSE)
}
message(
  "Installed SimpleITK ",
  as.character(utils::packageVersion("SimpleITK", lib.loc = lib)),
  " in ", lib
)
