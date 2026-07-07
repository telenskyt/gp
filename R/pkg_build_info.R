# custom-generated with ChatGPT: https://chatgpt.com/c/6a4cf4e4-e398-83ed-b8d8-e8785bf44076
#
# it should cover different build paths: RStudio package build, install_github, CRAN, ...
# see the files in ../tools/ that generate the ../inst/BUILDINFO

pkg_build_info <- function() 
{
	pkg <- this_package()

	read_buildinfo <- function() {
		path <- system.file("BUILDINFO", package = pkg, mustWork = FALSE)
		if (!nzchar(path)) {
			return(list())
		}
		as.list(read.dcf(path)[1, , drop = TRUE])
	}

	first_nonempty <- function(...) {
		xs <- list(...)
		for (x in xs) {
			if (!is.null(x) && length(x) > 0L && !is.na(x[[1]]) && nzchar(x[[1]])) {
				return(as.character(x[[1]]))
			}
		}
		NA_character_
	}

	as_bool <- function(x) {
		if (is.null(x) || is.na(x) || !nzchar(x)) {
			return(NA)
		}
		identical(toupper(x), "TRUE")
	}

	desc <- as.list(utils::packageDescription(pkg))
	bi <- read_buildinfo()

	sha <- first_nonempty(bi$GitSHA, desc$RemoteSha, desc$GithubSHA1)
	dirty <- as_bool(first_nonempty(bi$GitDirty, NA_character_))

	list(
		package = first_nonempty(bi$Package, desc$Package, pkg),
		package_version = first_nonempty(bi$Version, desc$Version),
		git_sha = sha,
		git_ref = first_nonempty(bi$GitRef, desc$RemoteRef, desc$GithubRef),
		git_dirty = dirty,
		build_info_source = first_nonempty(bi$BuildInfoSource, NA_character_),
		build_time_utc = first_nonempty(bi$BuildTimeUTC, NA_character_),
		r_version = first_nonempty(bi$RVersion, NA_character_),
		provenance_exact = !is.na(sha) && !isTRUE(dirty)
	)
}


this_package <- function() 
{
	pkg <- utils::packageName()

	if (is.null(pkg) || length(pkg) == 0L || is.na(pkg) || !nzchar(pkg)) {
		pkg <- "gp"
	}

	pkg
}
