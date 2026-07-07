# Write the exact build info to be stored in the model objects
# 
# https://chatgpt.com/c/6a4cf4e4-e398-83ed-b8d8-e8785bf44076

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
buildinfo_path <- file.path(root, "inst", "BUILDINFO")

empty_to_na <- function(x) {
	if (is.null(x) || length(x) == 0L || is.na(x[[1]]) || !nzchar(x[[1]])) {
		return(NA_character_)
	}
	as.character(x[[1]])
}

first_nonempty <- function(...) {
	xs <- list(...)
	for (x in xs) {
		x <- empty_to_na(x)
		if (!is.na(x)) {
			return(x)
		}
	}
	NA_character_
}

read_dcf1 <- function(path) {
	if (!file.exists(path)) {
		return(list())
	}
	as.list(read.dcf(path)[1, , drop = TRUE])
}

field <- function(x, name) {
	empty_to_na(x[[name]])
}

env <- function(name) {
	empty_to_na(Sys.getenv(name, unset = NA_character_))
}

git <- function(args) {
	out <- tryCatch(
		system2("git", c("-C", root, args), stdout = TRUE, stderr = FALSE),
		error = function(e) character()
	)
	if (!is.null(attr(out, "status"))) {
		return(NA_character_)
	}
	empty_to_na(out)
}

git_dirty <- function() {
	out <- tryCatch(
		system2("git", c("-C", root, "status", "--porcelain"), stdout = TRUE, stderr = FALSE),
		error = function(e) NA_character_
	)
	if (!is.null(attr(out, "status"))) {
		return(NA_character_)
	}
	if (length(out) == 0L) {
		return("FALSE")
	}
	as.character(any(nzchar(out)))
}

desc <- read_dcf1(file.path(root, "DESCRIPTION"))
old <- read_dcf1(buildinfo_path)

sha <- first_nonempty(
	env("PKG_GIT_SHA"),
	env("GITHUB_SHA"),
	env("CI_COMMIT_SHA"),
	env("BITBUCKET_COMMIT"),
	field(desc, "RemoteSha"),
	field(desc, "GithubSHA1"),
	git(c("rev-parse", "HEAD")),
	field(old, "GitSHA")
)

ref <- first_nonempty(
	env("PKG_GIT_REF"),
	env("GITHUB_REF_NAME"),
	env("GITHUB_HEAD_REF"),
	env("CI_COMMIT_REF_NAME"),
	env("BITBUCKET_BRANCH"),
	field(desc, "RemoteRef"),
	field(desc, "GithubRef"),
	git(c("branch", "--show-current")),
	field(old, "GitRef")
)

dirty <- first_nonempty(
	env("PKG_GIT_DIRTY"),
	git_dirty(),
	if (!is.na(field(desc, "RemoteSha"))) "FALSE" else NA_character_,
	field(old, "GitDirty")
)

source <- if (!is.na(env("PKG_GIT_SHA")) || !is.na(env("GITHUB_SHA")) || !is.na(env("CI_COMMIT_SHA"))) {
	"env"
} else if (!is.na(field(desc, "RemoteSha"))) {
	"remotes"
} else if (!is.na(git(c("rev-parse", "HEAD")))) {
	"git"
} else if (!is.na(field(old, "GitSHA"))) {
	"existing-buildinfo"
} else {
	"unknown"
}

info <- as.data.frame(
	as.list(c(
		Package = field(desc, "Package"),
		Version = field(desc, "Version"),
		GitSHA = sha,
		GitRef = ref,
		GitDirty = dirty,
		BuildInfoSource = source,
		BuildTimeUTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
		RVersion = as.character(getRversion())
	)),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

dir.create(dirname(buildinfo_path), showWarnings = FALSE, recursive = TRUE)
write.dcf(info, buildinfo_path)

