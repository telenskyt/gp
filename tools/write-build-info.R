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

has_git_checkout <- function() {
	file.exists(file.path(root, ".git"))
}

git <- function(args) {
	if (!has_git_checkout()) {
		return(NA_character_)
	}

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
	if (!has_git_checkout()) {
		return(NA_character_)
	}

	out <- tryCatch(
		system2("git", c("-C", root, "status", "--porcelain"), stdout = TRUE, stderr = FALSE),
		error = function(e) character()
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

env_sha <- first_nonempty(
	env("PKG_GIT_SHA"),
	env("GITHUB_SHA"),
	env("CI_COMMIT_SHA"),
	env("BITBUCKET_COMMIT")
)

env_ref <- first_nonempty(
	env("PKG_GIT_REF"),
	env("GITHUB_REF_NAME"),
	env("GITHUB_HEAD_REF"),
	env("CI_COMMIT_REF_NAME"),
	env("BITBUCKET_BRANCH")
)

remote_sha <- first_nonempty(
	field(desc, "RemoteSha"),
	field(desc, "GithubSHA1")
)

remote_ref <- first_nonempty(
	field(desc, "RemoteRef"),
	field(desc, "GithubRef")
)

git_sha <- NA_character_
git_ref <- NA_character_
git_is_dirty <- NA_character_

if (is.na(env_sha) && is.na(remote_sha) && has_git_checkout()) {
	git_sha <- git(c("rev-parse", "HEAD"))
	git_ref <- git(c("branch", "--show-current"))
	git_is_dirty <- git_dirty()
}

sha <- first_nonempty(
	env_sha,
	remote_sha,
	git_sha,
	field(old, "GitSHA")
)

ref <- first_nonempty(
	env_ref,
	remote_ref,
	git_ref,
	field(old, "GitRef")
)

dirty <- first_nonempty(
	env("PKG_GIT_DIRTY"),
	if (!is.na(env_sha)) "FALSE" else NA_character_,
	if (!is.na(remote_sha)) "FALSE" else NA_character_,
	git_is_dirty,
	field(old, "GitDirty")
)

build_info_source <- if (!is.na(env_sha)) {
	"env"
} else if (!is.na(remote_sha)) {
	"remotes"
} else if (!is.na(git_sha)) {
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
		BuildInfoSource = build_info_source,
		BuildTimeUTC = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
		RVersion = as.character(getRversion())
	)),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

dir.create(dirname(buildinfo_path), showWarnings = FALSE, recursive = TRUE)
write.dcf(info, buildinfo_path)
