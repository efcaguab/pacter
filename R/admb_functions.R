# Functions adapted from the code written by Ben Bolker for R2admb
compile_admb <- function (fn, safe = FALSE, re = FALSE, verbose = FALSE, admb_errors = c("stop", 
                                                                                         "warn", "ignore")) 
{
  admb_errors <- match.arg(admb_errors)
  if (!file.exists(fn2 <- paste(fn, "tpl", sep = "."))) 
    stop("can't find TPL file ", fn2)
  test <- try(system("admb", intern = TRUE), silent = TRUE)
  if (inherits(test, "try-error")) 
    stop("base admb command failed: run setup_admb(), or check ADMB installation")
  args <- ""
  if (re) 
    args <- "-r"
  if (safe) 
    args <- paste(args, "-s")
  if (verbose) 
    message("compiling with args: '", args, "' ...\n")
  res0 <- system(paste("admb", args, fn, " 2>", paste(fn, ".cout", sep = "")),
                 intern = TRUE)
  coutfile <- readLines(paste(fn, ".cout", sep = ""))
  if (verbose) {
    message("compile output:\n", res0, "\n")
    message("compile log:\n")
    message(coutfile, sep = "\n")
  }
  admb_warnerr_index <- grepl("warning|error", coutfile)
  csplit <- split(coutfile, head(c(0, cumsum(admb_warnerr_index)), 
                                 -1))
  wchunks <- which(sapply(lapply(csplit, grep, pattern = "warning"), 
                          length) > 0)
  echunks <- which(sapply(lapply(csplit, grep, pattern = "error"), 
                          length) > 0)
  if (length(wchunks) > 0) {
    if (!verbose) {
      admb_warnings <- paste("from ADMB:", unlist(csplit[wchunks]))
      sapply(admb_warnings, warning)
    }
    csplit <- csplit[-wchunks]
  }
  Sys.chmod(fn, mode = "0755")
  if (length(echunks) > 0) {
    comperrmsg <- "errors detected in compilation: run with verbose=TRUE to view"
    if (admb_errors == "stop") 
      stop(comperrmsg)
    else if (admb_errors == "warn") 
      warning(comperrmsg)
  }
}


run_admb <- function (fn, verbose = FALSE, mcmc = FALSE, mcmc.opts = mcmc.control(), 
                      profile = FALSE, extra.args = "", admb_errors = c("stop", 
                                                                        "warn", "ignore")) 
{
  admb_errors <- match.arg(admb_errors)
  args <- ""
  if (mcmc) {
    if (is.null(mcmc.opts$mcmcpars)) 
      stop("you must specify at least one parameter in 'mcmc.opts$mcmcpars' (see ?mcmc.control)")
    args <- paste(args, mcmc.args(mcmc.opts))
  }
  if (profile) 
    args <- paste(args, "-lprof")
  if (!missing(extra.args)) {
    args <- paste(args, extra.args)
  }
  if (verbose) {
    message("running compiled executable with args: '", args, 
            "'...\n")
    std_pipe <- "| tee"
  } else {
    std_pipe = ">"
  }
    
  outfn <- paste(fn, "out", sep = ".")
  if (.Platform$OS.type == "windows") {
    cmdname <- paste(fn, ".exe", sep = "")
    shellcmd <- shell
  }
  else {
    cmdname <- paste("./", fn, sep = "")
    shellcmd <- system
  }
  if (!file.exists(cmdname)) 
    stop("executable ", cmdname, " not found: did you forget to compile it?")
  res <- shellcmd(paste(cmdname, args, std_pipe, outfn), intern = FALSE)
  outfile <- readLines(paste(fn, ".out", sep = ""))
  if (mcmc) {
    mcinfofile <- file(paste(fn, "mcinfo", sep = "."), "w")
    mctab <- unlist(mapply(function(x, y) {
      c(paste("# ", x), if (is.null(y)) "" else paste(y, 
                                                      collapse = " "))
    }, names(mcmc.opts), mcmc.opts))
    writeLines(mctab, paste(fn, "mcinfo", sep = "."))
  }
  if (verbose) {
    message("Run output:\n", res, "\n", sep = "\n")
    # message(outfile, "\n", sep = "\n")
  }
  if (length(grep("^Error", outfile) > 0)) {
    runerrmsg <- "errors detected in ADMB run: run with verbose=TRUE to view"
    if (admb_errors == "stop") 
      stop(runerrmsg)
    else if (admb_errors == "warn") 
      warning(runerrmsg)
  }
  invisible(outfile)
}


read_admb <- function (fn, verbose = FALSE, profile = FALSE, mcmc = FALSE, 
                       mcmc.opts = NULL, admbOut = NULL, checkterm = TRUE) 
{
  tpldat <- read_tpl(fn)
  if (verbose) 
    message("reading output ...\n")
  parfn <- paste(fn, "par", sep = ".")
  if (!file.exists(parfn)) 
    stop("couldn't find parameter file ", parfn)
  L <- c(list(fn = fn, txt = admbOut), R2admb:::read_pars(fn))
  if (mcmc) {
    if (file.exists(paste(fn, "hst", sep = "."))) {
      L <- c(L, list(hist = R2admb:::read_hst(fn)))
    }
    if (is.null(mcmc.opts)) {
      mcinfofile <- paste(fn, "mcinfo", sep = ".")
      if (file.exists(mcinfofile)) {
        w <- readLines(mcinfofile)
        wnames <- gsub("^# +", "", w[seq(1, length(w), 
                                         by = 2)])
        wvals <- as.list(w[seq(2, length(w), by = 2)])
        wvals[c(1, 2, 5)] <- as.numeric(wvals[c(1, 2, 
                                                5)])
        wvals[3:4] <- as.logical(wvals[3:4])
        wvals[[6]] <- strsplit(wvals[[6]], " ")[[1]]
        names(wvals) <- wnames
        mcmc.opts <- wvals
      }
      else warning("having difficulty retrieving MCMC info, will try to continue anyway")
    }
    sdinfo <- read.table(paste(fn, "std", sep = "."), skip = 1)
    pnames <- R2admb:::rep_pars(sdinfo[, 2])
    sdreportvars <- as.character(tpldat$info$sdnums$vname)
    pnames <- setdiff(pnames, sdreportvars)
    if (is.null(mcmc.opts) || mcmc.opts[["mcsave"]] > 0) {
      mcpnames <- pnames[!grepl("^r_", pnames)]
      L$mcmc <- read_psv(fn, names = mcpnames)
      attr(L$mcmc, "mcpar") <- c(1, mcmc.opts[["mcmc"]], 
                                 mcmc.opts[["mcsave"]])
    }
  }
  if (profile) {
    if (!is.null(tpldat$info$raneff)) {
      stop("something's wrong -- profiling is not implemented for random effects models")
    }
    profpars <- tpldat$info$profparms$vname
    L$prof <- lapply(profpars, R2admb:::read_plt)
    names(L$prof) <- gsub("p_", "", profpars)
  }
  if (checkterm) {
    v <- with(L, vcov[seq(npar), seq(npar)])
    ev <- try(eigen(solve(v))$value, silent = TRUE)
    L$eratio <- if (inherits(ev, "try-error")) 
      NA
    else min(ev)/max(ev)
  }
  class(L) <- "admb"
  L
}

read_tpl <- function (f) 
{
  r <- readLines(paste(f, "tpl", sep = "."))
  if (!nzchar(r[1])) 
    r <- r[-(1:rle(!nzchar(r))$lengths[1])]
  secStart <- grep("^ *[A-Z]+(_[A-Z]+)+$", r)
  calcLines <- grep("_CALCS *$", r)
  secStart <- setdiff(secStart, calcLines)
  if (length(secStart) == 0) 
    stop("tpl file must contain at least one section (capitalized header)")
  if (secStart[1] != 1) {
    secStart <- c(1, secStart)
  }
  nsec <- length(secStart)
  L <- c(secStart[-1], length(r) + 1) - secStart
  sec <- rep(1:nsec, L)
  splsec <- split(r, sec)
  splnames <- sapply(splsec, "[", 1)
  names(splsec) <- gsub("_.+", "", splnames)
  splsec_proc <- lapply(splsec, R2admb:::drop_calcs)
  L1 <- L2 <- NULL
  pp <- splsec_proc$PARAMETER
  pp <- pp[!grepl("^ *!!", pp)]
  if (!is.null(pp)) {
    pp <- proc_var(pp, maxlen = 7)
    type <- 1
    if (!is.null(pp)) {
      L1 <- with(pp, list(inits = pp[grep("^init", type), 
                                     ], raneff = pp[grep("^random", type), ], sdnums = pp[grep("^sdreport_number", 
                                                                                               type), ], sdvecs = pp[grep("^sdreport_vector", 
                                                                                                                          type), ], other = pp[grep("^init|random|sdreport", 
                                                                                                                                                    type, invert = TRUE), ], profparms = pp[grep("^likeprof", 
                                                                                                                                                                                                 type), ]))
    }
  }
  pp <- splsec_proc$DATA
  if (!is.null(pp)) {
    pp <- proc_var(pp, maxlen = 7)
    L2 <- with(pp, list(data = pp[grep("^init", type), ]))
  }
  L <- c(L1, L2)
  L <- L[!sapply(L, is.null)]
  list(secs = splsec, info = L[sapply(L, nrow) > 0])
}

proc_var <- function (s, drop.first = TRUE, maxlen) 
{
  if (drop.first) 
    s <- s[-1]
  calclocs <- grep("_CALCS", s)
  if (length(calclocs) > 0) {
    droplines <- unlist(apply(matrix(-calclocs, ncol = 2, 
                                     byrow = TRUE), 1, function(x) seq(x[1], x[2])))
    s <- s[droplines]
  }
  s2 <- gsub("^[ \\\t]*", "", gsub("[;]*[ \\\t]*$", "", R2admb:::strip_comments(s)))
  s2 <- s2[nchar(s2) > 0]
  s2 <- s2[!grepl("+[ \\\t]*!!", s2)]
  words <- strsplit(s2, " ")
  words <- lapply(words, function(x) x[x != ""])
  type <- sapply(words, "[[", 1)
  rest <- sapply(words, "[[", 2)
  rest2 <- strsplit(gsub("[(),]", " ", rest), " ")
  vname <- sapply(rest2, "[[", 1)
  if (length(rest2) == 0) 
    ret <- NULL
  else {
    maxlen0 <- max(sapply(rest2, length))
    if (missing(maxlen)) 
      maxlen <- maxlen0
    else maxlen <- pmax(maxlen, maxlen0)
    opts <- t(sapply(rest2, function(w) {
      c(w[-1], rep(NA, maxlen + 1 - length(w)))
    }))
    ret <- data.frame(type, vname, opts, stringsAsFactors = FALSE)
  }
  ret
}
