
myMax <- function (..., na.rm = FALSE)
{
  elts <- list(...)
  if (length(elts) == 0L)
    stop("no arguments")
  if (all(vapply(elts, function(x) is.atomic(x) && !is.object(x),NA))) {
    mmm <- pmax(na.rm, ...)
  }
  else {
    mmm <- elts[[1L]]
    has.na <- FALSE
    for (e in seq_along(elts)[-1]) {
      each <- elts[[e]]
      l1 <- length(each)
      l2 <- length(mmm)
      if (l2 < l1) {
        if (l2 && l1%%l2)
          warning("an argument will be fractionally recycled")
        mmm <- rep(mmm, length.out = l1)
      }
      else if (l1 && l1 < l2) {
        if (l2%%l1)
          warning("an argument will be fractionally recycled")
        each <- rep(each, length.out = l2)
      }
      nammm <- is.na(mmm)
      naeach <- is.na(each)
      if (has.na || (has.na <- any(nammm*naeach))) {
        mmm[nammm] <- each[nammm]
        each[naeach] <- mmm[naeach]
      }
      change <- mmm < each
      change <- change & !is.na(change)
      mmm[change] <- each[change]
      if (has.na && !na.rm)
        mmm[nammm | naeach] <- NA
    }
  }
  mmm
}
