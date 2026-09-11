rxTest({
  test_that("a homogeneous event table draws its own residuals per subject", {
    skip_on_cran()

    ## A homogeneous event table keeps ONE representative record set and expands
    ## it to every id in the group at solve time.  The residual draw was sized
    ## from the un-expanded record count, so it came up short and every subject
    ## past the first group's worth reused the last drawn row.
    .m <- rxode2({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      cp <- linCmt()
      cp2 <- cp * (1 + prop.err)
    })
    .ev <- et(et(amt=100, id=1:6), seq(0, 24, by=8))
    .om <- lotri::lotri(eta.ka + eta.cl ~ c(0.1, 0.01, 0.1))
    .sg <- lotri::lotri(prop.err ~ 0.1)
    .p <- c(tka=0.45, tcl=1, tv=3.45)

    .s <- function(.e) {
      withr::with_seed(42, {
        rxSetSeed(1234)
        rxSolve(.m, .e, params=.p, omega=.om, sigma=.sg, nStud=3, dfSub=10)
      })
    }

    .hom <- .s(.ev)
    ## as.data.frame() expands the groups, so this is the same regimen without
    ## the homogeneous representation
    .exp <- .s(as.data.frame(.ev))

    expect_equal(as.data.frame(.hom), as.data.frame(.exp))

    ## 6 subjects x 4 observations x 3 studies, not 4 x 3
    expect_equal(nrow(attr(class(.hom), ".rxode2.env")$.sigma), 6 * 4 * 3)

    ## every subject gets its own residuals
    .r <- as.data.frame(.hom)
    .r <- .r[.r$cp > 0, ]
    .r$resid <- .r$cp2 / .r$cp - 1
    expect_equal(length(unique(.r$resid)), nrow(.r))
  })
})
