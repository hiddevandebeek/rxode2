rxTest({

  test_that(".rxIsReservedName tracks the parser's reserved names", {
    expect_true(all(.rxIsReservedName(c("t", "time", "tlast", "newind", "NEWIND",
                                        "rxFlag", "amt", "mixnum", "mixest",
                                        "mixunif", "M_PI", "M_E", "M_LN10",
                                        "pi", "NA", "NaN", "Inf"))))
    # the reserved variables the parser matches case-insensitively
    expect_true(all(.rxIsReservedName(c("Time", "TIME", "AMT"))))
    # ordinary model variables are not reserved
    expect_false(any(.rxIsReservedName(c("tka", "ka", "cl", "v", "eta.ka",
                                         "add.sd", "wt", "T"))))
    expect_equal(.rxIsReservedName(character(0)), logical(0))
    expect_equal(.rxIsReservedName(NA_character_), NA)
  })

  # a model with no reserved names in it, used as the piping base below
  .base <- function() {
    ini({
      tka <- log(1.5)
      tcl <- log(1)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl * center
      cp <- center
      cp ~ add(add.sd)
    })
  }

  test_that("reserved variables are not promoted when appending", {
    for (v in c("t", "time", "tlast", "newind", "rxFlag", "M_PI", "pi",
                "NA", "NaN", "Inf")) {
      ui <- .base()
      expect_error(
        ui <- do.call(model, list(ui, str2lang(paste0("cp2 <- cp * ", v)),
                                  append = quote(cp))),
        NA, info = v)
      expect_equal(ui$iniDf$name, c("tka", "tcl", "add.sd"), info = v)
      expect_false(v %in% ui$allCovs, info = v)
    }
  })

  test_that("reserved variables are not promoted when prepending", {
    ui <- .base()
    ui <- do.call(model, list(ui, str2lang("f <- t * 2"), append = FALSE))
    expect_equal(ui$iniDf$name, c("tka", "tcl", "add.sd"))
    expect_false("t" %in% ui$allCovs)
  })

  test_that("a piped model using t still solves", {
    ui <- .base()
    ui <- do.call(model, list(ui, quote(cp2 <- cp * exp(-t / 10)),
                              append = quote(cp)))
    s <- rxSolve(ui, et(amt = 100, ii = 24, addl = 2) |> et(seq(0, 72, by = 12)),
                 params = c(tka = log(1.5), tcl = log(1), add.sd = 0),
                 returnType = "data.frame")
    expect_true(all(is.finite(s$cp2)))
    expect_equal(s$cp2, s$cp * exp(-s$time / 10))
  })

  test_that("non-reserved variables are still promoted when appending", {
    ui <- .base()
    ui <- do.call(model, list(ui, quote(cp2 <- cp * exp(tf)), append = quote(cp)))
    expect_true("tf" %in% ui$iniDf$name)
    # a covariate-looking name still becomes a covariate rather than a parameter
    ui <- .base()
    ui <- do.call(model, list(ui, quote(cp2 <- cp * wt), append = quote(cp)))
    expect_false("wt" %in% ui$iniDf$name)
    expect_true("wt" %in% ui$allCovs)
  })

})
