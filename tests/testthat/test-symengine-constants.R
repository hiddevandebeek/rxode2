rxTest({
  # rxode2#1359: symengine's own constants (e, E, I, Catalan, GoldenRatio,
  # EulerGamma) are bound in the symengine environment, so a model variable of
  # one of those names used to read back as the constant and D() refused it.
  # NB: rxFromSE() poisons the next `$`/`[[` read of a symengine env, so every
  # Basic is captured BEFORE the first rxFromSE() in each test.
  .cnst <- c("e", "E", "I", "Catalan", "GoldenRatio", "EulerGamma")

  test_that("a parameter named like a symengine constant is not shadowed", {
    for (.v in .cnst) {
      .s <- rxS(rxModelVars(paste0(
        "cl=exp(tcl+", .v, ");\nd/dt(center)=-cl*center;\n")))
      .b <- .s[[.v]]
      .cl <- .s$cl
      expect_true(inherits(.b, "Basic"), info = .v)
      expect_equal(as.character(.b), paste0("rx_SymPy_Res_", .v), info = .v)
      expect_equal(rxFromSE(.b), .v, info = .v)
      # the model still reads the declared variable, not the constant
      expect_equal(rxFromSE(.cl), paste0("exp(", .v, "+tcl)"), info = .v)
    }
  })

  test_that("a state named like a symengine constant is not shadowed", {
    .s <- rxS(rxModelVars("d/dt(e)=-cl*e;\n"))
    .b <- .s$e
    expect_equal(as.character(.b), "rx_SymPy_Res_e")
    expect_equal(rxFromSE(.b), "e")
  })

  test_that("an lhs named like a symengine constant is not shadowed", {
    .s <- rxS(rxModelVars("e=exp(tcl);\nd/dt(center)=-e*center;\n"))
    .b <- .s$e
    expect_equal(rxFromSE(.b), "exp(tcl)")
  })

  test_that("a symengine constant is still bound when the model does not use it", {
    .s <- rxS(rxModelVars("cl=exp(tcl);\nd/dt(center)=-cl*center;\n"))
    expect_equal(as.character(.s$e), "2.71828182845905")
  })

  test_that("D() by a model variable named like a symengine constant works", {
    .s <- rxS(rxModelVars("cl=exp(tcl+e);\nd/dt(center)=-cl*center;\n"))
    .d <- with(.s, D(cl, e))
    expect_equal(rxFromSE(.d), "exp(e+tcl)")
  })

  test_that("sensitivities by a parameter named like a symengine constant", {
    .m <- rxode2({
      cl <- exp(tcl + e)
      d/dt(center) <- -cl * center
    }, calcSens = TRUE)
    .n <- rxNorm(.m)
    expect_true(grepl("rx__sens_center_BY_e__", .n, fixed = TRUE))
    expect_false(grepl("rx_SymPy_Res_", .n, fixed = TRUE))
    expect_true(grepl(paste0("d/dt(rx__sens_center_BY_e__)=-exp(e+tcl)*center-",
                             "exp(e+tcl)*rx__sens_center_BY_e__"),
                      .n, fixed = TRUE))
  })

  test_that(".rxSEres()/.rxSEsym() mangle only the reserved names", {
    expect_equal(.rxSEres(c("e", "cl", "I", "eta.cl")),
                 c("rx_SymPy_Res_e", "cl", "rx_SymPy_Res_I", "eta.cl"))
    expect_equal(as.character(.rxSEsym("e")), "rx_SymPy_Res_e")
    expect_equal(as.character(.rxSEsym("eta.cl")), "eta.cl")
  })
})
