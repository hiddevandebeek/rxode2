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

  test_that("a parameter named like a constant solves like any other name", {
    .mk <- function(v) {
      rxode2(sprintf(paste0("ka=exp(tka);\ncl=exp(tcl+%s);\n",
                            "d/dt(depot)=-ka*depot;\n",
                            "d/dt(center)=ka*depot-cl*center;\n"), v),
             calcSens = TRUE)
    }
    .ev <- et(amt = 100) %>% et(seq(0, 24, by = 2))
    .a <- rxSolve(.mk("e"), .ev, params = c(tka = 0.4, tcl = -0.1, e = 0.2),
                  returnType = "data.frame", atol = 1e-11, rtol = 1e-11)
    .b <- rxSolve(.mk("ee"), .ev, params = c(tka = 0.4, tcl = -0.1, ee = 0.2),
                  returnType = "data.frame", atol = 1e-11, rtol = 1e-11)
    names(.a) <- sub("_BY_e__", "_BY_ee__", names(.a), fixed = TRUE)
    expect_equal(sort(names(.a)), sort(names(.b)))
    expect_equal(as.matrix(.a[names(.b)]), as.matrix(.b))
  })

  test_that("jump event-sensitivities wrt a parameter named like a constant", {
    # .rxEventSensDExpr() tested the model-side name against symengine-side free
    # symbols, so the term was silently dropped rather than erroring
    .mk <- function(v) {
      sprintf(paste0("ka=exp(tka);\ncl=exp(tcl);\nf(depot)=expit(%s);\n",
                     "d/dt(depot)=-ka*depot;\n",
                     "d/dt(center)=ka*depot-cl*center;\n"), v)
    }
    .ev <- et(amt = 100) %>% et(seq(0, 24, by = 2))
    .p <- c(tka = 0.4, tcl = -0.1)
    .a <- rxSolve(rxode2(.mk("e"), calcSens = "e", eventSens = "jump"), .ev,
                  params = c(.p, e = 0.2), returnType = "data.frame",
                  atol = 1e-11, rtol = 1e-11)
    .b <- rxSolve(rxode2(.mk("ee"), calcSens = "ee", eventSens = "jump"), .ev,
                  params = c(.p, ee = 0.2), returnType = "data.frame",
                  atol = 1e-11, rtol = 1e-11)
    expect_equal(.a$rx__sens_center_BY_e__, .b$rx__sens_center_BY_ee__)
    # and it is the real derivative, not zero
    .mb <- rxode2(.mk("e"))
    .h <- 1e-6
    .fd <- (rxSolve(.mb, .ev, params = c(.p, e = 0.2 + .h),
                    returnType = "data.frame", atol = 1e-11, rtol = 1e-11)$center -
            rxSolve(.mb, .ev, params = c(.p, e = 0.2 - .h),
                    returnType = "data.frame", atol = 1e-11, rtol = 1e-11)$center) / (2 * .h)
    expect_lt(max(abs(.a$rx__sens_center_BY_e__ - .fd)), 1e-5)
    expect_gt(max(abs(.fd)), 1)
  })

  test_that("matExp()/indLin() rate constants keep a parameter named e", {
    # .multCollapse() re-parsed model-side text with symengine::S(), which read
    # the parameter `e` as Euler's number: k_p_q=exp(1) instead of k_p_q=e
    .n <- rxNorm(rxode2("d/dt(p)=-e*p;\nd/dt(q)=e*p-k*q;\n", indLin = TRUE))
    expect_true(grepl("k_p_q=e;", .n, fixed = TRUE))
    expect_false(grepl("exp(1)", .n, fixed = TRUE))
  })

  test_that("lag() of a variable named like a symengine constant round-trips", {
    # .rxToSELagOrLead()'s .vref() wraps the variable in symengine::S()
    expect_equal(rxNorm("b=lag(e,1);\nd/dt(center)=-b*center;\n"),
                 "b=lag(e,1);\nd/dt(center)=-b*center;\n")
    .s <- rxS(rxModelVars("b=lag(e,1);\nd/dt(center)=-b*center;\n"))
    expect_true(any(grepl("lag(e,1)", .s$..lhs, fixed = TRUE)))
  })

  test_that("the unshadowed name follows a later rxToSE() into the same env", {
    # doing the unshadowing once after rxS() loads the model would leave the
    # plain name stale as soon as the environment was extended
    .s <- rxS(rxModelVars("e=1;\nd/dt(center)=-e*center;\n"))
    expect_equal(as.character(.s$e), "1")
    invisible(rxToSE("e=2", envir = .s))
    expect_equal(as.character(.s$e), "2")
    expect_equal(as.character(.s$rx_SymPy_Res_e), "2")
  })

  test_that("adjoint sensitivities accept a parameter named like a constant", {
    .m <- rxS(rxGetModel("d/dt(depot)=-ka*depot;\nd/dt(center)=ka*depot-(e/v)*center;\n"),
              TRUE, promoteLinSens = FALSE)
    .v <- c("ka", "e", "v")
    invisible(.rxJacobian(.m, c(rxStateOde(.m), .v)))
    .adj <- .rxAdjoint(.m, .v, "center")
    expect_true(any(grepl("d/dt(rx__sens_center_BY_e__)=rx__adjLambda_center_center__*center/v",
                          .adj, fixed = TRUE)))
    expect_false(any(grepl("rx_SymPy_Res_", .adj, fixed = TRUE)))
    expect_false(any(grepl("2.718", .adj, fixed = TRUE)))
  })

  test_that("delay() terms resolve against a parameter named like a constant", {
    .m <- .rxode2({
      d/dt(cen) <- -e * cen + 0.1 * delay(cen, tau)
      tau <- 1.5
    })
    .t <- .rxDelayTerms(.m)
    expect_equal(.t$state, "cen")
    expect_equal(.t$tau, "tau")
    .s <- rxS(.m)
    .b <- .s$..ddt
    expect_equal(.b, "d/dt(cen)=-cen*e+0.1*delay(cen, 1.5)")
  })

  test_that("mu-referencing keeps a covariate parameter named like a constant", {
    .f <- function() {
      ini({
        tcl <- 1
        e <- 0.5
        add.sd <- 0.7
        eta.cl ~ 0.1
      })
      model({
        cl <- exp(tcl + e * WT + eta.cl)
        d/dt(center) <- -cl * center
        cp <- center
        cp ~ add(add.sd)
      })
    }
    .d <- .f()$muRefCovariateDataFrame
    expect_equal(.d$covariateParameter, "e")
    expect_equal(.d$covariate, "WT")
  })

  test_that(".rxSEres() mangles only the reserved names", {
    expect_equal(.rxSEres(c("e", "cl", "I", "eta.cl")),
                 c("rx_SymPy_Res_e", "cl", "rx_SymPy_Res_I", "eta.cl"))
    expect_equal(.rxSEres(character(0)), character(0))
    expect_equal(.rxSEres(names(.rxSEreserved)),
                 paste0("rx_SymPy_Res_", names(.rxSEreserved)))
  })
})
