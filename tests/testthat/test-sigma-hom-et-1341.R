rxTest({
  # A homogeneous event table stores one representative subject per group, so
  # the residual draw used to be sized from that single subject (issue #1341).
  .m1341 <- rxode2({
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    y  <- cp + err
  })
  .p1341 <- c(ka = 1.5, cl = 2.7, v = 31.5)
  .s1341 <- lotri(err ~ 0.25)
  .t1341 <- c(0.25, 1, 4, 12, 24)

  test_that("sigma is simulated for every subject that comes from et(id=) (#1341)", {
    .ev <- et(amt = 320) %>% et(.t1341) %>% et(id = 1:6)
    .b <- withr::with_seed(1, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, .ev, .p1341, sigma = .s1341,
                                             addDosing = FALSE)))
    })
    expect_equal(nrow(.b), 6L * length(.t1341))
    # one draw per observation row, not one subject's worth recycled
    expect_equal(length(unique(.b$y - .b$cp)), nrow(.b))

    # ... and it matches the nSub= path exactly
    .a <- withr::with_seed(1, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, et(amt = 320) %>% et(.t1341),
                                             .p1341, sigma = .s1341, nSub = 6,
                                             addDosing = FALSE)))
    })
    expect_equal(.b$y - .b$cp, .a$y - .a$cp)
  })

  test_that("sigma expands per group when the groups differ (#1341)", {
    .ev <- rbind(et(amt = 320) %>% et(.t1341) %>% et(id = 1:3),
                 et(amt = 100) %>% et(.t1341) %>% et(id = 4:5))
    .b <- withr::with_seed(4, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, .ev, .p1341, sigma = .s1341,
                                             addDosing = FALSE)))
    })
    expect_equal(length(unique(.b$y - .b$cp)), nrow(.b))
  })

  test_that("sigma expansion also covers the dosing rows kept by addDosing (#1341)", {
    .ev <- et(amt = 320) %>% et(.t1341) %>% et(id = 1:4)
    .b <- withr::with_seed(5, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, .ev, .p1341, sigma = .s1341,
                                             addDosing = TRUE)))
    })
    expect_equal(length(unique(.b$y - .b$cp)), nrow(.b))
  })

  # `curObs` -- the number of rows drawn -- comes from a different count for
  # each `addDosing`: `nobs2` (evid=0 only), `nobs` (observations, so evid=2
  # too) and `nall` (every record).  All three have to be expanded.
  test_that("every addDosing branch draws one sigma per output row (#1341)", {
    .ev <- et(amt = 320) %>% et(.t1341) %>% et(time = 2, evid = 2) %>% et(id = 1:4)
    for (.ad in list(NULL, FALSE, TRUE, NA)) {
      .b <- withr::with_seed(7, {
        suppressWarnings(as.data.frame(rxSolve(.m1341, .ev, .p1341, sigma = .s1341,
                                               addDosing = .ad)))
      })
      expect_equal(length(unique(.b$y - .b$cp)), nrow(.b),
                   label = paste0("addDosing=", if (is.null(.ad)) "NULL" else .ad))
    }
  })

  test_that("sigma expansion covers addl and steady-state doses (#1341)", {
    .addl <- et(time = 0, amt = 320, addl = 2, ii = 12) %>% et(.t1341) %>% et(id = 1:3)
    .b <- withr::with_seed(2, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, .addl, .p1341, sigma = .s1341,
                                             addDosing = TRUE)))
    })
    expect_equal(length(unique(.b$y - .b$cp)), nrow(.b))

    .ss <- et(time = 0, amt = 320, ii = 12, ss = 1) %>% et(.t1341) %>% et(id = 1:3)
    .b <- withr::with_seed(3, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, .ss, .p1341, sigma = .s1341,
                                             addDosing = FALSE)))
    })
    expect_equal(length(unique(.b$y - .b$cp)), nrow(.b))
  })

  test_that("sigma expansion discounts the evid=9 ini records (#1341)", {
    # no record at time 0, so etTrans adds one evid=9 ini record per subject;
    # those never take a residual draw and must not inflate the count
    .ev <- et(c(1, 2, 4, 8)) %>% et(id = 1:4)
    for (.ad in list(FALSE, TRUE)) {
      .b <- withr::with_seed(6, {
        suppressWarnings(as.data.frame(rxSolve(.m1341, .ev, .p1341, sigma = .s1341,
                                               addDosing = .ad)))
      })
      expect_equal(length(unique(.b$y - .b$cp)), nrow(.b),
                   label = paste0("addDosing=", .ad))
    }
  })

  test_that("a non-homogeneous multi-subject data set still draws one sigma per row (#1341)", {
    .ev <- rbind(data.frame(id = 1, time = c(0, 1, 2, 3), amt = c(320, NA, NA, NA),
                            evid = c(1, 0, 0, 0)),
                 data.frame(id = 2, time = c(0, 1, 2, 3, 8), amt = c(100, NA, NA, NA, NA),
                            evid = c(1, 0, 0, 0, 0)))
    .b <- withr::with_seed(3, {
      suppressWarnings(as.data.frame(rxSolve(.m1341, .ev, .p1341, sigma = .s1341,
                                             addDosing = FALSE)))
    })
    expect_equal(length(unique(.b$y - .b$cp)), nrow(.b))
  })
})
