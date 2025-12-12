# This file is part of the 'ralid' package for jamovi.

#' @rdname jamovi
#' @export
ralidClass <- R6::R6Class(
    "ralidClass",
    inherit = ralidBase,
    private = list(
        .run = function() {

            data <- self$data
            opts <- self$options
            res  <- self$results

            # Set visibility of optional tables/plots
            res$errorMetrics$setVisible(opts$showErrorMetrics)
            res$iccRetest$setVisible(opts$hasRetest)
            res$tostEquivPlot$setVisible(opts$doTOST)
            res$tostDecisionPlot$setVisible(opts$doTOST && opts$showTOSTDecisionCurve)
            res$withinSession$setVisible(opts$doWithinSession)

            vars <- opts$vars

            if (length(vars) < 2) {
                res$cccTable$setStatus('error')
                res$cccTable$setError('Select at least two numeric variables for Time 1.')
                return()
            }

            # Extract data for Time 1 vars and drop rows with any missing (for Time 1 analyses)
            sub <- data[vars]
            complete <- stats::complete.cases(sub)
            sub <- sub[complete, , drop = FALSE]

            if (nrow(sub) < 2) {
                res$cccTable$setStatus('error')
                res$cccTable$setError('Not enough complete cases for Time 1.')
                return()
            }

            # Helper: CCC and bootstrap CI ---------------------------------
            .ccc <- function(x, y) {
                mx <- mean(x)
                my <- mean(y)
                vx <- stats::var(x)
                vy <- stats::var(y)
                sxy <- stats::cov(x, y)

                r   <- sxy / sqrt(vx * vy)
                ccc <- (2 * sxy) / (vx + vy + (mx - my)^2)
                cb  <- if (isTRUE(all.equal(r, 0))) NA_real_ else ccc / r

                list(ccc = ccc, r = r, cb = cb)
            }

            .cccBootCI <- function(x, y, conf.level = 0.95, nboot = 1000) {
                n <- length(x)
                if (n < 3) {
                    return(c(NA_real_, NA_real_))
                }
                bvals <- numeric(nboot)
                for (b in seq_len(nboot)) {
                    idx <- sample.int(n, n, replace = TRUE)
                    bvals[b] <- .ccc(x[idx], y[idx])$ccc
                }
                alpha <- (1 - conf.level) / 2
                stats::quantile(bvals,
                                probs = c(alpha, 1 - alpha),
                                na.rm = TRUE,
                                names = FALSE)
            }

            # -----------------------------------------------------------------
            # CCC table (all unique pairs at Time 1)
            # -----------------------------------------------------------------
            nvars <- length(vars)
            for (i in seq_len(nvars - 1L)) {
                for (j in seq.int(i + 1L, nvars)) {
                    v1 <- vars[i]
                    v2 <- vars[j]
                    x  <- sub[[v1]]
                    y  <- sub[[v2]]

                    cc  <- .ccc(x, y)
                    ci  <- .cccBootCI(x, y,
                                      conf.level = opts$confLevel,
                                      nboot = opts$nBoot)

                    pairName <- paste0(v1, " vs ", v2)

                    res$cccTable$addRow(
                        rowKey = pairName,
                        values = list(
                            pair  = pairName,
                            ccc   = cc$ccc,
                            lower = ci[1],
                            upper = ci[2],
                            r     = cc$r,
                            cb    = cc$cb
                        )
                    )
                }
            }

            res$cccTable$setStatus('complete')

            # -----------------------------------------------------------------
            # If at least two Time 1 variables, compute BA, TOST, error metrics, plots (validation view)
            # -----------------------------------------------------------------
            if (nvars >= 2) {

                v1 <- vars[1]
                v2 <- vars[2]
                x  <- sub[[v1]]
                y  <- sub[[v2]]

                diff  <- x - y
                mean2 <- (x + y) / 2

                # Bland-Altman stats
                bias   <- mean(diff)
                sdDiff <- stats::sd(diff)
                loaLow <- bias - 1.96 * sdDiff
                loaUp  <- bias + 1.96 * sdDiff

                slope   <- NA_real_
                pSlope  <- NA_real_
                if (isTRUE(opts$baTrend) && length(diff) >= 2) {
                    fit <- try(stats::lm(diff ~ mean2), silent = TRUE)
                    if (!inherits(fit, "try-error")) {
                        coefs <- stats::coef(fit)
                        if (length(coefs) >= 2)
                            slope <- coefs[2]
                        smry <- summary(fit)$coefficients
                        if (!is.null(smry) &&
                            nrow(smry) >= 2 &&
                            "Pr(>|t|)" %in% colnames(smry)) {
                            pSlope <- smry[2, "Pr(>|t|)"]
                        }
                    }
                }

                res$baStats$setRow(
                    rowNo = 1,
                    values = list(
                        bias     = bias,
                        sdDiff   = sdDiff,
                        loaLower = loaLow,
                        loaUpper = loaUp,
                        slope    = slope,
                        pSlope   = pSlope
                    )
                )

                # Clinical limits comparison summary
                clinLow  <- opts$clinicalLower
                clinUp   <- opts$clinicalUpper
                includeClin <- isTRUE(opts$includeClinLimits)
                msg <- ""

                if (includeClin) {
                    if (!is.na(clinLow) && !is.na(clinUp) && clinLow < clinUp) {
                        if (!is.na(loaLow) && !is.na(loaUp)) {
                            withinClin <- (loaLow >= clinLow) && (loaUp <= clinUp)
                            if (withinClin) {
                                msg <- sprintf(
                                    "LoA [%.3f, %.3f] are within the clinical limits [%.3f, %.3f].",
                                    loaLow, loaUp, clinLow, clinUp
                                )
                            } else {
                                msg <- sprintf(
                                    "LoA [%.3f, %.3f] exceed the clinical limits [%.3f, %.3f].",
                                    loaLow, loaUp, clinLow, clinUp
                                )
                            }
                        } else {
                            msg <- "LoA could not be computed to compare with clinical limits."
                        }
                    } else {
                        msg <- "Clinical limits not set or invalid (lower must be < upper); no comparison with LoA made."
                    }
                } else {
                    msg <- "Clinical limits not included on the BA plot."
                }

                res$baSummary$setContent(msg)

                # ---------------- Common error metrics (for first pair) ------
                if (opts$showErrorMetrics) {
                    mae  <- mean(abs(diff))
                    rmse <- sqrt(mean(diff^2))

                    mn_y <- mean(y)
                    if (isTRUE(all.equal(mn_y, 0))) {
                        nrmse_mean <- NA_real_
                    } else {
                        nrmse_mean <- rmse / mn_y
                    }

                    # MAPE and sMAPE with protection against zero denominators
                    nonzero_y <- y != 0
                    if (any(nonzero_y)) {
                        ape  <- abs(diff[nonzero_y] / y[nonzero_y]) * 100
                        mape <- mean(ape)
                    } else {
                        mape <- NA_real_
                    }

                    denom <- (abs(x) + abs(y)) / 2
                    valid <- denom != 0
                    if (any(valid)) {
                        smape <- mean(abs(diff[valid]) / denom[valid] * 100)
                    } else {
                        smape <- NA_real_
                    }

                    TE   <- sdDiff / sqrt(2)
                    TEpc <- if (isTRUE(all.equal(mn_y, 0))) NA_real_ else TE / mn_y * 100

                    # Fill errorMetrics table
                    errTab <- res$errorMetrics
                    errTab$setRow(rowNo = 1, values = list(metric = "Bias (mean x - y)",        value = bias))
                    errTab$setRow(rowNo = 2, values = list(metric = "SD of differences",        value = sdDiff))
                    errTab$setRow(rowNo = 3, values = list(metric = "Mean absolute error (MAE)", value = mae))
                    errTab$setRow(rowNo = 4, values = list(metric = "Root mean square error (RMSE)", value = rmse))
                    errTab$setRow(rowNo = 5, values = list(metric = "Normalized RMSE (mean criterion)", value = nrmse_mean))
                    errTab$setRow(rowNo = 6, values = list(metric = "Mean absolute percent error (MAPE)", value = mape))
                    errTab$setRow(rowNo = 7, values = list(metric = "Symmetric MAPE (sMAPE)", value = smape))
                    errTab$setRow(rowNo = 8, values = list(metric = "Typical error (TE)", value = TE))
                    errTab$setRow(rowNo = 9, values = list(metric = "Typical error percent (TE%)", value = TEpc))
                }

                # ------------------------- TOST ------------------------------
                if (opts$doTOST) {
                    n  <- length(diff)
                    md <- bias
                    sd <- sdDiff
                    se <- sd / sqrt(n)
                    df_t <- n - 1

                    low <- opts$eqLower
                    up  <- opts$eqUpper

                    tLow  <- (md - low) / se
                    pLow  <- stats::pt(tLow, df = df_t, lower.tail = FALSE)
                    tUp   <- (md - up) / se
                    pUp   <- stats::pt(tUp, df = df_t, lower.tail = TRUE)
                    pTOST <- max(pLow, pUp)

                    res$tost$setRow(
                        rowNo = 1,
                        values = list(
                            meanDiff = md,
                            lower    = low,
                            upper    = up,
                            tLower   = tLow,
                            pLower   = pLow,
                            tUpper   = tUp,
                            pUpper   = pUp,
                            pTOST    = pTOST
                        )
                    )

                    # Equivalence region plot state
                    eqImage <- res$tostEquivPlot
                    eqImage$setState(list(
                        meanDiff = md,
                        se       = se,
                        df       = df_t,
                        lower    = low,
                        upper    = up
                    ))

                    # TOST decision curve state
                    if (opts$showTOSTDecisionCurve) {
                        decImage <- res$tostDecisionPlot
                        decImage$setState(list(
                            tLower = tLow,
                            tUpper = tUp,
                            df     = df_t
                        ))
                    }
                }

                # ------------------------- Plots -----------------------------
                if (opts$showBA) {
                    image <- res$baPlot
                    image$setState(list(
                        x        = x,
                        y        = y,
                        eqLower  = opts$eqLower,
                        eqUpper  = opts$eqUpper,
                        overlayEquivRegion = opts$overlayEquivRegion,
                        clinicalLower      = opts$clinicalLower,
                        clinicalUpper      = opts$clinicalUpper,
                        includeClinLimits  = opts$includeClinLimits
                    ))
                }

                if (opts$showMountain) {
                    image <- res$mountainPlot
                    image$setState(list(x = x, y = y))
                }
            }

            # -----------------------------------------------------------------
            # ICC table (all Time 1 vars as raters: agreement between measures at a time point)
            # -----------------------------------------------------------------
            if (opts$showICC && ncol(sub) >= 2) {

                mat <- as.matrix(sub)
                n   <- nrow(mat)
                k   <- ncol(mat)

                # One-way ANOVA (subjects only)
                dfLong <- data.frame(
                    subject = factor(rep(seq_len(n), times = k)),
                    rater   = factor(rep(seq_len(k), each = n)),
                    value   = as.vector(mat)
                )

                aov1 <- stats::aov(value ~ subject, data = dfLong)
                an1  <- summary(aov1)[[1]]
                MSB  <- an1["subject", "Mean Sq"]
                MSW  <- an1["Residuals", "Mean Sq"]

                ICC1  <- (MSB - MSW) / (MSB + (k - 1) * MSW)
                ICC1k <- (MSB - MSW) / MSB

                # Two-way ANOVA (subject + rater)
                aov2 <- stats::aov(value ~ subject + rater, data = dfLong)
                an2  <- summary(aov2)[[1]]
                MSR  <- an2["subject", "Mean Sq"]
                MSC  <- an2["rater", "Mean Sq"]
                MSE  <- an2["Residuals", "Mean Sq"]

                ICCC1 <- (MSR - MSE) / (MSR + (k - 1) * MSE)  # consistency, single
                ICCCk <- (MSR - MSE) / MSR                     # consistency, average

                ICCA1 <- (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n)  # agreement, single
                ICCAk <- (MSR - MSE) / (MSR + (MSC - MSE) / n)                      # agreement, average

                # Rows: Shrout & Fleiss style
                rows <- list(
                    list(model = "One-way random", type = "Agreement",   unit = "Single",  value = ICC1),
                    list(model = "One-way random", type = "Agreement",   unit = "Average", value = ICC1k),
                    list(model = "Two-way random", type = "Agreement",   unit = "Single",  value = ICCA1),
                    list(model = "Two-way random", type = "Agreement",   unit = "Average", value = ICCAk),
                    list(model = "Two-way mixed",  type = "Consistency", unit = "Single",  value = ICCC1),
                    list(model = "Two-way mixed",  type = "Consistency", unit = "Average", value = ICCCk)
                )

                for (i in seq_along(rows)) {
                    res$iccTable$setRow(rowNo = i, values = rows[[i]])
                }

            }

            # -----------------------------------------------------------------
            # Test-retest ICCs: consistency across time points
            # -----------------------------------------------------------------
            if (isTRUE(opts$hasRetest)) {

                trTestVars   <- opts$trTestVars
                trRetestVars <- opts$trRetestVars

                if (length(trTestVars) == 0 || length(trRetestVars) == 0) {
                    res$iccRetest$setStatus('error')
                    res$iccRetest$setError('Test-retest reliability is enabled, but no Test and/or Retest variables were supplied.')
                } else if (length(trTestVars) != length(trRetestVars)) {
                    res$iccRetest$setStatus('error')
                    res$iccRetest$setError('Number of test and retest variables must match (pair variables in order).')
                } else {

                    for (i in seq_along(trTestVars)) {

                        varTest   <- trTestVars[i]
                        varRetest <- trRetestVars[i]

                        xTest   <- jmvcore::toNumeric(data[[varTest]])
                        xRetest <- jmvcore::toNumeric(data[[varRetest]])

                        completePair <- stats::complete.cases(xTest, xRetest)
                        xTest   <- xTest[completePair]
                        xRetest <- xRetest[completePair]

                        n <- length(xTest)

                        if (n < 2) {
                            res$iccRetest$addRow(
                                rowKey = paste0("retest_", varTest),
                                values = list(
                                    measure = varTest,
                                    icc     = NA_real_,
                                    sem     = NA_real_,
                                    mdc95   = NA_real_,
                                    n       = n
                                )
                            )
                            next
                        }

                        dfLong_m <- data.frame(
                            subject = factor(rep(seq_len(n), times = 2)),
                            time    = factor(rep(c("T1", "T2"), each = n)),
                            value   = c(xTest, xRetest)
                        )

                        aov_m <- stats::aov(value ~ subject + time, data = dfLong_m)
                        anm   <- summary(aov_m)[[1]]
                        MSR_m <- anm["subject", "Mean Sq"]
                        MSE_m <- anm["Residuals", "Mean Sq"]

                        # Consistency ICC across time points (two-way mixed, single measure)
                        ICC_retest <- (MSR_m - MSE_m) / (MSR_m + (2 - 1) * MSE_m)

                        sd_pool <- stats::sd(c(xTest, xRetest), na.rm = TRUE)

                        if (!is.na(sd_pool) && !is.na(ICC_retest)) {
                            SEM   <- sd_pool * sqrt(1 - ICC_retest)
                            MDC95 <- SEM * 1.96 * sqrt(2)
                        } else {
                            SEM   <- NA_real_
                            MDC95 <- NA_real_
                        }

                        rowKey <- paste0("retest_", varTest)

                        res$iccRetest$addRow(
                            rowKey = rowKey,
                            values = list(
                                measure = varTest,
                                icc     = ICC_retest,
                                sem     = SEM,
                                mdc95   = MDC95,
                                n       = n
                            )
                        )
                    }
                }
            }

            # -----------------------------------------------------------------
            # Within-session repeatability (up to 5 measures)
            # -----------------------------------------------------------------
            if (isTRUE(opts$doWithinSession)) {

                wsTab <- res$withinSession
                rowKey <- 1L

                wsLists <- list(
                    opts$wsMeasure1,
                    opts$wsMeasure2,
                    opts$wsMeasure3,
                    opts$wsMeasure4,
                    opts$wsMeasure5
                )

                for (m in seq_along(wsLists)) {
                    varsWS <- wsLists[[m]]
                    if (length(varsWS) < 2)
                        next

                    subWS <- data[varsWS]
                    completeWS <- stats::complete.cases(subWS)
                    subWS <- subWS[completeWS, , drop = FALSE]

                    if (nrow(subWS) < 2)
                        next

                    matWS <- as.matrix(subWS)
                    n  <- nrow(matWS)
                    k  <- ncol(matWS)

                    dfLongWS <- data.frame(
                        subject = factor(rep(seq_len(n), times = k)),
                        trial   = factor(rep(seq_len(k), each = n)),
                        value   = as.vector(matWS)
                    )

                    aovWS <- stats::aov(value ~ subject, data = dfLongWS)
                    anWS  <- summary(aovWS)[[1]]
                    MSBws <- anWS["subject", "Mean Sq"]
                    MSWws <- anWS["Residuals", "Mean Sq"]

                    ICC_single  <- (MSBws - MSWws) / (MSBws + (k - 1) * MSWws)
                    ICC_average <- (MSBws - MSWws) / MSBws

                    TEws    <- sqrt(MSWws)
                    meanAll <- mean(dfLongWS$value, na.rm = TRUE)
                    if (isTRUE(all.equal(meanAll, 0))) {
                        CVws <- NA_real_
                    } else {
                        CVws <- TEws / meanAll * 100
                    }

                    SEMws   <- TEws
                    MDC95ws <- SEMws * 1.96 * sqrt(2)

                    measLabel <- if (length(varsWS) > 0) varsWS[0 + 1] else sprintf("Measure %i", m)

                    wsTab$addRow(
                        rowKey = as.character(rowKey),
                        values = list(
                            measure   = measLabel,
                            n         = n,
                            nTrials   = k,
                            iccSingle = ICC_single,
                            iccAverage= ICC_average,
                            sem       = SEMws,
                            cv        = CVws,
                            mdc95     = MDC95ws
                        )
                    )
                    rowKey <- rowKey + 1L
                }
            }
        },

        .baPlot = function(image, ...) {
            state <- image$state
            if (is.null(state))
                return()

            x <- state$x
            y <- state$y

            diff  <- x - y
            mean2 <- (x + y) / 2

            bias   <- mean(diff)
            sdDiff <- stats::sd(diff)
            loaLow <- bias - 1.96 * sdDiff
            loaUp  <- bias + 1.96 * sdDiff

            eqLower <- state$eqLower
            eqUpper <- state$eqUpper
            overlay <- isTRUE(state$overlayEquivRegion)

            clinLow <- state$clinicalLower
            clinUp  <- state$clinicalUpper
            includeClin <- isTRUE(state$includeClinLimits)

            df <- data.frame(
                mean = mean2,
                diff = diff
            )

            p <- ggplot2::ggplot(df, ggplot2::aes(x = mean, y = diff)) +
                 ggplot2::theme_classic(base_size = 16) +
                 ggplot2::theme(axis.line = element_line(linewidth = 1))

            # Equivalence region shading
            if (overlay && !is.null(eqLower) && !is.null(eqUpper) &&
                !is.na(eqLower) && !is.na(eqUpper) && eqLower < eqUpper) {
                p <- p +
                    ggplot2::annotate(
                        "rect",
                        xmin = -Inf, xmax = Inf,
                        ymin = eqLower, ymax = eqUpper,
                        fill = "blue",
                        alpha = 0.075
                    )
            }

            p <- p +
                ggplot2::geom_point(alpha = 0.6) +
                ggplot2::geom_hline(yintercept = bias, linetype = "solid", linewidth = 1.2, color = "blue") +
                ggplot2::geom_hline(yintercept = loaLow, linetype = "dashed", linewidth = 1, color = "red") +
                ggplot2::geom_hline(yintercept = loaUp,  linetype = "dashed", linewidth = 1, color = "red")

            if (isTRUE(self$options$baTrend)) {
                if (identical(self$options$baTrendMethod, "lm")) {
                    p <- p + ggplot2::geom_smooth(
                        method   = "lm",
                        se       = TRUE,
                        linetype = "solid",
                        #color    = "#24c724",
                        #fill      = "green",
                        color    = "#969696",
                        fill      = "gray",
                        alpha     = 0.2,
                        linewidth = 1
                    )
                } else if (identical(self$options$baTrendMethod, "loess")) {
                    p <- p + ggplot2::geom_smooth(
                        se       = FALSE,
                        linetype = "solid",
                        color    = "#969696",
                        linewidth = 1.2
                    )
                }
            }

            # Clinical limit lines
            if (includeClin && !is.null(clinLow) && !is.null(clinUp) &&
                !is.na(clinLow) && !is.na(clinUp) && clinLow < clinUp) {
                p <- p +
                    ggplot2::geom_hline(yintercept = clinLow, linetype = "dotdash", linewidth = 1, color = "#8947af") +
                    ggplot2::geom_hline(yintercept = clinUp,  linetype = "dotdash", linewidth = 1, color = "#8947af")
            }

            p <- p +
                ggplot2::labs(
                    x = "Mean of measures",
                    y = "Difference (x - y)"
                )

            print(p)
        },

        .mountainPlot = function(image, ...) {
            state <- image$state
            if (is.null(state))
                return()

            x <- state$x
            y <- state$y

            diff <- sort(x - y)
            n    <- length(diff)
            if (n < 2)
                return()

            pvals <- (seq_len(n) - 0.5) / n
            folded <- abs(pvals - 0.5)

            df <- data.frame(
                diff   = diff,
                folded = folded,
                mountain = 0.5 - folded   # peak = 0.5, base = 0
            )

            x.max <- max(abs(df$diff), na.rm = TRUE)

            g <- ggplot2::ggplot(df, ggplot2::aes(x = diff, y = mountain)) +
                ggplot2::theme_classic(base_size = 16) +
                ggplot2::theme(axis.line = element_line(linewidth = 1)) +
                ggplot2::geom_line(linewidth = 1.2, color = "blue") +
                ggplot2::labs(
                    x = "Difference (x - y)",
                    y = "Folded cumulative proportion"
                )

            g <- ggplot2::ggplot(df, ggplot2::aes(x = diff, y = mountain)) +
                ggplot2::theme_classic(base_size = 16) +
                ggplot2::theme(
                    axis.line  = ggplot2::element_line(linewidth = 1),
                    axis.title = ggplot2::element_text(margin = margin(t = 6, r = 6, b = 6, l = 6))
                ) +
                ggplot2::geom_vline(
                    xintercept = 0,
                    linetype   = "dashed",
                    linewidth  = 1,
                    color      = "red"
                ) +
                ggplot2::geom_area(
                    fill  = "gray",
                    alpha = 0.2
                ) +
                ggplot2::geom_line(
                    linewidth = 1.4,
                    color     = "blue"
                ) +
                ggplot2::scale_x_continuous(
                    limits = c(-x.max, x.max),         # symmetric around 0
                    expand = ggplot2::expansion(mult = 0.02)
                ) +
                ggplot2::scale_y_continuous(
                    limits = c(0, 0.5),
                    expand = ggplot2::expansion(mult = c(0.02, 0.05))
                ) +
                ggplot2::labs(
                    x = "Difference (x - y)",
                    y = "Folded cumulative proportion"
                )

            print(g)
        },

        .tostEquivPlot = function(image, ...) {
            state <- image$state
            if (is.null(state))
                return()

            md    <- state$meanDiff
            se    <- state$se
            df_t  <- state$df
            low   <- state$lower
            up    <- state$upper

            if (is.null(md) || is.null(se) || is.null(df_t) ||
                is.null(low) || is.null(up) ||
                any(is.na(c(md, se, df_t, low, up)))) {
                return()
            }

            # 90% CI for mean difference (TOST uses 1 - 2*alpha)
            tcrit <- stats::qt(0.95, df = df_t)
            ciLow <- md - tcrit * se
            ciUp  <- md + tcrit * se

            df_ci <- data.frame(
                y      = 1,
                est    = md,
                ciLow  = ciLow,
                ciUp   = ciUp
            )

            p <- ggplot2::ggplot(df_ci, ggplot2::aes(y = y)) +
                ggplot2::geom_errorbarh(
                    ggplot2::aes(xmin = ciLow, xmax = ciUp),
                    height = 0.1,
                    size = 1.3,
                    color = "blue"
                ) +
                ggplot2::geom_point(
                    ggplot2::aes(x = est),
                    shape = 24,
                    size = 3.5,
                    color = "blue",
                    fill = "white",
                    stroke = 2.5
                ) +
                ggplot2::geom_vline(xintercept = low, linetype = "dashed", linewidth = 1, color = "red") +
                ggplot2::geom_vline(xintercept = up,  linetype = "dashed", linewidth = 1, color = "red") +
                ggplot2::annotate(
                    "rect",
                    xmin = low, xmax = up,
                    ymin = -Inf, ymax = Inf,      # full height of the panel
                    #fill = "#579ceb",
                    alpha = 0.075
                ) +
                ggplot2::labs(
                    x = "Mean difference with 90% CI",
                    y = ""
                    #title = "TOST equivalence region"
                ) +
                ggplot2::theme_classic(base_size = 16) +
                ggplot2::theme(
                    axis.line = element_line(linewidth = 1),
                    axis.line.y = ggplot2::element_blank(),
                    axis.text.y = ggplot2::element_blank(),
                    axis.ticks.y = ggplot2::element_blank(),
                    aspect.ratio = 0.25
                )

                p <- p +
                    #ggplot2::coord_cartesian(ylim = c(-0.5, 1.5))
                    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.5))

            print(p)
        },

        .tostDecisionPlot = function(image, ...) {
            state <- image$state
            if (is.null(state))
                return()

            tLow <- state$tLower
            tUp  <- state$tUpper
            df_t <- state$df

            if (any(is.null(c(tLow, tUp, df_t))) ||
                any(is.na(c(tLow, tUp, df_t)))) {
                return()
            }

            tCrit <- stats::qt(0.95, df = df_t)
            df_tost <- data.frame(
                test   = c("Lower TOST", "Upper TOST"),
                tvalue = c(abs(tLow), abs(tUp))
            )

            p <- ggplot2::ggplot(df_tost, ggplot2::aes(x = test, y = tvalue)) +
                ggplot2::geom_col(alpha = 0.6) +
                ggplot2::geom_hline(yintercept = tCrit, linetype = "dashed") +
                ggplot2::labs(
                    x = "",
                    y = "Absolute t-value",
                    title = "TOST decision curve"
                ) +
                ggplot2::theme_minimal()

            print(p)
        }
    )
)

# Convenience wrapper when using as an R package
agreement <- function(data,
                      vars,
                      hasRetest = FALSE,
                      trTestVars = NULL,
                      trRetestVars = NULL,
                      eqLower = -0.1,
                      eqUpper = 0.1,
                      confLevel = 0.95,
                      nBoot = 1000,
                      clinicalLower = 0,
                      clinicalUpper = 0,
                      includeClinLimits = FALSE,
                      showBA = TRUE,
                      baTrend = FALSE,
                      baTrendMethod = "lm",
                      overlayEquivRegion = FALSE,
                      showMountain = TRUE,
                      doTOST = TRUE,
                      showTOSTDecisionCurve = FALSE,
                      showICC = TRUE,
                      showErrorMetrics = TRUE,
                      doWithinSession = FALSE,
                      wsMeasure1 = NULL,
                      wsMeasure2 = NULL,
                      wsMeasure3 = NULL,
                      wsMeasure4 = NULL,
                      wsMeasure5 = NULL) {

    if (!requireNamespace("jmvcore", quietly = TRUE))
        stop("The 'jmvcore' package is required for this function.")

    options <- ralidOptions$new(
        vars               = vars,
        hasRetest          = hasRetest,
        trTestVars         = trTestVars,
        trRetestVars       = trRetestVars,
        eqLower            = eqLower,
        eqUpper            = eqUpper,
        confLevel          = confLevel,
        nBoot              = nBoot,
        clinicalLower      = clinicalLower,
        clinicalUpper      = clinicalUpper,
        includeClinLimits  = includeClinLimits,
        showBA             = showBA,
        baTrend            = baTrend,
        baTrendMethod      = baTrendMethod,
        overlayEquivRegion = overlayEquivRegion,
        showMountain       = showMountain,
        doTOST             = doTOST,
        showTOSTDecisionCurve = showTOSTDecisionCurve,
        showICC            = showICC,
        showErrorMetrics   = showErrorMetrics,
        doWithinSession    = doWithinSession,
        wsMeasure1         = wsMeasure1,
        wsMeasure2         = wsMeasure2,
        wsMeasure3         = wsMeasure3,
        wsMeasure4         = wsMeasure4,
        wsMeasure5         = wsMeasure5
    )

    analysis <- ralidClass$new(
        options = options,
        data    = data
    )

    analysis$run()

    analysis$results
}
