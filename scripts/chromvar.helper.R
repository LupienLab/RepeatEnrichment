differentialDeviations2 <- function(object,
                                    groups,
                                    alternative = c("two.sided", "less",
                                                    "greater"),
                                    parametric = TRUE) {
  stopifnot(is(object,"chromVARDeviations"))
  if (length(groups) == 1 && groups %in% colnames(colData(object))) {
    groups <- colData(object)[[groups]]
  } else if (length(groups) != ncol(object)) {
    stop("invalid groups input, must be vector of lench ncol(object) or column",
         " name from colData(object)")
  }

  groups <- as.factor(groups)

  alternative <- match.arg(alternative)
  inputs <- deviations(object)
  inputs <- na.omit(inputs)
  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      p_val <- apply(inputs, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(inputs, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(inputs, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(inputs, 1, kw_helper, groups)
    }
  }

  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}


t_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(t.test(splitx[[1]],splitx[[2]],
                alternative = alternative,
                paired = FALSE,
                var.equal = FALSE)$p.value)
}

anova_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- oneway.test(devs ~ groups, tmpdf, var.equal = FALSE)
  return(res$p.value)
}

kw_helper <- function(x, groups) {
  tmpdf <- data.frame(groups = groups, devs = x)
  res <- kruskal.test(devs ~ groups, tmpdf)
  return(res$p.value)
}

wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
                     alternative = alternative,
                     paired = FALSE)$p.value)
}

differentialVariability2 <- function (object, groups, parametric = TRUE)
{
    stopifnot(is(object, "chromVARDeviations"))
    if (length(groups) == 1 && groups %in% colnames(colData(object))) {
        groups <- colData(object)[[groups]]
    }
    else if (length(groups) != ncol(object)) {
        stop("invalid groups input, must be vector of lench ncol(object) or column",
            " name from colData(object)")
    }
    groups <- as.factor(groups)
    inputs <- deviationScores(object)
    inputs <- na.omit(inputs)

    if (parametric) {
        p_val <- apply(inputs, 1, bf_var_test, groups)
    }
    else {
        p_val <- apply(inputs, 1, bf_kw_var_test, groups)
    }
    p_adj <- p.adjust(p_val, method = "BH")
    return(data.frame(p_value = p_val, p_value_adjusted = p_adj))
}



bf_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(anova(lm(median_diff ~ groups))[1, 5])
}

bf_kw_var_test <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  median_diff <- abs(x - unsplit(medians, groups))
  return(kruskal.test(median_diff ~ groups)$p.value)
}

## Run differential testing on means/medians of the zscore and not check variance

mediancalc <- function(x, groups) {
  medians <- aggregate(x, list(groups), median)$x
  return(medians)
}


differentialDeviations3 <- function(inputs,
                                    groups,
                                    alternative = c("two.sided", "less","greater"),  parametric = TRUE) {

  groups <- as.factor(groups)
  alternative <- match.arg(alternative)
  #inputs <- deviationScores(object)
  inputs <- na.omit(inputs)

  if (parametric) {
    if (nlevels(groups) == 2) {
      # t-test
      p_val <- apply(inputs, 1, t_helper, groups, alternative)
    } else {
      # anova
      p_val <- apply(inputs, 1, anova_helper, groups)
    }
  } else {
    if (nlevels(groups) == 2) {
      # wilcoxon
      p_val <- apply(inputs, 1, wilcoxon_helper, groups, alternative)
    } else {
      # kruskal-wallis
      p_val <- apply(inputs, 1, kw_helper, groups)
    }
  }
  medians <- apply(inputs, 1, function(e) aggregate(e, list(groups), median)$x)
  p_adj <- p.adjust(p_val, method = "BH")
  return(data.frame(p_value = p_val, p_value_adjusted = p_adj, median_grp1=medians[1,], median_grp2=medians[2,] ))
}
