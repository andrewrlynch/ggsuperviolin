StatSuperviolin <- ggproto(
  "StatSuperviolin", Stat,
  required_aes = c("x", "y", "rep"),
  default_aes  = aes(weight = 1),

  compute_panel = function(data, scales,
                           bw       = NULL,
                           n        = 128,
                           trim     = TRUE,
                           width    = 0.9,
                           preserve = c("pooled", "total", "rep"),
                           na.rm    = FALSE) {

    preserve <- match.arg(preserve)

    data <- data[stats::complete.cases(data[c("x", "y", "rep")]), , drop = FALSE]
    if (nrow(data) == 0L) return(data.frame())

    # normalise bw
    if (is.null(bw) || length(bw) == 0L) {
      use_bw <- NULL
    } else {
      if (length(bw) != 1L || !is.finite(bw))
        stop("`bw` must be a single finite numeric value or NULL.", call. = FALSE)
      use_bw <- bw
    }

    x_res   <- ggplot2::resolution(data$x, zero = FALSE)
    base_w  <- width * x_res

    x_vals  <- sort(unique(data$x))
    out     <- vector("list", length(x_vals))
    idx_out <- 1L
    group_id <- 1L

    for (xi in x_vals) {
      df_x <- data[data$x == xi, , drop = FALSE]
      reps <- unique(df_x$rep)
      n_rep <- length(reps)
      if (n_rep == 0L) next

      y_rng <- range(df_x$y, na.rm = TRUE)
      y_min <- y_rng[1]
      y_max <- y_rng[2]
      if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) next

      y_grid <- seq(y_min, y_max, length.out = n)

      dens_mat <- matrix(0, nrow = n, ncol = n_rep)
      colnames(dens_mat) <- as.character(reps)

      for (k in seq_along(reps)) {
        yk <- df_x$y[df_x$rep == reps[k]]
        if (length(unique(yk[is.finite(yk)])) < 2L) next

        if (is.null(use_bw)) {
          d <- stats::density(yk,
                              from = y_min,
                              to   = y_max,
                              n    = n,
                              na.rm = TRUE)
        } else {
          d <- stats::density(yk,
                              from = y_min,
                              to   = y_max,
                              n    = n,
                              bw   = use_bw,
                              na.rm = TRUE)
        }
        dens_mat[, k] <- d$y
      }

      # if everything is zero, skip
      if (all(rowSums(dens_mat) <= 0)) next

      polys <- vector("list", n_rep)

      if (preserve %in% c("pooled", "total")) {
        # stripe proportions come from replicate densities
        total_rep <- rowSums(dens_mat)

        # pooled KDE for the envelope when preserve == "pooled"
        pool <- rep(0, n)
        yk_all <- df_x$y
        if (length(unique(yk_all[is.finite(yk_all)])) >= 2L) {
          if (is.null(use_bw)) {
            d_pool <- stats::density(yk_all, from = y_min, to = y_max, n = n, na.rm = TRUE)
          } else {
            d_pool <- stats::density(yk_all, from = y_min, to = y_max, n = n, bw = use_bw, na.rm = TRUE)
          }
          pool <- d_pool$y
        }

        envelope <- if (preserve == "pooled") pool else total_rep
        max_env  <- max(envelope)
        total_scaled <- if (max_env > 0) envelope / max_env else envelope

        p_mat    <- matrix(0, nrow = n, ncol = n_rep)
        positive <- total_rep > 0
        p_mat[positive, ] <- dens_mat[positive, , drop = FALSE] / total_rep[positive]

        if (n_rep == 1L) {
          cum_mat   <- p_mat
          start_mat <- matrix(0, nrow = n, ncol = 1L)
        } else {
          cum_mat   <- t(apply(p_mat, 1L, cumsum))
          start_mat <- cbind(0, cum_mat[, -ncol(cum_mat), drop = FALSE])
        }

        for (k in seq_along(reps)) {
          rep_id  <- reps[k]
          start_f <- start_mat[, k]
          end_f   <- cum_mat[, k]
          if (all(end_f <= 0)) next

          total_w <- base_w * total_scaled

          left_x  <- xi - total_w / 2 + total_w * start_f
          right_x <- xi - total_w / 2 + total_w * end_f

          poly_x <- c(left_x, rev(right_x))
          poly_y <- c(y_grid, rev(y_grid))

          base_row <- df_x[df_x$rep == rep_id, , drop = FALSE][1, ]

          poly_df <- base_row[rep(1, length(poly_x)), , drop = FALSE]
          poly_df$x     <- poly_x
          poly_df$y     <- poly_y
          poly_df$rep   <- rep_id
          poly_df$group <- group_id

          polys[[k]] <- poly_df
          group_id <- group_id + 1L
        }

      } else {
        # ---- preserve == "rep": keep each replicate's own KDE shape ----
        slot_w <- base_w / n_rep

        for (k in seq_along(reps)) {
          rep_id <- reps[k]
          d_k <- dens_mat[, k]
          if (all(d_k <= 0)) next

          d_scaled <- d_k / max(d_k)  # 0–1 within this replicate

          # centre of this replicate's slot
          center_frac <- (k - 0.5) / n_rep - 0.5  # -0.5 .. 0.5
          x_center <- xi + base_w * center_frac

          half_w <- (slot_w * d_scaled) / 2

          left_x  <- x_center - half_w
          right_x <- x_center + half_w

          poly_x <- c(left_x, rev(right_x))
          poly_y <- c(y_grid, rev(y_grid))

          base_row <- df_x[df_x$rep == rep_id, , drop = FALSE][1, ]

          poly_df <- base_row[rep(1, length(poly_x)), , drop = FALSE]
          poly_df$x     <- poly_x
          poly_df$y     <- poly_y
          poly_df$rep   <- rep_id
          poly_df$group <- group_id

          polys[[k]] <- poly_df
          group_id <- group_id + 1L
        }
      }

      out[[idx_out]] <- do.call(rbind, polys[!vapply(polys, is.null, logical(1))])
      idx_out <- idx_out + 1L
    }

    out_df <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
    rownames(out_df) <- NULL
    out_df
  }
)

stat_superviolin <- function(mapping = NULL, data = NULL, geom = "superviolin",
                             position = "identity",
                             ...,
                             bw       = NULL,
                             n        = 128,
                             trim     = TRUE,
                             width    = 0.9,
                             preserve = c("pooled", "total", "rep"),
                             na.rm    = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {

  preserve <- match.arg(preserve)

  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatSuperviolin,
    geom        = geom,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      bw       = bw,
      n        = n,
      trim     = trim,
      width    = width,
      preserve = preserve,
      na.rm    = na.rm,
      ...
    )
  )
}

GeomSuperviolin <- ggproto(
  "GeomSuperviolin", GeomPolygon,
  default_aes = aes(
    colour   = NA,
    fill     = "grey70",
    linewidth = 0.2,
    linetype  = 1,
    alpha     = 1
  ),
  required_aes = c("x", "y", "group"),
  draw_key = draw_key_polygon
)

geom_superviolin <- function(mapping = NULL, data = NULL,
                             stat = "superviolin",
                             position = "identity",
                             ...,
                             bw       = NULL,
                             n        = 128,
                             trim     = TRUE,
                             width    = 0.9,
                             preserve = c("pooled", "total", "rep"),
                             na.rm    = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {

  preserve <- match.arg(preserve)

  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatSuperviolin,
    geom        = GeomSuperviolin,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      bw       = bw,
      n        = n,
      trim     = trim,
      width    = width,
      preserve = preserve,
      na.rm    = na.rm,
      ...
    )
  )
}

StatSupercentral <- ggproto(
  "StatSupercentral", Stat,
  required_aes = c("x", "y", "rep"),
  default_aes  = aes(),

  compute_panel = function(data, scales,
                           bw       = NULL,
                           n        = 128,
                           trim     = TRUE,
                           width    = 0.9,
                           preserve = c("pooled", "total", "rep"),
                           fun      = "median",
                           na.rm    = FALSE) {

    preserve <- match.arg(preserve)
    fun_y    <- match.fun(fun)

    data <- data[stats::complete.cases(data[c("x", "y", "rep")]), , drop = FALSE]
    if (nrow(data) == 0L) return(data.frame())

    # normalise bw
    if (is.null(bw) || length(bw) == 0L) {
      use_bw <- NULL
    } else {
      if (length(bw) != 1L || !is.finite(bw))
        stop("`bw` must be a single finite numeric value or NULL.", call. = FALSE)
      use_bw <- bw
    }

    x_res  <- ggplot2::resolution(data$x, zero = FALSE)
    base_w <- width * x_res

    x_vals  <- sort(unique(data$x))
    out     <- vector("list", length(x_vals))
    idx_out <- 1L

    for (xi in x_vals) {
      df_x <- data[data$x == xi, , drop = FALSE]
      reps <- unique(df_x$rep)
      n_rep <- length(reps)
      if (n_rep == 0L) next

      y_rng <- range(df_x$y, na.rm = TRUE)
      y_min <- y_rng[1]
      y_max <- y_rng[2]
      if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) next

      y_grid <- seq(y_min, y_max, length.out = n)

      dens_mat <- matrix(0, nrow = n, ncol = n_rep)
      colnames(dens_mat) <- as.character(reps)

      for (k in seq_along(reps)) {
        yk <- df_x$y[df_x$rep == reps[k]]
        if (length(unique(yk[is.finite(yk)])) < 2L) next

        if (is.null(use_bw)) {
          d <- stats::density(yk,
                              from = y_min,
                              to   = y_max,
                              n    = n,
                              na.rm = TRUE)
        } else {
          d <- stats::density(yk,
                              from = y_min,
                              to   = y_max,
                              n    = n,
                              bw   = use_bw,
                              na.rm = TRUE)
        }
        dens_mat[, k] <- d$y
      }

      if (all(rowSums(dens_mat) <= 0)) next

      pts <- vector("list", n_rep)

      if (preserve %in% c("pooled", "total")) {
        # stripe proportions come from replicate densities
        total_rep <- rowSums(dens_mat)

        # pooled KDE for the envelope when preserve == "pooled"
        pool <- rep(0, n)
        yk_all <- df_x$y
        if (length(unique(yk_all[is.finite(yk_all)])) >= 2L) {
          if (is.null(use_bw)) {
            d_pool <- stats::density(yk_all, from = y_min, to = y_max, n = n, na.rm = TRUE)
          } else {
            d_pool <- stats::density(yk_all, from = y_min, to = y_max, n = n, bw = use_bw, na.rm = TRUE)
          }
          pool <- d_pool$y
        }

        envelope <- if (preserve == "pooled") pool else total_rep
        max_env  <- max(envelope)
        total_scaled <- if (max_env > 0) envelope / max_env else envelope

        p_mat    <- matrix(0, nrow = n, ncol = n_rep)
        positive <- total_rep > 0
        p_mat[positive, ] <- dens_mat[positive, , drop = FALSE] / total_rep[positive]

        if (n_rep == 1L) {
          cum_mat   <- p_mat
          start_mat <- matrix(0, nrow = n, ncol = 1L)
        } else {
          cum_mat   <- t(apply(p_mat, 1L, cumsum))
          start_mat <- cbind(0, cum_mat[, -ncol(cum_mat), drop = FALSE])
        }

        total_w <- base_w * total_scaled

        for (k in seq_along(reps)) {
          rep_id <- reps[k]
          yk     <- df_x$y[df_x$rep == rep_id]
          if (!any(dens_mat[, k] > 0) || length(yk) == 0L) next

          y_t <- fun_y(yk)
          if (!is.finite(y_t)) next

          # clamp into grid range
          y_t <- min(max(y_t, y_min), y_max)
          j   <- which.min(abs(y_grid - y_t))

          start_f <- start_mat[j, k]
          end_f   <- cum_mat[j, k]

          # if stripe is effectively zero here, skip
          if (end_f <= start_f) next

          left_x  <- xi - total_w[j] / 2 + total_w[j] * start_f
          right_x <- xi - total_w[j] / 2 + total_w[j] * end_f
          x_mid   <- (left_x + right_x) / 2

          base_row <- df_x[df_x$rep == rep_id, , drop = FALSE][1, , drop = FALSE]
          pt <- base_row
          pt$x <- x_mid
          pt$y <- y_t

          pts[[k]] <- pt
        }

      } else {
        # per-replicate layout (preserve = "rep")
        slot_w <- base_w / n_rep

        for (k in seq_along(reps)) {
          rep_id <- reps[k]
          yk     <- df_x$y[df_x$rep == rep_id]
          d_k    <- dens_mat[, k]
          if (!any(d_k > 0) || length(yk) == 0L) next

          y_t <- fun_y(yk)
          if (!is.finite(y_t)) next

          # centre of this replicate's slot
          center_frac <- (k - 0.5) / n_rep - 0.5
          x_center <- xi + base_w * center_frac

          base_row <- df_x[df_x$rep == rep_id, , drop = FALSE][1, , drop = FALSE]
          pt <- base_row
          pt$x <- x_center
          pt$y <- y_t

          pts[[k]] <- pt
        }
      }

      out[[idx_out]] <- do.call(rbind, pts[!vapply(pts, is.null, logical(1))])
      idx_out <- idx_out + 1L
    }

    out_df <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
    rownames(out_df) <- NULL
    out_df
  }
)

stat_supercentral <- function(mapping = NULL, data = NULL,
                              geom = "point", position = "identity",
                              ...,
                              bw       = NULL,
                              n        = 128,
                              trim     = TRUE,
                              width    = 0.9,
                              preserve = c("pooled", "total", "rep"),
                              fun      = "median",
                              na.rm    = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {

  preserve <- match.arg(preserve)

  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatSupercentral,
    geom        = geom,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      bw       = bw,
      n        = n,
      trim     = trim,
      width    = width,
      preserve = preserve,
      fun      = fun,
      na.rm    = na.rm,
      ...
    )
  )
}

StatSuperpoint <- ggproto(
  "StatSuperpoint", Stat,
  required_aes = c("x", "y", "rep"),
  default_aes  = aes(),

  compute_panel = function(data, scales,
                           bw          = NULL,
                           n           = 128,
                           trim        = TRUE,
                           width       = 0.9,
                           preserve    = c("pooled", "total", "rep"),
                           point_inset = 0.9,
                           na.rm       = FALSE) {

    preserve <- match.arg(preserve)

    data <- data[stats::complete.cases(data[c("x", "y", "rep")]), , drop = FALSE]
    if (nrow(data) == 0L) return(data.frame())

    # normalise bw
    if (is.null(bw) || length(bw) == 0L) {
      use_bw <- NULL
    } else {
      if (length(bw) != 1L || !is.finite(bw))
        stop("`bw` must be a single finite numeric value or NULL.", call. = FALSE)
      use_bw <- bw
    }

    if (!is.numeric(point_inset) || length(point_inset) != 1L ||
        point_inset <= 0 || point_inset > 1)
      stop("`point_inset` must be in (0, 1].", call. = FALSE)

    x_res  <- ggplot2::resolution(data$x, zero = FALSE)
    base_w <- width * x_res

    x_vals  <- sort(unique(data$x))
    out     <- vector("list", length(x_vals))
    idx_out <- 1L

    for (xi in x_vals) {

      df_x <- data[data$x == xi, , drop = FALSE]
      reps <- unique(df_x$rep)
      n_rep <- length(reps)
      if (n_rep == 0L) next

      y_rng <- range(df_x$y, na.rm = TRUE)
      y_min <- y_rng[1]
      y_max <- y_rng[2]
      if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) next

      y_grid <- seq(y_min, y_max, length.out = n)

      dens_mat <- matrix(0, nrow = n, ncol = n_rep)
      colnames(dens_mat) <- as.character(reps)

      # KDE per replicate on common grid
      for (k in seq_along(reps)) {
        yk <- df_x$y[df_x$rep == reps[k]]
        if (length(unique(yk[is.finite(yk)])) < 2L) next

        if (is.null(use_bw)) {
          d <- stats::density(yk,
                              from = y_min,
                              to   = y_max,
                              n    = n,
                              na.rm = TRUE)
        } else {
          d <- stats::density(yk,
                              from = y_min,
                              to   = y_max,
                              n    = n,
                              bw   = use_bw,
                              na.rm = TRUE)
        }
        dens_mat[, k] <- d$y
      }

      if (all(rowSums(dens_mat) <= 0)) next

      pts <- vector("list", n_rep)

      if (preserve %in% c("pooled", "total")) {
        # stripe proportions come from replicate densities
        total_rep <- rowSums(dens_mat)

        # pooled KDE for the envelope when preserve == "pooled"
        pool <- rep(0, n)
        yk_all <- df_x$y
        if (length(unique(yk_all[is.finite(yk_all)])) >= 2L) {
          if (is.null(use_bw)) {
            d_pool <- stats::density(yk_all, from = y_min, to = y_max, n = n, na.rm = TRUE)
          } else {
            d_pool <- stats::density(yk_all, from = y_min, to = y_max, n = n, bw = use_bw, na.rm = TRUE)
          }
          pool <- d_pool$y
        }

        envelope <- if (preserve == "pooled") pool else total_rep
        max_env  <- max(envelope)
        total_scaled <- if (max_env > 0) envelope / max_env else envelope

        p_mat    <- matrix(0, nrow = n, ncol = n_rep)
        positive <- total_rep > 0
        p_mat[positive, ] <- dens_mat[positive, , drop = FALSE] / total_rep[positive]

        if (n_rep == 1L) {
          cum_mat   <- p_mat
          start_mat <- matrix(0, nrow = n, ncol = 1L)
        } else {
          cum_mat   <- t(apply(p_mat, 1L, cumsum))
          start_mat <- cbind(0, cum_mat[, -ncol(cum_mat), drop = FALSE])
        }

        total_w <- base_w * total_scaled

        for (k in seq_along(reps)) {
          rep_id <- reps[k]
          d_k <- dens_mat[, k]
          df_rep <- df_x[df_x$rep == rep_id, , drop = FALSE]
          if (!any(d_k > 0) || nrow(df_rep) == 0L) next

          x_vals_pt <- numeric(nrow(df_rep))

          for (i in seq_len(nrow(df_rep))) {
            y_i <- df_rep$y[i]
            if (!is.finite(y_i)) next

            # clamp into grid, find nearest row
            y_i <- min(max(y_i, y_min), y_max)
            j   <- which.min(abs(y_grid - y_i))

            start_f <- start_mat[j, k]
            end_f   <- cum_mat[j, k]

            if (end_f <= start_f || total_w[j] <= 0) {
              # fallback: tiny jitter around centre
              x_vals_pt[i] <- stats::runif(1, xi - base_w * 0.02, xi + base_w * 0.02)
            } else {
              left  <- xi - total_w[j] / 2 + total_w[j] * start_f
              right <- xi - total_w[j] / 2 + total_w[j] * end_f
              mid   <- (left + right) / 2
              left_i  <- mid + (left  - mid) * point_inset
              right_i <- mid + (right - mid) * point_inset
              x_vals_pt[i] <- stats::runif(1, left_i, right_i)
            }
          }

          df_rep$x <- x_vals_pt
          pts[[k]] <- df_rep
        }

      } else {
        # ---- per-replicate slots (preserve = "rep") ----
        slot_w <- base_w / n_rep

        for (k in seq_along(reps)) {
          rep_id <- reps[k]
          d_k <- dens_mat[, k]
          df_rep <- df_x[df_x$rep == rep_id, , drop = FALSE]
          if (!any(d_k > 0) || nrow(df_rep) == 0L) next

          d_scaled <- d_k / max(d_k)

          center_frac <- (k - 0.5) / n_rep - 0.5
          x_center <- xi + base_w * center_frac

          x_vals_pt <- numeric(nrow(df_rep))

          for (i in seq_len(nrow(df_rep))) {
            y_i <- df_rep$y[i]
            if (!is.finite(y_i)) next

            y_i <- min(max(y_i, y_min), y_max)
            j   <- which.min(abs(y_grid - y_i))

            hwidth <- (slot_w * d_scaled[j]) / 2
            if (hwidth <= 0) hwidth <- slot_w * 0.05

            left  <- x_center - hwidth
            right <- x_center + hwidth
            mid   <- (left + right) / 2
            left_i  <- mid + (left  - mid) * point_inset
            right_i <- mid + (right - mid) * point_inset
            x_vals_pt[i] <- stats::runif(1, left_i, right_i)
          }

          df_rep$x <- x_vals_pt
          pts[[k]] <- df_rep
        }
      }

      out[[idx_out]] <- do.call(rbind, pts[!vapply(pts, is.null, logical(1))])
      idx_out <- idx_out + 1L
    }

    out_df <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
    rownames(out_df) <- NULL
    out_df
  }
)

stat_superpoint <- function(mapping = NULL, data = NULL,
                            geom = "point", position = "identity",
                            ...,
                            bw          = NULL,
                            n           = 128,
                            trim        = TRUE,
                            width       = 0.9,
                            preserve    = c("pooled", "total", "rep"),
                            point_inset = 0.9,
                            na.rm       = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {

  preserve <- match.arg(preserve)

  layer(
    data        = data,
    mapping     = mapping,
    stat        = StatSuperpoint,
    geom        = geom,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      bw          = bw,
      n           = n,
      trim        = trim,
      width       = width,
      preserve    = preserve,
      point_inset = point_inset,
      na.rm       = na.rm,
      ...
    )
  )
}
