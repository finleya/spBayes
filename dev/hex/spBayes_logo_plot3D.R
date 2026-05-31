rm(list = ls())

library(scico)
library(plot3D)

set.seed(8)

# ---- Logo controls ---------------------------------------------------------

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = TRUE))
} else {
  getwd()
}

out_file <- file.path(script_dir, "spBayes_gp_surface_plot3D.png")
img_width <- 2200
img_height <- 1800
gp_grid_n <- 96          # Number of GP sample locations; changing this changes the realization.
surface_render_n <- 150  # Render/interpolation resolution; increase this to smooth facets.

color_palette <- "roma" # Any scico palette name, e.g., "romaO", "batlow", "vikO"
surface_alpha <- 1        # Surface opacity; background stays transparent.
surface_brightness <- 1 # Darkens the plot3D surface colors; 1 leaves palette unchanged.

# plot3D camera and surface controls.
# These affect spBayes_gp_surface_plot3D.png first; the hex sticker then embeds that PNG.
# The sticker stage uses a fresh temporary copy of this PNG so cached image paths
# do not make the hex sticker reuse an older surface.
plot3d_theta <- -38       # Rotate around the vertical axis. Try changing by 20-40 degrees.
plot3d_phi <- 42          # Tilt view: smaller tilts toward you; larger is more top-down.
plot3d_expand <- 0.78     # Additional plot3D z exaggeration; surface_height is usually more obvious.
plot3d_ltheta <- -35      # Light direction if plot3d_lighting is TRUE.
plot3d_lphi <- 55         # Light elevation if plot3d_lighting is TRUE.
plot3d_lighting <- FALSE  # FALSE gives a flatter/matte surface; TRUE uses controls below.
plot3d_ambient <- 0.32     # Base light level when plot3d_lighting is TRUE.
plot3d_diffuse <- 0.75     # Directional matte light strength.
plot3d_specular <- 0.05      # Highlight strength.
plot3d_exponent <- 0.18     # Highlight tightness; larger gives smaller/brighter highlights.
plot3d_sr <- 1            # Specular reflection mix; 0 is white shine, 1 tints shine by surface color.
plot3d_shade <- NA        # Use NA with lighting FALSE; set numeric values only if using lighting.

surface_height <- 1.55    # Main vertical GP scale before plotting; lower to flatten the surface.

mesh_every <- 8           # Draw every nth grid line over the surface.
mesh_color <- "#000000"
mesh_alpha <- 0.3
mesh_lwd <- 0.7

sticker_file <- file.path(script_dir, "spBayes_hex.png")
sticker_hex_fill <- "#000000"
sticker_hex_border <- "#000000"
sticker_hex_border_size <- 0.58

# Surface image placement inside the hex sticker.
# These do NOT change the 3D camera; they only place/scale the already-rendered PNG in the hex.
sticker_image_x <- 1.00
sticker_image_y <- 1.14
sticker_image_width <- 0.96  # Main control for surface size inside the hex.
sticker_image_height <- 0.74 # Less important for PNG subplots, but keep proportional.

sticker_text <- "spBayes"
sticker_text_x <- 0.97
sticker_text_y <- 0.48
sticker_text_left <- "sp"
sticker_text_right <- "Bayes"
sticker_text_color_offset <- 0.28 # Distance from palette center; must be between 0 and 0.5.
sticker_text_join_x <- NULL        # NULL centers the combined two-color word.
sticker_text_width <- 1.0          # Approximate full word width in sticker coordinates.
sticker_text_family <- "Latin Modern Roman"
sticker_text_face <- "plain"
sticker_text_size <- 45
sticker_panel_pad <- 0.008    # Extra plotting-room so the hex border is not clipped.
sticker_padding_px <- 5       # Transparent padding added evenly around final sticker.

font_files <- c(
  "Latin Modern Roman" = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-regular.otf",
  "Liberation Mono" = "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf"
)

if (requireNamespace("sysfonts", quietly = TRUE) &&
    requireNamespace("showtext", quietly = TRUE) &&
    sticker_text_family %in% names(font_files) &&
    file.exists(font_files[[sticker_text_family]])) {
  sysfonts::font_add(sticker_text_family, regular = font_files[[sticker_text_family]])
  showtext::showtext_auto(TRUE)
}

gp_grid <- function(n = 96, sigma.sq = 1, range = 0.26, nugget = 1e-8) {
  x <- seq(-1, 1, length.out = n)
  y <- seq(-1, 1, length.out = n)

  # A squared-exponential GP on a regular grid can be sampled efficiently
  # with separable covariance matrices, avoiding a huge dense 2-D covariance.
  cov_1d <- function(u) {
    D <- as.matrix(dist(u))
    sigma.sq * exp(-0.5 * (D / range)^2)
  }

  Lx <- chol(cov_1d(x) + diag(nugget, n))
  Ly <- chol(cov_1d(y) + diag(nugget, n))
  z <- t(Lx) %*% matrix(rnorm(n * n), n, n) %*% Ly

  z <- z - min(z)
  z <- z / max(z)
  z <- z^1.18
  z <- z - 0.45
  z <- z * surface_height

  list(x = x, y = y, z = z)
}

resample_surface <- function(gp, n = surface_render_n) {
  if (is.null(n) || n <= length(gp$x)) {
    return(gp)
  }

  x_new <- seq(min(gp$x), max(gp$x), length.out = n)
  y_new <- seq(min(gp$y), max(gp$y), length.out = n)

  # Interpolate the already sampled GP. This increases render resolution
  # without drawing a new random surface.
  z_x <- apply(
    gp$z,
    2,
    function(z_col) {
      approx(gp$x, z_col, xout = x_new, ties = "ordered")$y
    }
  )
  z_new <- t(apply(
    z_x,
    1,
    function(z_row) {
      approx(gp$y, z_row, xout = y_new, ties = "ordered")$y
    }
  ))

  list(x = x_new, y = y_new, z = z_new)
}

surface_cols <- function(z, n_col = 256) {
  pal <- scico(n_col, palette = color_palette)
  z_scaled <- (z - min(z)) / diff(range(z))
  pal[pmax(1, pmin(n_col, floor(z_scaled * (n_col - 1)) + 1))]
}

darken_cols <- function(cols, brightness = surface_brightness) {
  rgb <- matrix(grDevices::col2rgb(cols), nrow = 3) / 255
  rgb <- matrix(pmax(0, pmin(1, rgb * brightness)), nrow = 3)
  grDevices::rgb(rgb[1, ], rgb[2, ], rgb[3, ])
}

darkest_right_palette_color <- function(n_col = 256) {
  pal <- scico(n_col, palette = color_palette)
  right <- pal[seq(floor(n_col / 2), n_col)]
  rgb <- grDevices::col2rgb(right)
  luminance <- 0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  right[which.min(luminance)]
}

plot3d_lighting_control <- function() {
  if (!isTRUE(plot3d_lighting)) {
    return(plot3d_lighting)
  }

  list(
    type = "light",
    ambient = plot3d_ambient,
    diffuse = plot3d_diffuse,
    specular = plot3d_specular,
    exponent = plot3d_exponent,
    sr = plot3d_sr,
    alpha = surface_alpha
  )
}

sticker_text_colors <- function(n_col = 256) {
  offset <- max(0, min(0.5, sticker_text_color_offset))
  pal <- scico(n_col, palette = color_palette)
  positions <- c(0.5 - offset, 0.5 + offset)
  idx <- pmax(1, pmin(n_col, round(positions * (n_col - 1)) + 1))
  pal[idx]
}

draw_with_plot3D <- function(gp, file) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file, width = img_width, height = img_height,
                  units = "px", background = "transparent", res = 300)
  } else {
    png(file, width = img_width, height = img_height,
        bg = "transparent", res = 300)
  }
  on.exit(dev.off(), add = TRUE)

  par(mar = rep(0, 4), bg = NA)

  z <- gp$z
  pal <- darken_cols(scico(512, palette = color_palette))
  plot3D::persp3D(
    x = gp$x,
    y = gp$y,
    z = z,
    colvar = z,
    theta = plot3d_theta,
    phi = plot3d_phi,
    scale = FALSE,
    expand = plot3d_expand,
    col = pal,
    border = NA,
    shade = plot3d_shade,
    lighting = plot3d_lighting_control(),
    alpha = surface_alpha,
    ltheta = plot3d_ltheta,
    lphi = plot3d_lphi,
    bty = "n",
    axes = FALSE,
    ticktype = "simple",
    colkey = FALSE,
    xlab = NA,
    ylab = NA,
    zlab = NA
  )

  if (mesh_every > 0) {
    mesh_col <- grDevices::adjustcolor(mesh_color, alpha.f = mesh_alpha)
    mesh_idx_x <- seq(1, length(gp$x), by = mesh_every)
    mesh_idx_y <- seq(1, length(gp$y), by = mesh_every)

    for (j in mesh_idx_y) {
      plot3D::lines3D(
        x = gp$x,
        y = rep(gp$y[j], length(gp$x)),
        z = z[, j],
        add = TRUE,
        colvar = NULL,
        col = mesh_col,
        lwd = mesh_lwd,
        colkey = FALSE
      )
    }
    for (i in mesh_idx_x) {
      plot3D::lines3D(
        x = rep(gp$x[i], length(gp$y)),
        y = gp$y,
        z = z[i, ],
        add = TRUE,
        colvar = NULL,
        col = mesh_col,
        lwd = mesh_lwd,
        colkey = FALSE
      )
    }
  }

  invisible(TRUE)
}

save_sticker_plot <- function(sticker_plot, output_file, panel_pad) {
  center <- 1
  radius <- 1
  half_width <- sqrt(3) / 2 * radius
  border_room <- sticker_hex_border_size * 0.04
  x_pad <- half_width * panel_pad + border_room
  y_pad <- radius * panel_pad + border_room

  sticker_plot <- sticker_plot +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      limits = c(center - half_width - x_pad, center + half_width + x_pad)
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(center - radius - y_pad, center + radius + y_pad)
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, unit = "lines"))

  ggplot2::ggsave(
    filename = output_file,
    plot = sticker_plot,
    width = 43.9,
    height = 50.8,
    units = "mm",
    bg = "transparent",
    dpi = 600
  )
}

add_sticker_padding <- function(input_file, output_file, padding_px) {
  if (!requireNamespace("magick", quietly = TRUE) || padding_px <= 0) {
    file.copy(input_file, output_file, overwrite = TRUE)
    return(invisible(output_file))
  }

  img <- magick::image_read(input_file)
  info <- magick::image_info(img)
  geometry <- sprintf("%dx%d", info$width + 2 * padding_px, info$height + 2 * padding_px)
  img <- magick::image_extent(img, geometry = geometry, gravity = "center", color = "none")
  magick::image_write(img, path = output_file, format = "png")

  invisible(output_file)
}

add_sticker_text <- function(sticker_plot) {
  cols <- sticker_text_colors()
  split <- if (is.null(sticker_text_join_x)) {
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      widths <- systemfonts::string_width(
        c(sticker_text_left, sticker_text_right),
        family = sticker_text_family,
        size = sticker_text_size
      )
      widths[1] / sum(widths)
    } else {
      nchar(sticker_text_left) / nchar(paste0(sticker_text_left, sticker_text_right))
    }
  } else {
    0.5
  }
  join_x <- if (is.null(sticker_text_join_x)) {
    sticker_text_x + (split - 0.5) * sticker_text_width
  } else {
    sticker_text_join_x
  }
  sticker_plot +
    ggplot2::annotate(
      "text",
      x = join_x,
      y = sticker_text_y,
      label = sticker_text_left,
      hjust = 1,
      color = cols[1],
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    ) +
    ggplot2::annotate(
      "text",
      x = join_x,
      y = sticker_text_y,
      label = sticker_text_right,
      hjust = 0,
      color = cols[2],
      family = sticker_text_family,
      fontface = sticker_text_face,
      size = sticker_text_size
    )
}

gp <- resample_surface(gp_grid(gp_grid_n), surface_render_n)

if (is.null(sticker_hex_border)) {
  sticker_hex_border <- darkest_right_palette_color()
}

draw_with_plot3D(gp, out_file)

message("Wrote surface PNG: ", normalizePath(out_file, mustWork = FALSE))

sticker_surface_tmp <- tempfile(fileext = ".png")
sticker_tmp <- tempfile(fileext = ".png")
sticker_panel_tmp <- tempfile(fileext = ".png")

if (!file.copy(out_file, sticker_surface_tmp, overwrite = TRUE)) {
  stop("Could not copy the rendered surface PNG for the sticker.")
}

sticker_plot <- hexSticker::sticker(
  subplot = sticker_surface_tmp,
  s_x = sticker_image_x,
  s_y = sticker_image_y,
  s_width = sticker_image_width,
  s_height = sticker_image_height,
  package = "",
  p_x = sticker_text_x,
  p_y = sticker_text_y,
  p_color = "transparent",
  p_family = sticker_text_family,
  p_fontface = sticker_text_face,
  p_size = 0,
  h_fill = sticker_hex_fill,
  h_color = sticker_hex_border,
  h_size = sticker_hex_border_size,
  white_around_sticker = FALSE,
  filename = sticker_tmp,
  dpi = 600
)

sticker_plot <- add_sticker_text(sticker_plot)
save_sticker_plot(sticker_plot, sticker_panel_tmp, sticker_panel_pad)
add_sticker_padding(sticker_panel_tmp, sticker_file, sticker_padding_px)

unlink(out_file)
message("Wrote hex sticker PNG: ", normalizePath(sticker_file, mustWork = FALSE))
