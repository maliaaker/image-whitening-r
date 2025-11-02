
rm(list = ls())


library(OpenImageR)
library(plotly)
library(ggplot2)


setwd("~/Desktop/ProgrammingR/Projects")

getwd()
file.exists("ImageR.jpg")

# 1) read the image ------------------------------------------------------------
photo_name <- "ImageR.jpg"  
pic <- readImage("ImageR.jpg")
dim(pic)
imageShow(pic)

n_rows <- nrow(pic)
n_cols <- ncol(pic)

# split into RGB vectors (pixels x 3)
chan_r <- as.vector(pic[, , 1])
chan_g <- as.vector(pic[, , 2])
chan_b <- as.vector(pic[, , 3])
X_rgb  <- cbind(chan_r, chan_g, chan_b)

# 2) whitening by hand (like in class) ----------------------------------------
X_center <- scale(X_rgb, center = TRUE, scale = FALSE)   # center only
covX     <- cov(X_center)
eigX     <- eigen(covX)

eps      <- 1e-8  # small safety for divisions
D_inv2   <- diag(1 / sqrt(pmax(eigX$values, eps)))
W_whiten <- eigX$vectors %*% D_inv2
Z_white  <- X_center %*% W_whiten

# quick checks (should be ~0 means and ~I covariance)
apply(Z_white, 2, mean)
cov(Z_white)

# 3) small sample to search faster (create Z_sub) -----------------------------
set.seed(123)
n_take   <- max(1, round(0.03 * nrow(Z_white)))
pick_idx <- sample(1:nrow(Z_white), size = n_take)
Z_sub    <- Z_white[pick_idx, , drop = FALSE]

# 4) quick try: one direction + histogram -------------------------------------
v_try <- c(1, -1, 0); v_try <- v_try / sqrt(sum(v_try^2))
proj_try <- as.vector(Z_sub %*% v_try)

ggplot(data.frame(proj_try), aes(x = proj_try)) +
  geom_histogram(bins = 50) +
  labs(title = "Initial projection  v = (1, -1, 0)",
       x = "Projected value", y = "Count")

# 5) GRID SEARCH (theta = 1..360, phi = 1..180) -------------------------------
if (!exists("Z_sub")) {
  set.seed(123)
  n_take   <- max(1, round(0.03 * nrow(Z_white)))
  pick_idx <- sample(1:nrow(Z_white), size = n_take)
  Z_sub    <- Z_white[pick_idx, , drop = FALSE]
}

ang_theta <- 360   # azimuth range
ang_phi   <- 180   # polar range

grid_scores <- data.frame(
  az_deg = integer(), el_deg = integer(),
  vx = numeric(), vy = numeric(), vz = numeric(),
  Jval = numeric()
)

row_id <- 1
for (az in 1:ang_theta) {
  for (el in 1:ang_phi) {
    # radians
    az_r <- az * pi / 180
    el_r <- el * pi / 180
    
    # unit vector from spherical coords
    vx <- cos(az_r) * sin(el_r)
    vy <- sin(az_r) * sin(el_r)
    vz <- cos(el_r)
    
    # project subsample
    proj_vec <- as.vector(Z_sub %*% c(vx, vy, vz))
    
    # k-means (2 groups) + Fisher index
    km_fit <- kmeans(proj_vec, centers = 2, nstart = 10)
    grp1   <- proj_vec[km_fit$cluster == 1]
    grp2   <- proj_vec[km_fit$cluster == 2]
    
    m1 <- mean(grp1); m2 <- mean(grp2)
    v1 <- var(grp1);  v2 <- var(grp2)
    if (is.na(v1) || v1 == 0) v1 <- eps
    if (is.na(v2) || v2 == 0) v2 <- eps
    
    J_now <- (m1 - m2)^2 / (v1 + v2)
    
    grid_scores[row_id, ] <- list(az, el, vx, vy, vz, J_now)
    row_id <- row_id + 1
  }
}

# best direction
best_idx <- which.max(grid_scores$Jval)
best_row <- grid_scores[best_idx, ]
best_row  # shows angles and (vx, vy, vz)

v_best <- c(best_row$vx, best_row$vy, best_row$vz)
v_best <- v_best / sqrt(sum(v_best^2))

# optional: surface view -------------------------------------------------------
score_matrix <- matrix(grid_scores$Jval, nrow = ang_phi, ncol = ang_theta, byrow = TRUE)
p3d <- plot_ly(x = 1:ang_theta, y = 1:ang_phi, z = score_matrix)
p3d <- add_surface(p3d)
p3d <- layout(
  p3d,
  scene = list(
    xaxis = list(title = "theta (deg)"),
    yaxis = list(title = "phi (deg)"),
    zaxis = list(title = "Fisher J")
  ),
  title = "Fisher index surface (grid search)"
)
p3d

# 6) segment ALL pixels with v_best -------------------------------------------
proj_full <- as.vector(Z_white %*% v_best)
km_full   <- kmeans(proj_full, centers = 2, nstart = 10)

seg_map  <- matrix(km_full$cluster, nrow = n_rows, ncol = n_cols)
seg_view <- (seg_map - min(seg_map)) / (max(seg_map) - min(seg_map))  # 0..1
imageShow(seg_view)


