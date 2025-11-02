# Image Whitening & Segmentation (R)

**Author:** Malia Aker  
**Tools:** R | OpenImageR | ggplot2 | plotly | k-means | Linear Algebra  

---

## 🧠 Overview
This project performs **image whitening and segmentation** on a dermoscopy image using pixel-level RGB data.  
It applies whitening transformations, explores projections that maximize **Fisher’s discriminant index**, and segments the image using **k-means clustering**.

The work was inspired by independent component analysis (ICA) concepts used in pattern recognition and computer vision.

---

## ⚙️ Steps
1. **Image Loading & Preprocessing** – Read and split RGB channels into pixel vectors.  
2. **Whitening Transformation** – Apply centering and eigenvalue decomposition to decorrelate features.  
3. **Grid Search** – Explore 3D projection directions (θ, φ) to maximize Fisher index.  
4. **Visualization** – Render Fisher index surface using `plotly`.  
5. **Segmentation** – Cluster whitened pixels into two groups using k-means, producing a binary segmentation map.

---

## 📊 Example Outputs
| Visualization | Description |
|----------------|-------------|
| ![Fisher Surface](reports/figures/fisher_surface.png) | Fisher index surface showing discriminant strength across projection angles |
| ![Segmented Image](reports/figures/segmented_output.png) | Segmented binary image (lesion vs. background) after optimal projection |

*(Note: sample plots can be saved from your script using `ggsave()` or `export(p3d)`.)*

---

## 🧩 How to Run
```r
# Install dependencies
install.packages(c("OpenImageR", "ggplot2", "plotly"))

# Run the script
source("image_whitening_segmentation.R")
