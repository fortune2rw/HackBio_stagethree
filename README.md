# 🧬 Simulated Cell Type Viewer

This Shiny app lets you explore a simulated single-cell RNA‑seq dataset, visualize cells on a UMAP, and automatically identify a “best” marker gene for each cell type based on a simple specificity score.

---

## ▶️ How to run the app

1. Make sure you have R (and preferably RStudio) installed.
2. Install the required packages (only needed once):

   ```r
   install.packages(c("shiny", "dplyr"))
