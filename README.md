# 🧬 Simulated Cell Type Viewer

Interactive **Shiny app** for exploring simulated single-cell RNA-seq data, visualising UMAP embeddings, and identifying marker genes for selected cell types.

---

## 🚀 Running the App

### 1️⃣ Install Required Packages (one-time)
```r
install.packages(c("shiny", "dplyr"))
```

### 2️⃣ Save Script
Save the application code as:

```
app.R
```

### 3️⃣ Launch App
```r
shiny::runApp("app.R")
```

The browser interface includes:

- 📍 Cell type dropdown selector  
- 📊 Gene statistics table  
- 🖼️ UMAP plot  
- 📝 Marker gene summary  

---

## 🔢 Gene Specificity Score

**Definition**
```
diff = mean_in − mean_out
```

Where:

- `mean_in`  = average expression inside selected cell type  
- `mean_out` = average expression outside selected cell type  

### Interpretation

| diff value | Meaning |
|-----------|--------|
| > 0 | Good marker (higher in target type) |
| ≈ 0 | Not specific |
| < 0 | Marker for other cell types |

### Additional Statistics

- `det_in`   = % of cells in type expressing gene (>0)
- `det_out`  = % of cells outside type expressing gene
- `mean_in`  = mean expression inside type
- `mean_out` = mean expression outside type

---

## 🧬 Marker Gene Selection

Function: `pick_marker_gene()`

**Selection Strategy**

1. Sort genes by `diff` (descending) → maximize specificity  
2. Break ties using `det_in` (descending) → maximize detection robustness  

**Result:**  
Top gene with highest specificity and strongest detection consistency.

---

## 📈 App Workflow

```
Select cell type
      ↓
Recompute gene statistics
      ↓
Identify best marker gene
      ↓
Color UMAP by marker expression
      ↓
Display statistics table + summary
```

---

## 🎯 Purpose

This app demonstrates:

- Marker gene specificity logic
- Expression detection metrics
- Cluster-level vs gene-level visualisation
