# Spatial Exploratory Data Analysis: SOM, PCA & HC

R code for the paper **"Uncovering Spatial Patterns in Archaeological Data through the Integration of Self-Organizing Maps, Principal Component Analysis and Hierarchical Clustering."**

The analysis applies three complementary approaches to explore spatial patterning among early medieval settlement sites in the Upper Arlanza Basin (northern Spain):

- **Principal Component Analysis (PCA)** — dimensionality reduction and variable contribution analysis
- **Hierarchical Clustering (HC)** — Ward's method with Euclidean distance
- **Self-Organizing Maps (SOM)** — a 4×5 hexagonal Kohonen network

The three methods are compared against one another to assess how well each captures underlying structure in the site-location data, and results are exported for further cartographic analysis in GIS.

## Citation

If you use this code, please cite:

TBC

## Repository contents

```
├── code.R   # Full analysis script (PCA → HC → SOM)
├── data/
│   └── data.csv   # Site sampling data (replace with your own)
├── figs/
│   ├── pca/                # PCA figures (scree plot, biplot, contributions)
│   ├── dend/                # Hierarchical clustering figures (elbow, dendrogram, heatmap)
│   ├── som/                 # SOM figures (training, U-matrix, cluster map, codebook heatmaps)
│   └── 5_comparative_map/  # Exported GeoPackage for GIS cartography
└── README.md
```


## Data

The analysis expects a CSV file with one row per site and the following variables (column names as used in the script, followed by their description):

| Variable | Description |
|---|---|
| `HubDist` | Distance to nearest site |
| `dem05_med` | Altitude |
| `current_med` | Connectivity |
| `aspect05_med` | Aspect |
| `slope05_med` | Slope |
| `distRiver_med` | Distance to nearest river |
| `visInOut_med` | Visibility index (outward) |
| `visInIng_med` | Visibility index (inward) |
| `prom10k_med` | Prominence |
| `rough_med` | Roughness |
| `globRad_med` | Solar radiance |

Site coordinates (`x`, `y`) cannot be provided following the regulations of the archaeological inventory from Junta de Castilla y León (Spain).

The published analysis uses 57 sites and 11 standardized variables.

## Requirements

- R (≥ 4.0 recommended)
- R packages:
  `readr`, `ggplot2`, `ggfortify`, `FactoMineR`, `factoextra`, `dplyr`, `dendextend`, `ggdendro`, `patchwork`, `kohonen`, `RColorBrewer`, `tidyr`, `sf`

Install all dependencies:

```r
install.packages(c("readr", "ggplot2", "ggfortify", "FactoMineR", "factoextra",
                    "dplyr", "dendextend", "ggdendro", "patchwork",
                    "kohonen", "RColorBrewer", "tidyr", "sf"))
```

## Usage

1. Place your sampling data CSV in the working directory and update `setwd()` and the `read.csv()` call at the top of the script.
2. Update `vars_keep` and `var_labels` if your variable names differ.
3. Run the script section by section — it is organized into three self-contained parts (PCA, hierarchical clustering, SOM), each generating its own figures.
4. Figure-generating code is commented with the corresponding figure number from the paper (e.g. `# FIGURE 3`), making it straightforward to locate and reproduce any individual plot.
5. The final section exports a GeoPackage (`sites_classified.gpkg`) containing site coordinates, PCA scores, and cluster assignments from both HC and SOM, for further mapping in GIS software (e.g. QGIS).

## Analysis parameters

- SOM grid: 4×5 hexagonal, Gaussian neighbourhood function, 400 training iterations
- Clustering: four-cluster solution (Ward's method for HC; hierarchical clustering of SOM codebook vectors for SOM)
- Cluster number selected via elbow method (within-cluster sum of squares)

## License

GPL-3.0 license

## Contact

Guillem Domingo-Ribas
[Guillem.Domingo-Ribas@newcastle.ac.uk]
