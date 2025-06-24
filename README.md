# RNA-seq Expression Database Viewer

A Shiny web application for interactive exploration, filtering, and visualization of RNA-seq gene expression data.

## Features

- **Load Data:** Upload your RNA-seq `.rds` file and annotation file (`.csv`, `.tsv`, `.txt`, `.xlsx`, `.xls`).
- **Gene Expression:** Filter samples by material, controls, family, or individual. Search and visualize gene expression.
- **Housekeeper Expression:** Visualize expression of housekeeper genes across selected samples.
- **Download:** Export plots and updated data for downstream analysis.
- **Annotation Guidance:** Sidebar tabs provide clear instructions and annotation file structure requirements.

## Annotation File Requirements

Your annotation file should be a tab-separated (`.tsv`), comma-separated (`.csv`), or Excel file (`.xlsx`, `.xls`) with the following columns:

- **InternalID:** Unique identifier for each sample (must match column names in the RDS file).
- **IndividualID:** Identifier for the individual/sample.
- **Family:** Family or group identifier.
- **Material:** Sample material type (e.g., blood, fibroblast, muscle).
- **Salmon:** (Optional) Path to quantification files for new samples.
- **Affected:** (Optional) Status of the sample (e.g., affected/unaffected).  
  *Column name can be "Affected" or "affected".*
- Additional columns are allowed but not required.

## Getting Started

### 1. Clone the Repository

```sh
git clone https://github.com/YOURUSERNAME/GUI_RNAseq_expression_data.git
cd GUI_RNAseq_expression_data
```

### 2. Set Up the R Environment

This app uses [`renv`](https://rstudio.github.io/renv/) for reproducible environments.

- Open R or RStudio in this directory.
- Run:

```r
install.packages("renv") # if not already installed
renv::restore()
```

This will install all required R packages.

### 3. Launch the App

You can run the app with:

```r
shiny::runApp()
```

Or, if you have a launch script (e.g., `run_shiny_app.bat` or `run_shiny_app.command`), double-click it or run from the command line.

### 4. Using the App

- Go to the **Load data** tab and upload your `.rds` and annotation files.
- Use the **Gene Expression** and **Housekeeper Expression** tabs to explore and visualize your data.
- Download plots or updated data as needed.

## Project Structure

```
.
├── app.R
├── modules/
│   ├── helpers.R
│   ├── load_data.R
│   ├── gene_violinplot.R
│   ├── gene_selector_module.R
│   ├── filtering_annot_module.R
│   ├── filtering_data_module.R
│   ├── housekeeper_module.R
├── complementary_scripts/
│   ├── make_initial_rds.R
│   ├── make_tx2gene.R
├── data/
│   └── (your data files, e.g., tx2gene.rds)
├── renv/
├── renv.lock
├── .gitignore
└── README.md
```

## Notes

- The `data/` folder and any `.bat` files are excluded from version control via `.gitignore`.
- For reproducibility, always use the provided `renv.lock` and `renv` environment.
- For questions or issues, see the documentation or contact the app maintainer.

## License

[MIT License](LICENSE) (or specify your license here)

---

*Created by [Your Name or Lab/Group].*