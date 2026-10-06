# Hands-on Tutorial: Quality Control, Population Structure, and Multi-Model GWAS

This tutorial walks through the bovine GWAS workflow used for the workshop. Run the commands in order and keep all generated output files in your own `~/workshop` directory unless otherwise stated.  

<sub>*For shorter commands, it is better to type out each command manually rather than copying directly. This helps you get accustomed to the bash and R environments.*</sub>

## Workflow overview

1. Quality control filtering on call rates and MAF.
2. Hardy-Weinberg equilibrium (HWE) analysis with full-distribution and zoomed-tail plots, followed by data-informed HWE filtering.
3. PCA calculation.
4. Scree plot generation, metadata PCA plots, covariate testing, and covariate-file export.
5. Whole-genome representation including chromosome X.
6. In class, run **additive GWAS only** for the **unadjusted** and **PC1 + PC2** configurations.
7. Generate matching Manhattan and Q-Q plots and compare the effect of population-structure adjustment.
8. Use the flexible GWAS script later to run individual covariates, custom covariate combinations, or the full homework analysis across additive, dominant, and recessive models.


---

> ⚠️ **Before starting:** Make sure you completed the setup in [`01_Before_Class_Setup.md`](01_Before_Class_Setup.md), are connected to the WSU network or VPN, and can log in to the workshop server.

> **Working directory:** Before starting, make sure you are in your own workshop directory.

```bash
cd ~/workshop     # cd = change directory; directory means folder
pwd               # pwd = print working directory
ls                # ls = list contents
```

<details>
<summary><strong>Click to view some useful bash commands</strong></summary>

<br>

Below is a list of useful commands 😏

```bash
# --- Navigation & Path Inspection ---
pwd                                 # Print the absolute path of the current working directory
cd my_folder                        # Change directory into 'my_folder'
cd ..                               # Move up one directory level (parent directory)
cd ~                                # Jump directly to your user's home directory
cd -                                # Switch back to the previous directory you were in

# --- Listing Files & Folders ---
ls                                  # List names of files and folders in the current directory
ls -l                               # Detailed list view showing permissions, owner, file size, and modification date
ls -lh                              # Detailed list view with human-readable file sizes (K, M, G)
ls -la                              # List all files, including hidden files and dotfiles (names starting with .)
ls -lt                              # List files sorted by last modified time (newest first)

# --- Viewing File Content ---
cat file.txt                        # Print the entire contents of file.txt to the terminal
less file.txt                       # Open a scrollable, paginated viewer (press 'q' to quit, '/' to search)
head file.txt                       # Display the first 10 lines of file.txt
head -n 25 file.txt                 # Display the first 25 lines of file.txt
tail file.txt                       # Display the last 10 lines of file.txt
tail -n 20 file.txt                 # Display the last 20 lines of file.txt
tail -f process.log                 # Follow file updates in real-time as new lines are appended

# --- File & Directory Management ---
mkdir new_folder                    # Create a new directory named 'new_folder'
mkdir -p path/to/nested/folder      # Create nested parent directories automatically without error
touch new_file.txt                  # Create an empty file or update the timestamp of an existing file
cp source.txt backup.txt            # Copy 'source.txt' to a new file called 'backup.txt'
cp -r folder_a folder_backup        # Recursively copy an entire folder and its contents
mv old_name.txt new_name.txt        # Rename a file or directory
mv file.txt /path/to/destination/   # Move a file to another location
rm unwanted_file.txt                # Permanently delete a file (no trash bin, cannot undo)
rm -r unwanted_folder               # Recursively delete a directory and all files inside it

# --- Searching & Inspection ---
wc -l file.txt                      # Count the total number of lines in file.txt
grep "pattern" file.txt             # Search and print lines containing 'pattern' in file.txt
grep -i "pattern" file.txt          # Case-insensitive search for 'pattern' in file.txt
grep -c "pattern" file.txt          # Count the number of lines matching 'pattern' in file.txt
find . -name "*.bed"                # Search the current directory recursively for files ending with .bed
which plink                         # Locate the executable path of a program in your system's PATH

# --- Pipes & Redirection ---
command > output.txt                # Run command and write its output to a file, overwriting existing content
command >> output.txt               # Run command and append its output to the end of a file
command1 | command2                 # Pipe: pass the stdout output of command1 directly into command2 as input
cat file.txt | grep "Chr1" | wc -l  # Count how many lines in file.txt contain "Chr1"

# --- System, Memory & Process Monitoring ---
df -h                               # Show available and used disk space across mounted filesystems
du -sh folder_name                  # Show the total disk space consumed by folder_name
free -h                             # Show total, used, and available system RAM memory
top                                 # Open an interactive monitor of CPU and RAM usage by running processes
htop                                # An enhanced, color-coded interactive process monitor (if installed)
kill 12345                          # Terminate the process with Process ID (PID) 12345
history                             # Display a numbered list of previously executed commands
clear                               # Clear terminal window output (shortcut: Ctrl + L)
```

</details>

The main genotype input used in this tutorial is:

```text
/workshop/data/SRD_HFL_AI_50K.ped
/workshop/data/SRD_HFL_AI_50K.map
```

The metadata file used later is:

```text
/workshop/data/SRD_HFL_AI_50K_metadata.csv
```

Before jumping into PCA and GWAS, it is a good idea to first inspect the metadata structure so you know what variables are available and what may be useful as covariates.

For example:

```bash
head /workshop/data/SRD_HFL_AI_50K_metadata.csv
```

or if you want a little more:

```bash
head -n 20 /workshop/data/SRD_HFL_AI_50K_metadata.csv
```
If you prefer to view the metadata in a spreadsheet-like format, you can use Gnumeric:

```bash
gnumeric /workshop/data/SRD_HFL_AI_50K_metadata.csv
```


This lets you quickly check the column names, overall structure, possible covariates, and the way the metadata is stored.

### What does **metadata** mean?

**Metadata means “data about the data.”** The genotype files contain the SNP genotypes for each animal, while the metadata file contains additional information describing the animals or how the observations were collected. Examples can include the sample ID, sire, birth year, technician, breeding protocol, and other study information.

These variables are important because some may influence the phenotype independently of the SNP being tested. If that happens, they may act as **potential covariates** in the GWAS.

As you inspect the first few rows, ask yourself:

1. Which column identifies each animal?
2. Which variables look numerical?
3. Which variables represent categories or groups?
4. Which variables might plausibly be associated with the phenotype or population structure?
5. Which of those variables might therefore need to be considered as GWAS covariates?

> 💡 We will explore several metadata variables during the PCA/covariate section. The in-class GWAS will focus only on the **unadjusted** and **PC1 + PC2** additive models, but the same workflow can later be used with other covariates.

---

## 1. Initial Quality Control: Call Rate and MAF

This step applies the initial SNP- and animal-level QC filters and creates a binary PLINK dataset for downstream analyses.

First, check the version and help details for plink.

```bash
plink --version
plink --help
```

If the help output is too long, pipe the output to `less`, so you can scroll pages using the space bar, and use `g`/`G` to go to the first/last page respectively. Use `q` to exit when you have finished viewing the file.

```bash
plink --help | less
```

Now that you have verified the version and viewed the help menu, run the following command:

```bash
plink \
    --file /workshop/data/SRD_HFL_AI_50K \
    --allow-no-sex \
    --maf 0.05 \
    --geno 0.10 \
    --mind 0.10 \
    --memory 4000 \
    --make-bed \
    --out srd_qc
```

This command tells PLINK to read the PED/MAP dataset, keep SNPs with minor allele frequency of at least 5% (`--maf 0.05`), remove markers with more than 10% missing genotypes (`--geno 0.10`), remove animals with more than 10% missing genotypes (`--mind 0.10`), and write the output in binary PLINK format (`.bed`, `.bim`, `.fam`) using the prefix `srd_qc`.

`--allow-no-sex` is included because sex coding is not the focus here, and we do not want missing/ambiguous sex values to stop the analysis.

### Oops! What happened? 🤭

What error message are you getting and why? Read the error message carefully.

**Question:** Why is PLINK having a problem with the chromosome numbers in this dataset?

<details>
<summary><strong>Clue 🤔🧐</strong></summary>

<p align="left">
  <img src="images/holstein.jpg" width="300" alt="Tutorial Overview" />
  <br>
  <sub><i>Source: <a href="https://www.agdaily.com/livestock/facts-about-holstein-cattle-cows/">AGDAILY</a></i></sub>
</p>

</details>

<details>
<summary><strong>Click to reveal the solution 🤫</strong></summary>

<br>

By default, PLINK assumes that the dataset is **human**. Human autosomes are numbered 1–22, while cattle have **29 autosomes**.

We therefore need to tell PLINK that we are working with cattle by adding:

```bash
--cow
```

The corrected command is:

```bash
plink \
    --cow \
    --file /workshop/data/SRD_HFL_AI_50K \
    --allow-no-sex \
    --maf 0.05 \
    --geno 0.10 \
    --mind 0.10 \
    --memory 4000 \
    --make-bed \
    --out srd_qc
```

After adding `--cow`, PLINK recognizes the bovine chromosome set and the QC analysis can proceed. 😮‍💨

</details>

### Main output

```text
srd_qc.bed
srd_qc.bim
srd_qc.fam
srd_qc.log
```

---

## 2. Hardy-Weinberg Equilibrium and Post-HWE Filtering

First, calculate HWE statistics from the QC-filtered dataset.

```bash
plink \
    --bfile srd_qc \
    --hardy \
    --out plink_results_hwinb
```

`--hardy` asks PLINK to calculate Hardy-Weinberg equilibrium statistics for each marker. HWE is useful here as a quality-control check because markers showing extreme deviation from expected genotype proportions can sometimes reflect genotyping error, batch effects, or problematic loci.

If you got an error, click below.

<details>
<summary><strong>Click here if you got an error</strong></summary>

<br>

So you fell for this again eh? 🫣

You might want to think through the error first, before peeking below 🧑‍💻🧠
</details>

<details>
<summary><strong>Do you want to reveal the correct code? 😏</strong></summary>

<br>

OK, enough teasing, here is the correct code. Just add `--cow`. 🫩

```bash
plink \
    --cow \
    --bfile srd_qc \
    --hardy \
    --out plink_results_hwinb
```
</details>

The HWE results are written to:

```text
plink_results_hwinb.hwe
```

### Plot the HWE distribution in R

For a longer R block, **Geany** is preferred because it is more flexible and easier to edit visually. If you are already comfortable in the command line, `nano` is also very fast and works well.

Create a script with **Geany**:

```bash
geany hwe_plots.R
```

If you prefer `nano`, you can use:

```bash
nano hwe_plots.R
```

Paste the R code below into the file. In `nano`, save with **Ctrl+O**, press **Enter**, then exit with **Ctrl+X**. In Geany, you can use the normal menu or keyboard shortcuts to save.

Run the completed script with:

```bash
Rscript hwe_plots.R
```

Alternatively, start an interactive R session with `R` and paste the same code directly.

```r
library(tidyverse)
library(patchwork)
library(scales)

hwe_raw <- read.table("plink_results_hwinb.hwe", header = TRUE, stringsAsFactors = FALSE)

if ("TEST" %in% colnames(hwe_raw)) {
  hwe_clean <- hwe_raw %>% filter(TEST == "ALL")
} else {
  hwe_clean <- hwe_raw
}

hwe_clean <- hwe_clean %>%
  filter(!is.na(P) & P >= 0 & P <= 1) %>%
  mutate(log10_P = -log10(ifelse(P == 0, 1e-100, P)))

n_failed <- sum(hwe_clean$P < 3e-9)

p_hwe_all <- ggplot(hwe_clean, aes(x = log10_P)) +
  geom_histogram(bins = 60, fill = "steelblue", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = -log10(3e-9), color = "red", linetype = "dashed", linewidth = 0.9) +
  scale_y_continuous(
    trans = "pseudo_log",
    breaks = c(1, 10, 100, 1000, 10000, 50000),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  annotate(
    "text",
    x = -log10(3e-9) + 0.8,
    y = 5000,
    label = paste0("p < 3e-9 Filter\n(", n_failed, " SNPs removed)"),
    color = "red",
    hjust = 0,
    fontface = "bold",
    size = 3.8
  ) +
  labs(
    title = "HWE Spectrum (Full Distribution)",
    subtitle = "Log10 y-axis highlights outlier departures",
    x = expression(-log[10](italic(p))),
    y = expression(Variant ~ Count ~ (log[10] ~ scale))
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10)),
    panel.grid.minor = element_blank()
  )

hwe_zoom_data <- hwe_clean %>% filter(log10_P >= 2)

p_hwe_zoom <- ggplot(hwe_zoom_data, aes(x = log10_P)) +
  geom_histogram(binwidth = 0.5, fill = "darkorange", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = -log10(3e-9), color = "red", linetype = "dashed", linewidth = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  annotate(
    "rect",
    xmin = -log10(3e-9), xmax = Inf, ymin = 0, ymax = Inf,
    fill = "red", alpha = 0.1
  ) +
  annotate(
    "text",
    x = -log10(3e-9) + 0.5,
    y = max(table(cut_width(hwe_zoom_data$log10_P, 0.5))) * 0.85,
    label = "Excluded SNPs",
    color = "darkred",
    hjust = 0,
    fontface = "bold",
    size = 4
  ) +
  labs(
    title = expression(Zoomed ~ HWE ~ Tail ~ and ~ Empirical ~ Breakoff),
    subtitle = "Extreme HWE departures separate from the main distribution",
    x = expression(-log[10](italic(p))),
    y = "Variant Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, margin = margin(b = 10)),
    panel.grid.minor = element_blank()
  )

hwe_final_plot <- p_hwe_all + p_hwe_zoom
ggsave("hwe_distribution_plots.png", plot = hwe_final_plot, width = 13, height = 5.5, dpi = 300)
```

The plot is saved as:

```text
hwe_distribution_plots.png
```

To view the plot, use:

```bash
feh hwe_distribution_plots.png
```

If the image is too large, play around with the `--zoom` flag to fit your desired viewing %. Here, we use 50:

```bash
feh --zoom 50 hwe_distribution_plots.png
```

> 💡 If you have trouble viewing the plot, you can download the plot using `scp` or FileZilla.

Based on the plot, what threshold should be used for HWE?

### Choosing the HWE cutoff from the observed distribution

A commonly used HWE threshold such as `1e-6` can be useful as a general rule, but it should not be treated as a universal biological boundary between a “good” and “bad” SNP.

For this exercise, look at the **extreme tail of the HWE distribution** and identify the point where a relatively small group of SNPs begins to separate sharply from the majority of markers. We can think of this as an **empirical breakoff point** or **data-informed QC threshold**.

For this dataset, that breakoff occurs at approximately:

```text
-log10(p) ≈ 8.5
```

which corresponds approximately to:

```text
p ≈ 3 × 10^-9
```

We will therefore use `3e-9` for this tutorial.

> 💡 This does **not** mean that `3e-9` is a universal HWE cutoff. The purpose is to learn how to inspect the distribution and identify unusually extreme deviations rather than automatically applying the same threshold to every dataset.

### Apply the HWE filter

After reviewing the HWE distribution and identifying the breakoff point, apply the selected threshold:

```bash
plink \
    --cow \
    --bfile srd_qc \
    --hwe 3e-9 \
    --make-bed \
    --out srd_qc_hwe
```

`--hwe 3e-9` removes markers with very extreme deviation from Hardy-Weinberg equilibrium. We are not trying to remove every slight departure from HWE; we are mainly removing the tail of markers that look unusually problematic based on the HWE distribution.

### Main output

```text
srd_qc_hwe.bed
srd_qc_hwe.bim
srd_qc_hwe.fam
srd_qc_hwe.log
```

---

## 3. PCA Calculation

Calculate principal components from the post-HWE dataset.

```bash
plink \
    --cow \
    --bfile srd_qc_hwe \
    --pca \
    --out srd_pca
```

PCA is calculated to summarize major patterns of genetic similarity and population structure in the dataset. This is important because hidden structure can confound GWAS and produce misleading association signals if not accounted for.

### Main output

```text
srd_pca.eigenval
srd_pca.eigenvec
```

---

## 4. PCA Scree Plot, Metadata PCA, Covariate Testing, and Export

This R section:

- calculates the percentage of variance explained by each PC;
- creates a PCA scree plot;
- joins PCA results to the workshop metadata;
- creates PCA plots colored by selected metadata variables;
- tests candidate covariates against the binary phenotype;
- exports a covariate file for the GWAS models.

We use this step to decide which non-genetic and structure-related variables may need to be carried forward into the GWAS models.

We also focus mainly on **PC1 and PC2** because they usually explain the largest share of the structure, and in this tutorial they are selected based on the scree plot "elbow" idea — after the first few PCs, the additional variance explained starts to level off.

Because this is a long block, save it as a script.

Geany is preferred:

```bash
geany pca_covariates.R
```

If you prefer the command line:

```bash
nano pca_covariates.R
```

Paste the code below, save the file, and run it with:

```bash
Rscript pca_covariates.R
```

```r
library(tidyverse)

eigenvalues <- read.table("srd_pca.eigenval", header = FALSE)$V1
var_explained <- (eigenvalues / sum(eigenvalues)) * 100

scree_df <- data.frame(
  PC = factor(paste0("PC", seq_along(eigenvalues)), levels = paste0("PC", seq_along(eigenvalues))),
  PC_num = seq_along(eigenvalues),
  Variance = var_explained,
  Cumulative = cumsum(var_explained)
)

scree_plot <- ggplot(scree_df, aes(x = PC_num)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.85) +
  geom_line(aes(y = Variance), color = "darkblue", linewidth = 1) +
  geom_point(aes(y = Variance), color = "darkblue", size = 2) +
  scale_x_continuous(breaks = seq_along(eigenvalues)) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA Scree Plot",
    x = "Principal Component",
    y = "% Variance Explained"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("pca_scree_plot.png", plot = scree_plot, width = 8, height = 5, dpi = 300)

pca_data <- read.table("srd_pca.eigenvec", header = FALSE)
colnames(pca_data)[1:2] <- c("FID", "IID")
num_pcs <- ncol(pca_data) - 2
colnames(pca_data)[3:ncol(pca_data)] <- paste0("PC", 1:num_pcs)

metadata <- read.csv("/workshop/data/SRD_HFL_AI_50K_metadata.csv", stringsAsFactors = FALSE)

metadata <- metadata %>%
  mutate(
    Birth_Year_Label = ifelse(as.numeric(Birth_Year) >= 2020, "Post-2020", "Pre-2020"),
    Birth_Year_Group = ifelse(as.numeric(Birth_Year) >= 2020, 2, 1)
  )

fam_data <- read.table("srd_qc_hwe.fam", header = FALSE)[, c(1, 2, 6)]
colnames(fam_data) <- c("FID", "IID", "Phenotype")

plot_data <- inner_join(pca_data, metadata, by = c("IID" = "SampleID")) %>%
  inner_join(fam_data, by = c("FID", "IID")) %>%
  mutate(
    Pheno_Binary = case_when(
      Phenotype == 2 ~ 1,
      Phenotype == 1 ~ 0,
      TRUE ~ NA_real_
    )
  )

variables_to_plot <- list(
  "Sire"              = "Sire ID",
  "Birth_Year"        = "Birth Year",
  "Birth_Year_Label"  = "Birth Year Grouping (Pre vs Post 2020)",
  "Technician"        = "Technician ID",
  "Protocol"          = "Protocol"
)

pc1_var <- round(var_explained[1], 2)
pc2_var <- round(var_explained[2], 2)

for (var_name in names(variables_to_plot)) {
  legend_title <- variables_to_plot[[var_name]]
  plot_data[[var_name]] <- as.factor(plot_data[[var_name]])

  pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = .data[[var_name]])) +
    geom_point(alpha = 0.8, size = 2.5) +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Cattle Population PCA -", legend_title),
      x = paste0("PC1 (", pc1_var, "% variance)"),
      y = paste0("PC2 (", pc2_var, "% variance)"),
      color = legend_title
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )

  out_name <- ifelse(var_name == "Birth_Year_Label", "srd_pca_Birth_Year_Group.jpg", paste0("srd_pca_", var_name, ".jpg"))

  ggsave(
    filename = out_name,
    plot = pca_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
}

test_vars <- c("Birth_Year", "Birth_Year_Group", "Technician", "Sire", "Protocol", "PC1", "PC2", "PC3")

covar_results <- list()

for (v in test_vars) {
  sub_df <- plot_data %>% filter(!is.na(Pheno_Binary) & !is.na(.data[[v]]))
  n_unique <- length(unique(sub_df[[v]]))

  if (n_unique < 2) {
    covar_results[[v]] <- data.frame(
      Covariate = v,
      Type = "Constant/Single-Level",
      Df = 0,
      Deviance = 0,
      P_Value = NA_real_,
      Significant_05 = "NO (Invariable)"
    )
    next
  }

  is_cat <- is.factor(sub_df[[v]]) || is.character(sub_df[[v]])
  var_type <- ifelse(is_cat, paste0("Factor (", n_unique, " levels)"), "Numeric")

  null_mod <- glm(Pheno_Binary ~ 1, data = sub_df, family = binomial)
  cand_formula <- as.formula(paste("Pheno_Binary ~", v))
  cand_mod <- glm(cand_formula, data = sub_df, family = binomial)

  lrt <- anova(null_mod, cand_mod, test = "Chisq")
  p_val <- lrt$`Pr(>Chi)`[2]
  dev_diff <- round(lrt$Deviance[2], 3)
  df_diff <- lrt$Df[2]

  covar_results[[v]] <- data.frame(
    Covariate = v,
    Type = var_type,
    Df = df_diff,
    Deviance = dev_diff,
    P_Value = signif(p_val, 4),
    Significant_05 = ifelse(!is.na(p_val) & p_val < 0.05, "YES", "NO")
  )
}

covar_stats <- bind_rows(covar_results)
cat("\n=== Covariate Association with Phenotype (Likelihood Ratio Test) ===\n")
print(covar_stats)
write.csv(covar_stats, "covariate_statistical_tests.csv", row.names = FALSE)

sub_pcs <- plot_data %>% filter(!is.na(Pheno_Binary) & !is.na(PC1) & !is.na(PC2))
null_pcs <- glm(Pheno_Binary ~ 1, data = sub_pcs, family = binomial)
joint_pcs <- glm(Pheno_Binary ~ PC1 + PC2, data = sub_pcs, family = binomial)
lrt_pcs <- anova(null_pcs, joint_pcs, test = "Chisq")

cat("\n=== Joint Test: Phenotype ~ PC1 + PC2 ===\n")
cat("P-value:", signif(lrt_pcs$`Pr(>Chi)`[2], 4), "\n\n")

# ------------------------------------------------------------------
# Prepare the covariate file used by PLINK
# ------------------------------------------------------------------
# PC1, PC2, Birth Year, and Birth Year Group are retained as numeric
# covariates. Sire and Technician are categorical identifiers, so we
# convert them to n-1 dummy variables before exporting the file.

covar_base <- inner_join(
  pca_data[, c("FID", "IID", "PC1", "PC2")],
  plot_data[, c("IID", "Birth_Year", "Birth_Year_Group", "Technician", "Sire")],
  by = "IID"
) %>%
  mutate(
    Birth_Year = as.numeric(as.character(Birth_Year)),
    Birth_Year_Group = as.numeric(as.character(Birth_Year_Group))
  )

make_dummy_block <- function(x, prefix) {
  if (any(is.na(x))) {
    stop(paste("Missing values found in", prefix, "while creating dummy variables."))
  }

  f <- factor(x)

  if (nlevels(f) < 2) {
    return(data.frame())
  }

  mm <- model.matrix(~ f)
  mm <- mm[, -1, drop = FALSE]   # drop intercept / retain n-1 indicators

  colnames(mm) <- paste0(
    prefix,
    "_",
    make.names(levels(f)[-1], unique = TRUE)
  )

  as.data.frame(mm, check.names = FALSE)
}

sire_dummy <- make_dummy_block(covar_base$Sire, "Sire")
technician_dummy <- make_dummy_block(covar_base$Technician, "Technician")

covar_df <- bind_cols(
  covar_base %>%
    select(FID, IID, PC1, PC2, Birth_Year, Birth_Year_Group),
  sire_dummy,
  technician_dummy
)

write.table(
  covar_df,
  "covariates.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

cat("\nCovariate file written with", ncol(covar_df) - 2, "covariate columns.\n")
cat("Sire dummy variables:", ncol(sire_dummy), "\n")
cat("Technician dummy variables:", ncol(technician_dummy), "\n")
```

If you used an interactive R session instead of `Rscript`, exit without saving the workspace:

```r
q("no")    # You can also use Ctrl + d and when prompted to save workspace, type n.
```

### Main outputs

```text
pca_scree_plot.png
srd_pca_Sire.jpg
srd_pca_Birth_Year.jpg
srd_pca_Birth_Year_Group.jpg
srd_pca_Technician.jpg
srd_pca_Protocol.jpg
covariate_statistical_tests.csv
covariates.txt
```

You can inspect the covariate test results directly in the terminal:

```bash
head covariate_statistical_tests.csv
```
or open the results as a spreadsheet:

```bash
gnumeric covariate_statistical_tests.csv
```
💡 Gnumeric is particularly useful here because the covariates, test statistics, and p-values are easier to compare when displayed as rows and columns.

> 💡 **Why were Sire and Technician dummy-coded?** Their values are identifiers for categories, not continuous measurements. For example, a larger sire ID does not represent “more sire.” The exported `covariates.txt` therefore contains `n-1` indicator variables for Sire and Technician so that these effects can be modeled as categorical covariates in PLINK. PC1, PC2, Birth Year, and Birth Year Group remain available as numeric covariates.

You can inspect the exported covariate columns with:

```bash
head -n 1 covariates.txt | tr '\t' '\n'
```



---

## 5. Enable Whole-Genome Testing: Autosomes + Chromosome X

Create the dataset used for the association models.

```bash
plink \
    --bfile srd_qc_hwe \
    --autosome-num 30 \
    --allow-extra-chr \
    --make-bed \
    --out srd_qc_allchr
```

This step prepares the genotype data for whole-genome association testing. `--autosome-num 30` tells PLINK how to handle the chromosome numbering scheme in this dataset, and `--allow-extra-chr` helps PLINK tolerate nonstandard chromosome coding beyond the usual human defaults.

Wait, why did it work 😲? Something looks different here, what is it? 🤔

<details>
<summary><strong>Clue 🤔🧐</strong></summary>

<video src="https://github.com/user-attachments/assets/5ac3ef17-9b6f-49af-b706-ac3de99d2182" controls autoplay loop playsinline preload="auto" width="100%"></video>

<sub><i>Source: <a href="https://www.tiktok.com/t/ZTUFuSbBq">TikTok</a></i></sub>

</details>

### Main output

```text
srd_qc_allchr.bed
srd_qc_allchr.bim
srd_qc_allchr.fam
```

---

## 6. Flexible Association Testing

For the **in-class exercise**, we will keep the GWAS focused and run only the **additive model** for two configurations:

1. **Unadjusted**
2. **PC1 + PC2 adjusted**

This gives us a direct comparison between a GWAS with no covariate adjustment and a GWAS adjusted for the major population-structure axes selected from the PCA.

Later, the same script can be used to run additional covariates, custom combinations, and different inheritance models.

### Covariate menu

The script will present this menu:

```text
Choose covariate configuration(s):

1. Unadjusted
2. PC1
3. PC2
4. PC1 + PC2
5. Sire
6. Birth Year
7. Birth Year Group
8. Technician
9. Birth Year + Sire
10. ALL

Use a comma (,) to run configurations separately.
Example: 1,4 runs Unadjusted and PC1 + PC2 as two separate GWAS analyses.

Use a plus sign (+) to combine covariates in the SAME GWAS model.
Example: 2+5 fits PC1 + Sire together.
Example: 6+5 fits Birth Year + Sire together.

ALL runs each predefined configuration 1-9 separately.
This includes the predefined PC1 + PC2 and Birth Year + Sire models.

Enter selection:
```

The inheritance-model menu is:

```text
Choose inheritance model(s):

1. Additive
2. Dominant
3. Recessive
4. ALL

Enter selection:
```

### How the selection syntax works

The punctuation has an important meaning:

```text
,  = separate GWAS analyses
+  = covariates combined in the same GWAS analysis
```

Examples:

| Entry | Meaning |
|---|---|
| `1` | Unadjusted only |
| `1,2` | Unadjusted and PC1 as two separate analyses |
| `2+3` | PC1 + PC2 in the same model |
| `6+5` | Birth Year + Sire in the same model |
| `1,2+3,6+5` | Three analyses: Unadjusted; PC1 + PC2; Birth Year + Sire |
| `ALL` | Run predefined configurations 1-9 separately |

> ⚠️ **Unadjusted cannot be combined with another covariate.** For example, `1+2` does not make sense because once PC1 is added, the model is no longer unadjusted.

The predefined combinations are provided for convenience:

```text
4 = PC1 + PC2
9 = Birth Year + Sire
```

Therefore, `4` and `2+3` describe the same covariate adjustment, while `9` and `6+5` describe the same adjustment.

---

### Create the flexible GWAS script

Open a new script:

```bash
geany run_gwas.sh
```

or:

```bash
nano run_gwas.sh
```

Paste:

```bash
#!/usr/bin/env bash
set -euo pipefail

COVAR_FILE="covariates.txt"
MANIFEST="gwas_run_manifest.tsv"

if [ ! -f "$COVAR_FILE" ]; then
    echo "ERROR: $COVAR_FILE was not found."
    echo "Run the PCA/covariate section first."
    exit 1
fi

# ------------------------------------------------------------
# Locate the dummy-coded categorical covariates
# ------------------------------------------------------------
HEADER=$(head -n 1 "$COVAR_FILE")

SIRE_COLS=$(printf '%s\n' "$HEADER" | tr '\t' '\n' | grep '^Sire_' | paste -sd, - || true)
TECH_COLS=$(printf '%s\n' "$HEADER" | tr '\t' '\n' | grep '^Technician_' | paste -sd, - || true)

# ------------------------------------------------------------
# Covariate menu
# ------------------------------------------------------------
echo
echo "Choose covariate configuration(s):"
echo
echo "1. Unadjusted"
echo "2. PC1"
echo "3. PC2"
echo "4. PC1 + PC2"
echo "5. Sire"
echo "6. Birth Year"
echo "7. Birth Year Group"
echo "8. Technician"
echo "9. Birth Year + Sire"
echo "10. ALL"
echo
echo "Use , to run configurations separately."
echo "Example: 1,4 = Unadjusted and PC1 + PC2 as separate runs."
echo
echo "Use + to combine covariates in the SAME model."
echo "Example: 2+5 = PC1 + Sire."
echo "Example: 6+5 = Birth Year + Sire."
echo
echo "ALL runs predefined configurations 1-9 separately."
echo
read -r -p "Enter selection: " COV_INPUT

COV_INPUT=$(echo "$COV_INPUT" | tr -d ' ')
COV_UPPER=$(echo "$COV_INPUT" | tr '[:lower:]' '[:upper:]')

if [ "$COV_UPPER" = "ALL" ] || [ "$COV_INPUT" = "10" ]; then
    COV_INPUT="1,2,3,4,5,6,7,8,9"
elif [[ "$COV_UPPER" == *"ALL"* ]] || [[ ",$COV_INPUT," == *",10,"* ]]; then
    echo "ERROR: ALL/10 must be selected by itself."
    exit 1
fi

# ------------------------------------------------------------
# Inheritance-model menu
# ------------------------------------------------------------
echo
echo "Choose inheritance model(s):"
echo
echo "1. Additive"
echo "2. Dominant"
echo "3. Recessive"
echo "4. ALL"
echo
read -r -p "Enter selection: " MODEL_INPUT

MODEL_INPUT=$(echo "$MODEL_INPUT" | tr -d ' ')
MODEL_UPPER=$(echo "$MODEL_INPUT" | tr '[:lower:]' '[:upper:]')

if [ "$MODEL_UPPER" = "ALL" ] || [ "$MODEL_INPUT" = "4" ]; then
    MODEL_INPUT="1,2,3"
fi

# ------------------------------------------------------------
# Helper: expand one menu number into abstract covariate names
# ------------------------------------------------------------
expand_choice() {
    case "$1" in
        1) echo "UNADJUSTED" ;;
        2) echo "PC1" ;;
        3) echo "PC2" ;;
        4) echo "PC1 PC2" ;;
        5) echo "Sire" ;;
        6) echo "Birth_Year" ;;
        7) echo "Birth_Year_Group" ;;
        8) echo "Technician" ;;
        9) echo "Birth_Year Sire" ;;
        *)
            echo "ERROR: Unknown covariate option '$1'." >&2
            return 1
            ;;
    esac
}

# ------------------------------------------------------------
# Start a fresh manifest describing every GWAS actually run
# ------------------------------------------------------------
printf "PREFIX\tTAG\tTITLE\tMODEL\tTEST_ID\n" > "$MANIFEST"

IFS=',' read -ra COV_SPECS <<< "$COV_INPUT"
IFS=',' read -ra MODEL_SPECS <<< "$MODEL_INPUT"

for spec in "${COV_SPECS[@]}"; do

    IFS='+' read -ra PARTS <<< "$spec"

    abstract_covars=()
    has_unadjusted=0

    for part in "${PARTS[@]}"; do
        expanded=$(expand_choice "$part") || exit 1

        for item in $expanded; do
            if [ "$item" = "UNADJUSTED" ]; then
                has_unadjusted=1
            else
                abstract_covars+=("$item")
            fi
        done
    done

    if [ "$has_unadjusted" -eq 1 ] && [ "${#abstract_covars[@]}" -gt 0 ]; then
        echo "ERROR: Unadjusted (1) cannot be combined with another covariate using +."
        exit 1
    fi

    # --------------------------------------------------------
    # Canonicalize / de-duplicate covariates so 2+3 equals 4
    # and 6+5 equals 9.
    # --------------------------------------------------------
    if [ "$has_unadjusted" -eq 1 ]; then
        tag="unadjusted"
        title="Unadjusted"
        covar_csv=""
    else
        unique_covars=()

        for candidate in PC1 PC2 Birth_Year Birth_Year_Group Technician Sire; do
            for x in "${abstract_covars[@]}"; do
                if [ "$x" = "$candidate" ]; then
                    already=0
                    for y in "${unique_covars[@]:-}"; do
                        [ "$y" = "$candidate" ] && already=1
                    done
                    [ "$already" -eq 0 ] && unique_covars+=("$candidate")
                fi
            done
        done

        if [ "${#unique_covars[@]}" -eq 0 ]; then
            echo "ERROR: No covariates were selected."
            exit 1
        fi

        tag_parts=()
        title_parts=()
        plink_cols=()

        for cov in "${unique_covars[@]}"; do
            case "$cov" in
                PC1)
                    tag_parts+=("pc1")
                    title_parts+=("PC1")
                    plink_cols+=("PC1")
                    ;;
                PC2)
                    tag_parts+=("pc2")
                    title_parts+=("PC2")
                    plink_cols+=("PC2")
                    ;;
                Birth_Year)
                    tag_parts+=("birth_year")
                    title_parts+=("Birth Year")
                    plink_cols+=("Birth_Year")
                    ;;
                Birth_Year_Group)
                    tag_parts+=("birth_year_group")
                    title_parts+=("Birth Year Group")
                    plink_cols+=("Birth_Year_Group")
                    ;;
                Technician)
                    if [ -z "$TECH_COLS" ]; then
                        echo "ERROR: No Technician_* dummy columns were found in $COVAR_FILE."
                        exit 1
                    fi
                    tag_parts+=("technician")
                    title_parts+=("Technician")
                    IFS=',' read -ra temp_cols <<< "$TECH_COLS"
                    plink_cols+=("${temp_cols[@]}")
                    ;;
                Sire)
                    if [ -z "$SIRE_COLS" ]; then
                        echo "ERROR: No Sire_* dummy columns were found in $COVAR_FILE."
                        exit 1
                    fi
                    tag_parts+=("sire")
                    title_parts+=("Sire")
                    IFS=',' read -ra temp_cols <<< "$SIRE_COLS"
                    plink_cols+=("${temp_cols[@]}")
                    ;;
            esac
        done

        tag=$(IFS=_; echo "${tag_parts[*]}")

        title="${title_parts[0]}"
        for ((j = 1; j < ${#title_parts[@]}; j++)); do
            title+=" + ${title_parts[$j]}"
        done

        covar_csv=$(IFS=,; echo "${plink_cols[*]}")
    fi

    for model_choice in "${MODEL_SPECS[@]}"; do

        case "$model_choice" in
            1)
                model="ADD"
                test_id="ADD"
                model_label="Additive"
                model_args=(--logistic hide-covar)
                ;;
            2)
                model="DOM"
                test_id="DOM"
                model_label="Dominant"
                model_args=(--logistic dominant hide-covar)
                ;;
            3)
                model="REC"
                test_id="REC"
                model_label="Recessive"
                model_args=(--logistic recessive hide-covar)
                ;;
            *)
                echo "ERROR: Unknown inheritance-model option '$model_choice'."
                exit 1
                ;;
        esac

        prefix="gwas_${tag}_${model}"

        cmd=(
            plink
            --bfile srd_qc_allchr
            --autosome-num 30
            --allow-no-sex
        )

        if [ -n "$covar_csv" ]; then
            cmd+=(--covar "$COVAR_FILE" --covar-name "$covar_csv")
        fi

        cmd+=("${model_args[@]}" --out "$prefix")

        echo
        echo "============================================================"
        echo "Running: $title | $model_label"
        echo "Output:  $prefix"
        echo "============================================================"

        "${cmd[@]}"

        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$prefix" "$tag" "$title" "$model_label" "$test_id" \
            >> "$MANIFEST"
    done
done

echo
echo "Finished. Runs are recorded in $MANIFEST"

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$MANIFEST"
else
    cat "$MANIFEST"
fi
```

Save the script, then make it executable:

```bash
chmod +x run_gwas.sh
```

### In-class GWAS selection

For the classroom exercise, run:

```bash
./run_gwas.sh
```

At the covariate prompt enter:

```text
1,4
```

Remember that the comma means **two separate GWAS analyses**:

```text
1 = Unadjusted
4 = PC1 + PC2
```

At the inheritance-model prompt enter:

```text
1
```

which means:

```text
Additive only
```

Therefore, the in-class analysis produces only:

```text
gwas_unadjusted_ADD.assoc.logistic
gwas_pc1_pc2_ADD.assoc.logistic
```

This is intentional. The goal during class is to compare the **unadjusted additive GWAS** against the **PC1 + PC2 adjusted additive GWAS** without generating dozens of additional files.

### What is `gwas_run_manifest.tsv`?

Every time the script runs, it records the analyses that were actually performed in:

```text
gwas_run_manifest.tsv
```

For the classroom run, it should look approximately like:

```text
PREFIX                 TAG         TITLE        MODEL      TEST_ID
gwas_unadjusted_ADD     unadjusted  Unadjusted   Additive   ADD
gwas_pc1_pc2_ADD        pc1_pc2     PC1 + PC2    Additive   ADD
```

This manifest is important because the plotting script will read it automatically. You therefore do **not** need to edit the plotting script every time you choose a different covariate or covariate combination.

### Quickly inspect the two in-class GWAS results

```bash
head gwas_unadjusted_ADD.assoc.logistic
```

and:

```bash
head gwas_pc1_pc2_ADD.assoc.logistic
```

If you specifically want to view only the additive SNP-test rows:

```bash
awk 'NR==1 || $5=="ADD"' gwas_unadjusted_ADD.assoc.logistic | head
```

The association output contains marker-level results. Important columns include:

- `CHR` = chromosome
- `SNP` = marker name
- `BP` = base-pair position
- `A1` = coded allele
- `TEST` = test being reported
- `NMISS` = number of animals used for that marker
- `OR` = odds ratio
- `STAT` = test statistic
- `P` = p-value

---

### 🏠 Homework / later analysis

After you understand the in-class comparison, rerun:

```bash
./run_gwas.sh
```

For the covariate selection, enter:

```text
ALL
```

For the inheritance model, enter:

```text
ALL
```

This runs **each predefined covariate configuration 1-9 separately** under:

```text
Additive
Dominant
Recessive
```

The predefined configurations are:

```text
1. Unadjusted
2. PC1
3. PC2
4. PC1 + PC2
5. Sire
6. Birth Year
7. Birth Year Group
8. Technician
9. Birth Year + Sire
```

Thus, `ALL` includes both of the combined presets you will want to compare later:

```text
PC1 + PC2
Birth Year + Sire
```

> 💡 `ALL` does **not** generate every mathematically possible covariate combination. It runs the nine predefined configurations above as separate analyses. If you want a different custom combination, specify it explicitly with `+`, such as `2+5` for PC1 + Sire.


## 7. Generate Standalone Manhattan and Q-Q Plots

This section generates:

- FDR Manhattan plots with a red dashed line at `FDR < 0.05`;
- nominal Manhattan plots with a red dashed line at `p < 1e-5`;
- Q-Q plots with genomic inflation (`lambda`) displayed on each plot.

The plotting script is intentionally **not hard-coded to a fixed list of covariates**. Instead, it reads:

```text
gwas_run_manifest.tsv
```

and plots whatever GWAS analyses were actually generated by `run_gwas.sh`.

This means the same R script works for:

```text
Unadjusted
PC1
PC2
PC1 + PC2
Birth Year + Sire
PC1 + Sire
or another custom selection
```

without needing to rewrite the R code.

Save the plotting code as an R script.

Geany option:

```bash
geany gwas_plots.R
```

Command-line option:

```bash
nano gwas_plots.R
```

Run it with:

```bash
Rscript gwas_plots.R
```

Paste:

```r
library(tidyverse)
library(scales)

plot_single_manhattan <- function(df_model, model_name, run_title, out_png, mode = c("FDR", "Nominal")) {
  mode <- match.arg(mode)

  chr_info <- df_model %>%
    distinct(CHR, BP) %>%
    group_by(CHR) %>%
    summarise(chr_len = max(as.numeric(BP)), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(tot = lag(cumsum(as.numeric(chr_len)), default = 0))

  df_plot <- df_model %>%
    left_join(chr_info %>% select(CHR, tot), by = "CHR") %>%
    mutate(BP_cum = as.numeric(BP) + tot)

  axis_df <- df_plot %>%
    group_by(CHR_LABEL, CHR) %>%
    summarize(center = (max(BP_cum) + min(BP_cum)) / 2, .groups = "drop") %>%
    arrange(CHR)

  shade_rects <- chr_info %>%
    filter(CHR %% 2 == 0) %>%
    mutate(xmin = tot, xmax = tot + chr_len)

  N_samples  <- if ("NMISS" %in% colnames(df_model)) max(df_model$NMISS, na.rm = TRUE) else NA
  n_snps     <- nrow(df_model)
  chisq      <- qchisq(1 - df_model$P, df = 1)
  lambda_val <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1), 3)

  if (mode == "FDR") {
    df_plot$y_val <- df_plot$log10_FDR
    cutoff_y <- -log10(0.05)
    sig_count <- sum(df_plot$FDR < 0.05, na.rm = TRUE)
    subtitle_text <- "Red dashed line: FDR < 0.05"
    y_axis_label <- expression(-log[10](FDR))
    label_text <- paste0(
      "N = ", comma(N_samples),
      " | SNPs = ", comma(n_snps),
      " | lambda = ", lambda_val,
      " | FDR < 0.05: ", sig_count
    )
  } else {
    df_plot$y_val <- df_plot$log10_P
    cutoff_y <- -log10(1e-5)
    sig_count <- sum(df_plot$P < 1e-5, na.rm = TRUE)
    subtitle_text <- "Red dashed line: Nominal p < 1e-5"
    y_axis_label <- expression(-log[10](italic(p)))
    label_text <- paste0(
      "N = ", comma(N_samples),
      " | SNPs = ", comma(n_snps),
      " | lambda = ", lambda_val,
      " | p < 1e-5: ", sig_count
    )
  }

  p <- ggplot() +
    geom_rect(
      data = shade_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "grey93", alpha = 0.6, inherit.aes = FALSE
    ) +
    geom_point(
      data = df_plot,
      aes(x = BP_cum, y = y_val, color = as.factor(CHR %% 10)),
      alpha = 0.75, size = 1.3
    ) +
    geom_hline(yintercept = cutoff_y, color = "red", linetype = "dashed", linewidth = 0.7) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = label_text,
      hjust = -0.02, vjust = 1.6,
      size = 4, fontface = "bold"
    ) +
    scale_x_continuous(
      labels = axis_df$CHR_LABEL,
      breaks = axis_df$center,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    scale_color_brewer(palette = "Paired", guide = "none") +
    labs(
      title = paste0(run_title, " Manhattan Plot: ", model_name, " Model (", mode, ")"),
      subtitle = subtitle_text,
      x = "Chromosome",
      y = y_axis_label
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(b = 6)),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 8.5, vjust = 0.5)
    )

  ggsave(out_png, plot = p, width = 11, height = 5.5, dpi = 300)
}

plot_single_qq <- function(pvals, model_name, run_title, out_png) {
  pvals <- na.omit(pvals)
  pvals <- pvals[pvals > 0 & pvals <= 1]
  n <- length(pvals)

  obs <- -log10(sort(pvals, decreasing = FALSE))
  exp <- -log10(ppoints(n))

  chisq <- qchisq(1 - pvals, df = 1)
  lambda_val <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1), 3)

  qq_data <- data.frame(Observed = obs, Expected = exp)
  max_val <- ceiling(max(c(obs, exp)))

  p <- ggplot(qq_data, aes(x = Expected, y = Observed)) +
    geom_point(color = "steelblue", alpha = 0.6, size = 2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 0.8) +
    annotate(
      "text",
      x = max_val * 0.78,
      y = max_val * 0.95,
      label = paste0("lambda == ", lambda_val),
      parse = TRUE,
      size = 5.5,
      fontface = "bold"
    ) +
    coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
    labs(
      title = paste0(run_title, " Q-Q Plot: ", model_name, " Model"),
      x = expression(Expected ~ -log[10](italic(p))),
      y = expression(Observed ~ -log[10](italic(p)))
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  ggsave(out_png, plot = p, width = 6, height = 6, dpi = 300)
}

# ------------------------------------------------------------
# Read the GWAS runs that were actually performed
# ------------------------------------------------------------
manifest_file <- "gwas_run_manifest.tsv"

if (!file.exists(manifest_file)) {
  stop("gwas_run_manifest.tsv was not found. Run ./run_gwas.sh first.")
}

runs <- read.delim(
  manifest_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(runs) == 0) {
  stop("The GWAS manifest contains no completed runs.")
}

cat("\nGWAS runs found in manifest:\n")
print(runs)
cat("\n")

# ------------------------------------------------------------
# Plot every completed GWAS listed in the manifest
# ------------------------------------------------------------
for (i in seq_len(nrow(runs))) {

  r <- runs[i, ]
  file_path <- paste0(r$PREFIX, ".assoc.logistic")

  if (!file.exists(file_path)) {
    warning(paste("File missing:", file_path))
    next
  }

  cat("Plotting:", r$TITLE, "-", r$MODEL, "\n")

  df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE) %>%
    filter(TEST == r$TEST_ID) %>%
    filter(!is.na(P) & P > 0 & P <= 1) %>%
    mutate(
      CHR_LABEL = ifelse(CHR == 30, "X", as.character(CHR)),
      CHR = as.numeric(CHR),
      log10_P = -log10(P),
      FDR = p.adjust(P, method = "BH"),
      log10_FDR = -log10(FDR)
    ) %>%
    filter(CHR >= 1 & CHR <= 30)

  model_short <- tolower(substr(r$MODEL, 1, 3))

  plot_single_manhattan(
    df_model = df,
    model_name = r$MODEL,
    run_title = r$TITLE,
    out_png = paste0("manhattan_fdr_", r$TAG, "_", model_short, ".png"),
    mode = "FDR"
  )

  plot_single_manhattan(
    df_model = df,
    model_name = r$MODEL,
    run_title = r$TITLE,
    out_png = paste0("manhattan_nominal_", r$TAG, "_", model_short, ".png"),
    mode = "Nominal"
  )

  plot_single_qq(
    pvals = df$P,
    model_name = r$MODEL,
    run_title = r$TITLE,
    out_png = paste0("qq_", r$TAG, "_", model_short, ".png")
  )
}
```

### In-class plot outputs

Because the classroom GWAS contains only the two additive runs, the plotting script will initially create:

```text
manhattan_fdr_unadjusted_add.png
manhattan_nominal_unadjusted_add.png
qq_unadjusted_add.png

manhattan_fdr_pc1_pc2_add.png
manhattan_nominal_pc1_pc2_add.png
qq_pc1_pc2_add.png
```

View the Q-Q plots side by side conceptually and ask:

1. How does the unadjusted Q-Q plot differ from the PC1 + PC2 adjusted Q-Q plot?
2. What happens to genomic inflation (`lambda`)?
3. Do the strongest Manhattan-plot signals remain similar after PC adjustment?
4. Does adjustment appear to reduce broad inflation, or does it appear overly conservative?

### Later plots: PC1 alone, PC2 alone, and other combinations

No new R code is required.

For example, if you rerun `run_gwas.sh` and select:

```text
2,3,4
```

with:

```text
1
```

for the additive inheritance model, the manifest will contain PC1, PC2, and PC1 + PC2 as separate runs. Running:

```bash
Rscript gwas_plots.R
```

will then automatically generate Manhattan and Q-Q plots for all three.

Similarly, a custom selection such as:

```text
6+5
```

will generate a Birth Year + Sire GWAS and the plotting script will automatically produce the corresponding figures.

> 💡 The plotting script follows the manifest. If you change the GWAS selection, rerun `Rscript gwas_plots.R` after the GWAS has completed.

---

## 🧠 What have we done in this tutorial?

You have now worked through the major steps of a complete GWAS workflow:

```text
Raw PED/MAP genotype data
        ↓
Initial SNP and animal QC
        ↓
Hardy-Weinberg equilibrium assessment
        ↓
Data-informed HWE filtering
        ↓
Principal component analysis
        ↓
Metadata and potential-covariate exploration
        ↓
Covariate-file preparation
        ↓
Whole-genome GWAS dataset
        ↓
Additive GWAS: Unadjusted vs PC1 + PC2
        ↓
Manhattan and Q-Q plots
        ↓
Comparison of association signals and genomic inflation
```

More specifically, you learned how to:

- inspect genotype and metadata files before analysis;
- apply call-rate and minor-allele-frequency QC filters;
- calculate and visualize HWE statistics;
- use the observed HWE distribution to identify an extreme-deviation tail;
- calculate principal components to describe population structure;
- visualize metadata variables on the PCA;
- test potential covariates against the phenotype;
- prepare numerical and categorical covariates for GWAS;
- run logistic GWAS under an additive inheritance model;
- compare an unadjusted model against a PC-adjusted model;
- calculate FDR-adjusted p-values and genomic inflation (`lambda`);
- generate Manhattan and Q-Q plots;
- use a flexible GWAS script to construct additional covariate models without rewriting the analysis code.

The **in-class analysis deliberately stops at two additive models** so that the focus remains on understanding what covariate adjustment does. The homework extends the same workflow across the predefined covariate configurations and additive, dominant, and recessive inheritance models.

---


## Quick command summary

| Step | Main command/script | Main output |
|---|---|---|
| Initial QC | `plink --file ... --maf --geno --mind --make-bed` | `srd_qc.*` |
| HWE calculation | `plink --bfile srd_qc --hardy` | `plink_results_hwinb.hwe` |
| HWE plots | `Rscript hwe_plots.R` | `hwe_distribution_plots.png` |
| HWE filtering | `plink --bfile srd_qc --hwe 3e-9 --make-bed` | `srd_qc_hwe.*` |
| PCA | `plink --bfile srd_qc_hwe --pca` | `srd_pca.eigenval`, `srd_pca.eigenvec` |
| PCA/covariates | `Rscript pca_covariates.R` | PCA plots, `covariates.txt` |
| Whole-genome dataset | `plink ... --autosome-num 30 ...` | `srd_qc_allchr.*` |
| Flexible GWAS | `./run_gwas.sh` | `.assoc.logistic` files + `gwas_run_manifest.tsv` |
| In-class GWAS | Covariates `1,4`; model `1` | Unadjusted ADD + PC1/PC2 ADD |
| Homework GWAS | Covariates `ALL`; model `ALL` | Presets 1-9 × ADD/DOM/REC |
| Manhattan/Q-Q plots | `Rscript gwas_plots.R` | `.png` plots for runs in the manifest |

## Notes for students ✍️📖

- Run commands from your own `~/workshop` directory so your output files stay separate from other students' work.
- Read the PLINK `.log` file after every major PLINK command. It records how many animals and SNPs were loaded, removed, and retained.
- Do not delete intermediate files until the workflow is complete; later steps depend on several of them.
- Remember: a comma in the GWAS menu means **separate analyses**, while `+` means **covariates included together in one model**.
- `ALL` runs the nine predefined covariate configurations separately; it does not generate every possible combination.
- The classroom GWAS uses **additive only** for **Unadjusted** and **PC1 + PC2**.
- If you use **Geany**, it is a nice lightweight editor and easier for most people to navigate visually.
- If you use **nano**, it is usually faster if you are already comfortable in the terminal.
- If an R script stops with an error, read the **first** error message before rerunning the script. Later errors may simply be consequences of the first one.
- When comparing Q-Q plots, do not judge a model only by whether points are above or below the diagonal. Consider the overall pattern, genomic inflation (`lambda`), sample size, and whether covariate adjustment is biologically/statistically justified.

## The End! :grin: :clap:
<p align="center">
  <img src="images/celebrate.gif" width="1000" alt="Tutorial Overview" />
  <br>
  <sub><i>Source: <a href="https://www.pinterest.com/pin/55239532918424769/">Pinterest</a></i></sub>
</p>
