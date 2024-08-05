#dietary_indices_microbiota
##This repository includes the script to carry out Visualization and statistical techniques utilized included principal component analysis biplot, non-metric multidimensional scaling, causal mediation analysis, and correlation analysis.

#Non-metric Multidimensional Scaling (NMDS)

```{r}
library(phyloseq)
library(vegan)
library(phyloseq)
library(ggplot2)
library(compositions)
library(grid)
```

```{r}
physeq
physeq_MorI_infant <- subset_samples(physeq, MorI == "Infant")
```

#transform counts to abundances
```{r}
# Assuming your phyloseq object is named physeq
# Transform counts to relative abundances
physeq_MorI_infant <- transform_sample_counts(physeq_MorI_infant, function(x) x / sum(x))
physeq_ab <- physeq_MorI_infant
```

#if I do this I will eliminate all of the varaibility between speicies and families and keep only hte genus which is not the best approach.
```{r}
physeq_ab <- tax_glom(physeq_ab, taxrank = "Genus", NArm = TRUE)
```
#CLR normalitacion before NMDS (this time I am not doing otherswise I cannot do bray curtis distance)
```{r}
# Extract OTU table from physeq object
otu_table <- otu_table(physeq_ab)

# Ensure OTU table is in the right orientation (samples as rows, taxa as columns)
if(taxa_are_rows(physeq_ab)) {
  otu_table <- t(otu_table)
}

# Convert to a data frame
otu_df <- as.data.frame(otu_table)

# CLR transformation
clr_transformed <- clr(otu_df + 1)  # Adding 1 to avoid log of zero issues

# Convert back to matrix and then to phyloseq otu_table
clr_matrix <- as.matrix(clr_transformed)
otu_table(physeq_ab) <- otu_table(clr_matrix, taxa_are_rows = FALSE)

# Check the transformed OTU table
head(otu_table(physeq_ab))
```
#incorporating bfidobacterium categorized
```{r}
summary(infant_core_divers$g__Bifidobacterium)

# Calculate the median abundance of Bifidobacterium
global_median <- median(infant_core_divers$g__Bifidobacterium)

# Create a new categorical variable based on the median abundance
abundance_category <- ifelse(infant_core_divers$g__Bifidobacterium > global_median, "Above Median", "Below Median")

table(abundance_category)
# Add the new variable to the sample data in the physeq object
sample_data(physeq_ab)$Bifido_Abundance_Category <- abundance_category
table(abundance_category)
```


#preprocesing for NMDS

#if I am going to do bray curtis I will have to Remove empty rows and check dimensions of the data firsst
```{r}
# Extract OTU table from phyloseq object
otu_table_ab <- otu_table(physeq_ab)

# Convert to matrix if necessary
otu_matrix <- as(otu_table_ab, "matrix")

# Remove rows with all zeros (samples with no counts)
otu_matrix_clean <- otu_matrix[rowSums(otu_matrix) != 0, ]

# Remove columns with all zeros (OTUs with no counts)
otu_matrix_clean <- otu_matrix_clean[, colSums(otu_matrix_clean) != 0]

# Check the dimensions
if (nrow(otu_matrix_clean) == 0 || ncol(otu_matrix_clean) == 0) {
  stop("Your dataset has no non-zero rows or columns. Please provide a valid dataset.")
}

# Create a new OTU table
otu_table_clean <- otu_table(otu_matrix_clean, taxa_are_rows = FALSE)

# Create a new phyloseq object with the cleaned OTU table
physeq_ab_clean <- phyloseq(otu_table_clean, sample_data(physeq_ab), tax_table(physeq_ab))

# Check the cleaned phyloseq object
physeq_ab_clean


```

```{r}
# Convert phyloseq object to a community matrix
com <- otu_table(physeq_ab_clean)
m_com <- as.matrix(com)
dim(m_com)
```


```{r}
# Extract sample data
sample_data_df <- data.frame(sample_data(physeq_ab_clean))


# Ensure environmental variables are correctly formatted in sample data
env <- sample_data_df[, c("Total_DII_sum", "Total_MMDS_sum", "Total_DQI_sum", "Total_HEI_sum")]
dim(env)
env1 <- sample_data_df[, c("Total_DII_sum", "Total_MMDS_sum", "Total_DQI_sum", "Total_HEI_sum", "Parto", "Bifido_Abundance_Category")]
dim(env)
```



```{r}
# Perform NMDS
set.seed(123)
# Run NMDS
nmds <- metaMDS(m_com, distance = "euclidean", k = 3, trymax = 50)
print(nmds)
```

#The Shepard diagram helps you see how well the NMDS distances match the original dissimilarities. You can use the stressplot function from the vegan package.
```{r}
# Shepard Diagram
stressplot(nmds)
```

#The stress plot is useful for determining the appropriate number of dimensions for the NMDS. It plots the stress value against the number of dimensions.
```{r}
?metaMDS
```
# Creating a stress plot for different dimensions
```{r}
stress_values <- sapply(1:10, function(k) {
  nmds <- metaMDS(m_com, distance = "euclidean", k = k, trymax = 50)
  return(nmds$stress)
}) 
```

```{r}
# Plotting stress values
plot(1:10, stress_values, type = "b", pch = 19, 
     xlab = "Number of Dimensions", ylab = "Stress",
     main = "Stress Plot")
```


```{r}
# Run envfit
en <- envfit(nmds, env, permutations = 999, na.rm = TRUE)
print(en)

```

```{r}
# Extract NMDS scores for plotting
nmds_scores <- scores(nmds, display = "sites")
sample_data_df$NMDS1 <- nmds_scores[, 1]
sample_data_df$NMDS2 <- nmds_scores[, 2]

# Prepare coordinates for continuous and categorical environmental variables
arrow_mul <- ordiArrowMul(en)
en_coord_cont <- as.data.frame(scores(en, "vectors")) * arrow_mul
en_coord_cat <- as.data.frame(scores(en, "factors")) * arrow_mul
```

```{r}
gg <- ggplot(data = sample_data_df, aes(x = NMDS1, y = NMDS2)) +  
  geom_point(aes(colour = as.factor(Parto), shape = as.factor(Bifido_Abundance_Category)), 
             size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("orange", "steelblue")) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size = 1, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
            colour = "grey30", fontface = "bold", label = rownames(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Delivery Mode", shape = "Bifido Abundance Category")

# Print the plot
print(gg)

```

```{r}
print(gg)
```



```{r}
ggsave("NMDSeuclidean_vector.png")
```
# Significance of Environmental Factors (Vectors or Factors)
#A plot showing the significance (p-values) of the environmental vectors or factors can help in understanding which variables are significantly influencing the ordination.

#. Bi-plot with Environmental Variables
#This is similar to the ordination plot but includes both species and environmental variables. This can provide a more comprehensive view of the relationships.

#plot
```{r}
# Extract p-values and R values for vectors and factors
pvals_vectors <- en$vectors$pvals

#optional (r2)^2 to explain the correlation coefficient (option1)
r_vectors_r2_s <- sqrt(en$vectors$r)

#optional R2 to explain the explainable variance (option2)
r_vectors <- en$vectors$r

# Combine them into a single data frame
pvals_df <- data.frame(
  Variable = c(names(pvals_vectors)),
  P_values = c(pvals_vectors),
  R_Squared_squared = c(r_vectors_r2_s),
  R_Squared = c(r_vectors, r_factors)
)

print(pvals_df)
str(pvals_df)
```

#r values and p values combined 


```{r}
pvals_df <- data.frame(
  Variable = c(names(pvals_vectors)),
  P_values = c(pvals_vectors),
  R_Squared_squared = c(r_vectors_r2_s),
  R_Squared = c(r_vectors)
)

```
```{r}
# Create the plot
plot_combined <- ggplot(pvals_df, aes(x = reorder(Variable, R_Squared_squared), y = R_Squared_squared)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = sprintf("P: %.3f", P_values)), nudge_y = 0.01, hjust = -0.1, size = 3, color = "black") +
  geom_text(aes(label = sprintf("R²: %.3f", R_Squared)), nudge_y = -0.01, hjust = 1.1, size = 3, color = "black") +
  coord_flip(clip = "off") +
  theme_minimal() +
  labs(title = "R² and P-values of Environmental Variables",
       x = "Environmental Variable",
       y = "Suqare Root of R² (correlation coefficient)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.margin = margin(t = 10, r = 50, b = 8, l = 8))

# Display the plot
print(plot_combined)
```

```{r}
# Add a new column for custom labels
pvals_df$Custom_Variable <- c("Dietary Inflammatory Index", 
                              "Mediterranean Diet Score", 
                              "Diet Quality Index", 
                              "Healthy Eating Index")

# Create the plot using the custom labels
plot_combined <- ggplot(pvals_df, aes(x = reorder(Custom_Variable, R_Squared_squared), y = R_Squared_squared)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_text(aes(label = sprintf("P: %.3f", P_values)), nudge_y = 0.01, hjust = -0.001, size = 3, color = "black") +
  geom_text(aes(label = sprintf("R²: %.3f", R_Squared)), nudge_y = -0.02, hjust = 0.8, size = 3, color = "black") +
  coord_flip(clip = "off") +
  theme_minimal() +
  labs(title = "R² and Significance of Environmental Variables",
       x = "Environmental Variable",
       y = "Square Root of R² (correlation coefficient)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.margin = margin(t = 10, r = 60, b = 10, l = 10))

# Display the plot
print(plot_combined)
```

```{r}
ggsave("R² and Significance of Environmental Variables_vectors.png")
```
