
# function for cleaning corncob DA differentialTest output
# retrieves only the mu (abundance) scores and intercepts for each significant model
cleanDA <- function(da, varname){
  df <- tibble()
  n <- da$significant_models[[1]]$np.mu
  for(i in 1:length(da$significant_models)){
    df <- rbind(df, da$significant_models[[i]]$coefficients[1:n,] %>% as.data.frame() %>% rownames_to_column(var="variable"))
  }
  df$variable <- df$variable %>% gsub(x=., pattern=str_c("mu.", varname), replace="") %>% gsub(x=., pattern="mu.\\(Intercept)", "Intercept") %>% factor()
  tax <- da$data %>% tax_table() %>% as("matrix") %>% data.frame() %>% rownames_to_column(var="ASV") %>%
    filter(ASV %in% da$significant_taxa) %>%  slice(rep(1:n(), each = n))
  df <- bind_cols(df, tax) %>% as_tibble() %>% rename(StdE = "Std. Error", t= "t value", p = "Pr(>|t|)")
  return(df)
}

# small shortcut function that I add whenever I subset samples. Removes taxa with zero abundance
prune_zeros <- function(ps){
  prune_taxa(taxa_sums(ps)>0, ps)
}

# function makes a null matrix from an otu matrix by reshuffling observations for each ASV.
# first is the otu table and then second parameter is the number of matrices to make
make_null <- function(otu_table,x){
  randomize <- function(otu_table){
    nm <- sapply(1:nrow(otu_table), 
                 function(x){sample(as.numeric(otu_table[x,]))}) %>% t()
    rownames(nm) <- rownames(otu_table)
    colnames(nm) <- colnames(otu_table)
    return(nm)
  }
  replicate(x, randomize(otu_table))
}

# borrowed these functions from https://rdrr.io/github/mworkentine/mattsUtils/src/R/microbiome_helpers.R
#' Add taxonomy label
#'
#' add a column to the taxonomy table of a phyloseq object that lists the
#' lowest rank taxonomy assigned to that OTU along with a prefix indicating
#' the taxonomic rank.
#'
#' Example: g:Pseudomonas
#'
#' @param physeq a valid phyloseq object that contains a taxonomy table
#' @param num_species the number of species to retain if more than one are identified
#' @return a phyloseq object with an additional column on the taxonomy
#'         table called "Taxonomy"
#' @export

split_species = function(string, n = 2) {
  splits = str_split(string, "/", n + 1)
  res = map_if(splits, ~length(.x) > 2, ~.x[1:n]) %>%
    map_chr(str_c, collapse = "/")
  return(res)
}

add_taxonomy_column = function(physeq, num_species = 2) {
  tax_df = as.data.frame(tax_table(physeq)@.Data) %>%
    rownames_to_column("OTU") %>%
    mutate(Species = split_species(Species, n = num_species)) %>%
    mutate(Taxonomy =
             case_when(
               is.na(Class)  ~ str_c("p:", Phylum),
               is.na(Order)  ~ str_c("c:", Class),
               is.na(Family)  ~ str_c("o:", Order),
               is.na(Genus)   ~ str_c("f:", Family),
               is.na(Species) ~ str_c("g:", Genus),
               TRUE ~ str_c(Genus, " ", Species)
             )
    )
  
  tax = as.matrix(tax_df[, -1])
  rownames(tax) = tax_df$OTU
  tax_table(physeq) = tax_table(tax)
  
  return(physeq)
}