library(readxl)
library(dplyr)
library(ggplot2)
library(lazyeval)

# filter the passed data frame(df) by the passed column name and passed value,
# which is expected to be vector and may be NULL; if is_numeric is TRUE, converts
# the value to numeric type before filtering
filterByIn <- function(df, col_name, value, is_numeric = FALSE) {
  if (!is.null(value)) {
    #value <- as.character(value)
    #value <- strsplit(value, ',')[[1]]
    if (is_numeric)
      value <- as.numeric(value)
    in_criteria <- interp(~y %in% x, .values=list(y = as.name(col_name), x = value))
    df <- filter_(df, in_criteria)
  }
  return(df)
}

# filter the passed data frame(df) by the passed column name and passed value,
# which is expected to be a number and may be NA
filterByGT <- function(df, col_name, value) {
  if (!is.null(value)) {
    gt_criteria <- interp(~y > x, .values=list(y = as.name(col_name), x = value))
    df <- filter_(df, gt_criteria)
  }
  return(df)
}

# calculates the weighted average of %Identity and total coverage from the passed vectors;
# the vectors contain %Coverage and %Identity respectively for a gene, at a given node and
# for a given isolate; both vectors are expected to be of the same length; the weighted average
# is the sum of pairwise products of the vectors; the total coverage is just the sum of coverages
# in the first vector; calculations are done only if the vectors have multiple elements;
# otherwise just returns the %Coverage and %Identity from the respective vector
calcCovIden <- function(covs, idens) {
  cov <- covs[1]
  iden <- idens[1]*100
  len <- length(idens)
  if (len > 1) {
    cov <- sum(covs)
    iden <- sum(covs*idens)
  }
  c(cov, iden)
}

updateProg <- function(updateProgress = NULL, text) {
  # If we were passed a progress update function, call it
  if (is.function(updateProgress))
    updateProgress(detail = text)
}

# link various sheets in the excel file to get the final list of fields
# (including SEQUENCE, which will be removed later); some fields are renamed 
# to make the linking work
getInitialData <- function(excel_file, updateProgress = NULL) {
  updateProg(updateProgress, "Extracting fields from pipeline run")
  vfdb <- read_excel(excel_file,sheet='VFDB') %>%
                                 select(FILE, GENE, SEQUENCE, `%COVERAGE`, `%IDENTITY`)
  updateProg(updateProgress, "Renaming fields and linking sheets")
  linking <- read_excel(excel_file,sheet='Linking') %>% rename(FILE = File)
  
  join_by_file <- inner_join(vfdb, linking, by = "FILE")
  
  meta <- read_excel(excel_file,sheet='Meta')
  
  join_by_id <- inner_join(join_by_file, meta, by = "ID#") %>% rename(ID = 'ID#')
  
  updateProg(updateProgress, "Extracting final fields")
  fields_with_seq <- select(join_by_id, ID, GENE, SEQUENCE, Species, Serogroup, `Sequence Type`, `Collected Date`, `%COVERAGE`,
                            `%IDENTITY`, Facility, Room) %>%
                                  arrange(ID, GENE, SEQUENCE)
  return(fields_with_seq)
}

getFilteredData <- function(fields_with_seq, species = NULL, serogroup = NULL, seq_type = NULL, genes = NULL, 
                            coverage = 0.0, identity = 0.0, updateProgress = NULL) {
  # filter by species, serogroup, sequence type, genes 
  # if one or more of them have been specified
  updateProg(updateProgress, "Filtering by species")
  fields_with_seq <- filterByIn(fields_with_seq, "Species", species)
  updateProg(updateProgress, "Filtering by serogroups")
  fields_with_seq <- filterByIn(fields_with_seq, "Serogroup", serogroup, TRUE)
  updateProg(updateProgress, "Filtering by sequence types")
  fields_with_seq <- filterByIn(fields_with_seq, "Sequence Type", seq_type, TRUE)
  updateProg(updateProgress, "Filtering by genes")
  fields_with_seq <- filterByIn(fields_with_seq, "GENE", genes)
  
  updateProg(updateProgress, "Generating gene fields")
  aggr_by_gene <- tibble()
  if (nrow(fields_with_seq) > 0) {
    prev_id <- ""
    prev_seq <- ""
    prev_gene <- ""
    max_cov <- -1.0
    tot_cov <- 0.0
    max_tot_cov <- -1.0
    gene_max_tot_cov <- -1.0
    max_iden <- -1.0
    wgt_iden <- 0.0
    max_wgt_iden <- -1.0
    gene_max_wgt_iden <- -1.0
    covs <- c()
    idens <- c()
    j <- 0
    new_row <- FALSE
    new_gene <- FALSE
    prev_gene_first_index <- -1
    gene_first_index <- 1
    last_index = nrow(fields_with_seq)
  
    # each row in the output shows specimen id, gene name, sequence (node) name,
    # the maximum raw %COVERAGE and %IDENTITY for a node, the maximum calculated 
    # %COVERAGE and %IDENTITY for a node (sequence), as explained in calcCovIden(),
    # and the maximum %COVERAGE and %IDENTITY for a gene across all nodes; for a specimen,
    # a gene has as many rows as the numbers of nodes that have that gene; each node
    # only has one row per gene
    for (i in 1:last_index) {
      id <- fields_with_seq[[i,'ID']]
      gene <- fields_with_seq[[i,'GENE']]
      seq <- fields_with_seq[[i,'SEQUENCE']]
      cov <- fields_with_seq[[i,'%COVERAGE']]
      iden <- fields_with_seq[[i,'%IDENTITY']]
      if (prev_id == id) {
        if (prev_gene == gene) {
          if (prev_seq == seq) {
            covs <- c(covs, cov)
            idens <- c(idens, iden/100)
          } else {
            new_row <- TRUE
          }
          new_gene <- FALSE
        } else {
          new_row <- TRUE
          new_gene <- TRUE
        }
      } else {
        new_row <- TRUE
        new_gene <- TRUE
      }
      
      # if it's a new row then it's a new node; so calculate the %coverage and %identity 
      # for the previous node, so that its (previous) row can be updated
      if (new_row || i == last_index) {
        if (length(covs) > 0) {
          calc <- calcCovIden(covs, idens)
          tot_cov <- calc[1]
          wgt_iden <- calc[2]
          # may need only one of tot_cov and max_tot_cov and only one of wgt_iden and max_wgt_iden
          if (tot_cov > max_tot_cov)
            max_tot_cov <- tot_cov
          if (wgt_iden > max_wgt_iden)
            max_wgt_iden <- wgt_iden
        }
        # get the maximum %coverage and %identity for a gene across all nodes
        if (gene_max_tot_cov < max_tot_cov)
          gene_max_tot_cov <- max_tot_cov
        if (gene_max_wgt_iden < max_wgt_iden)
          gene_max_wgt_iden <- max_wgt_iden
        if (length(covs) > 0)
          covs <- c()
        if (length(idens) > 0)
          idens <- c()
        covs <- c(covs, cov)
        idens <- c(idens, iden/100)
        if (new_row) {
          j <- j + 1
          max_cov <- -1.0
          max_iden <- -1.0
        }
      }
      
      if (cov > max_cov)
        max_cov <- cov
      if (iden > max_iden)
        max_iden <- iden
      # update/add a row with the data obtained above 
      aggr_by_gene[j,1] <- id
      aggr_by_gene[j,2] <- gene
      aggr_by_gene[j,3] <- seq
      aggr_by_gene[j,4] <- max_cov
      aggr_by_gene[j,5] <- max_iden
      if (new_row || i == last_index) {
        if (j > 1) {
          k <- ifelse(i == last_index, j, j - 1)
          aggr_by_gene[k,6] <- max_tot_cov
          aggr_by_gene[k,7] <- max_wgt_iden
          # if it's a new gene, update the previous gene's maximum %coverage and %identity
          # in all the rows that are showing it; also set this row's index as the first
          # index for the current gene
          if (new_gene || i == last_index) {
            aggr_by_gene[gene_first_index : k,8] <- gene_max_tot_cov
            aggr_by_gene[gene_first_index : k,9] <- gene_max_wgt_iden
            gene_first_index <- j
            gene_max_tot_cov <- -1.0
            gene_max_wgt_iden <- -1.0
          }
        }
        max_tot_cov <- 1.0
        max_wgt_iden <- -1.0
        new_row <- FALSE
      }
      prev_id <-id
      prev_gene <- gene
      prev_seq <- seq
    }
    colnames(aggr_by_gene) <- c('ID', 'GENE', 'SEQUENCE', 
                                'SEQ MAX %COVERAGE', 'SEQ MAX %IDENTITY', 
                                'SEQ TOTAL %COVERAGE',  'SEQ WEIGHTED %IDENTITY',
                                'GENE MAX TOTAL %COVERAGE', 'GENE MAX WEIGHTED %IDENTITY')
  }
  #filter by %coverage and/or %identity
  updateProg(updateProgress, "Filtering by %coverage")
  aggr_by_gene <- filterByGT(aggr_by_gene, "SEQ TOTAL %COVERAGE", coverage)
  updateProg(updateProgress, "Filtering by %weighted identity")
  aggr_by_gene <- filterByGT(aggr_by_gene, "SEQ WEIGHTED %IDENTITY", identity)
  
  fields_with_seq <- select(fields_with_seq, -SEQUENCE)
  list(all_fields = fields_with_seq, gene_fields = aggr_by_gene)
}