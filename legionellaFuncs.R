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

# since the rows shown on Gene Fields tab have sequence (node) in them, there can be many repeated rows for the
# same combination of gene, species, serogroup and sequence type; the tibble returned here makes sure that
# there's only one row per unique combination of those fields; it's used for the Virulence Hits tab
getVirulenceRows <- function(gene_rows) {
  if (nrow(gene_rows) > 0) {
    vir_rows <- tibble()
    prev_id <- ""
    prev_gene <- ""
    prev_spec <- ""
    prev_sero <- ""
    prev_seqt <- ""
    j <- 0
    for (i in 1:nrow(gene_rows)) {
      id <- gene_rows[[i,'ID']]
      gene <- gene_rows[[i,'GENE']]
      spec <- gene_rows[[i,'Species']]
      sero <- gene_rows[[i,'Serogroup']]
      seqt <- gene_rows[[i,'Sequence Type']]
      if (prev_id != id || prev_gene != gene || prev_spec != spec || prev_sero != sero || prev_seqt != seqt) {
        j <- j + 1
        vir_rows[j,1] <- id
        vir_rows[j,2] <- gene
        vir_rows[j,3] <- spec
        vir_rows[j,4] <- sero
        vir_rows[j,5] <- seqt
      }
      prev_id <- id
      prev_gene <- gene
      prev_spec <- spec
      prev_sero <- sero
      prev_seqt <- seqt
    }
    colnames(vir_rows) <- c('ID','GENE', 'Species', 'Serogroup', 'Sequence Type')
    return(vir_rows)
  }
  gene_rows
}

#Returns the virulence hits for the passed genes in the passed data frame
getGeneHits <- function(gene_rows, meta_rows, genes) {
  hits <- c()
  tot_num <- nrow(meta_rows)
  hits <- c(hits, toString(tot_num))
  for (g in genes) {
    gene_num <- nrow(filter(gene_rows, GENE == g))
    perc <- (gene_num/tot_num)*100
    hits <- c(hits, paste(gene_num, "(", format(perc, digits = 2), ")", sep = ""))
  }
  hits
}

# Returns a vector for showing virulence hits for a species, it'll be a row in the table shown on virulence hits tab;
# adds empty values in the columns for serogroup and/or sequence type if they are shown
getSpeciesHitsRow <- function(spec_name, spec_hits, show_sero, show_seq_type) {
  hits_row  <- c(spec_name)
  if (show_sero) {
    hits_row <- c(hits_row, "")
    if (show_seq_type) {
      hits_row <- c(hits_row, "")
    }
  }
  hits_row <- c(hits_row, spec_hits)
}

# Returns a vector for showing virulence hits for a serogroup, it'll be a row in the table shown on virulence hits tab;
# adds empty values in the columns for species and sequence type
getSerogroupHitsRow <- function(sero_name, sero_hits, show_seq_type) {
  hits_row  <- c("", sero_name)
  if (show_seq_type) {
    hits_row <- c(hits_row, "")
  }
  hits_row <- c(hits_row, sero_hits)
}

# Returns a vector for showing virulence hits for a sequence type, it'll be a row in the table shown on virulence hits tab;
# adds empty values in the columns for species and serogroup
getSeqTypeHitsRow <- function(seq_name, seq_hits) {
  hits_row <- c("", "", seq_name, seq_hits)
}

#Returns the data to be shown in the "Virulence Hits" tab 
getVirulenceHits <- function(gene_rows, meta_rows, show_spec, show_sero, show_seq_type) {
  vir_hits <- data.frame()
  if (!show_spec)
    return(vir_hits)
  vir_rows <- getVirulenceRows(gene_rows)
  # get the names of the species, serogroups, sequence types and genes from the data
  species <- unlist(distinct(vir_rows, Species), use.names = FALSE)
  genes <- unlist(distinct(vir_rows, GENE), use.names = FALSE)
  # set the genes as the column names and the species, serogroups etc. as the row names;
  # also set default values in the rows
  options(stringsAsFactors = FALSE)
  # set the number tested for each species, serogroup etc. (rows) for each gene (columns)
  for (s in species) {
    spec_gene_rows <- filter(vir_rows, Species == s)
    spec_meta_rows <- filter(meta_rows, Species == s)
    # gene hits per species
    spec_hits <- getGeneHits(spec_gene_rows, spec_meta_rows, genes)
    spec_hits <- getSpeciesHitsRow(s, spec_hits, show_sero, show_seq_type)
    vir_hits <- rbind(vir_hits, spec_hits)
    if (show_sero) {
      # gene hits per serogroup
      sero_groups <- unlist(distinct(spec_gene_rows, Serogroup), use.names = FALSE)
      for (sg in sero_groups) {
        sero_gene_rows <- filter(spec_gene_rows, Serogroup == sg)
        sero_meta_rows <- filter(meta_rows, Serogroup == sg)
        sero_hits <- getGeneHits(sero_gene_rows, sero_meta_rows, genes)
        sero_hits <- getSerogroupHitsRow(sg, sero_hits, show_seq_type)
        vir_hits <- rbind(vir_hits, sero_hits)
        if (show_seq_type) {
          # gene hits per sequence type
          seq_types <- unlist(distinct(sero_gene_rows, `Sequence Type`), use.names = FALSE)
          for (st in seq_types) {
            st_gene_rows <- filter(sero_gene_rows, `Sequence Type` == st)
            st_meta_rows <- filter(meta_rows, `Sequence Type` == st)
            seq_hits <- getGeneHits(st_gene_rows, st_meta_rows, genes)
            seq_hits <- getSeqTypeHitsRow(st, seq_hits)
            vir_hits <- rbind(vir_hits, seq_hits)
          }
        }
      }
    }
  }
  col_names <- c("Species")
  if (show_sero) {
    col_names <- c(col_names,"Serogroup")
    if (show_seq_type) {
      col_names <- c(col_names, "Sequence Type")
    }
  }
  col_names <- c(col_names, "# Tested", genes)
  if (nrow(vir_hits) > 0) {
    colnames(vir_hits) <- col_names
  }
  vir_hits
}

# Returns the data to be shown in the "Statistical analysis" tab; the data is pairwise chi-square tests between
# species, serogroups, sequence types, depending upon what's present in the passed virulence hits data frame
getStatAnalysis <- function(vir_hits) {
  stat_analysis <- list()
  #get the rows just for species data
  spec_rows <- filter(vir_hits, Species != "") %>% select(-Serogroup, -'Sequence Type')
  spec_names <- spec_rows[1]
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
  list(fields = fields_with_seq, meta = meta)
}

getFilteredData <- function(fields_with_seq, meta, species = NULL, serogroup = NULL, seq_type = NULL, genes = NULL, 
                            coverage = 0.0, identity = 0.0, show_spec, show_sero, show_seq_type, updateProgress = NULL) {
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
      spec <- fields_with_seq[[i,'Species']]
      sero <- fields_with_seq[[i,'Serogroup']]
      seqt <- fields_with_seq[[i,'Sequence Type']]
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
      aggr_by_gene[j,3] <- spec
      aggr_by_gene[j,4] <- sero
      aggr_by_gene[j,5] <- seqt
      aggr_by_gene[j,6] <- seq
      aggr_by_gene[j,7] <- max_cov
      aggr_by_gene[j,8] <- max_iden
      if (new_row || i == last_index) {
        if (j > 1) {
          k <- ifelse(i == last_index, j, j - 1)
          aggr_by_gene[k,9] <- max_tot_cov
          aggr_by_gene[k,10] <- max_wgt_iden
          # if it's a new gene, update the previous gene's maximum %coverage and %identity
          # in all the rows that are showing it; also set this row's index as the first
          # index for the current gene
          if (new_gene || i == last_index) {
            aggr_by_gene[gene_first_index : k,11] <- gene_max_tot_cov
            aggr_by_gene[gene_first_index : k,12] <- gene_max_wgt_iden
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
    colnames(aggr_by_gene) <- c('ID', 'GENE', 'Species', 'Serogroup', 'Sequence Type',
                                'SEQUENCE', 'SEQ MAX %COVERAGE', 'SEQ MAX %IDENTITY', 
                                'SEQ TOTAL %COVERAGE',  'SEQ WEIGHTED %IDENTITY',
                                'GENE MAX TOTAL %COVERAGE', 'GENE MAX WEIGHTED %IDENTITY')
  }
  #filter by %coverage and/or %identity
  updateProg(updateProgress, "Filtering by %coverage")
  aggr_by_gene <- filterByGT(aggr_by_gene, "SEQ TOTAL %COVERAGE", coverage)
  updateProg(updateProgress, "Filtering by %weighted identity")
  aggr_by_gene <- filterByGT(aggr_by_gene, "SEQ WEIGHTED %IDENTITY", identity)
  
  fields_with_seq <- select(fields_with_seq, -SEQUENCE)
  vir_hits = getVirulenceHits(aggr_by_gene, meta, show_spec, show_sero, show_seq_type)
  list(all_fields = fields_with_seq, gene_fields = aggr_by_gene, virulence_hits = vir_hits)
}
