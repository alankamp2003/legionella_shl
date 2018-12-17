library(readxl)
library(dplyr)
library(ggplot2)
library(lazyeval)

# filter the passed data frame(df) by the passed column name and passed value,
# which is expected to be vector and may be NULL; if is_numeric is TRUE, converts
# the value to numeric type before filtering
filterByIn <- function(df, col_name, value, is_numeric = FALSE) {
  if (!is.null(value)) {
    if (is_numeric)
      value <- as.numeric(value)
    in_criteria <- interp(~y %in% x, .values=list(y = as.name(col_name), x = value))
    df <- filter_(df, in_criteria)
  }
  return(df)
}

# filter the passed data frame(df) by the passed column name and passed value,
# which is expected to be a number and may be null
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

# Sets the passed coverage and identity in the row at index "row" in the passed data frame;
# The coverage and identity are set at "col" and "col + 1" respectively 
setCovAndIden <- function(cov, iden, row, col, aggr_by_gene) {
  aggr_by_gene[row, col] <- cov
  aggr_by_gene[row, col+1] <- iden
}

setCovAndIdenMultiRow <- function(cov, iden, row_start, row_end, col, aggr_by_gene) {
  #print(paste("row_start:", row_start, "row_end:", row_end, "col:", col))
  aggr_by_gene[row_start : row_end, col] <- cov
  aggr_by_gene[row_start : row_end, col+1] <- iden
}

updateProg <- function(updateProgress = NULL, text) {
  # If we were passed a progress update function, call it
  if (is.function(updateProgress))
    updateProgress(detail = text)
}

# Since the rows shown on Gene Fields tab have sequence (node) in them, there can be many repeated rows for the
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
      fac <- gene_rows[[i,'Facility']]
      cov <- gene_rows[[i,'SEQ TOTAL %COVERAGE']]
      iden <- gene_rows[[i,'SEQ WEIGHTED %IDENTITY']]
      if (prev_id != id || prev_gene != gene || prev_spec != spec || prev_sero != sero || prev_seqt != seqt) {
        j <- j + 1
        vir_rows[j,1] <- id
        vir_rows[j,2] <- gene
        vir_rows[j,3] <- spec
        vir_rows[j,4] <- sero
        vir_rows[j,5] <- seqt
        vir_rows[j,6] <- fac
        vir_rows[j,7] <- cov
        vir_rows[j,8] <- iden
      }
      prev_id <- id
      prev_gene <- gene
      prev_spec <- spec
      prev_sero <- sero
      prev_seqt <- seqt
    }
    colnames(vir_rows) <- c('ID','GENE', 'Species', 'Serogroup', 'Sequence Type', 'Facility', 
                            'SEQ TOTAL %COVERAGE', 'SEQ WEIGHTED %IDENTITY')
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

# Returns the data to be shown in the "Virulence Hits" tab and to be used for the scatter plots shown
# under the "Visualization" tab
getVirulenceData <- function(gene_rows, meta_rows, show_spec, show_sero, show_seq_type) {
  vir_hits <- data.frame()
  vir_rows <- NULL
  if (!show_spec)
    return(list(vir_hits = vir_hits, vir_rows = vir_rows))
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
  return(list(vir_hits = vir_hits, vir_rows = vir_rows))
}

# The passed virulence hit is supposed to be in the format "AB(X)", where A and B are integers
# whereas X can be a double; based on the passed pattern, gets rid of a particular part from the hit and
# returns the rest converted to an integer
getInteger <- function(pattern, vir_hit) {
  int_val <- gsub(pattern = pattern, "", vir_hit)
  int_val <- as.integer(int_val)
  int_val
}

# The passed virulence hit is supposed to be in the format "AB(X)", where A and B are integers
# whereas X can be a double; gets rid of the parenthses and X; returns AB converted to an integer
getCount <- function(vir_hit) {
  getInteger("\\(.+)", vir_hit)
}

# The passed virulence hit is supposed to be in the format "AB(X)", where A and B are integers
# whereas X can be a double; gets rid of the parenthses and X; returns AB converted to an integer
getPercentage <- function(vir_hit) {
  getInteger("^.+\\(|\\)", vir_hit)
}

# Returns a string containing zeros repeated the passed number of times
repeatZeros <- function(num) {
  paste(rep("0", num), collapse = "")
}

# Truncates the passed number value to the specified number of decimal digits and returns it 
# as a double; doesn't truncate if the value contains the exponential symbol "e"
getTruncValue <- function(number_val, digits) {
  #print(paste("start", number_val))
  val <- as.character(number_val)
  if (!is.null(number_val) && !is.nan(number_val)) {
     val <- format(number_val, digits = digits)
     #print(val)
     # find out if number_val contains "e"; if it does, return number_val as is;
     # otherwise truncate the digits after the decimal point to the passed number of digits ("digits");
     # appends zeroes at the end if number_val doesn't contain enough digits as specified by "digits" 
     val_vec <- unlist(grep("e", val, ignore.case = TRUE))
     #print(paste("val_vec1", val_vec))
     if (length(val_vec) == 0) {
       val_vec <- unlist(strsplit(val, "\\."))
       #print(paste("val_vec2", val_vec))
       if (length(val_vec) > 0) {
          if (length(val_vec) > 1) {
              #print(paste("length val_vec[2]", length(val_vec[2])))
              num_chars <- nchar(val_vec[2])
              if (num_chars >= digits) {
                 suffix <- substr(val_vec[2], 1, digits)
              } else {
                 suffix <- paste(val_vec[2], repeatZeros(digits-num_chars), sep = "")
              }
          } else {
              suffix <- repeatZeros(digits)
          }
          #val <- paste(val_vec[1], ".", substr(val_vec[2], 1, digits), sep = "")
          #print(paste("suffix", suffix))
          val <- paste(val_vec[1], ".", suffix, sep = "")
          #print(paste("val1", val))
       }
     }
     #print(paste("val2", val))
     #number_val <- as.double(val)
     #print(paste("number_val", number_val))
  }
  #print(paste("end", val))
  #number_val
  val
} 

# Returns a vector containing the number of hits and no hits for the gene at the passed column
# in the row at the passed index of the passed data frame
getHitsVector <- function(vir_hits, row, col) {
  hit <- vir_hits[row,col]
  no_hit <- vir_hits[row,1] - hit
  c(hit, no_hit)
}

# The passed data frame contains two rows of virulence hits; the rows are for hits between
# two entities at the same level i.e. species, serogroup or sequence type; the first column
# is for an entity's total hits across all genes; the subsequent columns are for specific 
# genes' hits; computes Fisher's exact p-values for all genes and returns the values as a vector 
getFisherExactVector <- function(vir_hits) {
  p_vals <- character(ncol(vir_hits)-1)
  j = 1
  for (i in 2:ncol(vir_hits)) {
    r1 <- getHitsVector(vir_hits, 1, i)
    #print(r1)
    r2 <- getHitsVector(vir_hits, 2, i)
    #print(r2)
    matrixa <- matrix(c(r1,r2), 2, byrow=TRUE)
    fisher_test <- fisher.test(matrixa)
    p_vals[j] <- getTruncValue(fisher_test$p.value, 4)
    j <- j + 1
  }
  p_vals
}

# The passed data frame (vir_hits) contains virulence hits in all "numeric" columns at a particular level
# e.g. species; "fisher_exact_vals" contains the Fisher's exact test p-values computed so far at other levels;
# calculates p-values at this level, adds them to fisher_exact_vals and returns that; the first column in the
# returned data frame shows the names of the entities (e.g. species) being compared; adds "prefix" (e.g. "SG"),
# before those names, if specified
getFisherExactVals <- function(vir_hits, fisher_exact_vals, prefix="") {
  #print(spec_hits)
  if (is.null(fisher_exact_vals))
    fisher_exact_vals <- tibble()
  i = 1
  num_rows <- nrow(vir_hits)
  # calculate p-values for each pair of species, serogroups etc. i.e. each pair of rows
  # num_vals only contains the columns with "numeric" data i.e. # tested and gene hits e.g. 20(80)
  num_vals <- vir_hits[, 2:ncol(vir_hits)]
  while (i < num_rows) {
    j <- i+1
    while (j <= num_rows) {
      int_vals <- tibble()
      comp <- c(paste(prefix, vir_hits[i,1], sep = ""), paste(prefix, vir_hits[j,1], sep = ""))
      int_vals <- rbind(int_vals, sapply(num_vals[i,], getCount))
      int_vals <- rbind(int_vals, sapply(num_vals[j,], getCount))
      p_vals <- getFisherExactVector(int_vals)
      p_vals <- c(paste(comp, collapse = " vs. "), p_vals)
      fisher_exact_vals <- rbind(fisher_exact_vals, p_vals)
      j <- j+1
    }
    i <- i+1
  }
  fisher_exact_vals
  #print(int_hits)
}

# Returns the column names used for the table showing Fisher's exact p-values; the first column is for showing 
# which entities e.g. species are being compared; the rest of the columns show gene names
getFisherExactColNames <- function(col_names, spec_hits) {
  if (is.null(col_names)) {
    col_names <- c("Comparison", colnames(spec_hits)[3:ncol(spec_hits)])
    #print(cols)
  }
  col_names
}

# Returns the column names and Fisher's exact p-values to be shown at the species level
getSpecFisherExactVals <- function(vir_hits, col_names, fisher_exact_vals) {
  spec_hits <- filter(vir_hits, Species != "")
  if (nrow(spec_hits) > 1) {
    num_col_start <- 2
    # find out where the "numeric" (e.g "20(80)") columns start; only keep columns that either
    # have the names of the species or are "numeric"
    if ("Serogroup" %in% colnames(spec_hits))
      num_col_start <- num_col_start + 1
    if ("Sequence Type" %in% colnames(spec_hits))
      num_col_start <- num_col_start + 1
    spec_hits <- spec_hits[, c(1, num_col_start:ncol(spec_hits))]
    col_names <- getFisherExactColNames(col_names, spec_hits)
    fisher_exact_vals <- getFisherExactVals(spec_hits, fisher_exact_vals)
  }
  list(col_names = col_names, fisher_exact_vals = fisher_exact_vals)
}

# Returns the column names and Fisher's exact p-values to be shown at the serogroup level
getSeroFisherExactVals <- function(vir_hits, col_names, fisher_exact_vals) {
  if ("Serogroup" %in% colnames(vir_hits)) {
    sero_hits <- filter(vir_hits, Serogroup != "")
    if (nrow(sero_hits) > 1) {
      num_col_start <- 3
      # find out where the "numeric" (e.g "20(80)") columns start; only keep columns that either
      # have the names of the serogroup or are "numeric"
      if ("Sequence Type" %in% colnames(sero_hits))
        num_col_start <- num_col_start + 1
      sero_hits <- sero_hits[, c(2, num_col_start:ncol(sero_hits))]
      col_names <- getFisherExactColNames(col_names, sero_hits)
      fisher_exact_vals <- getFisherExactVals(sero_hits, fisher_exact_vals, prefix = "SG")
    }
  }
  list(col_names = col_names, fisher_exact_vals = fisher_exact_vals)
}

# Returns the column names and Fisher's exact p-values to be shown at the sequence type level
getSeqTypeFisherExactVals <- function(vir_hits, col_names, fisher_exact_vals) {
  if ("Sequence Type" %in% colnames(vir_hits)) {
    seq_type_hits <- filter(vir_hits, `Sequence Type` != "") 
    if (nrow(seq_type_hits) > 1) {
      num_col_start <- 4
      # only keep columns that either have the names of the sequence or are "numeric"
      seq_type_hits <- seq_type_hits[, c(3, num_col_start:ncol(seq_type_hits))]
      col_names <- getFisherExactColNames(col_names, seq_type_hits)
      fisher_exact_vals <- getFisherExactVals(seq_type_hits, fisher_exact_vals, prefix = "ST")
    }
  }
  list(col_names = col_names, fisher_exact_vals = fisher_exact_vals)
}

# Goes through the passed data frame containing virulence hits at multiple levels e.g. species, serogroup, sequence type and returns a data frame
# containing Fisher's exact p-values at all of those levels
getAllFisherExactVals <- function(vir_hits) {
  col_names <- NULL
  fisher_exact_vals <- NULL

  spec_vals <- getSpecFisherExactVals(vir_hits, col_names, fisher_exact_vals)
  sero_vals <- getSeroFisherExactVals(vir_hits, spec_vals$col_names, spec_vals$fisher_exact_vals)
  seq_type_vals <- getSeqTypeFisherExactVals(vir_hits, sero_vals$col_names, sero_vals$fisher_exact_vals)
  col_names <- seq_type_vals$col_names
  fisher_exact_vals <- seq_type_vals$fisher_exact_vals
  #print(fisher_exact_vals)
  if (!is.null(fisher_exact_vals)) {
    colnames(fisher_exact_vals) <- col_names
  }
  fisher_exact_vals
}

# Returns a scatter plot to be shown under the "Visualization" tab; the passed data fram provides the data for the plot;
# the passed indexes specify the columns for the x axis, y axis, shape and color respectively
getPlot <- function(vir_rows, x_col, y_col, shape_col, color_col, x_lab, y_lab) {
  plot <- ggplot(vir_rows, aes(x = vir_rows[x_col][[1]], y = vir_rows[y_col][[1]], shape = as.factor(vir_rows[shape_col][[1]]), 
                               color = vir_rows[color_col][[1]])) + geom_point(size = 5) +
                               xlab(x_lab) + ylab(y_lab) + scale_shape_discrete(shape_col) + scale_color_discrete(color_col) +
                               guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
                               theme(
                                     legend.box = "horizontal", 
                                     legend.position = "bottom", 
                                     axis.title.x = element_text(color="#993333", size=14, face="bold"),
                                     axis.title.y = element_text(color="#993333", size=14, face="bold")
                                    )
                              
  #geom_density()
  plot
}

# Returns the index of the column with name column from the passed vector of column names
# returns zero if the column was not found in the vector
getColIndex <- function(cols, col) {
  inds <- which(cols == col)
  ind <- ifelse(length(inds) > 0, inds[1], 0)
  ind
}

# Returns a list containing the scatter plots to be shown under the "Visualization" tab;
# scatter_noise is used to specify the noise introduced in %identity and %coverage to visualize
# the clusters more clearly
getPlots <- function(vir_rows, scatter_noise) {
  plots <- list()
  if (!is.null(vir_rows)) {
    cols <- colnames(vir_rows)
    x_col <- "IDENTITY"
    y_col <- "COVERAGE"
    color_col <- "Facility"
    num_rows <- nrow(vir_rows)
    #print(paste("num vir_rows", num_rows))
    x_lab <- paste("Seq Weighted % Identity +/-", scatter_noise)
    y_lab <- paste("Seq Total \n % Coverage \n +/-", scatter_noise)
    vir_rows <- mutate(vir_rows, IDENTITY = `SEQ WEIGHTED %IDENTITY` + runif(num_rows, -scatter_noise, scatter_noise)) %>%
                mutate(COVERAGE = `SEQ TOTAL %COVERAGE` + runif(num_rows, -scatter_noise, scatter_noise))
    #print(paste("num vir_rows mute", num_rows))
    plots$spec_plot <- getPlot(vir_rows, x_col, y_col, "Species", color_col, x_lab, y_lab)
    plots$sero_plot <- getPlot(vir_rows, x_col, y_col, "Serogroup", color_col, x_lab, y_lab)
    plots$seq_type_plot <- getPlot(vir_rows, x_col, y_col, "Sequence Type", color_col, x_lab, y_lab)
  }
  plots
}

# Creates and returns a list of charts from the passed data frame; the first column is the name of the level e.g. "Species" and shows the 
# entities at that level; the second column onwards show the percentage of virulence hits for each gene for a given entity;
# the charts show bars depicting the percentage for each gene for a given entity
getChartList <- function(chart_hits) {
  cols <- colnames(chart_hits)
  chart_list <- list()
  for (i in 2:ncol(chart_hits)) {
    local({
      i <- i
      #print(paste("i:", i, cols[i]))
      #print(paste("chart_list[,i]", chart_list[,i]))
      #ylabel <- ifelse(i == 2, "Percentage", "")
      gplot <- ggplot(chart_hits, aes(x = chart_hits[,1], y = chart_hits[,i], fill = chart_hits[,1], ymax=100, ymin=0)) +
                          geom_bar(stat="identity", position = "dodge", width = 0.2, show.legend = FALSE) +
                          #xlab(cols[1]) + ylab("Percentage") +  ggtitle(cols[i]) +
                          xlab("") +ylab("Percentage") +  ggtitle(cols[i]) + coord_cartesian(ylim=c(0,100)) +
                          # + theme_bw() +                                                                                              
                          theme(
                            plot.title = element_text(color="red", size=16, face="bold.italic", hjust = 0.5),
                            axis.title.x = element_text(color="#993333", size=14, face="bold"),
                            axis.title.y = element_text(color="#993333", size=14, face="bold")
                          )
      chart_list[[i-1]] <<- gplot
    })
  }
  #print(length(chart_list))
  chart_list
}

# Returns the species level charts to be shown in the "Visualization" tab
getSpecChartList <- function(vir_hits) {
  spec_charts <- list()
  cols <- colnames(vir_hits)
  #print(cols)
  if ("Species" %in% cols) {
    spec_hits <- filter(vir_hits, Species != "")
    if (nrow(spec_hits) > 0) {
      gene_col_start <- 3
      # find out where the columns showing virulence hits for genes (e.g "20(80)") start;
      # only keep columns that either have the names of the species or virulence hits
      if ("Serogroup" %in% cols)
        gene_col_start <- gene_col_start + 1
      if ("Sequence Type" %in% cols)
        gene_col_start <- gene_col_start + 1
      spec_hits <- spec_hits[, c(1, gene_col_start:ncol(spec_hits))]
      num_cols <- ncol(spec_hits)
      spec_hits[2:num_cols] <- lapply(spec_hits[2:num_cols], getPercentage)
      #print(spec_hits)
      spec_charts <- getChartList(spec_hits)
    }
  }
  spec_charts
}

# Returns the serogroup level charts to be shown in "Visualization" tab
getSeroChartList <- function(vir_hits) {
  sero_charts <- list()
  cols <- colnames(vir_hits)
  if ("Serogroup" %in% cols) {
    sero_hits <- filter(vir_hits, Serogroup != "")
    if (nrow(sero_hits) > 0) {
      gene_col_start <- 4
      # find out where the columns showing virulence hits for genes (e.g "20(80)") start;
      # only keep columns that either have the names of the species or virulence hits
      if ("Sequence Type" %in% cols)
        gene_col_start <- gene_col_start + 1
      sero_hits <- sero_hits[, c(2, gene_col_start:ncol(sero_hits))]
      num_cols <- ncol(sero_hits)
      sero_hits[2:num_cols] <- lapply(sero_hits[2:num_cols], getPercentage)
      #print(sero_hits)
      sero_charts <- getChartList(sero_hits)
    }
  }
  sero_charts
}

# Returns the sequence type level charts to be shown in "Visualization" tab
getSeqTypeChartList <- function(vir_hits) {
  seq_type_charts <- list()
  cols <- colnames(vir_hits)
  if ("Sequence Type" %in% cols) {
    seq_type_hits <- filter(vir_hits, `Sequence Type` != "")
    if (nrow(seq_type_hits) > 0) {
      gene_col_start <- 5
      # find out where the columns showing virulence hits for genes (e.g "20(80)") start;
      # only keep columns that either have the names of the species or virulence hits
      seq_type_hits <- seq_type_hits[, c(3, gene_col_start:ncol(seq_type_hits))]
      num_cols <- ncol(seq_type_hits)
      seq_type_hits[2:num_cols] <- lapply(seq_type_hits[2:num_cols], getPercentage)
      #print(seq_type_hits)
      seq_type_charts <- getChartList(seq_type_hits)
    }
  }
  seq_type_charts
}

# Returns the list containing charts to be shown in "Visualization" tab
getCharts <- function(vir_hits) {
  spec_charts <- getSpecChartList(vir_hits)
  sero_charts <- getSeroChartList(vir_hits)
  seq_type_charts <- getSeqTypeChartList(vir_hits)
  list(spec_charts = spec_charts, sero_charts = sero_charts, seq_type_charts = seq_type_charts)
}

# link various sheets in the excel file to get the final list of fields
# (including SEQUENCE, which will be removed later); some fields are renamed 
# to make the linking work
getInitialData <- function(db_file, link_file, meta_file, updateProgress = NULL) {
  updateProg(updateProgress, "Extracting fields from pipeline run")
  #vfdb <- read_excel(db_file,sheet='VFDB') %>%
  vfdb <- read_excel(db_file) %>%
    select(FILE, GENE, SEQUENCE, `%COVERAGE`, `%IDENTITY`)
  updateProg(updateProgress, "Renaming fields and linking sheets")
  #linking <- read_excel(db_file,sheet='Linking') %>% rename(FILE = File)
  linking <- read_excel(link_file) %>% rename(FILE = File)
  
  join_by_file <- inner_join(vfdb, linking, by = "FILE")
  
  #meta <- read_excel(db_file,sheet='Meta')
  meta <- read_excel(meta_file)
  
  join_by_id <- inner_join(join_by_file, meta, by = "ID#") %>% rename(ID = 'ID#')
  
  updateProg(updateProgress, "Extracting final fields")
  fields_with_seq <- select(join_by_id, ID, GENE, SEQUENCE, Species, Serogroup, `Sequence Type`, `Collected Date`, `%COVERAGE`,
                            `%IDENTITY`, Facility, Room, Source) %>%
    arrange(ID, GENE, SEQUENCE)
  list(fields = fields_with_seq, meta = meta)
}

getFilteredData <- function(fields_with_seq, meta, species = NULL, serogroup = NULL, seq_type = NULL, genes = NULL, 
                            coverage = 0.0, identity = 0.0, show_spec, show_sero, show_seq_type, 
                            scatter_noise = 0.0, updateProgress = NULL) {
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
    prev_spec <- ""
    prev_sero <- ""
    prev_seqt <- ""
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
    l <- 0
    new_row <- FALSE
    new_gene <- FALSE
    prev_gene_first_index <- -1
    gene_first_index <- 1
    last_index = nrow(fields_with_seq)
    #vir_rows <- tibble()
    
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
      fac <- fields_with_seq[[i,'Facility']]
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
      aggr_by_gene[j,6] <- fac
      aggr_by_gene[j,7] <- seq
      aggr_by_gene[j,8] <- max_cov
      aggr_by_gene[j,9] <- max_iden
      if (new_row || i == last_index) {
        if (j > 1) {
          #k <- ifelse(i == last_index, j, j - 1)
          k <- j - 1
          aggr_by_gene[k,10] <- max_tot_cov
          aggr_by_gene[k,11] <- max_wgt_iden
          if (i == last_index) {
            k <- j
            aggr_by_gene[k,10] <- max_tot_cov
            aggr_by_gene[k,11] <- max_wgt_iden
          }
          # if it's a new gene, update the previous gene's maximum %coverage and %identity
          # in all the rows that are showing it; also set this row's index as the first
          # index for the current gene
          if (new_gene || i == last_index) {
            k <- j - 1
            aggr_by_gene[gene_first_index : k,12] <- gene_max_tot_cov
            aggr_by_gene[gene_first_index : k,13] <- gene_max_wgt_iden
            if (i == last_index) {
              k <- j
              aggr_by_gene[gene_first_index : k,12] <- gene_max_tot_cov
              aggr_by_gene[gene_first_index : k,13] <- gene_max_wgt_iden
            }
            gene_first_index <- j
            gene_max_tot_cov <- -1.0
            gene_max_wgt_iden <- -1.0
          }
        }
        # if (prev_id != id || prev_gene != gene || prev_spec != spec || prev_sero != sero || prev_seqt != seqt) {
        #   l <- l + 1
        #   vir_rows[l,1] <- id
        #   vir_rows[l,2] <- gene
        #   vir_rows[l,3] <- spec
        #   vir_rows[l,4] <- sero
        #   vir_rows[l,5] <- seqt
        #   vir_rows[l,6] <- fac
        #   vir_rows[l,7] <- max_tot_cov
        #   vir_rows[l,8] <- max_wgt_iden
        # }
        max_tot_cov <- 1.0
        max_wgt_iden <- -1.0
        new_row <- FALSE
      }
      prev_id <-id
      prev_gene <- gene
      prev_seq <- seq
      prev_spec <- spec
      prev_sero <- sero
      prev_seqt <- seqt
    }
    colnames(aggr_by_gene) <- c('ID', 'GENE', 'Species', 'Serogroup', 'Sequence Type',
                                'Facility', 'SEQUENCE', 'SEQ MAX %COVERAGE', 'SEQ MAX %IDENTITY',
                                'SEQ TOTAL %COVERAGE',  'SEQ WEIGHTED %IDENTITY',
                                'GENE MAX TOTAL %COVERAGE', 'GENE MAX WEIGHTED %IDENTITY')
    # colnames(vir_rows) <- c('ID','GENE', 'Species', 'Serogroup', 'Sequence Type', 'Facility', 
    #                         'SEQ TOTAL %COVERAGE', 'SEQ WEIGHTED %IDENTITY')
  }
  #filter by %coverage and/or %identity
  updateProg(updateProgress, "Filtering by %coverage")
  aggr_by_gene <- filterByGT(aggr_by_gene, "SEQ TOTAL %COVERAGE", coverage)
  updateProg(updateProgress, "Filtering by %weighted identity")
  aggr_by_gene <- filterByGT(aggr_by_gene, "SEQ WEIGHTED %IDENTITY", identity)
  
  fields_with_seq <- select(fields_with_seq, -SEQUENCE)
  vir_data <- getVirulenceData(aggr_by_gene, meta, show_spec, show_sero, show_seq_type)
  #vir_data <- getVirulenceData(vir_rows, meta, show_spec, show_sero, show_seq_type)
  fisher_exact_vals <- getAllFisherExactVals(vir_data$vir_hits)
  charts <- getCharts(vir_data$vir_hits)
  plots <- getPlots(vir_data$vir_rows, scatter_noise)
  list(all_fields = fields_with_seq, gene_fields = aggr_by_gene, virulence_hits = vir_data$vir_hits, 
                    stats_analysis = fisher_exact_vals, plots = plots, charts = charts)
}
