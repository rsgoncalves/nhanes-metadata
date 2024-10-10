install.packages("devtools")
devtools::install_github("cjendres1/nhanes")
library(nhanesA)
library(progress)
library(rvest)
library(dplyr)
library(xml2)


TABLES_WITHOUT_CODEBOOKS <- c("DOC_2000", "DRXFCD_I", "DRXFCD_J", "DRXFMT", "FOODLK_C", "FOODLK_D", 
                              "LA_DEMO", "P_DRXFCD", "PAX80_G_R", "PAXLUX_G_R", "POOLTF_D", "POOLTF_E", 
                              "RXQ_DRUG", "SSUIFG_R", "VARLK_C", "VARLK_D", "VID_2_00", "YDQ", "DRXFMT_B")


# Get variable-level metadata from cached NHANES-snapshot at the specified location
get_variables_metadata <- function(nhanes_snapshot_folder) {
  tables <- list.files(path = nhanes_snapshot_folder) 
  # Write a table containing metadata about all variables
  all_variables <- data.frame(matrix(ncol=7, nrow = 0))
  colnames(all_variables) <-c("Variable", "Table", "SASLabel", "EnglishText", "EnglishInstructions", "Target", "UseConstraints")
  variables_metadata_file <- paste(output_folder, "nhanes_variables.tsv", sep="")
  print(paste("Extracting variables metadata to", variables_metadata_file))
  write.table(all_variables, file=variables_metadata_file, row.names=FALSE, sep="\t", na="")
  
  # Write a table containing the codebooks of all variables
  all_codebooks <- data.frame(matrix(ncol=7, nrow = 0))
  variables_codebooks_file <- paste(output_folder, "nhanes_variables_codebooks.tsv", sep="")
  colnames(all_codebooks) <-c("Variable", "Table", "CodeOrValue", "ValueDescription", "Count", "Cumulative", "SkipToItem")
  write.table(all_codebooks, file=variables_codebooks_file, row.names=FALSE, na="", sep="\t")
  pb <- progress_bar$new(total = length(tables))
  for (i in 1:length(tables)) {
    pb$tick()
    table_name <- sub("\\.html$", "", tables[[i]])
    variable_details <- get_variables_in_table(table_name=table_name, nhanes_snapshot_folder=nhanes_snapshot_folder)
    if(length(variable_details) > 0) {
      table_variables <- variable_details[[1]]
      write.table(table_variables, file=variables_metadata_file, sep="\t", na="", append=TRUE, col.names=FALSE, row.names=FALSE)
      table_codebooks <- variable_details[[2]]
      tryCatch({
        write.table(table_codebooks, file=variables_codebooks_file, sep="\t", na="", append=TRUE, col.names=FALSE, row.names=FALSE)
      }, 
      error=function(e) {
        print(paste("Error writing out table: ", table_name, ":", conditionMessage(e)))
      })
    }
  }
}

# Get metadata about the variables in the given NHANES table (specified by the table name) and data group.
get_variables_in_table <- function(table_name, nhanes_snapshot_folder) {
  all_variables <- data.frame(matrix(ncol=7, nrow=0))
  all_variable_codebooks <- data.frame(matrix(ncol=7, nrow=0))
  table_codebooks = {}
  tryCatch({
    codebook_url <- paste0(nhanes_snapshot_folder, table_name, ".html")
    table_codebooks <- nhanesCodebookFromURL(url=codebook_url)
    doc_html_file <- read_html(codebook_url)
  }, error=function(msg) {
    log_missing_codebook("Error in nhanesA::nhanesCodebook()", conditionMessage(msg), table_name=table_name, variable_name="")
  }, warning=function(msg) {
    log_missing_codebook("Warning in nhanesA::nhanesCodebook()", conditionMessage(msg), table_name=table_name, variable_name="")
  })
  
  if(is.null(table_codebooks)) {
    if(!table_name %in% TABLES_WITHOUT_CODEBOOKS) {
      print(paste("nhanesA::nhanesCodebook() returns NULL for table ", table_name))
    }
    return(list())
  }
  else {
    for (i in 1:length(table_codebooks)) {
      variable <- table_codebooks[[i]]
      variable_name <- variable['Variable Name:'][[1]]
      variable_details <- table_codebooks[variable_name][[1]]
      if (!is.null(variable_details)) {
        label <- clean(variable_details['SAS Label:'][[1]])
        text <- clean(variable_details['English Text:'][[1]])
        instructions <- variable_details['English Instructions:'][[1]]
        if(is.null(instructions) || all(is.na(instructions))) {
          instructions <- ""
        } else {
          instructions <- clean(instructions)
        }
        targets <- ""
        for (i in seq_along(variable_details)) {
          if ('Target:' %in% names(variable_details[i])) {
            target <- clean(variable_details[[i]])
            if (targets == "") {
              targets <- target
            } else {
              targets <- paste(targets, target, sep = " *AND* ")
            }
          }
        }
        use_constraints <- has_use_constraints(html_file=doc_html_file)
        variable_details_vector <- c(toupper(variable_name), toupper(table_name), label, text, instructions, targets, use_constraints)
        all_variables[nrow(all_variables) + 1, ] <- variable_details_vector
        variable_codebook <- variable_details[variable_name][[1]]
        if (!is.null(variable_codebook)) {
          if (!all(is.na(variable_codebook))) {
            if (variable_name != "SEQN" && variable_name != "SAMPLEID") {
              if (length(variable_codebook) == 0) {
                print(paste("Variable codebook list for variable", variable_name, "in table", table_name, " is empty"))
              }
              if ("Value Description" %in% names(variable_codebook) && !is.null(variable_codebook[["Value Description"]])) {
                variable_codebook[["Value Description"]] <- clean(variable_codebook[["Value Description"]])
              } else {
                variable_codebook[["Value Description"]] <- ""
              }
              names(variable_codebook)[names(variable_codebook) == 'Value Description'] <- 'ValueDescription'
              # Make sure NA values are included in the output table
              variable_codebook$`Code or Value`[is.na(variable_codebook$`Code or Value`)] <- "NA"
              variable_codebook$ValueDescription[is.na(variable_codebook$ValueDescription)] <- "NA"
              variable_codebook <- cbind('Table'=table_name, variable_codebook)
              variable_codebook <- cbind('Variable'=variable_name, variable_codebook)
              variable_codebook$Table <- toupper(variable_codebook$Table)  # Ensure Variable and Table IDs are capitalized (#20)
              variable_codebook$Variable <- toupper(variable_codebook$Variable)
              all_variable_codebooks <- rbind(all_variable_codebooks, variable_codebook)
            }
          }
        } else {
          if (variable_name != "SEQN" && variable_name != "SAMPLEID") {
            log_missing_codebook("Warning", "Null or NA codebook for this variable", table_name = table_name, variable_name = variable_name)
          }
        }
      } else {
        if (!is.null(variable_name) && variable_name != "SEQN" && variable_name != "SAMPLEID") {
          log_missing_codebook("Warning", "nhanesA::nhanesCodebook() does not return a codebook for this variable", 
                               table_name = table_name, variable_name = variable_name)
        }
      }
    }
    return(list(all_variables, all_variable_codebooks))
  }
}

has_use_constraints <- function(html_file) {
  # Find all h4 elements
  h4_elements <- xml_find_all(html_file, "//h4")
  
  # Check if any h4 element contains the specific text
  return(any(grepl("RDC Only", xml_text(h4_elements))))
}

log_missing_codebook <- function(exception_type, exception_msg, table_name, variable_name) {
  if( !(table_name %in% TABLES_WITHOUT_CODEBOOKS) ) {
    log_msg <- paste(exception_type, ": ", exception_msg, " (Table-Variable: ", 
                     table_name, "-", variable_name, ")", sep="")
    cat(log_msg, "\n", file=log_file, append=TRUE)
    cat(table_name, "\t", variable_name, "\n", file=missing_codebooks_file, append=TRUE)
    print(log_msg)
  }
}

# Remove carriage return, new line, comma, backslash and quote characters
clean <- function(text) {
  text <- gsub("[\r\n,\\\"]", "", text)
  return(text)
}


# Prepare outputs
output_folder <- "/Users/rsgoncalves/Documents/Workspace/nhanes-metadata/metadata/"
if (!dir.exists(output_folder)) {
  dir.create(output_folder, showWarnings=FALSE, recursive=TRUE)
}
log_file <- paste(output_folder, "log.txt", sep="")
missing_codebooks_file <- paste(output_folder, "missing_codebooks.tsv", sep="")

# Get metadata about all variables in cached NHANES tables
start_vars = proc.time()
get_variables_metadata(nhanes_snapshot_folder = "/Users/rsgoncalves/Documents/Workspace/nhanes-snapshot/docs/")
diff <- proc.time() - start_vars
print(paste("Variables metadata acquired in ", diff[3][["elapsed"]], " seconds"))

# Write out to log file the date when metadata finished downloading
cat(paste("Downloaded on:", format(Sys.time(), format="%m-%d-%YT%H:%M:%S")), file=log_file, append=TRUE)
print("finished")