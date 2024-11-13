# please install package httr
library(httr)

# const values
url <- "http://ophid.utoronto.ca/pathDIP/Http_API"

searchOnUniprot_IDs <- function(IDs, component, sources, pathwayTypes, cutoffQ) {
    
    parameters <- list(
        typeChoice = "Uniprot ID",
        IDs = IDs,
        TableName = component,
        DataSet = sources,
        PathwayTypes = pathwayTypes,
        CutoffQ = cutoffQ
    )
    
    # ... send http POST
    res <- POST(url, body = parameters, encode = "form", verbose())
}

searchOnGenesymbols <- function(IDs, component, sources, pathwayTypes, cutoffQ) {
    
    parameters <- list(
        typeChoice = "Gene Symbol",
        IDs = IDs,
        TableName = component,
        DataSet = sources,
        PathwayTypes = pathwayTypes,
        CutoffQ = cutoffQ
    )
    
    # ... send http POST
    res <- POST(url, body = parameters, encode = "form", verbose())
}

searchOnEntrez_IDs <- function(IDs, component, sources, pathwayTypes, cutoffQ) {
    parameters <- list(
        typeChoice = "Egid",
        IDs = IDs,
        TableName = component,
        DataSet = sources,
        PathwayTypes = pathwayTypes,
        CutoffQ = cutoffQ
    )
    
    # ... send http POST
    res <- POST(url, body = parameters, encode = "form", verbose())
}


# make results-map as keyword - value
makeMap <- function(res) {

    ENTRY_DEL = "\001"
    KEY_DEL = "\002"

    response = content(res, "text")

    arr = unlist(strsplit(response, ENTRY_DEL, fixed = TRUE))

    list_map <- list("")
    vec_map_names <- c("");

    for (str in arr) {
        arrKeyValue = unlist(strsplit(str, KEY_DEL, fixed = TRUE));

        if (length(arrKeyValue) > 1) {
            list_map[length(list_map) + 1] <- arrKeyValue[2]
            vec_map_names[length(vec_map_names) + 1] <- arrKeyValue[1]
        }
    }

    names(list_map) <- vec_map_names

    list_map
}


# Pick-up right sample for your search
# Note: Adjust 'ophid.utoronto.ca' if you are using our development server.

if (TRUE) { # change to FALSE to comment this sample

# First, install and load necessary packages if they are not already installed
# install.packages("httr") # Uncomment if not already installed
# install.packages("jsonlite") # Uncomment if not already installed

# Load necessary packages
library(httr)
library(jsonlite)

# Function to handle command line arguments and extract IDs
getIDsFromArgs <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  # Check if the help flag is present
  if ('-h' %in% args || '--help' %in% args) {
    cat("Usage: Rscript Pathway.R <GeneSymbols>\n")
    cat("GeneSymbols: Comma-separated list of gene symbols to query\n")
    cat("Note: The script will only be executed successfully when you provide multiple genes\n")
    cat("Example: Rscript Pathway.R \"ATF2,FOS,JUN,MAPK1\"\n")
    cat("Help: Rscript Pathway.R -h\n")
    
    quit(status = 0) # Exit the script after displaying help
  }
  if (length(args) == 0) {
    stop("No gene symbols were provided. Please pass gene symbols as arguments when calling the script.")
  }
  return(args[1])
}

# Extract IDs from command line arguments
IDs <- getIDsFromArgs()

# Rest of the inputs
component <- "Literature curated (core) pathway memberships"
sources <- "ACSN2,BioCarta,HumanCyc,KEGG,Panther_Pathway,PathBank,PharmGKB,REACTOME,MetabolicAtlas,SIGNOR 3.0,UniProt_Pathways,WikiPathways"
pathwayTypes <- "Cellular processes and organization,Diseases,Drugs and vitamins,Environmental information processing,Genetic information processing,Metabolism,Organismal systems"
significance <- "0.05" # Recommended

# Here we would define the 'searchOnGenesymbols' function
# Since it's not provided, I'll assume it's defined elsewhere in your environment

# Search on gene symbols using the provided function
# The searchOnGenesymbols function is assumed to be previously defined

# Search on gene symbols using the provided function
# The searchOnGenesymbols function is assumed to be previously defined
res <- searchOnGenesymbols(IDs, component, sources, pathwayTypes, significance)

# Check if the 'res' is not NULL and has a successful status code
if (!is.null(res) && status_code(res) == 200) {
    list_map <- makeMap(res)

    # Output results
    cat("\nSearch on Gene Symbols:\n")
    cat("Generated at: ", unlist(list_map["GeneratedAt"]), "\n")
    cat("IDs: ", unlist(list_map["IDs"]), "\n")
    cat("DataComponent: ", unlist(list_map["TableName"]), "\n")
    cat("Sources: ", unlist(list_map["DataSet"]), "\n")
    cat("PathwayTypes: ", unlist(list_map["PathwayTypes"]), "\n")
    cat("q-value (FDR: BH-method): ", unlist(list_map["CutoffQ"]), "\n")

    cat("\nNote: Pathways, where only one ID from the input list was found, were excluded from results.\n")

    cat("\nSummary size: ", unlist(list_map["SummarySize"]), "\n")
    cat("Summary: \n", unlist(list_map["Summary"]), "\n") # formatted as tab - delimited spreadsheet

    cat("\nDetails size: ", unlist(list_map["DetailsSize"]), "\n")
    cat("Details: \n", unlist(list_map["Details"]), "\n") # formatted as tab - delimited spreadsheet

} else {
    responseCode <- ifelse(is.null(res), "NULL", status_code(res))
    cat("Error: Response Code is ", responseCode, ", or the 'res' object is NULL.\n")
}

cat("Script execution finished.\n")}
