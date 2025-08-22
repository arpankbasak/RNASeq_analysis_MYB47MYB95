# List of paths
data.files <- "./data/"
fig.path <- "./figures/"
stat.path <- "./statistics/"
out.path <- "./output/"
gtf.path <- "/klaster/scratch/CAMDA/CAMDA2019/MetagenomicForensics/REPOSITORY/RNAseq/TAIR10db/ref_seq/TAIR.10/gtf/"
bam.path <- paste(data.files, "bam_files/", sep = "")
salmon.path <- "/klaster/scratch/CAMDA/CAMDA2019/MetagenomicForensics/REPOSITORY/RNAseq/aligned_salmon"
# Parameter file for plotting and indexing
# General parameters for plotting
genotype <- data.frame(name = c("Col-0", "myb47/95", "myc2/3/4"),
                 short = c("WT", "myb", "myc.tKO"),
                 col.idx = c("black", "dodgerblue", "indianred"),
                 lty.idx = c(1, 5, 5),
                 pch.idx = c(15, 16, 17), 
                 stringsAsFactors = FALSE)

treatment <- data.frame(name = c("DW", "JA"),
                        col.idx = c("grey", "black"),
                        lty.idx = c(1, 5),
                        pch.idx = c(1, 16),
                        alpha.idx = c(0.4, 1),
                        stringsAsFactors = FALSE)

label.idx <- c("Col-0", 
               expression(italic("myb47/95")),
               expression(italic("myc3/4/5")))

# Parameters for statistics
p.adj.method <- "fdr"
FC_threshold <- 1
alpha <- 0.01
threshold <- .05
noise <- 0.0001
seed <- 1
set.seed(seed)

# Theme set for plotting canvas inspired from RTN
require(ggplot2)
theme_AKB <- theme_bw() +
    theme(text = element_text(size = 16, hjust = 0.3, vjust = 0.3),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          legend.key=element_blank(),
          legend.text.align=0.5,
          legend.position="top",
          strip.text=element_text(face="bold", size = 16),
          axis.line=element_line(),
          axis.line.x=element_line(size = 1),
          axis.line.y=element_line(size = 1),
          panel.background=element_rect(fill="transparent", colour=NA),
          plot.background=element_rect(fill="transparent", colour=NA),
          strip.background=element_rect(fill="transparent", colour=NA),
          strip.placement="outside")

# Colour palette
gradient.low = "darkmagenta"
gradient.high = "darkgreen"
gradient.mid = "whitesmoke"
gradient.sig.high = "darkred"
gradient.sig.low = "darkblue"
gradient.wide = "palegoldenrod"
gradient.na = "lightgrey"
col.pal = c("lightblue", "indianred", 
            "darkblue", "darkred")

# Saturate function for plotting
saturate <- function(x){
    max <- quantile(x, .99)
    min <- quantile(x, .01)
    
    idx <- (x > max)
    x[idx] <- max
    
    idx <- x < min
    x[idx] <- min
    return(x)
}

# Biomart
call_biomart <- function(x, columns=NULL, mart=NULL, filter=NULL, distinct=TRUE){
    
    options(fill = TRUE)
    if(is.null(mart)){
        mart_obj <- biomaRt::useMart("plants_mart", 
                             host = "plants.ensembl.org", 
                             dataset = "athaliana_eg_gene")
    }
    if(is.null(columns)){
        columns <- c("ensembl_gene_id", "external_gene_name", "gene_biotype",
                     "chromosome_name", "start_position", "end_position", "strand",
                     "go_id", "name_1006", "definition_1006", 
                     "po_id", "po_name_1006", "po_definition_1006")
        
    }
    if(is.null(filter)){
        filter = "ensembl_gene_id"
    }

    bm_res <- biomaRt::getBM(attributes = columns, 
                             values = x,
                             filters = filter,
                             mart = mart_obj) 
    
    #if(distinct != FALSE){
     #      bm_res <- dplyr::distinct(bm_res, paste(filter, sep = ""), .keep_all = T)
      #  }
    return(bm_res)
    
}

# Gene Ontology
call_GO <- function(x, columns=NULL, by = "PATH"){
    
    require(org.At.tair.db)
    require(ensembldb)
    if(is.null(columns)){
        flag = 0
        columns = c("GO",
                    "GENENAME",
                    "PATH",  
                    "SYMBOL",  
                    "ENZYME",   
                    "ONTOLOGY",  
                    "ARACYC")
    }
    qr <- ensembldb::select(org.At.tair.db, 
                        keys = as.character(x), 
                        keytype = by, 
                        columns = columns)
    if(by == "PATH"){

        go_obj <- qr %>% group_by(PATH) %>% 
                summarise(symbol = paste(unique(SYMBOL), collapse = "; "),
                          enzyme = paste(unique(ENZYME), collapse = "; "),
                          aracyc = paste(unique(ARACYC), collapse = "; "),
                          gene_name = paste(unique(GENENAME), collapse = "; "),
                          path_name = paste(unique(PATH), collapse = "; "),
                          go_id = paste(unique(GO), collapse = "; ")) %>% 
        data.frame(.)
    }else {
        
        go_obj <- qr %>% data.frame(.)
    
    }
    
    
    return(go_obj)

}

# Function to tidy concatenated annotations
flow <- function(x){
    temp <- NULL
    for(i in 1:nrow(x)){
        
        dev <- unlist(str_split(x[i,2]), use.names = FALSE)        
        for(j in 1:length(dev))
            temp[j,"col2"] <- dev[j]
            temp[j, "col1"] <- x[i,1]
        
        temp <- rbind(temp, temp)
        dev <- NULL
    }
    return(temp)
}

# Optimum clusters -- @RTN
kmeansAIC <- function(fit){
    m <- ncol(fit$centers)
    n <- length(fit$cluster)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(c(
        AIC=(D + 2*m*k),
        BIC=(D + log(n)*m*k)
        )
    )
}

repZ <- function(fit, Z, x) {
    ids <- names(fit$cluster[fit$cluster==x])
    
    idx <- rownames(Z) %in% ids
    Z_temp <- Z[idx, ]

    if(length(ids)==1){ medZ <- Z_temp } else { medZ <- apply(Z_temp, 2, median) }

    return(medZ)
}

spccAIC <- function(fit){
    m <- ncol(fit$centers)
    n <- length(fit$cluster)
    k <- nrow(fit$centers)
    D <- fit$tot.withinss
    return(c(
        AIC=(D + 2*m*k),
        BIC=(D + log(n)*m*k)
        )
    )
}

# Plot dendrogram
feature_dendrogram <- function(hc){

    require(ggdendro)
    # hc <- clust_obj
    #d <- as.dendrogram(hc)
    plot_obj <- ggdendrogram(hc, 
        rotate = 270, 
        size = 1,
        theme_dendro = FALSE) + 
    theme_void()

    return(plot_obj)

}

# Gene list
## -- JA Response query
target = c("AT1G52410", "AT5G24780", "AT3G15950", 
            "AT4G31500", "AT4G13770",
            "AT4G39950", "AT2G22330", "AT1G16410", "AT1G16400", 
            "AT3G15950", "AT1G52400", "AT5G24770", "AT4G29260")

# gls_pathway

# Indexing for pathway terms
path.idx <- list(jasmonate = c("JA", "jasmonate", "jaz", "jasm"), 
                 salicylate = c("SA", "salicylate", "salicylic"), 
                 auxin = c("IAA", "indole acetic acid"), 
                 glucosinolate = c("aliphatic glucosinolate", 
                                   "glucosinolate", 
                                   "indole glucosinolate", 
                                   "aromatic glucosinolate"), 
                 ethylene = c("ethylene", "ripening", "senescence"), 
                 flavonoids = c("flavanoids", "flav", "flavanols"),
                 gibberllin = c("GA", "gibberllin"),
                 cytokinin = c("cytokinin"),
                 
                 metabolism = c("metabolism", "metabolites"), 
                 biosynthesis = c("biosynthesis"), 
                 tfs = c("transcription factors", "transcription"), 
                 biogenesis = c("biogenesis"),
                 regulation = c("regulatory, regulation"),
                 aba = c("absiscic acid"),
                 defense = c("defense"),
                 coregulation = c("co-regulation"),
                 promotor = c("promoter"))

# Glucosinolate Pathway

# Transcription regulator

