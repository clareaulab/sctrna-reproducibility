library(stringdist)
library(reshape2)
library(dplyr)
library(readxl)
library(data.table)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename

# helpful function to generate reverse complement
reverse_complement <- function(seq) {
  # define conversions
  complement_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  # split the sequence, find complements, reverse, and combine back to a string
  rev_comp <- rev(complement_map[unlist(strsplit(seq, ""))])
  paste(rev_comp, collapse = "")
}

# get input reference
full_ref <- fread("../output/trna_reference.csv")
full_ref$V1 <- c()

# reference for exons
exon_ref <- fread("../output/exon_reference.csv")
exon_ref$V1 <- c()
exon_ref <- exon_ref %>%
  mutate(name_long = paste0(exon_ref$name_long, " exon")) %>% 
  mutate(name_short = paste0(exon_ref$name_short, "-exon")) #%>%

ref <- rbind(full_ref, exon_ref)

# new full reference containing exon sequences for housekeeping
full_reference <- ref %>% select(name_short, rhs_probe, lhs_probe)

# create pairwise distances LHS probes
lhs_dist_mat <- data.matrix(stringdistmatrix(ref$lhs_probe, method = "lv"))
rownames(lhs_dist_mat) <- ref$name_short; colnames(lhs_dist_mat) <- ref$name_short
melted_lhs_df <- reshape2::melt(lhs_dist_mat) %>% 
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  filter(Var1 >= Var2)

# do same for RHS probes
rhs_dist_mat <- data.matrix(stringdistmatrix(ref$rhs_probe, method = "lv"))
rownames(rhs_dist_mat) <- ref$name_short; colnames(rhs_dist_mat) <- ref$name_short
melted_rhs_df <- reshape2::melt(rhs_dist_mat) %>% 
  mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
  filter(Var1 >= Var2)

ref$designed_probe_lhs <- "missing"
while(any(ref$designed_probe_lhs == "missing")){
  print(sum(ref$designed_probe_lhs == "missing"))
  
  # extract missing sequences
  need <- ref %>%
    filter(designed_probe_lhs == "missing") %>%
    pull(name_short)
  
  # identify pairs of values where both are needed and are within a distance tolerance of each other
  melted_lhs_df_filt <- melted_lhs_df %>% filter(value < 2) %>%
    filter(Var1 %in% need & Var2 %in% need)
  
  # select the probe with the lowest number of occurrences (maximize diversity)
  chosen <- names(sort(table(c(melted_lhs_df_filt$Var1, melted_lhs_df_filt$Var2)))[1])
  aligned_with_chosen <- melted_lhs_df_filt %>%
    filter(Var1 %in% chosen | Var2 %in% chosen) 
  
  # put all of them in a "bag"
  bag <- c(aligned_with_chosen  %>% pull(Var1), 
           aligned_with_chosen  %>% pull(Var2)) %>% unique()
  
  ref$designed_probe_lhs <- ifelse(ref$name_short %in% bag, chosen, ref$designed_probe_lhs)
}

# look at unique probes
unique_probes_lhs <- unique(ref$designed_probe_lhs) 
length(unique_probes_lhs) #145
lhs_dist_mat[unique_probes_lhs, unique_probes_lhs] %>% reshape2::melt() %>%
  filter(Var1 != Var2 ) %>% arrange(value) %>% head()

#----------------------------
# now do the same for the RHS
#----------------------------
ref$designed_probe_rhs <- "missing"
while(any(ref$designed_probe_rhs == "missing")){
  print(sum(ref$designed_probe_rhs == "missing"))
  
  # extract missing sequences
  need <- ref %>%
    filter(designed_probe_rhs == "missing") %>%
    pull(name_short)
  
  # identify pairs of values where both are needed and are within a distance tolerance of each other
  melted_rhs_df_filt <- melted_rhs_df %>% filter(value < 2) %>%
    filter(Var1 %in% need & Var2 %in% need)
  
  # select the probe with the lowest number of occurrences (maximize diversity)
  chosen <- names(sort(table(c(melted_rhs_df_filt$Var1, melted_rhs_df_filt$Var2)))[1])
  aligned_with_chosen <- melted_rhs_df_filt %>%
    filter(Var1 %in% chosen | Var2 %in% chosen) 
  
  # put all of them in a "bag"
  bag <- c(aligned_with_chosen  %>% pull(Var1), 
           aligned_with_chosen  %>% pull(Var2)) %>% unique()
  
  ref$designed_probe_rhs <- ifelse(ref$name_short %in% bag, chosen, ref$designed_probe_rhs)
}

# look at unique probes
unique_probes_rhs <- unique(ref$designed_probe_rhs) 
length(unique_probes_rhs) #89
rhs_dist_mat_old[unique_probes_rhs, unique_probes_rhs] %>% reshape2::melt() %>% ## error out here
  filter(Var1 != Var2 ) %>% arrange(value) %>% head()

# see summary statistics
dim(ref) #460 ; 12
length(unique(paste0(ref$designed_probe_lhs, ref$designed_probe_rhs)))
# 181 after adding exons

length(unique(paste0(ref$lhs_probe, ref$rhs_probe))) # 235 after adding exons

# do some house keeping
lhs_seq <- full_reference$lhs_probe; names(lhs_seq) <- full_reference$name_short
rhs_seq <- full_reference$rhs_probe; names(rhs_seq) <- full_reference$name_short

# simplified design of everything
order_df  <- data.frame(
  ref[,c(2,3,4,5)],
  lhs_probe_id = paste0("lhs_probe_",as.numeric(as.factor(ref$designed_probe_lhs))),
  rhs_probe_id = paste0("rhs_probe_",as.numeric(as.factor(ref$designed_probe_rhs))),
  seq1 = lhs_seq[as.character(ref$designed_probe_lhs)],
  seq2 = rhs_seq[as.character(ref$designed_probe_rhs)]
)


# input for giftwrap
unique_pairs_df <- data.frame(
  order_df[,c(5,6,7,8)]
) %>% unique() %>%
  mutate(probePairID = paste0("PP", 1:n()))

extra_utility <- lapply(1:dim(unique_pairs_df)[1], function(i){
  matches <- order_df %>%
    filter(lhs_probe_id  == unique_pairs_df[i,"lhs_probe_id"] &
             rhs_probe_id  == unique_pairs_df[i,"rhs_probe_id"] ) %>% pull(name_short)
  
  data.frame(
    matches = paste0(matches, collapse = "_"),
    amino_acid = unique(substr(matches, 1,3)),
    expected_anticodons = paste0(unique(substr(matches, 5,7)),collapse = "_")
  )
}) %>% rbindlist() %>% data.frame()
unique_pairs_df <- cbind(unique_pairs_df, extra_utility)

# manual fix for iMet 
unique_pairs_df$amino_acid[unique_pairs_df$amino_acid == "iMe"] <- "iMet"
unique_pairs_df$expected_anticodons[unique_pairs_df$expected_anticodons == "-CA"] <- "CAT"

write.csv(unique_pairs_df, file = "../output/tRNA_probe_design.csv", row.names = TRUE)

# simplified design of everything
order_df <- data.frame(
  ref[,c(2,3,4,5)],
  lhs_probe_id = paste0("hs_trnav1_lhs_probe_",as.numeric(as.factor(ref$designed_probe_lhs))),
  rhs_probe_id = paste0("hs_trnav1_rhs_probe_",as.numeric(as.factor(ref$designed_probe_rhs))),
  seq1 = lhs_seq[as.character(ref$designed_probe_lhs)],
  seq2 = rhs_seq[as.character(ref$designed_probe_rhs)]
)

# assemble probes
Tru_LHS <- c("GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
Tru_RHS <- c("ACGCGGTTAGCACGTANNACTTTAGGCGGTCCTAGCAA")
phos_RHS <- c("/5Phos/")

order_df$RHS_probe_Tru <- paste0(phos_RHS, order_df$seq2, Tru_RHS)
order_df$LHS_probe_Tru <- paste0(Tru_LHS, order_df$seq1)

length(unique(order_df$RHS_probe_Tru))
#89

length(unique(order_df$LHS_probe_Tru))
#145

write.csv(order_df, file = "../output/tRNA_probe_order.csv", row.names = TRUE)



# design additional RHS probes for multiplexing

# BC001: ACTTTAGG
# BC002: AACGGGAA
Tru_RHS_bc2 <- c("ACGCGGTTAGCACGTANNAACGGGAACGGTCCTAGCAA")
# BC003: AGTAGGCT
Tru_RHS_bc3 <- c("ACGCGGTTAGCACGTANNAGTAGGCTCGGTCCTAGCAA")
# BC004: ATGTTGAC
Tru_RHS_bc4 <- c("ACGCGGTTAGCACGTANNATGTTGACCGGTCCTAGCAA")
phos_RHS <- c("/5Phos/")

order_df$RHS_probe_Tru_bc2 <- paste0(phos_RHS, order_df$seq2, Tru_RHS_bc2)
order_df$RHS_probe_Tru_bc3 <- paste0(phos_RHS, order_df$seq2, Tru_RHS_bc3)
order_df$RHS_probe_Tru_bc4 <- paste0(phos_RHS, order_df$seq2, Tru_RHS_bc4)

unique(order_df$RHS_probe_Tru_bc2)
#89
unique(order_df$RHS_probe_Tru_bc3)
#89
unique(order_df$RHS_probe_Tru_bc4)
#89

write.csv(order_df, file = "../output/tRNA_probes_4barcodes.csv", row.names = TRUE)
