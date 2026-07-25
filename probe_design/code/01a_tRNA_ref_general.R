library(dplyr)
library(stringr)
library(tidyr)
library(Biostrings)
library(data.table)

revComp <- function(invec) {
  chartr(old = "AaTtGgCc", new = "TTAACCGG", stringi::stri_reverse(invec))
}

dss <- readDNAStringSet("../data/human_tRNA_dB.txt")
sequences <- unname(as.character(dss))
names <- names(dss)
trna_orig <- data.frame(
  sequences,
  str_split_fixed(names, "-", 4)[,c(2,3)]
)

trna_orig$name_long <- names(dss)
trna_orig <- trna_orig %>%
  mutate(name_short = str_extract(name_long, "(?<=Homo_sapiens_tRNA-)[^ ]+"))

# define starting position to find codon
starting_position <- sapply( 1:dim(trna_orig)[1], function(idx){ 
  start_vec <- str_locate_all(substr(trna_orig[idx, "sequences"],31, 40), pattern = trna_orig[idx,"X2"])[[1]][,"start"] + 31
  if(length(start_vec) > 1){
    if(max(start_vec) > 37) { # Catch Tyr / GTA
      min(start_vec)
    } else {
      max(start_vec)
    } 
  }else{
    start_vec[1]
  }
}) %>% unname()
table(starting_position)

trna_orig <- trna_orig %>% mutate(idx = 1:n())
trna_orig$starting_position <- starting_position
trna_orig$seq_length <- nchar(trna_orig$sequences)

# manually change Glu CTC start from 33 to 35
trna_orig <- trna_orig %>%
  mutate(starting_position = ifelse(X1 == "Glu" & X2 == "CTC", 35, starting_position))

# make function to split sequence according to codon in X2
split_seq_i <- function(row) {
  s <- row[["sequences"]]
  starting_position <- row[["starting_position"]]
  RHS <- substr(s, 1, starting_position - 2)
  LHS <- substr(s, starting_position + 2, nchar(s))  # +2 to skip codon
  return(list(RHS = RHS, LHS = LHS))
}

result_list_trna <- lapply(seq_len(nrow(trna_orig)), function(i) split_seq_i(trna_orig[i, ]))
result_df_trna <- do.call(rbind, lapply(result_list_trna, as.data.frame))
trna_orig <- cbind(trna_orig, result_df_trna)

# measure length of initial RHS and LHS
trna_orig$RHS_nchar <- nchar(trna_orig$RHS)
trna_orig$LHS_nchar <- nchar(trna_orig$LHS)

trna_df <- trna_orig
#trna_df <- trna_df %>% rename("aa"="X1") %>% rename("ac"="X2")

# get 25mers around anticodon
trna_df$rhs_probe <- sapply(trna_df$RHS, function(seq) {
  dna_string <- DNAString(str_sub(seq, -25, -1))
  as.character(revComp(dna_string))
})
trna_df$lhs_probe <- sapply(trna_df$LHS, function(seq) {
  dna_string <- DNAString(str_sub(seq, 1, 25))
  as.character(revComp(dna_string))
})
trna_df$gapfill <- sapply(trna_df$X2, function(seq) {
  dna_string <- DNAString(seq)
  as.character(revComp(dna_string))
})

trna_ref <- trna_df %>% 
  select(X = idx, name_long, name_short, AA = X1, codon = X2, RHS, LHS,
         rhs_probe, lhs_probe, gapfill)
write.csv(trna_ref, file = "../output/trna_reference.csv", row.names = TRUE)
