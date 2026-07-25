library(dplyr)
library(stringr)
library(tidyr)
library(Biostrings)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename

revComp <- function(invec) {
  chartr(old = "AaTtGgCc", new = "TTAACCGG", stringi::stri_reverse(invec))
}

dss <- readDNAStringSet("../data/human_tRNA_dB.txt")
sequences <- unname(as.character(dss))
names <- names(dss)
trna_df <- data.frame(
  sequences,
  str_split_fixed(names, "-", 4)[,c(2,3)]
)

seq_intron <- c("Arg-TCT-1-1",
                "Arg-TCT-2-1",
                "Arg-TCT-3-1",
                "Arg-TCT-3-2",
                "Arg-TCT-5-1",
                "Tyr-GTA-1-1",
                "Tyr-GTA-2-1",
                "Tyr-GTA-3-1",
                "Tyr-GTA-4-1",
                "Tyr-GTA-5-1",
                "Tyr-GTA-5-2",
                "Tyr-GTA-5-3",
                "Tyr-GTA-5-4",
                "Tyr-GTA-5-5",
                "Tyr-GTA-6-1",
                "Tyr-GTA-7-1",
                "Tyr-GTA-8-1",
                "Tyr-GTA-9-1",
                "Ile-TAT-1-1",
                "Ile-TAT-2-1",
                "Ile-TAT-2-2",
                "Ile-TAT-2-3",
                "Ile-TAT-3-1",
                "Leu-CAA-1-1",
                "Leu-CAA-1-2",
                "Leu-CAA-2-1",
                "Leu-CAA-3-1",
                "Leu-CAA-4-1")

trna_df$names <- names(dss)
trna_df <- trna_df %>%
  mutate(name = str_extract(names, "(?<=Homo_sapiens_tRNA-)[^ ]+"))

trna_intron <- trna_df[trna_df$name %in% seq_intron, ]
trna_TAA <- trna_df[trna_df$X2 == "TAA", ]

# define starting position to find codon
starting_position <- sapply( 1:dim(trna_intron)[1], function(idx){ 
  start_vec <- str_locate_all(substr(trna_intron[idx, "sequences"],31, 40), pattern = trna_intron[idx,"X2"])[[1]][,"start"] + 31
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

trna_intron <- trna_intron %>% dplyr::mutate(idx = 1:n())
trna_intron$starting_position <- starting_position

trna_intron$seq_length <- nchar(trna_intron$sequences)

trna_intron_dupl <- trna_intron %>%
  group_by(sequences) %>%
  filter(n() > 1) %>%
  ungroup()

# length of LHS without intron
trna_intron$LHS_ai <- c(rep(36, 10), 45, 45, rep(46, 3), rep(36, 12), 35)

# find intron
split_seq_intron <- function(row) {
  s <- row[["sequences"]]
  starting_position <- row[["starting_position"]]
  LHS_ai <- row[["LHS_ai"]]
  intron <- substr(s, starting_position + 3, nchar(s) - LHS_ai) # +3 to skip codon +1
  return(list(intron = intron))
}

result_list_trna_intron2 <- lapply(seq_len(nrow(trna_intron)), function(i) split_seq_intron(trna_intron[i, ]))
result_df_trna_intron2 <- do.call(rbind, lapply(result_list_trna_intron2, as.data.frame))
trna_exon <- cbind(trna_intron, result_df_trna_intron2)
trna_exon <- trna_exon %>%
  mutate(exon = str_replace(sequences, intron, ""))

# make function to split sequence according to codon in X2
split_seq_e <- function(row){
  s <- row[["exon"]]
  starting_position <- row[["starting_position"]]
  RHS <- substr(s, 1, starting_position - 2)
  LHS_exon <- substr(s, starting_position + 3, nchar(s))
  return(list(RHS = RHS, LHS_exon = LHS_exon))
}

result_list_trna_exon <- lapply(seq_len(nrow(trna_exon)), function(i) split_seq_e(trna_exon[i, ]))
result_df_trna_exon <- do.call(rbind, lapply(result_list_trna_exon, as.data.frame))
trna_exon <- cbind(trna_exon, result_df_trna_exon)

# get 25mers around anticodon
trna_exon$rhs_probe <- sapply(trna_exon$RHS, function(seq) {
  dna_string <- DNAString(str_sub(seq, -25, -1))
  as.character(revComp(dna_string))
})
trna_exon$lhs_probe <- sapply(trna_exon$LHS_exon, function(seq) {
  dna_string <- DNAString(str_sub(seq, 1, 25))
  as.character(revComp(dna_string))
})
trna_exon$gapfill <- sapply(trna_exon$X2, function(seq) {
  dna_string <- DNAString(seq)
  as.character(revComp(dna_string))
})

exon_ref <- trna_exon %>% 
  select(X = idx, name_long = names, name_short = name, AA = X1, codon = X2, RHS, LHS = LHS_exon,
         rhs_probe, lhs_probe, gapfill)
write.csv(exon_ref, file = "../output/exon_reference.csv", row.names = TRUE)
