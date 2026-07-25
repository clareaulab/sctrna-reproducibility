library(tidyverse)
library(ggsankey)

trna_pp <- read.csv("../data/tRNA_genomic_location_PP.csv")

length(unique(trna_pp$name_long)) # 432 == nrow(trna_pp)

# assign individual genes
trna_pp$locus_id = paste0("locusid_", 1:nrow(trna_pp))

# A34->I34 editable codons
edited_anticodons <- c("AGC", "AGG", "AGT", "AGA", "AAG", "AAT", "AAC", "ACG")
trna_pp <- 
  trna_pp %>%
  mutate(AC_read_by_RT = AC) %>%
  dplyr::mutate(is_wobble_editable = AC %in% edited_anticodons) %>%
  mutate(encoded = "encoded") %>%
  dplyr::rename(probepair_id = probePairID)

# Verifying all anticodons in the wobble_editable sequences live between base 33 and 37
# trna_pp %>%
#   filter(is_wobble_editable) %>%
#   mutate(substr_3337 = substr(sequences, 33, 37),
#          is_in_there = str_detect(substr_3337, AC)) %>%
#   dplyr::select(AC, substr_3337, is_in_there, strand) %>%
#   View

trna_pp_wobble_edited = trna_pp %>%
  filter(is_wobble_editable)
  
for (i in 1:nrow(trna_pp_wobble_edited)){
  this_sequence = trna_pp_wobble_edited[[i,"sequences"]]
  this_AC_old = trna_pp_wobble_edited[[i,"AC"]]
  this_AC_read_by_RT_old = trna_pp_wobble_edited[[i,"AC_read_by_RT"]]
  
  this_AC_new = gsub("^A", "I", this_AC_old)
  this_AC_read_by_RT_new = gsub("^A", "G", this_AC_read_by_RT_old)
  
  stringr::str_sub(this_sequence, 33, 37) = gsub(this_AC_read_by_RT_old, this_AC_read_by_RT_new, stringr::str_sub(this_sequence, 33, 37)) 
  
  trna_pp_wobble_edited[[i,"sequences"]] = this_sequence
  trna_pp_wobble_edited[[i,"AC"]] = this_AC_new
  trna_pp_wobble_edited[[i,"AC_read_by_RT"]] = this_AC_read_by_RT_new
  trna_pp_wobble_edited[[i,"encoded"]] = "non-encoded"
}

trna_pp = bind_rows(trna_pp, trna_pp_wobble_edited)

# individual sequences
trna_pp <- trna_pp %>%
  dplyr::mutate(sequence_id = paste0("seqid_", as.character(dplyr::dense_rank(sequences))))
length(unique(trna_pp$sequences)) # 318
length(unique(trna_pp$sequence_id)) # 318

# number 'loci' within each sequence, number one is the 'representative locus', all others represent 'redundand loci'
trna_pp <- trna_pp %>%
  dplyr::group_by(sequence_id) %>%
  dplyr::mutate(sequence_number= dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  mutate(locus_redundance = ifelse(sequence_number == 1, "unique sequence", "redundant loci"))
table(trna_pp$locus_redundance) # 318 'unique sequences', 208 redundant loci

# probe pair id
length(unique(trna_pp$probepair_id)) # 172

# readout id (= probe pair + gapfill anticodon)
trna_pp$readout_id = paste0(trna_pp$probepair_id, "_", trna_pp$AC_read_by_RT)
length(unique(trna_pp$readout_id)) # 217

# readout_id (= probe pair + gapfill anticodon), within readoutgroup numbering
trna_pp = trna_pp %>% 
  group_by(readout_id) %>% 
  dplyr::mutate(readout_number = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  mutate(readout_redundance = ifelse(readout_number == 1, "unique readout", "redundant readouts")) %>%
  mutate(readout_redundance = ifelse(as.character(locus_redundance) == "redundant loci", "redundant loci, redundant readouts", as.character(readout_redundance)))

table(trna_pp$encoded)
table(trna_pp$locus_redundance)
table(trna_pp$readout_redundance)

trna_pp_long <- trna_pp %>%
  ggsankey::make_long(encoded, locus_redundance, readout_redundance) 
node_order = rev(c("encoded", "non-encoded",
                   "unique sequence", "redundant loci",
                   "unique readout", "redundant readouts", "redundant loci, redundant readouts"))
trna_pp_long = trna_pp_long %>%
  mutate(node = factor(node, levels = node_order))
ggplot(trna_pp_long, aes(x = x, 
                            next_x = next_x, 
                            node = node, 
                            next_node = next_node,
                            fill = factor(node),
                            label = node
)) +
  geom_sankey(alpha = 0.8) +
  geom_sankey_text(size = 1) +
  scale_fill_viridis_d() +
  theme_sankey(base_size = 8) +
  theme(legend.position = "none", axis.text.x = element_text(color = "black"))


cowplot::ggsave2("../plots/tRNA-sankey.pdf", height = 2, width = 3)

