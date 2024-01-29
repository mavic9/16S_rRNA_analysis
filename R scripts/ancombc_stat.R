#by Viktor Mamontov, 2022
pseq <- carbom

pseq_genus <- phyloseq::tax_glom(pseq, taxrank = "Genus")
pseq_genus


out = ancombc(
  phyloseq = pseq_genus, 
  formula = "Kits_type", 
  p_adj_method = "fdr", 
  #zero_cut = 0.90, # by default prevalence filter of 10% is applied
  lib_cut = 0, 
  group = "Kits_type", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)

res <- out$res
head(res)
