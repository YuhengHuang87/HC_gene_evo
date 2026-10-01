suppressMessages({
  library(ggtree); library(treeio); library(ggplot2); library(phangorn); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
no_caption <- "--no-caption" %in% args; args <- args[args != "--no-caption"] # flag (any position): omit the caption below the plot
treefile <- args[1]; outprefix <- args[2]; # optional 3rd arg: circular (default), rectangular, slanted, or clado (circular cladogram)
mode <- if (length(args) > 2) args[3] else "circular"
clado <- mode == "clado"
layout <- if (mode %in% c("rectangular", "slanted")) mode else "circular"
circ <- layout == "circular"
model_txt <- if (length(args) > 3) args[4] else "Q.PLANT+I+G4, MACSE codon-aware protein alignment"

tr <- read.newick(treefile, node.label = "label")
outgroup <- grep("Scaptodrosophila", tr$tip.label, value = TRUE)
tr <- if (length(outgroup) == 1) {
  ape::root(tr, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
} else {
  midpoint(tr, node.labels = "label")
}
tr <- ladderize(tr)

# copy name -> figure label; colors follow cerv_tree_copies.png
copy_lab <- c(
  cerv_X = "cerv (X)", qjt_3L = "qjt (3L)", CG42299_X = "CG42299 (X)", CG42300_X = "CG42300 (X)",
  cerv89B_3R = "cerv-L1 (3R)", cerv100B_3R = "cerv-L2 (3R)",
  cerv36B_2L = "cerv-L3 (2L)", cerv36B = "cerv-L3 (2L)",
  cerv83B_3R = "cerv-L4 (3R)", cerv90F_3R = "cerv-L5 (3R)",
  cerv86C_3R = "cerv-L6 (3R)", cerv62F_3L = "cerv-L7 (3L)",
  cerv21B_2L = "cerv-L8 (2L)", `cerv64E-753L` = "cerv-L9 (3L)",
  cervunk_2L = "cerv-L10 (X)", nerv = "cerv (outgroup)", cerv = "cerv (outgroup)")
pal <- c(
  "cerv (X)" = "#F5A623", "qjt (3L)" = "#4DB6AC", "CG42299 (X)" = "#E6C619",
  "CG42300 (X)" = "#C9B458", "cerv-L1 (3R)" = "#8BC34A", "cerv-L2 (3R)" = "#F08A80",
  "cerv-L3 (2L)" = "#EC7FB0", "cerv-L4 (3R)" = "#5DADE2", "cerv-L5 (3R)" = "#1A9A8A",
  "cerv-L6 (3R)" = "#1F5FA8", "cerv-L7 (3L)" = "#4E9A2E", "cerv-L8 (2L)" = "#E0463C",
  "cerv-L9 (3L)" = "#B5176B", "cerv-L10 (X)" = "#7FCFD0", "cerv (outgroup)" = "#555555")
# 3-letter species codes; species sharing a prefix get distinct codes
sp_code <- function(sp) {
  special <- c(Scaptodrosophila_lebanonensis = "Sleb", pseudoobscura = "pse", pseudoananassae = "psa", pseudotakahashii = "pst",
               subobscura = "sbo", subpulchrella = "spu")
  ifelse(sp %in% names(special), special[sp], substr(sp, 1, 3))
}

tips <- data.frame(label = tr$tip.label) %>%
  mutate(species = sub("__.*", "", label),
         copy = sub(".*__", "", label),
         paralog = factor(copy_lab[copy], levels = names(pal)),
         display = paste(sp_code(species), paralog))
stopifnot(!any(is.na(tips$paralog)), !anyDuplicated(tips$display))

p <- ggtree(tr, layout = layout, size = 0.4, branch.length = if (clado) "none" else "branch.length") %<+% tips

# UFBoot from IQ-TREE "SH-aLRT/UFBoot" node labels
nd <- p$data %>% filter(!isTip, grepl("/", label)) %>%
  mutate(ufb = as.numeric(sub(".*/", "", label))) %>% filter(ufb >= 95)

p <- p +
  geom_point(data = nd, aes(x, y), size = 1, colour = "grey20") +
  geom_tippoint(aes(colour = paralog), shape = 15, size = 2.6) +
  geom_tiplab(aes(label = display),
              size = 2.8, offset = if (clado) 0.3 else 0.02,
              align = circ, linesize = 0.15, linetype = if (circ) "dotted" else NA) +
  scale_colour_manual(values = pal, name = "cerv paralog", drop = TRUE) +
  (if (!clado && circ) geom_treescale(x = 0, y = -3, width = 0.2, fontsize = 2.8, linesize = 0.5, offset = 1.2)) +
  # rect/slanted: manual scale bar at the bottom left (tip 1 is at the bottom; root/outgroup ladderized there)
  (if (!circ) annotate("segment", x = 0, xend = 0.2, y = -2.5, yend = -2.5, linewidth = 0.5)) +
  (if (!circ) annotate("text", x = 0.1, y = -1.3, label = "0.2", size = 2.8)) +
  xlim(0, max(p$data$x) * if (circ) 1.75 else 1.32) +
  guides(colour = guide_legend(override.aes = list(size = 4), ncol = 1)) +
  labs(caption = if (no_caption) NULL else paste0("Maximum-likelihood tree (IQ-TREE; ", model_txt, "); rooted on Scaptodrosophila lebanonensis cerv.\nDots: UFBoot >= 95.\nSpecies codes: first 3 letters, except pse = pseudoobscura, psa = pseudoananassae, pst = pseudotakahashii, sbo = subobscura, spu = subpulchrella.")) +
  theme(legend.position = "right", legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        plot.caption = element_text(size = 8, hjust = 0))

w <- if (circ) 12 else 10; h <- if (circ) 10 else 12
ggsave(paste0(outprefix, ".pdf"), p, width = w, height = h)
ggsave(paste0(outprefix, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
cat("tips:", Ntip(tr), "\n")
