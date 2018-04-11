### A version of cmapR's GCTX rank function with 
### an option to rank ascending or descending

library(cmapR)

rank.gctx <- function(gctx.path, out.path) {
  gctx <- parse.gctx(gctx.path)
  gctx@mat <- (apply(gctx@mat, 2, function(x) rank(-x)))
  write.gctx(gctx, out.path)
}

gctx.input.path = "/mnt/oncogxA/Alex/l1k/10x_ilincs_sigs_top500.gct"
gctx.output.path ="/mnt/oncogxA/Alex/l1k/10x_ilincs_sigs_top500_custom-rank"
rank.gctx(gctx.input.path, gctx.output.path)
