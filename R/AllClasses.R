setClass("adj2",representation(G1 = "matrix",
                               G2 = "matrix",
                               N1 = "numeric",
                               N2 = "numeric",
                               L1 = "matrix",
                               L2 = "matrix",
                               method = "character"),
         contains=c("matrix","list"))
