#' Provide common color palettes for amino acids
#'
#' @name name of palette to retrieve
#' @keywords color protein amino-acid
#' @export
#' @examples
#' amino_acid_pal('shapely')

#' #' @export
amino_acid_pal <- function(name) {
    clustal <- c(rep("orange", 4),
                 rep("red", 3),
                 rep("blue", 3),
                 rep("green", 4))
    names(clustal) <- c("G", "P", "S", "T",
                        "H", "K", "R",
                        "F", "W", "Y",
                        "I", "L", "M", "V")

    lesk <- c(rep("orange", 4),
              rep("green", 9),
              rep("magenta", 3),
              rep("red", 2),
              rep("blue", 2))
    names(lesk) <- c("G", "A", "S", "T",
                     "C", "V", "I", "L", "P", "Y", "F", "Y", "M", "W",
                     "N", "Q", "H",
                     "D", "E",
                     "K", "R")

    maeditor <- c(rep("lightgreen", 4),
                  rep("green", 9),
                  rep("darkgreen", 3),
                  rep("blue", 2),
                  rep("lilac", 2),
                  "darkblue",
                  rep("orange", 2),
                  "pink",
                  rep("red", 2))
    names(maeditor) <- c("A", "G",
                         "C",
                         "S", "T",
                         "D", "E", "N", "Q",
                         "I", "L", "M", "V",
                         "F", "W", "Y",
                         "H",
                         "K", "R",
                         "P",
                         "S", "T")

    cinema <- c(rep("blue", 2),
                rep("red", 2),
                rep("green", 4),
                rep("white", 5),
                rep("magenta", 3),
                rep("brown", 2),
                "yellow")
    names(cinema) <- c("H", "K", "R",
                       "D", "E",
                       "S", "T", "N", "Q",
                       "A", "V", "L", "I", "M",
                       "F", "W", "Y",
                       "P", "G",
                       "C")

    shapely <- c(rep("E60A0A", 2),
                 rep("E6E600", 2),
                 rep("145AFF", 2),
                 rep("FA9600", 2),
                 rep("3232AA", 2),
                 "00DCDC",
                 rep("0F820F", 3),
                 "C8C8C8",
                 "B45AB4",
                 "8282D2",
                 "DC9682")
    names(shapely) <- c("D", "E",
                        "C", "M",
                        "K", "R",
                        "S", "T",
                        "F", "Y",
                        "N", "Q",
                        "G",
                        "L", "V", "I",
                        "A",
                        "W",
                        "H",
                        "P")

    heliquest <- c("gray", "yellow",
                   "red", "red",
                   "yellow", "gray", "lightblue", "yellow", "blue",
                   "yellow", "yellow",
                   "pink", "green", "pink",
                   "blue", "purple", "purple",
                   "yellow", "yellow", "yellow")
    names(heliquest) <- c("A", "C",
                          "D", "E",
                          "F", "G", "H", "I", "K",
                          "L", "M",
                          "N", "P", "Q",
                          "R", "S", "T",
                          "V", "W", "Y")

    switch(names,
           clustal = clustal,
           lesk = lesk,
           maeditor = maeditor,
           heliquest = heliquest,
           shapely)
}
