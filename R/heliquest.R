#' Draw a helical wheel diagram
#'
#' This function allows you to draw a helical wheel diagram
#' given the sequence and various optional parameters
#' @param helix_seq The sequence for which to draw a helical wheel plot
#' @param Ang The angle by which to rotate the raw plot by
#' @param Mom The hydrophobic moment of the helix
#' @param helix_type The type of helix
#' @param window_size The window over which to calculate Ang and Mom
#' @keywords helical-wheel bioinformatics
#' @examples
#' draw_helical_wheel('ACDEFGHIKLMNPQRSTVWY')
#' @export
draw_helical_wheel <- function(helix_seq, Ang, Mom,
                               helix_type = "alpha",
                               FactC = 0.05, FONT1 = 3, FONT2 = 5,
                               CEXFT = 1, CEXTEXT = 0.8, FlFH = 0, ANGT = 100,
                               NBMIN = 18, NBM2 = 36, NBMAX = 54,
                               circle_size = 3, tail = circle_size, ...) {


    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mfrow = c(1, 1), pty = "s")
    par(mar = c(4, 4, 3, 3))

    # residue color lookup table
    colorlookup <- c(4, 1, 3, 3, 1, 4, 8, 1, 2, 1, 1, 6, 7, 6, 2, 5, 5, 1, 1, 1)
    names(colorlookup) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

    # residues
    data <- strsplit(helix_seq, "")[[1]]

    # colors
    seq <- colorlookup[data]

    # color non-residue data as grey
    seq[is.na(seq)] <- 4

    # circle size
    tail <- rep(circle_size, length(data))

    color <- c("yellow", "blue", "red", "gray", "purple", "pink", "green", "lightblue")
    colorTxt <- c("black", "white", "white", "black", "white", "black", "white", "black")

    # If angle or moment are not specified, use the online server to get them
    if(missing(Ang) | missing(Mom)) {
        # if sequence, is too short, handle differently
        web_params <- get_params(helix_seq,
                                 helix_type = helix_type)
        if(missing(Ang)) Ang <- web_params[1, 'Val_angleM']
        if(missing(Mom)) Mom <- web_params[1, 'Hyd.Moment']
    }

    # For long helicies, limit to 3 turns

    # texte de base modif xlim et ylim a cause des cercles passage de 1.5 a 1.8
    if (length(data) <= NBMIN) {
        Xlimit <- c(-1.2, 1.2)
        Ylimit <- c(-1.2, 1.2)
    } else {
        if ((length(data) > NBMIN) & (length(data) <= NBM2)) {
            Xlimit <- c(-1.5, 1.5)
            Ylimit <- c(-1.5, 1.5)
        } else {
            Xlimit <- c(-1.73, 1.73)
            Ylimit <- c(-1.73, 1.73)
        }
    }

    x <- 1
    y <- 0
    ang <- 0
    Flag <- 0
    j <- 0
    k <- 0
    l <- 0
    # Outline for the first round of helicies
    for (i in 1:length(data)) {
        xsv <- x
        ysv <- y
        if (i > NBMIN) {
            Flag <- 1
        } else {
            Flag <- 0
            # flag pour recaler les positions suivant face hydrophobe
            # flag to reset next hydrophobic face positions
            if (FlFH == 0) {
                newang <- 270 - (Ang * (180/pi))
                x <- cos((ang + newang) * (pi/180))
                y <- sin((ang + newang) * (pi/180))
            } else {
                x <- cos(ang * (pi/180))
                y <- sin(ang * (pi/180))
            }
        }
        if (Flag == 0) {
            # pb de trait au debut quand recalage face
            if ((FlFH == 0) & (i == 1)) {
                xsv <- x
                ysv <- y
                plot(c(xsv, x), c(ysv, y),
                     type = "l", xlim = Xlimit, xlab = "",
                     ylab = "", ylim = Ylimit,
                     xaxt = "n", yaxt = "n",
                     cex = 0.7, lwd = 0.8)
                par(new = T)
            } else {
                plot(c(xsv, x), c(ysv, y),
                     type = "l", xlim = Xlimit, xlab = "",
                     ylab = "", ylim = Ylimit,
                     xaxt = "n", yaxt = "n",
                     cex = 0.7, lwd = 0.8)
                par(new = T)
            }
        }
        ang <- ang - ANGT
    }
    # hydrophobic moment arrow
    if (FlFH == 0) {
        yM <- (1 * Mom) * sin(-pi/2)
        xM <- (1 * Mom) * cos(-pi/2)
    } else {
        yM <- (1 * Mom) * sin(Ang)
        xM <- (1 * Mom) * cos(Ang)
    }

    # Set length/direction of the moment
    if (abs(xM) < 1e-04) {
        xM <- ifelse(xM == 0, 0.001, sign(xM) * 0.001)
    } else if (abs(xM) > 1) {
        xM <- sign(xM) * 0.001
    }

    if (abs(yM) < 1e-04) {
        yM <- ifelse(yM == 0, 0.001, sign(yM) * 0.001)
    } else if (abs(yM) > 1) {
        yM <- sign(yM) * 0.001
    }

    # if ((xM != 0.001) && (xM != -0.001) && (yM != 0.001) && (yM != -0.001)){
    arrows(0, 0, xM, yM, cex = 2, lwd = 2)
    par(new = T)
    # }

    # trace des residus reinitialisation param
    x <- 1
    y <- 0
    ang <- 0
    j <- 0
    k <- 0
    l <- 0
    m <- 0
    n <- 0
    for (i in 1:length(data)) {
        xsv <- x
        ysv <- y
        if ((i > NBMIN) & (i <= NBM2)) {
            if (FlFH == 0) {
                newang <- 270 - (Ang * (180/pi))
                x <- cos((ang + newang) * (pi/180))
                y <- sin((ang + newang) * (pi/180))
            } else {
                x <- cos(ang * (pi/180))
                y <- sin(ang * (pi/180))
            }
            # x<-x*1.27
            x <- x * (1 + (FactC * tail[i - NBMIN] + FactC * tail[i]))
            # y<-y*1.27
            y <- y * (1 + (FactC * tail[i - NBMIN] + FactC * tail[i]))
            if ((i == NBM2) & (NBM2 == 20)) {
                ang <- 0 - ANGT
            }
        } else {
            if (i > NBM2) {
                if (FlFH == 0) {
                  newang <- 270 - (Ang * (180/pi))
                  x <- cos((ang + newang) * (pi/180))
                  y <- sin((ang + newang) * (pi/180))
                } else {
                  x <- cos(ang * (pi/180))
                  y <- sin(ang * (pi/180))
                }
                # x<-x*1.47
                x <- x * (1 + (FactC * tail[i - NBM2] + FactC * tail[i - NBMIN] + FactC * tail[i - NBMIN] + FactC * tail[i]))
                # y<-y*1.47
                y <- y * (1 + (FactC * tail[i - NBM2] + FactC * tail[i - NBMIN] + FactC * tail[i - NBMIN] + FactC * tail[i]))
            } else {
                if (FlFH == 0) {
                  newang <- 270 - (Ang * (180/pi))
                  x <- cos((ang + newang) * (pi/180))
                  y <- sin((ang + newang) * (pi/180))
                } else {
                  x <- cos(ang * (pi/180))
                  y <- sin(ang * (pi/180))
                }
                if ((i == NBMIN) & (NBMIN == 10)) {
                  ang <- 0
                }
            }
        }
        if (i == 1) {
            symbols(x, y,
                    circles = FactC * tail[i],
                    fg = "black", bg = color[seq[i]],
                    xlim = Xlimit, ylim = Ylimit,
                    xlab = "", ylab = "", xaxt = "n",
                yaxt = "n", inches = FALSE)
            par(new = T)
            plot(x, y,
                 xlim = Xlimit, ylim = Ylimit,
                 xlab = "", ylab = "",
                 type = "p",
                 xaxt = "n", yaxt = "n",
                 pch = data[i], cex = tail[i], lwd = 2,
                 col = colorTxt[seq[i]])
            par(new = T)
            Delt <- (FactC * tail[i])/5
            # modif pour mettre N et C a la place de 1 et fin
            plot(x + (sqrt(2)/2) * (FactC * tail[i]) + Delt,
                 y - (sqrt(2)/2) * (FactC * tail[i]) - Delt,
                 xlim = Xlimit, ylim = Ylimit,
                 xlab = "", ylab = "",
                 type = "p",
                 xaxt = "n", yaxt = "n",
                 pch = "N", col = "red",
                 cex = 2, lwd = 1)
        } else {
            if (i <= NBMAX) {
                symbols(x, y,
                        circles = FactC * tail[i],
                        fg = "black", bg = color[seq[i]],
                        xlim = Xlimit, ylim = Ylimit,
                        xlab = "", ylab = "", xaxt = "n",
                        yaxt = "n", inches = FALSE)
                par(new = T)
                plot(x, y,
                     xlim = Xlimit, ylim = Ylimit,
                     xlab = "", ylab = "",
                     type = "p",
                     xaxt = "n", yaxt = "n",
                     pch = data[i], cex = tail[i], lwd = 1,
                  col = colorTxt[seq[i]])
            }
        }
        # modif pour mettre N et C a la place de 1 et fin
        if ((i == length(data)) | (i == NBMAX)) {
            par(new = T)
            plot(x + (sqrt(2)/2) * (FactC * tail[i]) + Delt,
                 y - (sqrt(2)/2) * (FactC * tail[i]) - Delt,
                 xlim = Xlimit, ylim = Ylimit, xlab = "",
                 ylab = "", type = "p",
                 xaxt = "n", yaxt = "n",
                 pch = "C", col = "red", cex = 2, lwd = 1)
        }
        if (i != length(data)) {
            par(new = T)
        }
        ang <- ang - ANGT
    }
    # if ((xM != 0.001) && (xM != -0.001) && (yM != 0.001) && (yM != -0.001)){ par(new=T) arrows(0,0,xM,yM, cex=2, lwd=2) }
    par(new = F)
}

get_params <- function(sequence,
                       helix_type = c("alpha", "3-10", "3-11", "pi"),
                       window_size = c("FULL", "1_TURN", 11, 12, 14, 16, 18, 20, 22, 25, 30, 36)) {
    require(httr)

    # Interpret options to pass to web server
    FHTYPE <- setNames(c("0", 1, 2, 3),
                       c("alpha", "3-10", "3-11", "pi"))[[match.arg(helix_type)]]
    Taille <- match.arg(window_size)

    # For short segments, pad with AA to reach length required
    if(nchar(sequence) < 9 ) {
        warning("Sequence is less than 9 residues. Hydrophobic moment may not be a valid metric.")
        adj_sequence <- paste0("AAAA", sequence, "AAAA")
    } else {
        adj_sequence <- sequence
    }

    # Calculate helical wheel parameters using Heliquest website
    # clef - Not certain what this is for (using value on webpage)

    # PPLOT and FHPLOT specify graphical options:
    # PPLOT - Whether the residue symbol is drawn proportional to its volume
    # FHPLOT - Whether to rotate the helix to vertically align the <µH> vector
    heliquest_request <- POST(
        url = "http://heliquest.ipmc.cnrs.fr/cgi-bin/ComputParamsV3.py",
        body = list(clef = "1",
                    FHTYPE = FHTYPE,
                    Taille = Taille,
                    sequence = adj_sequence,
                    PPLOT = "1",
                    FHPLOT = "0"))

    # Find the URL for the data text file
    page_text <- content(heliquest_request, "text")

    data_url <- regmatches(page_text,
                           regexpr('<a href=\\"(.*)\\s*">\\s*Data.txt<\\/a>',
                                   page_text,
                                   perl = TRUE))

    data_url <- regmatches(data_url,
                           regexpr('tmp.*\\/.*\\/Data.txt',
                                   data_url,
                                   perl = TRUE))

    # Read and return data file
    params_df <- read.delim(paste0("http://heliquest.ipmc.cnrs.fr/", data_url))
    if (nchar(sequence) < 9) {
        params_df$Val_angleM <- params_df$Val_angleM + 20 * 5
    }
    params_df
}
