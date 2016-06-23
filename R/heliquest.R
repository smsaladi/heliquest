#' Draw a helical wheel diagram
#'
#' This function allows you to draw a helical wheel diagram
#' given the sequence and various optional parameters
#' @param helix_seq Sequence draw the helical wheel plot for
#' @keywords helical-wheel bioinformatics
#' @export
#' @examples
#' draw_helical_wheel('ACDEFGHIKLMNPQRSTVWY')

#' @export
draw_helical_wheel <- function(helix_seq, FactC = 0.05, FONT1 = 3, FONT2 = 5,
                               CEXFT = 1, CEXTEXT = 0.8, Ang = 3.83972435439,
                               Mom = 0.046225448735, FlFH = 0, ANGT = 100,
                               NBMIN = 18, NBM2 = 36, NBMAX = 54,
                               circle_size = 3, tail = circle_size, ...) {
    par(mfrow = c(1, 1), pty = "s")
    par(mar = c(4, 4, 3, 3))

    # residue color lookup table
    colorlookup <- c(4, 1, 3, 3, 1, 4, 8, 1, 2, 1, 1, 6, 7, 6, 2, 5, 5, 1, 1, 1)
    names(colorlookup) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

    # residues
    data <- strsplit(helix_seq, "")[[1]]

    # colors
    seq <- colorlookup[data]

    # circle size
    tail <- rep(circle_size, length(data))

    # FOR several helices VM limite a 3 tours
    color <- c("yellow", "blue", "red", "gray", "purple", "pink", "green", "lightblue")
    colorTxt <- c("black", "white", "white", "black", "white", "black", "white", "black")

    # texte de base modif xlim et ylim a cause des cercles passage de 1.5 a 1.8
    if (length(data) <= NBMIN) {
        Xlimit = c(-1.2, 1.2)
        Ylimit = c(-1.2, 1.2)
    } else {
        if ((length(data) > NBMIN) & (length(data) <= NBM2)) {
            Xlimit = c(-1.5, 1.5)
            Ylimit = c(-1.5, 1.5)
        } else {
            Xlimit = c(-1.73, 1.73)
            Ylimit = c(-1.73, 1.73)
        }
    }

    x = 1
    y = 0
    ang = 0
    Flag = 0
    j = 0
    k = 0
    l = 0
    # trace du premier tour d'helice
    for (i in 1:length(data)) {
        xsv = x
        ysv = y
        if (i > NBMIN) {
            Flag = 1
        } else {
            Flag = 0
            # flag pour recaler les positions suivant face hydrophobe
            if (FlFH == 0) {
                newang = 270 - (Ang * (180/pi))
                x = cos((ang + newang) * (pi/180))
                y = sin((ang + newang) * (pi/180))
            } else {
                x = cos(ang * (pi/180))
                y = sin(ang * (pi/180))
            }
        }
        if (Flag == 0) {
            # pb de trait au debut quand recalage face
            if ((FlFH == 0) & (i == 1)) {
                xsv = x
                ysv = y
                plot(c(xsv, x), c(ysv, y), type = "l", xlim = Xlimit, xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = Ylimit, cex = 0.7, lwd = 0.8)
                par(new = T)
            } else {
                plot(c(xsv, x), c(ysv, y), type = "l", xlim = Xlimit, xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = Ylimit, cex = 0.7, lwd = 0.8)
                par(new = T)
            }
        }
        ang = ang - ANGT
    }
    # plot du moment
    if (FlFH == 0) {
        yM = (1 * Mom) * sin(-pi/2)
        xM = (1 * Mom) * cos(-pi/2)
    } else {
        yM = (1 * Mom) * sin(Ang)
        xM = (1 * Mom) * cos(Ang)
    }
    if (xM > 1) {
        xM = 1
    }
    if (xM < -1) {
        xM = -1
    }
    if (yM > 1) {
        yM = 1
    }
    if (yM < -1) {
        yM = -1
    }
    if (yM >= 0 && yM < 1e-04) {
        yM = 0.001
    }
    if (yM < 0 && yM > -1e-04) {
        yM = -0.001
    }
    if (xM >= 0 && xM < 1e-04) {
        xM = 0.001
    }
    if (xM < 0 && xM > -1e-04) {
        xM = -0.001
    }
    # if ((xM != 0.001) && (xM != -0.001) && (yM != 0.001) && (yM != -0.001)){
    arrows(0, 0, xM, yM, cex = 2, lwd = 2)
    par(new = T)
    # }

    # trace des residus reinitialisation param
    x = 1
    y = 0
    ang = 0
    j = 0
    k = 0
    l = 0
    m = 0
    n = 0
    for (i in 1:length(data)) {
        xsv = x
        ysv = y
        if ((i > NBMIN) & (i <= NBM2)) {
            if (FlFH == 0) {
                newang = 270 - (Ang * (180/pi))
                x = cos((ang + newang) * (pi/180))
                y = sin((ang + newang) * (pi/180))
            } else {
                x = cos(ang * (pi/180))
                y = sin(ang * (pi/180))
            }
            # x=x*1.27
            x = x * (1 + (FactC * tail[i - NBMIN] + FactC * tail[i]))
            # y=y*1.27
            y = y * (1 + (FactC * tail[i - NBMIN] + FactC * tail[i]))
            if ((i == NBM2) & (NBM2 == 20)) {
                ang = 0 - ANGT
            }
        } else {
            if (i > NBM2) {
                if (FlFH == 0) {
                  newang = 270 - (Ang * (180/pi))
                  x = cos((ang + newang) * (pi/180))
                  y = sin((ang + newang) * (pi/180))
                } else {
                  x = cos(ang * (pi/180))
                  y = sin(ang * (pi/180))
                }
                # x=x*1.47
                x = x * (1 + (FactC * tail[i - NBM2] + FactC * tail[i - NBMIN] + FactC * tail[i - NBMIN] + FactC * tail[i]))
                # y=y*1.47
                y = y * (1 + (FactC * tail[i - NBM2] + FactC * tail[i - NBMIN] + FactC * tail[i - NBMIN] + FactC * tail[i]))
            } else {
                if (FlFH == 0) {
                  newang = 270 - (Ang * (180/pi))
                  x = cos((ang + newang) * (pi/180))
                  y = sin((ang + newang) * (pi/180))
                } else {
                  x = cos(ang * (pi/180))
                  y = sin(ang * (pi/180))
                }
                if ((i == NBMIN) & (NBMIN == 10)) {
                  ang = 0
                }
            }
        }
        if (i == 1) {
            symbols(x, y, circles = FactC * tail[i], fg = "black", bg = color[seq[i]], xlim = Xlimit, ylim = Ylimit, xlab = "", ylab = "", xaxt = "n",
                yaxt = "n", inches = FALSE)
            par(new = T)
            plot(x, y, xlim = Xlimit, ylim = Ylimit, xlab = "", ylab = "", type = "p", xaxt = "n", yaxt = "n", pch = data[i], cex = tail[i], lwd = 2,
                col = colorTxt[seq[i]])
            par(new = T)
            Delt = (FactC * tail[i])/5
            # modif pour mettre N et C a la place de 1 et fin
            plot(x + (sqrt(2)/2) * (FactC * tail[i]) + Delt, y - (sqrt(2)/2) * (FactC * tail[i]) - Delt, xlim = Xlimit, ylim = Ylimit, xlab = "",
                ylab = "", type = "p", xaxt = "n", yaxt = "n", pch = "N", col = "red", cex = 2, lwd = 1)
        } else {
            if (i <= NBMAX) {
                symbols(x, y, circles = FactC * tail[i], fg = "black", bg = color[seq[i]], xlim = Xlimit, ylim = Ylimit, xlab = "", ylab = "", xaxt = "n",
                  yaxt = "n", inches = FALSE)
                par(new = T)
                plot(x, y, xlim = Xlimit, ylim = Ylimit, xlab = "", ylab = "", type = "p", xaxt = "n", yaxt = "n", pch = data[i], cex = tail[i], lwd = 1,
                  col = colorTxt[seq[i]])
            }
        }
        # modif pour mettre N et C a la place de 1 et fin
        if ((i == length(data)) | (i == NBMAX)) {
            par(new = T)
            plot(x + (sqrt(2)/2) * (FactC * tail[i]) + Delt, y - (sqrt(2)/2) * (FactC * tail[i]) - Delt, xlim = Xlimit, ylim = Ylimit, xlab = "",
                ylab = "", type = "p", xaxt = "n", yaxt = "n", pch = "C", col = "red", cex = 2, lwd = 1)
        }
        if (i != length(data)) {
            par(new = T)
        }
        ang = ang - ANGT
    }
    # if ((xM != 0.001) && (xM != -0.001) && (yM != 0.001) && (yM != -0.001)){ par(new=T) arrows(0,0,xM,yM, cex=2, lwd=2) }
    par(new = F)
}
