div.abiotic2 <- div.abiotic[-c(186:190),]
filled.contour(interp(div.abiotic2$S, div.abiotic2$Theta, div.abiotic2$ShannonWiener, duplicate="error"), xlab="Sal", ylab="Theta", color=function(x)rev(rainbow(x)), main="Shannon Color Intensity Plot", key.title=title("     Shannon #"))
filled.contour(interp(div.abiotic2$S, div.abiotic2$Theta, div.abiotic2$Simpson, duplicate="error"), xlab="Sal", ylab="Theta", color=function(x)rev(rainbow(x)), main="Simpson Color Intensity Plot", key.title=title("     Simpson #"))
filled.contour(interp(div.abiotic2$S, div.abiotic2$Theta, div.abiotic2$Evenness.SW, duplicate="mean"), xlab="Sal", ylab="Theta", color=function(x)rev(rainbow(x)), main="Evenness SW Color Intensity Plot", key.title=title("        Evenness #"))
filled.contour(interp(div.abiotic2$S, div.abiotic2$Theta, div.abiotic2$Eve, duplicate="mean"), xlab="Sal", ylab="Theta", color=function(x)rev(rainbow(x)), main="Evenness Sim CI Plot", key.title=title("        Evenness #"))
div.abiotic4n <- div.abiotic2[-c(141, 150, 152, 170),]
filled.contour(interp(div.abiotic4n$NO3..uM., div.abiotic4n$sigPoDen, div.abiotic4n$ShannonWiener, duplicate="mean"), xlab="N.uM", ylab="sigPoDen", color=function(x)rev(rainbow(x)), main="Shannon C-Intensity Plot with N", key.title=title("     Shannon #"))
ucol <- colormap(div.abiotic2$ShannonWiener, col = oceColors9B)
imagep(interp(div.abiotic2$S, div.abiotic2$sigPoDen, div.abiotic2$ShannonWiener), colormap = ucol, filledContour=TRUE, xlab="S", ylab="sigPoDen")
persp(interp(div.abiotic2$sigPoDen, div.abiotic2$S, div.abiotic2$Simpson, duplicate="mean"), phi=45, theta=45)

ucol <- colormap(div.abiotic2$Eveness, col = oceColors9B)
filled.contour(interp(div.abiotic4n$S, div.abiotic4n$sigma_t, div.abiotic4n$Eveness, duplicate="mean"), xlab="S", ylab="sigma_T", color=function(x)rev(rainbow(x)), main="Eveness Intensity Plot in S and sigT space", key.title=title("    Eveness #"))


myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(dat$grp)
colScale <- scale_colour_manual(name = "grp",values = myColors)
ggplot(div.abiotic2, aes(x = S, y = Theta, colour=div.abiotic2$ShannonWiener)) + geom_point(alpha=.7, label= div.abiotic2$ShannonWiener)
