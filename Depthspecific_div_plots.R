Upperdepth <- Adiv.abiotic2[which(Adiv.abiotic2[, 7]<51),]
SaltyShal <- Upperdepth[which(Upperdepth[, 17 ]>34.0),]
SaltyShal <- SaltyShal[which(SaltyShal[, 17] < 34.4),]
Surface <- Adiv.abiotic2[which(Adiv.abiotic2[, 7]<1),]
Middle <- Upperdepth[which(Upperdepth[,7] >0),]
Low <- Adiv.abiotic2[which(Adiv.abiotic2[, 7]>50),]
filled.contour(interp(A$latitude, A$depth, A$sigPoDen, duplicate="mean"), xlab="latitude", ylab="depth", color=function(x)rev(rainbow(x)), main="Transect A", key.title=title("Potential Density"))
filled.contour(interp(Middle$S, Middle$sigma_t, Middle$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="sigmaT", color=function(x)rev(rainbow(x)), main="Middle Simpson CI Plot", key.title=title("     Shannon #"))
filled.contour(interp(Low$sigPoDen, Low$S, Low$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="SigmaT", color=function(x)rev(rainbow(x)), main="Low Eveness CI Plot", key.title=title("     Shannon #"))
