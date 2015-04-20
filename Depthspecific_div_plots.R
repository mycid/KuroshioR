Upperdepth <- Adiv.abiotic2[which(Adiv.abiotic2[, 6]<31),]
Surface <- Adiv.abiotic2[which(Adiv.abiotic2[, 6]<1),]
Middle <- Upperdepth[which(Upperdepth[,6] >0),]
Low <- Adiv.abiotic2[which(Adiv.abiotic2[, 6]>50),]
filled.contour(interp(Surface$S, Surface$sigma_t, Surface$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="sigmaT", color=function(x)rev(rainbow(x)), main=" Surface Shannon CI Plot", key.title=title("     Shannon #"))
filled.contour(interp(Middle$S, Middle$sigma_t, Middle$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="sigmaT", color=function(x)rev(rainbow(x)), main="Middle Simpson CI Plot", key.title=title("     Shannon #"))
filled.contour(interp(Low$sigPoDen, Low$S, Low$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="SigmaT", color=function(x)rev(rainbow(x)), main="Low Eveness CI Plot", key.title=title("     Shannon #"))
