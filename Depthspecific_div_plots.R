souUpperdepth <- div.abiotic4n[which(div.abiotic4n[, 3]<51),]
Surface <- div.abio2[which(div.abio2[, 3]<1),] 
Middle <- Upperdepth[which(Upperdepth[,3] >0),]
Low <- Upperdepth[which(Upperdepth[, 3]>31),]
filled.contour(interp(Surface$S, Surface$sigma_t, Surface$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="sigmaT", color=function(x)rev(rainbow(x)), main=" Surface Shannon CI Plot", key.title=title("     Shannon #"))
filled.contour(interp(Middle$S, Middle$sigma_t, Middle$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="sigmaT", color=function(x)rev(rainbow(x)), main="Middle Simpson CI Plot", key.title=title("     Shannon #"))
filled.contour(interp(Low$sigPoDen, Low$S, Low$ShannonWiener, duplicate="mean"), xlab="Sal", ylab="SigmaT", color=function(x)rev(rainbow(x)), main="Low Eveness CI Plot", key.title=title("     Shannon #"))
