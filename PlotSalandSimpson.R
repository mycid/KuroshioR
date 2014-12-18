y <- SurfC$Simpson
x <- SurfC$S-34.2
y1 <- MiddleC$Simpson
x1 <- MiddleC$S-34.2
x2 <- LowC$S-34.2
y2 <- LowC$Simpson
par(mar = rep(2, 4))
plot(x, y, xlab=expression("Sal"), ylab=expression("Simpson"), main="Simpson D. Compared to Sa (-34.2)", xlim=c(-1.2,.5), ylim=c(.2,1), col=1, pch=1) 
points(x1,y1,col=2,pch=2)
points(x2,y2,col=4,pch=5)
legend(.175,1, c("Surface","20-30m","40-50m"), col=c(1,2,4),pch=c(1,2,5))
abline(lm(c(y,y1,y2)~c(x,x1,x2)))
t.test(c(y,y1,y2), c(x,x1,x2))