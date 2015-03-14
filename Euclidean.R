#Making a Euclidean distance plot for IntDiverse
#mydata <- na.omit(phytoInt)
#mydata <- scale(phytoInt)
Euclid <- cbind(DIp$STATION, phytoInt)
d <- dist(Euclid, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) # display dendogram
groups <- cutree(fit, k=10) # cut tree into 5 clusters

# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=10, border="red")
