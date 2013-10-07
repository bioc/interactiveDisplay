rainbow <- function(object){
    my_mat <- object
    
    melted <- melt(my_mat)
    names(melted) <- c("Var1","Var2","value")
    
    #melted$Var1 <- factor(melted$Var1, rownames(my_mat)[hclust(dist(my_mat))$order])
    #melted$Var2 <- factor(melted$Var2, colnames(my_mat)[hclust(dist(t(my_mat)))$order])
    
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    
    gp <- ggplot(melted, aes(x = Var2, y = Var1, fill = value))
    gp <- gp + geom_tile()
    gp <- gp + coord_fixed()
    gp <- gp + scale_fill_gradientn(colours = myPalette(100))
    gp <- gp + scale_x_discrete(expand = c(0, 0))
    gp <- gp + scale_y_discrete(expand = c(0, 0))
    gp <- gp + coord_equal()
    gp <- gp + theme_bw()
    gp <- gp + theme(axis.text.x=element_text(angle = -45, hjust = 0))
    gp <- gp + xlab("Samples")
    gp <- gp + ylab("Region")
    gp <- gp + ggtitle("Heatmap") + theme(plot.title = element_text(lineheight=.8, face="bold", vjust = 2))
    gp
}
