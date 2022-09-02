#' Plotting a xy domain with custom boundaries, colour paletts and optional world map 
#'
#' @description Plotting a xy domain with custom boundaries, colour paletts and optional world map 
#' @param x array containing x-axis values (e.g. longitude)
#' @param y array containing y-axis values (e.g. latitude)
#' @param fld field (which should be plotted) with dimensions (x,y)
#' @param levels levels for colour bar
#' @param worldmap should the world map contours be plotted (default TRUE)
#' @param main character containing main title of the plot
#' @param legend_loc location of legend
#' @param legend_title character containing legend title
#' @param legend_only logical TRUE only legend should be pltted, or FALSE everything should be plotted (default)
#' @param Lab lab palette from type colorRampPalette
#' @param ... additional graphic parameters
#' @return no return
#' @export
#' @importFrom graphics .filled.contour legend par plot
#' @importFrom purrr map
#' @importFrom grDevices colorRampPalette
fill_horiz = function(x,y,fld,levels=1:100,main='',worldmap=TRUE,legend_loc='topright', legend_title='',legend_only=FALSE, Lab=NULL,...) {
	n_lev=length(levels)
	if (length(x)!=dim(fld)[1] | length(y)!=dim(fld)[2]) {
		message('Dimensions mismatch.')
	}
	if (is.null(Lab)) {
		Lab=colorRampPalette(c("blue","lightblue","yellow","orange","red","darkred"), space = "Lab")
	}
	if (legend_only==FALSE) {
		plot(NA,xlim=range(x),ylim=range(y),main=main,...)
		.filled.contour(x,y,fld,col=Lab(n_lev),levels=levels)
		if (worldmap==TRUE) {
			map(add=TRUE,col='black',lwd=2.5)
		}
	}
	if (legend_only==TRUE) {
		oldpar=par(no.readonly=TRUE)
		on.exit(par(oldpar))
		par(new=TRUE)
		plot(NA,xlim=range(x),ylim=range(y),main=main,xlab='',ylab='',xaxt='n',yaxt='n',...)
	}
	levels2=levels[seq(2,n_lev-1,length.out=5)]
	colors=character(length(levels2))
	for (i in 1:length(levels2)) {
		j=which(levels==levels2[i])
		colors[i]=Lab(n_lev)[j]
	}
	legend(legend_loc,legend=levels2,fill=colors,bg='white',cex=1.8,title=legend_title)
}



