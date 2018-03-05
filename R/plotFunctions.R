globalVariables(c("ggplot", "aes", "Cutoff", "Counts", "geom_line", "Methods", "xlab", "ylab", "theme_bw", "theme",
				"element_blank", "element_text", "unit", "scale_x_continuous", "scale_y_continuous", "guides", 
				"guide_legend", "Value", "..count..", "geom_freqpoly", "%>%", "set", "discrete_scale"));
change_colours <- function(p, palette, type) { 
	n <- nlevels(p$data[[deparse(p$mapping$group)]])
	tryCatch(as.character(palette), error=function(e) stop('palette should be a vector of colours', call.=FALSE))
	if(n > length(palette)) stop('Not enough colours in palette.')
	if(missing(type)) type <- grep('colour|fill', names(p$layers[[1]]$mapping), value=TRUE)[1]
	pal <- function(n) palette[seq_len(n)]
	p + discrete_scale(type, 'foo', pal)
}

plotHC <- function(data, method=c("spearman", "pearson", "kendall"), mar=c(9,1,0,20)){
	if(!is.data.frame(data)) data <- data.frame(data);
	method <- match.arg(method);
	hc <- hclust(as.dist(1-cor(data, method=method)));
	dend <- as.dendrogram(hc);
	dend <- dend %>% set("labels_cex", 6.5) %>% set("branches_lwd", 6.5);
	par(mar=mar, mgp=c(10,5,0), cex.axis=6);
	plot(dend, horiz = TRUE);
	axis(side=1,lwd=8);
}

plotCVs <- function(data, methods=c("None", "HG7", "ERCC", "TN", "TC", "CR", "NR", "DESeq", "UQ", "TMM", "TU"), legend.position=c(.85, .48)){
	if(!is.data.frame(data)) data <- data.frame(data);
	if(is.factor(data$Cutoff)) data$Cutoff <- as.numeric(as.character(data$Cutoff));
	if(is.factor(data$Counts)) data$Counts <- as.numeric(as.character(data$Counts));
	data$Methods <- factor(data$Methods, levels=methods, labels=methods);
	change_colours(ggplot(data = data, aes(x=Cutoff, y=Counts)) + geom_line(aes(group=Methods, color=Methods),size=3)
				+ xlab("Normalized CV cutoff") + ylab("Number of uniform genes")+theme_bw()
				+theme(panel.grid.minor = element_blank(),axis.title.x = element_text(size = 48),axis.title.y = element_text(size =48),
				axis.text.x = element_text(size = 38), axis.text.y = element_text(size = 38),
				legend.text = element_text(size=39),legend.title = element_text(size=43),
				legend.position=legend.position,legend.background = element_blank(),legend.key = element_blank(),
				legend.key.height=unit(1.8,"cm"),plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
				+scale_x_continuous(breaks=seq(0,1,0.2))+ scale_y_continuous()
				+ guides(color=guide_legend(title=NULL)),
				c("olivedrab","blue", "red",  "violet", "orange", "yellow", "magenta", "peru", "black", "maroon", "lightblue", "darkslateblue", 
				"seashell4", "tan2", "darkgreen", "springgreen"));
}

plotCors <- function(data, methods=c("None", "HG7", "ERCC", "TN", "TC", "CR", "NR", "DESeq", "UQ", "TMM", "TU"), legend.position=c(.15, .56)){
	if(!is.data.frame(data)) data <- data.frame(data);
	if(is.factor(data$Value)) data$Value <- as.numeric(as.character(data$Value));
	data$Methods <- factor(data$Methods, levels=methods, labels=methods);
	change_colours(ggplot(data = data,aes(x=Value, y=..count../sum(..count..))) + geom_freqpoly(aes(group=Methods, color=Methods),size=3,bins=50) 
				+ xlab("Spearman correlation") + ylab("Fraction of gene pairs")+theme_bw()
				+theme(panel.grid.minor = element_blank(),axis.title.x = element_text(size = 48),axis.title.y = element_text(size =48),
				axis.text.x = element_text(size = 38),axis.text.y = element_text(size = 38),legend.text = element_text(size=39),legend.title = element_text(size=43),
				legend.position=legend.position, legend.background = element_blank(), legend.key = element_blank(),
				legend.key.height=unit(1.8,"cm"),plot.margin=unit(c(0.5,1,0.5,0.5),"cm"))
				+scale_x_continuous(expand=c(0.01,0.01),breaks=round(seq(-1,1,0.25),2))+ scale_y_continuous(expand=c(0.01,0))
				+ guides(color=guide_legend(title=NULL)),
				c("olivedrab","blue", "red",  "violet", "orange", "yellow", "magenta", "peru", "black", "maroon", "lightblue", "darkslateblue", 
				"seashell4", "tan2", "darkgreen", "springgreen"));
}
