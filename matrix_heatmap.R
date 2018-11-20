# Script to create heatmap from data matrix.
# Packages used: gplots, RColorBrewer
# Arguments: matrix_to_plot = matrix from which to create heatmap,
#           x_labels = Array containing labels for first set of rays (Default = NULL),
#           y_labels = Array containing labels for second set of rays (Default = NULL),
#	    filenameprefix = Name to be assigned to output svg file (Default = "matrix_heatmap_results")
#           main_label_setting = Main label for heatmap figure (Default = "Matrix heat map")
#           margins_setting = Size of figure margins (Default = c(12,9))
#           dendogram_setting = Draw dendogram for "both","row","column", or "none" (Default = "both")
#           colors_setting = Colors for heatmap (Default = c("yellow","orange","red"))
#           png_width = Pixel width of png output file (Default = 1500)
#           png_height = Pixel height of png output file (Default = 1500)
#           png_point_size = Point size for characters on png output file (Default = 8)
#           png_resolution = Resolution of png output file (Default = 300 dpi
#           cluster_columns_setting = Cluster column variables (TRUE/FALSE)? (Default = TRUE)
#           cluster_rows_setting = Cluster rows variables (TRUE/FALSE)? (Default = TRUE)
#           density_info_settings = Draw density info in key ("histogram","density","none") (Default = "none")
#           trace_settings = Trace lines at "column","row","both", or "none" level (Default = "none"),
#           lhei_setting = Column height (Default = NULL),
#           lwid_setting = Column width (Default = NULL)

matrix_heatmap <- function(matrix_to_plot,
			  x_labels = NULL,
			  y_labels = NULL,
			  main_label_setting = "Matrix heat map",
                          margins_setting = c(12,9),
                          dendogram_setting = "both",
                          colors_setting = c("yellow","orange","red"),
			  png_width = 1500,
			  png_height = 1500,
			  png_point_size = 8,
                          png_resolution = 300,
                          cluster_columns_setting=TRUE, #Cluster columns by default
                          cluster_rows_setting=TRUE, #Cluster rows by default
                          density_info_settings="none",
                          trace_settings="none",
                          filenameprefix="matrix_heatmap_results",
                          cexCol_setting=1,
                          cexRow_setting=1,
                          lhei_setting=NULL,
                          lwid_setting=NULL) {

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
}

#If no row labels passed in, default to row names of X arrays
if (is.null(x_labels)) {
   array1_labels = row.names(matrix_to_plot)
}

if (is.null(y_labels)) {
   array2_labels = colnames(matrix_to_plot)
}

return_list <- list()

return_list$matrix <- matrix_to_plot


# creates a own color palette from red to green
my_palette <- colorRampPalette(colors_setting)(n = 299)


# creates a 5 x 5 inch svg file
png(paste(filenameprefix,"_png.png"),    # create PNG for the heat map        
  width = png_width,        # height and length of PNG file
  height = png_height,
  res = png_resolution,
  pointsize = png_point_size)        #font size

return_list$heatmap_object <- heatmap.2(return_list$matrix,
  		       main = main_label_setting, # heat map title
  		       density.info=density_info_settings,  # turns on/off density plot inside color legend
  		       trace=trace_settings,         # turns off/on trace lines inside the heat map
  		       margins=margins_setting,     # widens margins around plot
  		       col=my_palette,       # use on color palette defined earlier
  		       dendrogram=dendogram_setting,     # only draw a row dendrogram
  		       Colv=cluster_columns_setting,		# turn on/off column clustering
                       Rowv=cluster_rows_setting,
  		       labRow = array1_labels, # Set row names
  		       labCol = array2_labels, # Set column names 
                       cexCol = cexCol_setting,
                       cexRow = cexRow_setting,
                       lhei = lhei_setting,
                       lwid = lwid_setting)    
                                

dev.off()               # close the PNG device

return(return_list)

}