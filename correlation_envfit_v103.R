# SCRIPT: correlation_envfit_v102.R
# VERSION: 1.02
# DATE: 7/8/2016
# AUTHOR: Matthew R. Urschel
# DESCRIPTION: Script to create PCA and envfit plots for species matrix vs
# geochemical variables matrix (or any 2 matrices), calculate correlations
# between significant envfit variables and all other environmental variables, 
# and plot additonal envfit of significantly correlated variables. Uses 'vegan' and
# 'Hmisc' packages.
# ARGUMENTS: spec_matrix = Matrix containing species distribution data for PCA,
#            env_matrix = Matrix containing environmental data for envfit,
#            filenameprefix = Prefix for results file names (optional, default is 'Analysis_pipeline_results_'),
#            min_r_squared = Minimum r-squared value for plotting of correlated environmental variables (optional, default is 0),
#            plot_type = Options for biplot of PCA results (e.g. 'text','points', or 'n') (optional, default is c('text',text')),
#            x_limit = Limits of X axis (optional, default is c(-8,8)),
#            y_limit = Limits of Y axis (optional, default is c(-3,3)),
#            plot_scaling = Scaling of biplot (optional, default is 'sites'),
#            p_max = maximum acceptable p value for primary envfit variable plotting (optional, default is 0.1),
#            primary_envfit_color = plotting color for primary envfit (optional, default is 'blue'),
#            secondary_envfit_color = plotting color for secondary envfit (optional, default is 'darkgoldenrod1'),
#            tertiary_envfit_color = plotting color for tertiary envfit (optional, default is 'darkmagenta'),
#            envfit_permutations = number of permutations for envfits (optional, default is 3000),
#            svg_width = width of svg plot in inches (optional, default = 15),
#            svg_height = height of svg plot in inches (optional, default = 14),
#            svg_point_size = point size for svg plot (optional, default = 12),
#            run_secondary_envfit = run/don't run flag for secondary correlations envfit (optional, default is TRUE),
#            run_tertiary_envfit = run/don't run flag for tertiary correlations envfit (optional, default is TRUE)
# MODIFICATIONS: 1. Added code to check if secondary or tertiary correlation sets are empty before running envfit (v1.01)
#                2. Added code to output primary_envfit_significant_vector (v1.02)
#                3. Changed code for building primary_envfit_significant_vector to include factors with p<=0.1 rather than p<0.1 (v1.02)



correlation_envfit_v102 <- function(spec_matrix,
                               env_matrix,
                               filenameprefix="Correlation_envfit_results",
                               min_r_squared = 0,
                               plot_type = c("text","text"),
                               spec_pca_scale = TRUE,
                               x_limit = c(-8,8),
                               y_limit = c(-3,3),
                               plot_scaling = "sites",
                               p_max = 0.1,
                               primary_envfit_color = "blue",
                               secondary_envfit_color = "darkgoldenrod1",
                               tertiary_envfit_color = "darkmagenta",
                               envfit_permutations = 3000,
                               run_secondary_envfit = TRUE,
                               run_tertiary_envfit = TRUE,
			       run_correlations = TRUE,
                               svg_width = 15,
                               svg_height = 14,
                               svg_point_size = 12) {
 

  if(require("vegan")){
    print("Vegan is loaded correctly.")
  } else {
    print("Trying to install vegan.")
    install.packages("vegan")
    if(require("vegan")){
        print("Vegan installed and loaded.")
    } else {
        stop("Could not install vegan!")
    }
  } 

  if(require("Hmisc")){
    print("Hmisc is loaded correctly.")
  } else {
    print("Trying to install Hmisc.")
    install.packages("Hmisc")
    if(require("Hmisc")){
        print("Hmisc installed and loaded.")
    } else {
        stop("Could not install Hmisc!")
    }
  }

  return_list <- list() 

  # Calculate species distribution PCA
  return_list$spec_matrix_PCA <- rda(spec_matrix,scale=spec_pca_scale)

  # Compute eigenvalue percent variance explained and incorporate into PCA axis labels
  eig_vals_spec <- eigenvals(return_list$spec_matrix_PCA)
  eig_vals_spec_percent <- round((eig_vals_spec/(sum(eig_vals_spec)))*100)
  xlabel_spec = paste("PC1 (",eig_vals_spec_percent[1],"%)",sep="")
  ylabel_spec = paste("PC2 (",eig_vals_spec_percent[2],"%)",sep="")

  # Print PCA
  sink(paste(filenameprefix,"_spec_matrix_PCA.txt"))
  print(return_list$spec_matrix_PCA)
  sink()


  # Calculate primary PCA envfit
  return_list$spec_matrix_PCA_envfit <- envfit(return_list$spec_matrix_PCA,env_matrix,permutations=envfit_permutations)

  # Print PCA envfit
  sink(paste(filenameprefix,"_spec_matrix_PCA_envfit.txt"))
  print(return_list$spec_matrix_PCA_envfit)
  sink()


  primary_envfit_significant_vector <- vector()

  # Create list of significant primary envfit variables
  for (i in 1:NROW(return_list$spec_matrix_PCA_envfit$vectors$pvals)) {
      if (return_list$spec_matrix_PCA_envfit$vectors$pvals[i]<=0.1) {
         primary_envfit_significant_vector <- append(primary_envfit_significant_vector,i)
      }
  }

  return_list$primary_envfit_significant_vector <- primary_envfit_significant_vector

  #Print significant_vector_secondary
   sink(paste(filenameprefix,"primary_envfit_significant_vector.txt"))
   print(return_list$primary_envfit_significant_vector)
   sink()

  #print("return_list$spec_matrix_PCA_envfit$vectors$pvals")
  #print(return_list$spec_matrix_PCA_envfit$vectors$pvals)
  #print("primary_envfit_significant_vector")
  #print(primary_envfit_significant_vector)

  if (run_correlations) {
      # Calculate correlations between environmental variables
      return_list$rcorr_object <- rcorr(as.matrix(env_matrix))
  

      # Write p-values and r-squared values to csv files
      write.csv(return_list$rcorr_object$P, file = paste(filenameprefix,"_env_matrix_correlation_P_values.csv"))
      write.csv(return_list$rcorr_object$r, file = paste(filenameprefix,"_env_matrix_correlation_r_squared_values.csv"))
  } else {
      print("WARNING: Cannot run secondary or tertiary envfit without correlations. Secondary and tertiary envfit will not be run.")
      run_secondary_envfit = FALSE
      run_tertiary_envfit = FALSE
  }

  if (run_secondary_envfit) {

     significant_vector_secondary <- vector()

     #print(NROW(primary_envfit_significant_vector))

     # Create vector of environmental variables significantly correlated with primary envfit significant variables
     for (i in 1:NROW(primary_envfit_significant_vector)) {
         for (j in 1:NCOL(return_list$rcorr_object$P)) {
             if (primary_envfit_significant_vector[i]!=j) {
                if (return_list$rcorr_object$P[primary_envfit_significant_vector[i],j]<0.1) {
                   if (abs(return_list$rcorr_object$r[primary_envfit_significant_vector[i],j])>=min_r_squared) {
                      if (!j %in% primary_envfit_significant_vector) {
                         if (!j %in% significant_vector_secondary) {
                            significant_vector_secondary <- append(significant_vector_secondary,j)
                         }
                      }
                   }
                }
             }
         }
     }

     if (NROW(significant_vector_secondary) > 0) {

        return_list$significant_vector_secondary <- significant_vector_secondary

        #Print significant_vector_secondary
        sink(paste(filenameprefix,"significant_vector_secondary.txt"))
        print(return_list$significant_vector_secondary)
        sink()

        # Calculate secondary envfit (significant correlations)
        return_list$significant_secondary_correlations_envfit <- envfit(return_list$spec_matrix_PCA,env_matrix[,significant_vector_secondary],permutations=envfit_permutations)

        #Print secondary envfit (significant correlations)
        sink(paste(filenameprefix,"_significant_secondary_correlations_envfit.txt"))
        print(return_list$significant_secondary_correlations_envfit)
        sink()
     } else {
       print(paste("WARNING: No significant secondary correlations with minimum r-squared value of ",min_r_squared, "!"))
       run_tertiary_envfit = FALSE
       run_secondary_envfit = FALSE
     }

  }

  if (run_tertiary_envfit & run_secondary_envfit) {

     significant_vector_tertiary <- vector()


     # Create vector of environmental variables significantly correlated with primary envfit significant variables
     for (i in 1:NROW(significant_vector_secondary)) {
         for (j in 1:NCOL(return_list$rcorr_object$P)) {
             if (significant_vector_secondary[i]!=j) {
                if (return_list$rcorr_object$P[significant_vector_secondary[i],j]<0.1) {
                   if (abs(return_list$rcorr_object$r[significant_vector_secondary[i],j])>=min_r_squared) {
                      if (!j %in% primary_envfit_significant_vector) {
                         if (!j %in% significant_vector_secondary) {
                            if (!j %in% significant_vector_tertiary) {
                               significant_vector_tertiary <- append(significant_vector_tertiary,j)
                            }
                         }
                      }
                   }
                }
             }
         }
     }

     if (NROW(significant_vector_tertiary) > 0) {

        return_list$significant_vector_tertiary <- significant_vector_tertiary

        #Print significant_vector_tertiary
        sink(paste(filenameprefix,"_significant_vector_tertiary.txt"))
        print(return_list$significant_vector_tertiary)
        sink()

        # Calculate tertiary envfit (significant correlations)
        return_list$significant_vector_tertiary_correlations_envfit <- envfit(return_list$spec_matrix_PCA,env_matrix[,significant_vector_tertiary],permutations=envfit_permutations)

        #Print tertiary envfit (significant correlations)
        sink(paste(filenameprefix,"_significant_vector_tertiary_correlations_envfit.txt"))
        print(return_list$significant_vector_tertiary_correlations_envfit)
        sink()
     } else {
       print(paste("WARNING: No significant tertiary correlations with minimum r-squared value of ",min_r_squared, "!"))
       run_tertiary_envfit = FALSE
     }
  }


  # Open svg results file
  svg(filename=paste(filenameprefix,"_svg.svg"), 
      width=svg_width, 
      height=svg_height, 
      pointsize=svg_point_size)

  # Plot PCA, with primary, secondary and tertiary envfit results to svg
  biplot(return_list$spec_matrix_PCA,scaling=plot_scaling,xlim=x_limit,ylim=y_limit,type=plot_type,xlab=xlabel_spec,ylab=ylabel_spec)
  plot(return_list$spec_matrix_PCA_envfit,p.max=p_max,col=primary_envfit_color)

  if (run_secondary_envfit) {
     plot(return_list$significant_secondary_correlations_envfit,col=secondary_envfit_color)
  }

  if (run_tertiary_envfit & run_secondary_envfit) {
     plot(return_list$significant_vector_tertiary_correlations_envfit,col=tertiary_envfit_color)
  }
     
  dev.off()

  return(return_list)
 
}