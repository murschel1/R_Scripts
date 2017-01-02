# Script to create PCA, distance/dissimilarity, envfit and CCA plots
# for species matrix vs geochemical variables matrix (or any 2 matrices).
# Arguments: spec_matrix = matrix with species distribution data,
#            env_matrix = matrix with environmental measurement (geochemistry) data,
#            filenameprefix = name of pdf file to save results into (optional, default is 'Analysis_pipeline_results'),
#            env_pca_scale = whether to scale environmental pca results (true/false) (optional, default is TRUE),
#            spec_pca_scale = whether to scale species pca results (true/false) (optional, default is TRUE),
#            spec_matrix_pca_name = name of species matrix PCA to be used as plot title (optional, default is 'Species matrix PCA'),
#            env_matrix_pca_name = name of environmental matrix PCA to be used as plot title (optional, default is 'Environmental matrix PCA'),
#            envfit_name = name of envfit correlation of environmental data to species matrix PCA to be used as plot title (optional, default is 'Envfit environmental data to species PCA'),
#            envfit_permutations = number of permutations for envfit function (optional, default is 3000),
#            environmental_matrix_dendogram_name  = name of environmental matrix euclidean distance dendogram to be used as plot title (optional, default is 'Environmental matrix euclidean distance dendogram'),
#            species_distribution_matrix_dendogram_name = name of Bray-Curtis species dissimilarity dendogram to be used as plot title (optional, default is 'Bray-Curtis species dissimilarity dendogram'),
#            environmental_matrix_dissimilarity_method = dissimilarity analysis method for environmental data (optional, default is "euclidean"),
#            species_distribution_dissimilarity_method = dissimilarity analysis method for species distribution data (optional, default is "euclidean"),
#            presence_absence = whether species distribution data is binary presence/absence data (true/false) (optional, default is TRUE),
#            CCA_name = name of CCA to be used as plot title (optional, default is 'CCA significant environmental parameters vs species'),
#            acceptable_p_value = Highest acceptable p-value resulting from envfit for environmental factors (optional, default is 0.05),
#            spec_matrix_species_plot_type = Plot type for species matrix PCA 'species' data (optional, default is "text"),
#            spec_matrix_site_plot_type = Plot type for species matrix PCA 'sites' data (optional, default is "text"),
#            env_matrix_species_plot_type = Plot type for environmental matrix PCA 'species' data (optional, default is "text"),
#            env_matrix_site_plot_type = Plot type for environmental matrix PCA 'sites' data (optional, default is "text"),
#            spec_matrix_species_color = Color for environmental matrix PCA 'species' data (optional, default is "red"),
#            spec_matrix_site_color = Color for species matrix PCA 'sites' data (optional, default is "black"),
#            env_matrix_species_color = Color for environmental matrix PCA 'species' data (optional, default is "red"),
#            env_matrix_site_color = Color for environmental matrix PCA 'sites' data (optional, default is "black")


AnalysisPipeline_rda_v3_1 <- function(spec_matrix,
                             env_matrix,
                             filenameprefix='Analysis_pipeline_results',
                             env_pca_scale=TRUE,
                             spec_pca_scale=TRUE,
                             spec_matrix_pca_name='Species distribution PCA',
                             env_matrix_pca_name='Environmental factor PCA',
                             envfit_name='Envfit environmental factor data to species distribution PCA',
                             envfit_permutations=3000,
                             environmental_matrix_dendogram_name ='Environmental matrix dissimilarity dendogram',
                             environmental_matrix_dissimilarity_method="euclidean",
                             species_distribution_dendogram_name='Species distribution dissimilarity dendogram',
                             species_distribution_dissimilarity_method="jaccard",
                             presence_absence=TRUE,
                             CCA_name='CCA significant environmental parameters vs species',
                             RDA_name='RDA significant environmental parameters vs species',
                             acceptable_p_value=0.05,
                             run_env_pca=TRUE,
                             run_spec_pca=TRUE,
                             run_env_dist=TRUE,
                             run_spec_dist=TRUE,
                             run_envfit=TRUE,
                             run_cca_rda=TRUE,
                             spec_matrix_species_plot_type="text",
                             spec_matrix_site_plot_type="text",
                             env_matrix_species_plot_type="text",
                             env_matrix_site_plot_type="text",
                             spec_matrix_species_color="red",
                             spec_matrix_site_color="black",
                             env_matrix_species_color="red",
                             env_matrix_site_color="black") {
 

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

  return_list <- list() 

  ResizeVector_v3 <- function(xarray, 
                            yarray,
                            array_labels, 
                            percent_of_previous=25, 
                            increase=FALSE) {


     pointlist <- list()
     pointlist$x <- xarray
     pointlist$y <- yarray
     pointlist$angles <- array()
     pointlist$angle_labels <- array()
     temp_angles <- array()
     temp_angle_labels <- array()

     percent_of_previous = percent_of_previous / 100

   
     for (i in 1:NROW(pointlist$x)) {
     

        old_hyp = sqrt((pointlist$x[i]^2)+(pointlist$y[i]^2))

        ang=atan2(pointlist$y[i],pointlist$x[i])

        if (increase) {
           new_hyp = old_hyp + (old_hyp * percent_of_previous)
        } else {
           new_hyp = old_hyp - (old_hyp * percent_of_previous)
        }

        x2new = cos(ang) * new_hyp
        y2new = sin(ang) * new_hyp

        pointlist$x[i] <- x2new
        pointlist$y[i] <- y2new

        if (ang<=0) {
           temp_angles[i] <- round(as.double(abs(ang*57.29)+90))
        } else {
         
           temp_angles[i] <- round(as.double(270+(180-(ang*57.29))))

           if (temp_angles[i]>=360) {
              temp_angles[i] <- round(as.double(temp_angles[i])-360)
           }
        
        }

        temp_angle_labels[i] <- array_labels[i]
        

     }

     ang_counter = 1
     while (length(temp_angles)>0) {
        pointlist$angles[ang_counter] <- min(temp_angles)
        pointlist$angle_labels[ang_counter] <- temp_angle_labels[match(min(temp_angles),temp_angles)]
        temp_angle_labels <- temp_angle_labels[-1*match(min(temp_angles),temp_angles)]
        temp_angles <- temp_angles[-1*match(min(temp_angles),temp_angles)]
        ang_counter = ang_counter + 1
     }

     return(pointlist)

  }


  #Set pdf file to save plots in
  pdf(file=paste(filenameprefix,'.pdf')) 

  if (run_spec_pca) {
     #Compute PCA species matrix
     spec_matrix_PCA <- rda(spec_matrix,scale=spec_pca_scale)

     return_list$spec_matrix_PCA <- spec_matrix_PCA

     #Compute eigenvalue percent variance explained and incorporate into axis labels
     eig_vals_spec <- eigenvals(spec_matrix_PCA)
     eig_vals_spec_percent <- round((eig_vals_spec/(sum(eig_vals_spec)))*100)
     xlabel_spec = paste("PC1 (",eig_vals_spec_percent[1],"%)",sep="")
     ylabel_spec = paste("PC2 (",eig_vals_spec_percent[2],"%)",sep="")

     #Plot PCA species matrix
     biplot(spec_matrix_PCA, main=spec_matrix_pca_name, scaling="sites", display = c("species","sites"), type = c(spec_matrix_species_plot_type,spec_matrix_site_plot_type),xlab=xlabel_spec,ylab=ylabel_spec,col = c(spec_matrix_site_color,spec_matrix_species_color))

     #Print eigenvalues for species distribution PCA to results file
     sink(paste(filenameprefix,"_spec_PCA_PC_percent_variance.txt"))
     print(eig_vals_spec_percent)
     sink()

     return_list$eig_vals_spec_percent <- eig_vals_spec_percent
  } else {
     if (run_envfit) {
        print("WARNING: Cannot run envfit without species matrix PCA. Envfit will not be run.")
        run_envfit=FALSE
     }
  }

  if (run_env_pca) {
     #Compute PCA environmental matrix
     env_matrix_PCA <- rda(env_matrix,scale=env_pca_scale)

     return_list$env_matrix_PCA <- env_matrix_PCA

     #Compute eigenvalue percent variance explained and incorporate into axis labels
     eig_vals_env <- eigenvals(env_matrix_PCA)
     eig_vals_env_percent <- round((eig_vals_env/(sum(eig_vals_env)))*100)
     xlabel_env = paste("PC1 (",eig_vals_env_percent[1],"%)",sep="")
     ylabel_env = paste("PC2 (",eig_vals_env_percent[2],"%)",sep="")

     #Plot PCA environmental matrix
     biplot(env_matrix_PCA, main=env_matrix_pca_name, scaling="sites", display = c("species","sites"), type = c(env_matrix_species_plot_type,env_matrix_site_plot_type),xlab=xlabel_env,ylab=ylabel_env,col = c(env_matrix_site_color,env_matrix_species_color))

     #Print eigenvalues for environmental matrix PCA to results file
     sink(paste(filenameprefix,"_env_PCA_PC_percent_variance.txt"))
     print(eig_vals_env_percent)
     sink()

     return_list$eig_vals_env_percent <- eig_vals_env_percent
  } 

  
  if (run_env_dist) {
     #Calculate dissimilarity matrix for environmental data
     Env_dissim_matrix <- vegdist(env_matrix,method=environmental_matrix_dissimilarity_method)

     return_list$Env_dissim_matrix <- Env_dissim_matrix

     #Print dissimilarity matrix for environmental data to results file
     sink(paste(filenameprefix,"_environmental_matrix_dissimilarity.txt"))
     print(Env_dissim_matrix)
     sink()
  
     #Cluster dissimilarity matrix for environmental data
     Env_dissim_matrix_cluster <- hclust(Env_dissim_matrix)

     return_list$Env_dissim_matrix_cluster <- Env_dissim_matrix_cluster

     #Plot dendogram for dissimilarity matrix for environmental data
     plot(Env_dissim_matrix_cluster,main=environmental_matrix_dendogram_name )
  }

  if (run_spec_dist) {
     #Calculate dissimilarity matrix for species distribution data
     Spec_dist_dissim_matrix <- vegdist(spec_matrix,method=species_distribution_dissimilarity_method,binary=presence_absence)

     return_list$Spec_dist_dissim_matrix <- Spec_dist_dissim_matrix

     #Print dissimilarity matrix for species distribution data to results file
     sink(paste(filenameprefix,"_species_distribution_matrix_dissimilarity.txt"))
     print(Spec_dist_dissim_matrix)
     sink()

     #Cluster dissimilarity matrix for species distribution data
     Spec_dist_dissim_matrix_cluster <- hclust(Spec_dist_dissim_matrix)

     return_list$Spec_dist_dissim_matrix_cluster <- Spec_dist_dissim_matrix_cluster

     #Plot dendogram for dissimilarity matrix for species distribution data
     plot(Spec_dist_dissim_matrix_cluster,main=species_distribution_dendogram_name)
  }

  if (run_envfit) {

     #Calculate envfit of environmental data to species matrix PCA
     spec_matrix_PCA_env_matrix_envfit <- envfit(spec_matrix_PCA,env_matrix,permutations=envfit_permutations)

     #Print envfit vector r-squared and p-values to console
     print(spec_matrix_PCA_env_matrix_envfit)

     #Print envfit vector r-squared and p-values to results file
     sink(paste(filenameprefix,"_envfit_vectors.txt"))
     print(spec_matrix_PCA_env_matrix_envfit)
     sink()

     return_list$spec_matrix_PCA_env_matrix_envfit <- spec_matrix_PCA_env_matrix_envfit

     
     #Define vector to hold indexes environmental factors w/significant p-values
     significant_factor_vector_indexes <- c()
     significant_factor_vector_pvals <- c()
 
     #Populate vector to hold indexes and p-values for environmental factors w/significant p-values
     for (i in 1:length(spec_matrix_PCA_env_matrix_envfit$vectors$pvals)) {
         if (spec_matrix_PCA_env_matrix_envfit$vectors$pvals[i]<=acceptable_p_value) {
            significant_factor_vector_pvals <- append(significant_factor_vector_pvals,spec_matrix_PCA_env_matrix_envfit$vectors$pvals[i])
            significant_factor_vector_indexes <- append(significant_factor_vector_indexes,i)
         }
     }

      
     #If there are at least 2 significant factors, plot envfit results and do CCA on species matrix using significant factors as environmental matrix.
     #Otherwise throw error message.
     if (length(significant_factor_vector_indexes)>1) {


        
        #Calculate and sort envfit vector angles from smallest to largest
        envfit_vector_angles <- ResizeVector_v3(spec_matrix_PCA_env_matrix_envfit$vectors$arrows[significant_factor_vector_indexes,1],spec_matrix_PCA_env_matrix_envfit$vectors$arrows[significant_factor_vector_indexes,2],labels(spec_matrix_PCA_env_matrix_envfit$vectors$arrows[significant_factor_vector_indexes,1]))
                
        return_list$envfit_vector_angles <- envfit_vector_angles$angles
        return_list$envfit_vector_angle_labels <- envfit_vector_angles$angle_labels

        #Build vector of colors for enfit vectors based on angles        
        envfit_vector_colors <- vector()
        color_counter=2
        envfit_vector_colors[1] <- color_counter

        for (i in 1:NROW(envfit_vector_angles$angles)) {
           if (i>1) {
              if (abs(envfit_vector_angles$angles[i]-envfit_vector_angles$angles[i-1]) <= 15) {
                 envfit_vector_colors[i] <- color_counter
              } else {
                 color_counter = color_counter + 1
                 envfit_vector_colors[i] <- color_counter
              }
           } 
        }

        return_list$envfit_vector_colors <- envfit_vector_colors

        envfit_labels <- labels(spec_matrix_PCA_env_matrix_envfit$vectors$arrows[significant_factor_vector_indexes,1])

        for (i in NROW(envfit_labels)) {
        }

        #Print sorted envfit vector angles, labels and colors to file
        sink(paste(filenameprefix,"_sorted_envfit_vector_angles.txt"))
        print(envfit_vector_angles$angle_labels)
        print(envfit_vector_angles$angles)
        print(envfit_vector_colors)
        sink()

        #Plot PCA species matrix
        biplot(spec_matrix_PCA, main=envfit_name, scaling="sites", display = c("species","sites"), type = c("text","text"),xlab=xlabel_spec,ylab=ylabel_spec)

        #Plot envfit of environmental data to species matrix PCA
        plot(spec_matrix_PCA_env_matrix_envfit, p.max=acceptable_p_value,col="blue")

        if (run_cca_rda) {
           #If the number of significant factors is less than the number of observations, run CCA and RDA on all significant factors
           #otherwise run CCA and RDA on only first n-1 significant factors with lowest p-values
           if (length(significant_factor_vector_indexes)<nrow(spec_matrix)) {

              

              sig_factor_CCA <- cca(spec_matrix,env_matrix[,significant_factor_vector_indexes],scale=spec_pca_scale)
              sig_factor_RDA <- rda(spec_matrix,env_matrix[,significant_factor_vector_indexes],scale=spec_pca_scale)

              #Compute eigenvalue percent variance explained and incorporate into axis labels
              eig_vals_cca <- eigenvals(sig_factor_CCA)
              eig_vals_cca_percent <- round((eig_vals_cca/(sum(eig_vals_cca)))*100)
              xlabel_cca = paste("CCA1 (",eig_vals_cca_percent[1],"%)",sep="")
              ylabel_cca = paste("CCA2 (",eig_vals_cca_percent[2],"%)",sep="")

              #Compute eigenvalue percent variance explained and incorporate into axis labels
              eig_vals_rda <- eigenvals(sig_factor_RDA)
              eig_vals_rda_percent <- round((eig_vals_rda/(sum(eig_vals_rda)))*100)
              xlabel_rda = paste("RDA1 (",eig_vals_rda_percent[1],"%)",sep="")
              ylabel_rda = paste("RDA2 (",eig_vals_rda_percent[2],"%)",sep="")

           } else {

             print("WARNING: There are greater than n-1 significant factors, so only n-1 factors with lowest p-values will be used for CCA.")

             significant_factor_vector_indexes_trimmed <- c()
             significant_factor_vector_indexes <- c()

             #Sort p-values
             significant_factor_vector_pvals <- sort(significant_factor_vector_pvals,decreasing=FALSE)


             #Populate factor index vector in same order as sorted p-values (smallest to largest)
             for (i in 1:length(significant_factor_vector_pvals)) {
                 for (j in 1:length(spec_matrix_PCA_env_matrix_envfit$vectors$pvals)) {
               
                     if (significant_factor_vector_pvals[i]==spec_matrix_PCA_env_matrix_envfit$vectors$pvals[j]) {
             
                        significant_factor_vector_indexes <- append(significant_factor_vector_indexes,j)
                     }
                 }
             }


             #New vector containing n-1 significant factors with lowest p-value
             for (i in 1:(nrow(spec_matrix)-1)) {
                 significant_factor_vector_indexes_trimmed <- append(significant_factor_vector_indexes_trimmed,significant_factor_vector_indexes[i])
             }

   
             #CCA and RDA with trimmed vector
             sig_factor_CCA <- cca(spec_matrix,env_matrix[,significant_factor_vector_indexes_trimmed])
             sig_factor_RDA <- rda(spec_matrix,env_matrix[,significant_factor_vector_indexes])

             #Compute eigenvalue percent variance explained and incorporate into CCA axis labels
             eig_vals_cca <- eigenvals(sig_factor_CCA)
             eig_vals_cca_percent <- round((eig_vals_cca/(sum(eig_vals_cca)))*100)
             xlabel_cca = paste("CCA1 (",eig_vals_cca_percent[1],"%)",sep="")
             ylabel_cca = paste("CCA2 (",eig_vals_cca_percent[2],"%)",sep="")

             #Compute eigenvalue percent variance explained and incorporate into RDA  axis labels
             eig_vals_rda <- eigenvals(sig_factor_RDA)
             eig_vals_rda_percent <- round((eig_vals_rda/(sum(eig_vals_rda)))*100)
             xlabel_rda = paste("RDA1 (",eig_vals_rda_percent[1],"%)",sep="")
             ylabel_rda = paste("RDA2 (",eig_vals_rda_percent[2],"%)",sep="")
           }

           #Plot CCA
           plot(sig_factor_CCA,main=CCA_name,display = c("species","sites","bp"),type = c("text"),xlab=xlabel_cca,ylab=ylabel_cca)

           #Print eigenvalues for CCA to results file
           sink(paste(filenameprefix,"_CCA_percent_variance.txt"))
           print(eig_vals_cca_percent)
           sink()

           return_list$sig_factor_CCA <- sig_factor_CCA
           return_list$eig_vals_cca_percent <- eig_vals_cca_percent 

           #Plot RDA
           plot(sig_factor_RDA,main=RDA_name,display = c("species","sites","bp"),type = c("text"),xlab=xlabel_rda,ylab=ylabel_rda)

           #Print eigenvalues for RDA to results file
           sink(paste(filenameprefix,"_RDA_percent_variance.txt"))
           print(eig_vals_rda_percent)
           sink()

           return_list$sig_factor_RDA <- sig_factor_RDA
           return_list$eig_vals_rda_percent <- eig_vals_rda_percent 

     
        }  
     } else {
        factor_num_err_message <- paste("WARNING: Less than 2 significant environmental factors (at", acceptable_p_value, "level). Cannot do CCA or RDA.",sep=" ")
        print(factor_num_err_message)
     }
  }     
  dev.off()

  return(return_list)
 
}