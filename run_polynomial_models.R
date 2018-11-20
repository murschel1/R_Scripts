# SCRIPT: run_polynomial_models.R
# VERSION: 1.0
# DATE: 3/8/2018
# AUTHOR: Matthew R. Urschel
# DESCRIPTION: Script to generate polynomial models from a given data matrix using "leave one out cross-validation" (LOOCV) 
#              as a resampling based performance measure for model validation and selection. Uses package 'caret'. 
# Arguments: data_matrix = matrix with independant and dependant variables as X and Y axis, respectively,
#            filenameprefix = name of pdf file to save results into (optional, default is 'Analysis_pipeline_results'),
#            col_num_predictor = matrix column of predictor variable (optional, default = 1),
#            col_num_predicted = matrix column of predicted variable (optional, default = 1),
#            prediction_increment = value by which to increment predictor variable when generating model predictions (optional, default = 0.1),
#            prediction_interval_lowb = lower bound of prediction interval (optional, default = 0),
#            prediction_interval_upb = upper bound of prediction interval (optional, default = 300),
#            degrees = highest order polynomial model to construct (optional, default = 10, maximum of 10),
run_polynomial_models <- function(data_matrix,
                                  filenameprefix="polynomial_model_results",
                                  col_num_predictor=1,
                                  col_num_predicted=2,
                                  prediction_increment=0.1,
                                  prediction_interval_lowb=0,
                                  prediction_interval_upb=300,
                                  degrees=10) {
 

  if (degrees>10) {
     print("Maximum degree is 10. Setting degree to 10.")
     degrees=10
  }

  if (degrees>(NROW(data_matrix)-1)) {
     print(paste("Not enough degrees of freedom! Setting degree to ",NROW(data_matrix)-1),".")
     degrees=NROW(data_matrix)-1
  }

  if(require("caret")){
    print("caret is loaded correctly.")
  } else {
    print("Trying to install caret.")
    install.packages("caret")
    if(require("caret")){
        print("caret installed and loaded.")
    } else {
        stop("Could not install caret!")
    }
  }

  x.pred <- seq(prediction_interval_lowb, prediction_interval_upb, prediction_increment)

  x <- data_matrix[col_num_predictor]
  y <- data_matrix[col_num_predicted]
  datframe <- data.frame(x,y)
  names(datframe) <- c('x','y')
  sink(paste(filenameprefix,".txt"))

  for (i in 2:degrees) {
      if (i==2){
         model.loocv <- train(y ~ x + I(x^2), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI2.csv"))}
      if (i==3){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI3.csv"))}
      if (i==4){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI4.csv"))}
      if (i==5){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4,'I(x^5)'=x.pred^5),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI5.csv"))}
      if (i==6){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4,'I(x^5)'=x.pred^5,'I(x^6)'=x.pred^6),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI6.csv"))}
      if (i==7){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4,'I(x^5)'=x.pred^5,'I(x^6)'=x.pred^6,'I(x^7)'=x.pred^7),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI7.csv"))}
      if (i==8){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4,'I(x^5)'=x.pred^5,'I(x^6)'=x.pred^6,'I(x^7)'=x.pred^7,'I(x^8)'=x.pred^8),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI8.csv"))}
      if (i==9){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4,'I(x^5)'=x.pred^5,'I(x^6)'=x.pred^6,'I(x^7)'=x.pred^7,'I(x^8)'=x.pred^8,'I(x^9)'=x.pred^9),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI9.csv"))}
      if (i==10){
         model.loocv <- train(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10), method = "lm", data = datframe, trControl = trainControl(method = "LOOCV"))
         model.loocv.ci <- predict(model.loocv$finalModel,list(x=x.pred,'I(x^2)'=x.pred^2,'I(x^3)'=x.pred^3,'I(x^4)'=x.pred^4,'I(x^5)'=x.pred^5,'I(x^6)'=x.pred^6,'I(x^7)'=x.pred^7,'I(x^8)'=x.pred^8,'I(x^9)'=x.pred^9,'I(x^10)'=x.pred^10),interval="prediction")
         write.csv(model.loocv.ci,paste(filenameprefix,"_PI10.csv"))}
   
      print('---------------------------------------------------------------')
      print(paste('Degree - ',i))
      print('---------------------------------------------------------------')
      print('')

      print('LOOCV Stats')
      print('')
      print(model.loocv)
      print('')
      print('Best model summary')
      print(' ')
      print(summary(model.loocv$finalModel))

      

  }

  
  sink()



}