#Script to calculate new x,y coordinates for endpoint of resized vector
#(e.g. arrows drawn on a biplot using the arrows() function). Assumes vector starts at origin (0,0).
#Also returns array of angles between vectors and positive y axis (positive y axis is 0 degrees, with angle increasing in clockwise direction).
#Arguments: xarray = Array containing x coordinates of arrow endpoints,
#           yarray = Array containing y coordinates of arrow endpoints,
#           percent_of_previous = Percentage of original arrow length by which arrow length will be decreased or increased (optional, default = 10%),
#           increase = Boolean indicator controlling whether arrow length will be increased (if TRUE) or decreased (if FALSE) (optional, default = FALSE)

ResizeVector_v3 <- function(xarray, 
                            yarray, 
                            percent_of_previous=25, 
                            increase=FALSE) {


   pointlist <- list()
   pointlist$x <- xarray
   pointlist$y <- yarray
   pointlist$angles <- array()

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
         pointlist$angles[i] <- abs(ang*57.29)+90
      } else {
         
         pointlist$angles[i] <- 270+(180-(ang*57.29))

         if (pointlist$angles[i]>=360) {
            pointlist$angles[i]=pointlist$angles[i]-360
         }
        
      }

   }

   return(pointlist)

}