#' logreg.chisq.gof(my.glm)
#' 
#' This function performs a chi-squared goodness of fit test for the 
#' overall fit of a logistic regression model.  A large p-value indicates 
#' no evidence of poor fit; a small p-value indicates evidence of poor fit.
#' 
#' @param my.glm A glm model object.
#' 
#' @return A list containing the following components:
#' @return fit = Fitted values from glm model
#' @return r = Resdiausl from glm model
#' @return ssR = sum of squared residuals
#' @return df = degrees of freedom (n-p)
#' @return pval = p-value for chi-squared gof test
#' 
#' @examples
#' my.model <- glm(y ~ x1 + x2, data=mydata, family='binomial')
#' logistic.chisq.gof(my.model)
#' @export

logreg.chisq.gof = function(my.glm){
   #-------------------------------------------------------------
   # FUNCTION: Performs logistic regression chisq gof test.
   #-------------------------------------------------------------
   # INPUTS:	my.glm         = a glm model object
   #           response       = vector of observed responses
   #-------------------------------------------------------------	
   # OUTPUTS:	fit            = fitted values
   #           r              = residuals
   #           ssR            = sum of squared residuals
   #           df             = degrees of freedom (n-p)
   #           pval           = p-value for gof test
   #-------------------------------------------------------------
   #
   ## Chisq Goodness of fit test (deviance, chi-sq) for the model vs the null model.
   
   # Save fitted values.
   fit = fitted(my.glm)
   
   # Calculate residuals
   response = my.glm$y
   r = (response - fit) / sqrt(fit*(1-fit))
   
   #Sum of squares of these residuals follows a chi-square dist with n-p df.
   df = my.glm$df.residual
   ssResid = sum(r^2)
   
   # Calculate p-value.
   pval = 1 - pchisq(ssResid, df=df)
   
   # Return function output.
   return(list('fit'  = fit,
               'r'    = r,
               'ssR'  = ssResid,
               'df'   = df,
               'pval' = pval))
}
