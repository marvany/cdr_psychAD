# https://raw.githubusercontent.com/cran/BaylorEdPsych/master/R/PseudoR2.R
# Nagelkerke, N. J. D. 1991. A note on the general definition of the coefficient of determination. Biometrika, 78:3, 691-692.

## Example
## data(mtcars)
## fit <- glm(vs ~ cyl + hp, data = mtcars, family=binomial(link="logit"))
## PseudoR2(fit)

#' Calculate various pseudo R-squared values for a given glm model
#'
#' @param glmModel A fitted generalized linear model object of class "glm"
#' @return A named vector of pseudo R-squared values and AIC values
#' @examples
#' data(mtcars)
#' fit <- glm(vs ~ cyl + hp, data = mtcars, family=binomial(link="logit"))
#' PseudoR2(fit)
#'
#' Calculate Various Pseudo R-squared Values for a GLM Model
#'
#' This function computes several pseudo R-squared values for a given generalized linear model (GLM).
#'
#' @param glmModel A fitted GLM model object of class \code{glm}.
#'
#' @return A named vector containing the following pseudo R-squared values and model fit statistics:
#' \itemize{
#'   \item \code{McFadden}: McFadden's R-squared
#'   \item \code{Adj.McFadden}: Adjusted McFadden's R-squared
#'   \item \code{Cox.Snell}: Cox/Snell R-squared
#'   \item \code{Nagelkerke}: Nagelkerke/Cragg-Uhler R-squared
#'   \item \code{McKelvey.Zavoina}: McKelvey and Zavoina pseudo R-squared
#'   \item \code{Effron}: Effron R-squared
#'   \item \code{Count}: Count R-squared (proportion correctly classified)
#'   \item \code{Adj.Count}: Adjusted Count R-squared (proportion correctly classified)
#'   \item \code{AIC}: Akaike Information Criterion
#'   \item \code{Corrected.AIC}: Corrected Akaike Information Criterion for finite sample size
#' }
#'
#' @details
#' The function calculates the following pseudo R-squared values:
#' \itemize{
#'   \item \strong{McFadden's R-squared}: Based on the log-likelihood of the null and full models.
#'   \item \strong{Adjusted McFadden's R-squared}: Adjusted version of McFadden's R-squared.
#'   \item \strong{Cox/Snell R-squared}: Based on the likelihood ratio statistic.
#'   \item \strong{Nagelkerke/Cragg-Uhler R-squared}: Adjusted version of Cox/Snell R-squared.
#'   \item \strong{McKelvey and Zavoina pseudo R-squared}: Based on the variance of the predicted values.
#'   \item \strong{Effron R-squared}: Based on the sum of squared errors.
#'   \item \strong{Count R-squared}: Proportion of correctly classified observations.
#'   \item \strong{Adjusted Count R-squared}: Adjusted version of Count R-squared.
#' }
#'
#' @examples
#' \dontrun{
#'   model <- glm(y ~ x1 + x2, family = binomial, data = mydata)
#'   pseudoR2 <- PseudoR2(model)
#'   print(pseudoR2)
#' }
#'
#' @export
PseudoR2<-function(glmModel){
#if (length(class(object)) > 1)
#if (class(glmModel)!="glm" | class(glmModel)!="lm"){
#stop("Object not of class 'glm'")
#}
#else {
#	glmModel<-glmModel
#	}

# Calculate log likelihood for the null model
logLikN <- glmModel$null / -2  

# Calculate log likelihood for the full model
logLikF <- glmModel$dev / -2  

# Calculate the likelihood ratio statistic
G2 <- glmModel$null - glmModel$deviance

# Number of observations
n <- length(glmModel$y)	

# Predicted values from the model
ystar <- predict(glmModel, type="response") 

# Classify predictions as 1 if predicted value > 0.5, else 0
class1 <- ifelse(ystar > .5, 1, 0) 

# Create a classification table
classtab <- table(class1, glmModel$y, dnn = c("Predicted", "Actual")) 

# Maximum count of observations in any class
maxOut <- max(margin.table(classtab, 2))

# Number of predictors
p <- glmModel$rank

# Calculate penalty for AIC correction
penaltyN <- 2 * (p) * (p + 1)  
penaltyD <- n - p - 1  
penalty <- penaltyN / penaltyD

# Predicted values for McKelvey and Zavoina pseudo R-squared
MZystar <- predict(glmModel) 

# Sum of squared errors for McKelvey and Zavoina pseudo R-squared
sse <- sum((MZystar - mean(MZystar))^2) 

# Variance for the link function
s2 <- switch(glmModel$family$link, "probit" = 1, "logit" = pi^2 / 3, NA)  

# Sum of squared errors for Effron R-squared
Enum <- sum((glmModel$y - ystar)^2) 
Edenom <- sum((glmModel$y - mean(glmModel$y))^2) 

# Calculate McFadden's R-squared
r2McF <- 1 - logLikF / logLikN  

# Calculate adjusted McFadden's R-squared
r2McFA <- 1 - (logLikF - p - 1) / logLikN 

# Calculate Cox/Snell R-squared
r2CS <- 1 - exp(-G2 / n) 

# Calculate Nagelkerke/Cragg-Uhler R-squared
r2N <- (1 - exp((glmModel$dev - glmModel$null) / n)) / (1 - exp(-glmModel$null / n))

# Calculate McKelvey and Zavoina pseudo R-squared
r2MZ <- sse / (n * s2 + sse)  

# Calculate Effron R-squared
r2E <- 1 - (Enum / Edenom) 

# Calculate Count R-squared (proportion correctly classified)
r2C <- (classtab[1] + classtab[4]) / n 

# Calculate adjusted Count R-squared (proportion correctly classified)
r2CA <- (classtab[1] + classtab[4] - maxOut) / (n - maxOut) 

# Calculate Akaike Information Criterion (AIC)
aic <- 2 * (p) + glmModel$dev 

# Calculate corrected AIC for finite sample size
Caic <- aic + penalty 
 
results<-c(McFadden=r2McF, Adj.McFadden=r2McFA, Cox.Snell=r2CS, Nagelkerke=r2N, McKelvey.Zavoina=r2MZ, Effron=r2E, Count=r2C, Adj.Count=r2CA, AIC=aic, Corrected.AIC=Caic)
return(results)

}
