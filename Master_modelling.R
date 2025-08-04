###########################
########## TO DO ##########
###########################

#review GAM assumptions
#ind of random effects

###########################
###########################
##### Data exploration ####
###########################
###########################

#NAs
#outliers:cleveland dotplot for the response variable; for the covariates
#collinearity: continuous v continuous; continuous v. categorical
#relationships: continuous vars v response; categorical vars v. response
#dependency: time, random effect

#!!! for data exploration, must cite:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

#NAs
#colSums(is.na(dat))

#cont.var <- c("var1","var2")
#Outliers
Mydotplot <- function(DataSelected){
  P <- dotplot(as.matrix(as.matrix(DataSelected)),
               groups=FALSE,
               strip = strip.custom(bg = 'white',
                                    par.strip.text = list(cex = 1.2)),
               scales = list(x = list(relation = "free", draw = TRUE),
                             y = list(relation = "free", draw = FALSE)),
               col=1, cex  = 0.5, pch = 16,
               xlab = list(label = "Value of the variable", cex = 1.5),
               ylab = list(label = "Order of the data from text file", cex = 1.5))
  print(P)  
}
#Mydotplot(dat[,cont.var])

#Collinearity
#Continuous variables
Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) points(x, y, 
                                             pch = 16, cex = 0.8, 
                                             col = gray(0.1)))
  print(P)
}

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
    cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex * r)
}

#Mypairs(dat[, cont.var])
#manual alternative for matrix of correlations 
#cont.var.cor <- as.matrix(cor(na.omit((plot[, cont.var])))) 
#manual alternative for plot:
#pairs(na.omit((plot[, cont.var])))

#Multipanel boxplots
#Continuous v categorical variables
Mybwplot <- function(Z, MyVar, TargetVar){
  AllY <- as.vector(as.matrix(Z[,MyVar]))
  AllX <- rep(Z[,TargetVar], length(MyVar))
  ID <- rep(MyVar, each = nrow(Z))
  P <- bwplot(AllY ~ factor(AllX) | ID, horizontal = FALSE,
              ylab = "", xlab = "",
              scales = list(alternating = TRUE,cex.lab = 1.5,
                            x = list(relation = "same",rot =90, abbreviate = TRUE, cex = 1.5),
                            y = list(relation = "free", draw = FALSE)),
              strip = strip.custom(bg = 'white',
                                   par.strip.text = list(cex = 1.2)),
              cex = .5,
              par.settings = list(
                box.rectangle = list(col = 1),
                box.umbrella  = list(col = 1),
                plot.symbol   = list(cex = .5, col = 1)))
  print(P)
}
#Mybwplot(dat, cont.var, "Sex")
#manual alternative: 
#ggplot(data = plot, aes(x = x, y = y)) + geom_boxplot()

#Multipanel scatterplots
Myxyplot <- function(Z, MyV, NameY1, MyXlab = "", MyYlab="") {
  AllX  <- as.vector(as.matrix(Z[,MyV]))
  AllY  <- rep(Z[,NameY1] , length(MyV))
  AllID <- rep(MyV, each = nrow(Z))
  library(mgcv)
  library(lattice)
  P <- xyplot(AllY ~ AllX|factor(AllID), col = 1,
              xlab = list(MyXlab, cex = 1.5),
              #ylab = list("Response variable", cex = 1.5),
              #ylab = list("Pearson residuals", cex = 1.5),
              ylab = list(MyYlab, cex = 1.5),
              #layout = c(2,2),   #Modify
              strip = function(bg='white', ...)
                strip.default(bg='white', ...),
              scales = list(alternating = TRUE,
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              panel=function(x, y){
                panel.grid(h=-1, v= 2)
                panel.points(x, y, col = 1)
                panel.loess(x, y, span = 0.8,col = 1, lwd = 2)
              })
  print(P)
}
#Myxyplot(dat, var, "Wingcrd")
#manual alternative: 
#ggplot(data = plot, aes(x = x, y = y)) + geom_point() 

###########################
###########################
####### Model types #######
###########################
###########################

###########################
####### Count data ########
###########################

#Poisson
  #link = log
  #can be overdispersed

#correct under/overdispersion

#Negative binomial
  #link = log
  #cannot be overdispersed

#correct over and underdispersion
#Generalized poisson
  #cannot be overdispersed

#Bernoilli
  #link = logit
  #yes/no response
  #can be overdispersed

###########################
##### Proportion data #####
###########################

#Binomial - cbind(success,failure)
  #link = logit
  #can be overdispersed
#Beta - bounded between 0 and 1 (not inclusive) 
  #link = logit
  #cannot be overdispersed
  #can transform: (Y*(N-1)+0.5))/N

###########################
#### Continuous data ######
###########################

#Gamma - positive data
#Tweedie - positive data with many zeros; cannot be overdispersed

###########################
##### Zero-inflated #######
###########################

#glmmTMB package
#for intercept only for Bernoilli part
#ziformula = ~1

#binomial formula
#family = nbinom2

#tweedie
#family = tweedie()

###########################
###########################
#### Model validation #####
###########################
###########################

#These can be done manually or with DHARMa

#LM
#homogeneity: residuals v fitted
#independence: residuals v covariates in/out of model
#normality: of residuals
#influential observations: cooks distance

#GLM
#homogeneity: residuals v fitted (not informative for binomial)
#independence: residuals v covariates in/out of model (including random effects)
#influential observations: cooks distance
#overdispersion: >1 = overdispersed (except NB, GP, CMP, beta, gamma; these built in correct for dispersion)
#zero-inflation
#spatial/temporal autocorrelation: variogram

#GAM
#homogeneity: residuals v fitted (not informative for binomial)
#independence: residuals v covariates in/out of model
#normality: of residuals
#??? overdispersion: unsure
#spatial/temporal autocorrelation: variogram

#homogeneity (residuals v fitted) function
resid_fit_plot <- function(mod,dat) {
  E2 <- resid(mod, type = "pearson")
  F2 <- fitted(mod)
  plot(x = F2, 
       y = E2, 
       xlab = "Fitted values",
       ylab = "Residuals")
  abline(h = 0) 
}
#resid_fit_plot(mod,dat)

#independence function continuous
indep_cont_plot <- function(mod,dat,var) {
  E2 <- resid(mod, type = "pearson")
  plot(x = var, 
       y = E2, 
       xlab = "var",
       ylab = "Pearson residuals", 
       pch = 16)
  abline(h = 0) 
}
#indep_cont_plot(mod,dat,df$var)

#independence function categorical
indep_cat_plot <- function(mod,dat,var) {
  E2 <- resid(mod, type = "pearson")
  boxplot(E2 ~ var, 
          data = dat, 
          cex.lab = 1.5,
          xlab = "var",
          ylab = "Pearson residuals")
  abline(h = 0) 
}
#indep_cat_plot(mod,dat,df$var)

#independence of random effects

#ai <- ranef(mod)$fSite$`(Intercept)` 

#obtain standard errors of these random effects
#ai.se <- se.ranef(mod)$fSite

#plot estimated random intercepts (dots) and their 95% confidence intervals (bars)
#MyData3 <- data.frame(Plot = 1:35,
#                     ai   = ranef(mod)$fSite$`(Intercept)`,
#                      aiSE = as.numeric(se.ranef(mod)$fSite))

#MyData3$SeUp <- MyData3$ai + 1.96 * MyData3$aiSE
#MyData3$SeLo <- MyData3$ai - 1.96 * MyData3$aiSE

#p3 <- ggplot()
#p3 <- p3 + geom_point(data = MyData3, 
#                      aes(y = Plot, x = ai),
#                      shape = 16, 
#                      size = 2)
#p3 <- p3 + xlab("Random intercepts") + ylab("Site")
#p3 <- p3 + theme(text = element_text(size=15))
#p3 <- p3 + geom_vline(xintercept = 0, 
#                      lty = 2,
#                      colour = "black")
#p3 <- p3 + geom_errorbar(data = MyData3,
#                         aes(y = Plot, 
#                             xmax = SeUp, 
#                             xmin = SeLo), 
#                         width = 0.1)
#p3

#cooks distance function
cooks_dist_plot <- function(mod) {
  plot(cooks.distance(mod), 
       type = "h",
       ylim = c(0,1),
       xlab = "Observations",
       ylab = "Cook distance values")
  abline(h = 0, lwd = 2, lty = 2) 
}
#cooks_dist_plot(mod)

#normality function
resid_hist <- function(mod) {
  E2 <- resid(mod, type = "pearson")
  hist(E2, main = "", 
       xlab = "Residuals")
}
#resid_hist(mod)

#qqplots
resid_qqplot <- function(mod) {
  E2 <- rstandard(mod)
  qqnorm(E2, main = "")
  qqline(E2)
}
#resid_qqplot(mod)
#this doesn't work for certain model types

#overdispersion function (poisson; var = mean)
disp_check <- function(mod,dat) {
  N <- nrow(dat)
  p <- length(coef(mod))
  E1 <- resid(mod, type = "pearson")
  Dispersion <- sum(E1^2)/ (N-p)
  return(Dispersion)
}
#disp_check(mod,dat)

###########################
#### DHARMa validation ####
###########################

#set.seed(556)
#testDispersion(mod) #test for over/under dispersion, tends to lean toward no overdispersion
#testDispersion(mod, type = "PearsonChisq") #no simluation version
#testZeroInflation(mod) #test zero inflation

#simulationOutput <- simulateResiduals(mod, plot = F) #calculate residuals (randomized quantile residuals)
#par(mfrow = c(2, 2)) #set panel arrangement
#plotResiduals(simulationOutput) #homogeneity (residuals v fitted)
#plotResiduals(simulationOutput, form = dat$var1) #independence continuous
#plotResiduals(simulationOutput, form = dat$var2) #independence categorical
#plotQQunif(simulationOutput) #qqplot

###########################
###########################
### Spatial Correlation ###
###########################
###########################

#variogram of residuals
#spat.dat <- data.frame(E2 = resid(mod, type = "pearson"), 
#                      Xkm = RF$longitude, 
#                      Ykm = RF$latitude)

#coordinates(spat.dat)<-c("Xkm", "Ykm")

#V1 <- variogram(E2 ~ 1, 
#                spat.dat, 
#                cressie = TRUE)

#p <- ggplot()
#p <- p + geom_point(data = V1,
#                    aes(x = dist, 
#                        y = gamma))
#p <- p + geom_smooth(data = V1,
#                     se = FALSE,
#                     span = 0.9,
#                     aes(x = dist, 
#                         y = gamma))
#p <- p + xlab("Distance") + ylab("Semi-variogram")
#p <- p + theme(text = element_text(size = 15)) 
#p <- p + theme(legend.position="none") 
#p

#distinct lat/lon
#important: latitude and longitude must be named correctly
Spat.cor <- function(mod,dat,dist) {
  coords <- cbind(dat$longitude, dat$latitude)
  #matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  #matrix.dist[1:10, 1:10]
  #matrix.dist.inv <- 1/matrix.dist
  #matrix.dist.inv[1:10, 1:10]
  #diag(matrix.dist.inv) <- 0
  #matrix.dist.inv[1:10, 1:10]
  myDist = dist
  # calculate residuals autocovariate (RAC)
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T)
  return(rac)
}
#Spat.cor(mod,dat,dist)

#repeat lat/lon; this moves them slightly around
Spat.cor.rep <- function(mod,dat,dist) {
  coords <- cbind(dat$longitude, dat$latitude) + matrix(runif(2*nrow(dat), 0, 0.00001), nrow = nrow(dat), ncol = 2)
  #matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  #matrix.dist[1:10, 1:10]
  #matrix.dist.inv <- 1/matrix.dist
  #matrix.dist.inv[1:10, 1:10]
  #diag(matrix.dist.inv) <- 0
  #matrix.dist.inv[1:10, 1:10]
  myDist = dist
  # calculate residuals autocovariate (RAC)
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T)
  return(rac)
}
#Spat.cor.rep(mod,dat,dist)

###########################
###########################
##### Model selection #####
###########################
###########################

#approaches to model selection
#1. leave model as is.
#2. AIC selecion (step() function below)
#3. information theoretic approach
#4. p-value based drop one var at a time

#automated backward selection
#step(mod) 

#test whether vars matter
#drop1(mod)

#GLM generalised R2
#(null deviance - residual deviance)/ null deviance

#cond v marginial R2
#performance(mod)

#compare multiple models
#AIC(mod1, mod2) 
#anova(mod1, mod2, test = "F")

###########################
###########################
####### Predictions #######
###########################
###########################

#this code is to aboid negative CIs from the predict function

# Make a function that defines the log-link
MyLog <- function(x) {exp(x) }

#Pred1  <- predict(mod1, newdata = dat, se = TRUE, link = "predictor")
#dat$mu   <- MyLog(Pred1$fit)
#dat$SeUp <- MyLog(Pred1$fit + 1.96 * Pred1$se.fit)
#dat$SeLo <- MyLog(Pred1$fit - 1.96 * Pred1$se.fit)

# Define a function to calculate the logistic link function
MyLogit <- function(x) {exp(x) / (1 + exp(x))}

#Pred1 <- predict(mod1, newdata = dat, se = TRUE, link = "predictor")
#dat$Pi   <- MyLogit(Pred1$fit)
#dat$SeUp <- MyLogit(Pred1$fit + 1.96 * Pred1$se.fit)
#dat$SeLo <- MyLogit(Pred1$fit - 1.96 * Pred1$se.fit)

###########################
###########################
####### Miscelanious ######
###########################
###########################

#Standardize the continuous covariates
MyStd <- function(x) { (x - mean(x)) / sd(x)}

#count zeros in dataset
#100 * sum(dat$var== 0) / nrow(dat) 