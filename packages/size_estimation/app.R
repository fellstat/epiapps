#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)
library(MCMCpack)
library(msm)
library(ggplot2)
library(GGally)
library(SizeEstimation)

customPlot <- function(data, mapping, ...)
{
  ggplot(data = data, mapping = mapping) + geom_point(..., alpha = .01) + geom_smooth()
}

sizeestima <- function( DATA,size, burnin,thin) {
  
  data = DATA[which(DATA[,4]>-1 | DATA[,8]>-1),]
  N = data[,3]
  N.sum = sum(N)
  x = data[,5:7]
  r = data[,5]+data[,6]-data[,7]
  y = data[,4]
  z = data[,8]
  index1 = setdiff(which(data[,5]>-1), which(r>-1))	# BSS but not capture-recapture
  index2 = setdiff(which(data[,6]>-1), which(r>-1))	# NEP but not capture-recapture
  index3 = which(r>-1)						# capture recapture
  index4 = which(z>-1)						# RSA
  n.min = apply(cbind(y,x[,1:2],r,rep(-1,length(y))), 1, max, na.rm = T)+1
  y.sum = n.sum = rep(0,3)
  for (i in 1:3)	y.sum[i] = sum(data[which(data[,3+i]>0),i+3])
  option.mu = 1                   # allow heterogeneous mu?
  option.ht = 1                   # allow heterogeneous mu?
  rho = 1                         # pij = pi x pj x rho, so rho=1 means independence
  
  ## 		MCMC
  
  
  n.total.q = NULL
  #size = 500 * 1000
  a = b = matrix(NA, size, 4)
  n = phi = matrix(NA, size, length(N))
  p1 = p2 = theta = matrix(NA, size, length(N))
  mu = sig2 = rep(NA, size)
  zl = length(index4)
  tau02 = sig02 = (log(10)/2)^2
  
  
  
  # this contains the functions called by MCMC.R
  
  ##		Initial
  if (rho==1){
    a[1,] = 2
    b[1,1] = 2000
    b[1,2:4] = 2
    phi[1,] = rbeta(length(N), a[1,1], b[1,1])
    n[1,] = rbinom(length(N), N, phi[1,])
    index = which(n[1,]<n.min)
    n[1,index] = rnbinom(length(index), n.min[index], 0.9) + n.min[index]
    n.index = list(which(data[,4]>-1), which(data[,5]>-1), which(data[,6]>-1))
    s2 = sum((log(z[index4]/n[1,index4])-mean(log(z[index4]/n[1,index4])))^2)/(length(index4)-1)
    mu[1] = mean(log(z[index4]/n[1,index4]))
    sig2[1] = s2
    l1 = length(which(data[,4]>-1))
    l2 = length(which(data[,5]>-1))
    l3 = length(which(data[,6]>-1))
    l4 = length(N)
    if (option.ht==0){
      p1 = p2 = theta = rep(NA, size)
      S = rep(NA,3)
    }
  }
  
  if (rho!=1){
    a[1,] = apply(a[keep,], 2, mean)
    b[1,] = apply(b[keep,], 2, mean)
    phi[1,] = apply(phi[keep,], 2, mean)
    n[1,] = round(apply(n[keep,], 2, mean))
    s2 = sum((log(z[index4]/n[1,index4])-mean(log(z[index4]/n[1,index4])))^2)/(length(index4)-1)
    mu[1] = mean(log(z[index4]/n[1,index4]))
    sig2[1] = s2
  }
  
  reject = rep(0,8+length(N))
  check = 0
  set.seed(1000)
  ptm <- proc.time()
  for (t in 2:size){
    check = check+1
    # update p1, p2, theta, phi
    if (option.ht==1){
      p1[t-1,n.index[[2]]] = rbeta(l2, a[t-1,2]+x[n.index[[2]],1], b[t-1,2]+n[t-1,n.index[[2]]]-x[n.index[[2]],1])
      p2[t-1,n.index[[3]]] = rbeta(l3, a[t-1,3]+x[n.index[[3]],2], b[t-1,3]+n[t-1,n.index[[3]]]-x[n.index[[3]],2])
      theta[t-1,n.index[[1]]] = rbeta(l1, a[t-1,4]+y[n.index[[1]]], b[t-1,4]+n[t-1,n.index[[1]]]-y[n.index[[1]]])
    }
    if (option.ht==0){
      for (i in 1:3)	S[i] = sum(n[t-1,n.index[[i]]])
      theta[t-1] = rbeta(1, 1+y.sum[1], 1+S[1]-y.sum[1])
      p1[t-1] = rbeta(1, 1+y.sum[2], 1+S[2]-y.sum[2])
      p2[t-1] = rbeta(1, 1+y.sum[3], 1+S[3]-y.sum[3])
    }
    phi[t-1,] = rbeta(l4, a[t-1,1]+n[t-1,], b[t-1,1]+N-n[t-1,])
    
    # update n by H-M
    for (j in 1:length(N)){
      n[t,j] = rpois(1, n[t-1,j])
      if (n[t,j]<n.min[j])	n[t,j]=n[t-1,j]
      if (n[t,j]>=n.min[j]){
        if (option.ht==1){
          diff = dbinom(n[t,j], N[j], phi[t-1,j], log=T) - dpois(n[t,j], n[t-1,j], log=T) -
            dbinom(n[t-1,j], N[j], phi[t-1,j], log=T) + dpois(n[t-1,j], n[t,j], log=T)
          if (is.element(j, n.index[[1]]))
            diff = diff + dbinom(y[j], n[t,j], theta[t-1,j], log=T) - dbinom(y[j], n[t-1,j], theta[t-1,j], log=T)
          if (is.element(j, n.index[[2]]))
            diff = diff + dbinom(x[j,1], n[t,j], p1[t-1,j], log=T) - dbinom(x[j,1], n[t-1,j], p1[t-1,j], log=T)
          if (is.element(j, n.index[[3]]))
            diff = diff + dbinom(x[j,2], n[t,j], p2[t-1,j], log=T) - dbinom(x[j,2], n[t-1,j], p2[t-1,j], log=T)
          if (is.element(j, index3))
            #	diff = diff + dnbinom(n[t,j]-r[j], r[j]-1, 1-(1-p1[t-1,j])*(1-p2[t-1,j]), log=T) -
            #		 dnbinom(n[t-1,j]-r[j], r[j]-1, 1-(1-p1[t-1,j])*(1-p2[t-1,j]), log=T)
            diff = diff + dnbinom(n[t,j]-r[j], r[j]-1, p1[t-1,j]+p2[t-1,j]-p1[t-1,j]*p2[t-1,j]*rho, log=T) -
              dnbinom(n[t-1,j]-r[j], r[j]-1, p1[t-1,j]+p2[t-1,j]-p1[t-1,j]*p2[t-1,j]*rho, log=T)
          
          if (is.element(j, index4) & option.mu==1)
            diff = diff + dnorm(log(z[j]/n[t,j]), mu[t-1], sqrt(sig2[t-1]), log=T) -
              dnorm(log(z[j]/n[t-1,j]), mu[t-1], sqrt(sig2[t-1]), log=T)
        }
        if (option.ht==0){
          diff = dbinom(n[t,j], N[j], phi[t-1,j], log=T) - dpois(n[t,j], n[t-1,j], log=T) -
            dbinom(n[t-1,j], N[j], phi[t-1,j], log=T) + dpois(n[t-    input$start1,j], n[t,j], log=T)
          if (is.element(j, n.index[[1]]))
            diff = diff + dbinom(y[j], n[t,j], theta[t-1], log=T) - dbinom(y[j], n[t-1,j], theta[t-1], log=T)
          if (is.element(j, n.index[[2]]))
            diff = diff + dbinom(x[j,1], n[t,j], p1[t-1], log=T) - dbinom(x[j,1], n[t-1,j], p1[t-1], log=T)
          if (is.element(j, n.index[[3]]))
            diff = diff + dbinom(x[j,2], n[t,j], p2[t-1], log=T) - dbinom(x[j,2], n[t-1,j], p2[t-1], log=T)
          if (is.element(j, index3))
            diff = diff + dnbinom(n[t,j]-r[j], r[j]-1, 1-(1-p1[t-1])*(1-p2[t-1]), log=T) -
              dnbinom(n[t-1,j]-r[j], r[j]-1, 1-(1-p1[t-1])*(1-p2[t-1]), log=T)
          if (is.element(j, index4) & option.mu==1)
            diff = diff + dnorm(log(z[j]/n[t,j]), mu[t-1], sqrt(sig2[t-1]), log=T) -
              dnorm(log(z[j]/n[t-1,j]), mu[t-1], sqrt(sig2[t-1]), log=T)
        }
        if (diff<log(runif(1,0,1)))	n[t,j] = n[t-1,j]
      }
    }
    reject[which(n[t,]==n[t-1,])] = reject[which(n[t,]==n[t-1,])]+1
    
    # update (a,b) by H-M
    for (m in 1:(1+3*option.ht)){
      lower = log(1/(a[t-1,m]+b[t-1,m]))
      upper = log(1-1/(a[t-1,m]+b[t-1,m]))
      if (m==1) sd = 0.25
      if (m==2) sd = 1
      if (m==3) sd = 0.1
      if (m==4) sd = 0.5
      ratio.ab = exp( rtnorm(1, log(a[t-1,m]/(a[t-1,m]+b[t-1,m])), sd, lower, upper) )
      a[t,m] = (a[t-1,m]+b[t-1,m])*ratio.ab
      b[t,m] = a[t-1,m]+b[t-1,m]-a[t,m]
      if (m==1)
        diff = sum(dbeta(phi[t-1,], a[t,m], b[t,m], log=T))-sum(dbeta(phi[t-1,], a[t-1,m], b[t-1,m], log=T))
      if (m==2)
        diff = sum(dbeta(p1[t-1,n.index[[2]]], a[t,m], b[t,m], log=T))-sum(dbeta(p1[t-1,n.index[[2]]], a[t-1,m], b[t-1,m], log=T))
      if (m==3)
        diff = sum(dbeta(p2[t-1,n.index[[3]]], a[t,m], b[t,m], log=T))-sum(dbeta(p2[t-1,n.index[[3]]], a[t-1,m], b[t-1,m], log=T))
      if (m==4)
        diff = sum(dbeta(theta[t-1,n.index[[1]]], a[t,m], b[t,m], log=T))-sum(dbeta(theta[t-1,n.index[[1]]], a[t-1,m], b[t-1,m], log=T))
      diff = diff + log(a[t,m]) - dtnorm(log(ratio.ab), log(a[t-1,m]/(a[t-1,m]+b[t-1,m])), sd, lower, upper, log=T) -
        log(a[t-1,m]) + dtnorm(log(a[t-1,m]/(a[t-1,m]+b[t-1,m])), log(ratio.ab), sd, lower, upper, log=T)
      if (diff>log(runif(1,0,1))){
        a[t-1,m] = a[t,m]
        b[t-1,m] = b[t,m]
      }
      if (a[t,m]!=a[t-1,m])	reject[length(N)+m] = reject[length(N)+m]+1
      
      lower = log(max((a[t-1,m]+b[t-1,m])/a[t-1,m], (a[t-1,m]+b[t-1,m])/b[t-1,m]))
      upper = log(10^6)
      if (m==1) sd = 0.5
      if (m==2) sd = 1
      if (m==3) sd = 0.5
      if (m==4) sd = 0.5
      sum.ab = exp( rtnorm(1, log(a[t-1,m]+b[t-1,m]), sd, lower, upper) )
      a[t,m] = a[t-1,m]/(a[t-1,m]+b[t-1,m])*sum.ab
      b[t,m] = b[t-1,m]/(a[t-1,m]+b[t-1,m])*sum.ab
      if (m==1)
        diff = sum(dbeta(phi[t-1,], a[t,m], b[t,m], log=T))-sum(dbeta(phi[t-1,], a[t-1,m], b[t-1,m], log=T))
      if (m==2)
        diff = sum(dbeta(p1[t-1,n.index[[2]]], a[t,m], b[t,m], log=T))-sum(dbeta(p1[t-1,n.index[[2]]], a[t-1,m], b[t-1,m], log=T))
      if (m==3)
        diff = sum(dbeta(p2[t-1,n.index[[3]]], a[t,m], b[t,m], log=T))-sum(dbeta(p2[t-1,n.index[[3]]], a[t-1,m], b[t-1,m], log=T))
      if (m==4)
        diff = sum(dbeta(theta[t-1,n.index[[1]]], a[t,m], b[t,m], log=T))-sum(dbeta(theta[t-1,n.index[[1]]], a[t-1,m], b[t-1,m], log=T))
      diff = diff - dtnorm(log(sum.ab), log(a[t-1,m]+b[t-1,m]), sd, lower, upper, log=T) +
        dtnorm(log(a[t-1,m]+b[t-1,m]), log(sum.ab), sd, lower, upper, log=T)
      if (diff<log(runif(1,0,1))){
        a[t,m] = a[t-1,m]
        b[t,m] = b[t-1,m]
        reject[length(N)+4+m] = reject[length(N)+4+m]+1
      }
    }
    # update (mu,sig2) by Gibbs
    if (option.mu==1){
      tau.l = 1/(1/tau02+zl/sig2[t-1])
      mu.l = sum(log(z[index4]/n[t,index4]))/sig2[t-1]*tau.l
      mu[t] = rnorm(1, mu.l, tau.l)
      s2 = sum((log(z[index4]/n[t,index4])-mu[t])^2)
      sig2[t] = rinvgamma(1, (zl+1)/2, (sig02+s2)/2)
    }
    if ((100*t)%%size==0)	
    {
      print(paste(100*t/size, "percents has been done"))
      incProgress(1/100)
    }
  }
  print((proc.time() - ptm)/60)
  
  
  update.n <- function(option.ht, option.mu){
    if (option.ht==1){
      diff = dbinom(n[t,j], N[j], phi[t-1,j], log=T) - dpois(n[t,j], n[t-1,j], log=T) -
        dbinom(n[t-1,j], N[j], phi[t-1,j], log=T) + dpois(n[t-1,j], n[t,j], log=T)
      if (is.element(j, n.index[[1]]))
        diff = diff + dbinom(y[j], n[t,j], theta[t-1,j], log=T) - dbinom(y[j], n[t-1,j], theta[t-1,j], log=T)
      if (is.element(j, n.index[[2]]))
        diff = diff + dbinom(x[j,1], n[t,j], p1[t-1,j], log=T) - dbinom(x[j,1], n[t-1,j], p1[t-1,j], log=T)
      if (is.element(j, n.index[[3]]))
        diff = diff + dbinom(x[j,2], n[t,j], p2[t-1,j], log=T) - dbinom(x[j,2], n[t-1,j], p2[t-1,j], log=T)
      if (is.element(j, index3))
        diff = diff + dnbinom(n[t,j]-r[j], r[j]-1, 1-(1-p1[t-1,j])*(1-p2[t-1,j]), log=T) -
          dnbinom(n[t-1,j]-r[j], r[j]-1, 1-(1-p1[t-1,j])*(1-p2[t-1,j]), log=T)
      if (is.element(j, index4) & option.mu==1)
        diff = diff + dnorm(log(z[j]/n[t,j]), mu[t-1], sqrt(sig2[t-1]), log=T) -
          dnorm(log(z[j]/n[t-1,j]), mu[t-1], sqrt(sig2[t-1]), log=T)
    }
    if (option.ht==0){
      diff = dbinom(n[t,j], N[j], phi[t-1,j], log=T) - dpois(n[t,j], n[t-1,j], log=T) -
        dbinom(n[t-1,j], N[j], phi[t-1,j], log=T) + dpois(n[t-1,j], n[t,j], log=T)
      if (is.element(j, n.index[[1]]))
        diff = diff + dbinom(y[j], n[t,j], theta[t-1], log=T) - dbinom(y[j], n[t-1,j], theta[t-1], log=T)
      if (is.element(j, n.index[[2]]))
        diff = diff + dbinom(x[j,1], n[t,j], p1[t-1], log=T) - dbinom(x[j,1], n[t-1,j], p1[t-1], log=T)
      if (is.element(j, n.index[[3]]))
        diff = diff + dbinom(x[j,2], n[t,j], p2[t-1], log=T) - dbinom(x[j,2], n[t-1,j], p2[t-1], log=T)
      if (is.element(j, index3))
        diff = diff + dnbinom(n[t,j]-r[j], r[j]-1, 1-(1-p1[t-1])*(1-p2[t-1]), log=T) -
          dnbinom(n[t-1,j]-r[j], r[j]-1, 1-(1-p1[t-1])*(1-p2[t-1]), log=T)
      if (is.element(j, index4) & option.mu==1)
        diff = diff + dnorm(log(z[j]/n[t,j]), mu[t-1], sqrt(sig2[t-1]), log=T) -
          dnorm(log(z[j]/n[t-1,j]), mu[t-1], sqrt(sig2[t-1]), log=T)
    }
    return(diff)
  }
  
  
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks=20)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
  }
  panel.pearson <- function(x, y, ...) {
    horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
    vertical <- (par("usr")[3] + par("usr")[4]) / 2;
    text(horizontal, vertical, format(cor(x,y), digits=2), cex=3)
  }
  
  keep = seq(burnin, size, thin)
  N.r = DATA$PopSize[-which(DATA$NASROB>-1 | DATA$RSA>-1)]
  phi.r = n.r = NULL
  for (i in 1:length(N.r)){
    phi.i = rbeta(length(keep), a[keep,1],b[keep,1])
    n.i = rbinom(length(keep), N.r[i], phi.i)
    phi.r = cbind(phi.r, phi.i)
    n.r = cbind(n.r, n.i)
  }
  
  chosen = which(rowSums(is.na(DATA[,4:8])) < 5)
  notchosen = which(rowSums(is.na(DATA[,4:8])) == 5)
  
  temp = n[keep,]
  
  colnames(temp) =DATA[chosen,1]
  colnames(n.r) = DATA[notchosen,1]
  
  n = cbind(temp, n.r)
  
  colnames(phi) =  DATA[chosen,1]
  
  colnames(theta) =  DATA[chosen,1]
  
  colnames(p1) =  DATA[chosen,1]
  
  colnames(p2) = DATA[chosen,1]
  
  returnObj = list(a = a[keep,],b = b[keep,], mu = mu[keep], n
                   , phi = phi[keep,], theta = theta[keep,], p1 = p1[keep,], p2 = p2[keep,])
  
  
  return(returnObj)
  
}

# Define User Interface 
ui <- fluidPage(
  
  # Application title
  titlePanel("Size Estimation Shiny App"),
  
  # Create a sidbar with places to enter inputs for sizeestima
  sidebarLayout(
    sidebarPanel(
      
      numericInput("iterations", "Number of MCMC Iterations", 100, min = 1, max = NA, step = 1,
                   width = NULL),
      numericInput("burnin", "Number Burn-In Iterations", 1, min = 1, max = NA, step = 1,
                   width = NULL),
      numericInput("thin", "Keep Every _th Sample", 1, min = 1, max = NA, step = 1,
                   width = NULL),
      fileInput("data", "Upload Dataset", accept = c('text/csv', 
                                                     'text/comma-separated-values,text/plain', 
                                                     '.csv')),
      bsButton("start", "Sample Posterior"),
      downloadButton("template.csv", "Download Data Template"),
      downloadButton("posterior.csv", "Download Posterior Sample")
     
    ),
    
    # Show the selected plot
    mainPanel(
      tabsetPanel(
        tabPanel("Plots", selectInput("type", "", 
                                      choices = c("Total Population Size", "Prevalence by District",
                                                  "Participation Rate For Incomplete Count", "Participation Rate for First Listing",
                                                  "Participation Rate for Second Listing","Matrix Plot of Key Parameters"
                                      )), plotOutput("distPlot")),
        
        tabPanel("Tables", selectizeInput("tabletype", "", choices = c("Key Parameter Estimates", "Prevalence by District",
                                                                       "Participation Rate For Incomplete Count", "Participation Rate for First Listing",
                                                                       "Participation Rate for Second Listing")),
                 tableOutput("table")
                 ),
        tabPanel("Help", 
                HTML("<br><center><b>Frequently Asked Questions </b> </center><br>"),
                actionLink("uploadHelp","How do I upload my own data?"),
                HTML("<br>"),
                actionLink("uploadHelp2","What do the column labels in the data template refer to?"),
                HTML("<br>"),
                actionLink("posteriorHelp","I want a plot or table that is not provided by the app. How can I obtain these?"),
                HTML("<br>"),
                actionLink("posteriorHelp2","What do the column labels refer to in the posterior sample?"),
                HTML("<br>"),
                actionLink("parameterHelp", "What are the key parameters?"),
                HTML("<br>")
        ),
        
        tabPanel("References",
                 
                 HTML("<br>[1] Bao, L., Raftery, A. E., & Reddy, A. (2015). Estimating the sizes of populations at risk of HIV infection from multiple data sources using a Bayesian hierarchical model. Statistics and its Interface, 8(2), 125-136.")
                 
        )
        
      )
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
  
  addPopover(session, "iterations", "", "It is recommended to do a test run with a small 
             number of samples (around 100) to make sure your data is working. When you 
             are ready for a full run, use a larger number (probably > 100,000). This is the final count after taking
             other options into account.", placement = "bottom",
             trigger = "hover")
  addPopover(session, "burnin", "", "It takes a while for the chain to reach what is called
             'equilibrium', where the distribution that matches the true posterior. Throwing
             away early samples helps mitigate the problems this causes. This should be small
             percentage of the total number of samples (probably < 10%)", placement = "bottom"
             , trigger = "hover")
  addPopover(session, "thin", "", "A larger number here reduces dependence between consecutive 
             iterations. It also means that it takes longer to get a large sample. 
             If you intend to download the resulting sample to use elsewhere this can 
             decrease download time and the intensity of further calculations. Keep this 
             number relatively small.", placement = "bottom", trigger = "hover")
  
  addPopover(session, "uploadHelp", "How do I upload my own data?", 
             "By default, the data from Bao et al. is being used by the app. However, 
             it is possible to use your own data. To use a data set, the first step is to 
             convert it into a '.csv' file having a specific format. A template can be found by
             clicking 'Download Data Template' in the side bar. The template is editable by most commonly used spreadsheet programs
             and simply contains a single row providing labels to the column. Copy your data into the 
             appropriate columns and upload it using the menu found on the side bar", placement = "bottom", trigger = "hover")
  
  addPopover(session, "uploadHelp2", "What do the column labels in the template refer to?",
             " <b>District:</b> The names of the areas in your data. 
<br><b>PopSize:</b> The size of an appropriate reference class for target population (e.g. adult males for men who have sex with men).
<br><b>Incomplete Count:</b> An incomplete count of the target population from an intervention, survey, etc... 
<br><b>Listing 1, Listing 2, and Overlap:</b> Capture-recapture data or two listings with known overlap. 
<br><b>Estimate</b> An estimate or guestimate for the size of the popultion from another source.
", placement = "bottom", trigger = "hover")
  
  addPopover(session, "posteriorHelp", "I want a plot or table that is not provided by the app. How can I obtain these?",
             "Click 'Download Posterior Sample' on the side bar and save the '.csv' file. You can then load this
             into your statistics software of choice and create any additional plots you need.", placement = "bottom", trigger = "hover")
 
  addPopover(session, "posteriorHelp2", "What do the column labels refer to in the posterior sample?",
             "<b>a.i, b.i:</b> a.i / (a.i + b.i) is the mean of phi if i = 1,  p1 if i = 2,  p2 if i = 3, and theta if i = 4.  <br>
             <b>mu:</b> The bias of the estimate given in the data. <br>
             <b>n.district_name:</b> The size of the target population in the named district. <br>
             <b>phi.district_name:</b> The prevalence of the target population in the named district <br>
             <b>theta.district_name:</b> The participation rate for the incomplete count in the named district  <br>
             <b>pi.district_name:</b> The participation rate for the ith listing in the named district  <br>  ", placement = "bottom", trigger = "hover")
  
   addPopover(session, "parameterHelp", "What are the key parameters?", 
             "<b>n.total:</b> The size of the target population. <br> 
              <b>E(phi):</b> The expected prevalence of the target population.<br> 
              <b>mu:</b> The bias of the estimate given in the data. <br> 
              <b>E(p1):</b> The expected participation rate in the first listing<br> 
              <b>E(p2):</b> The expected participation rate in the second listing<br> 
              <b>E(theta):</b> The expected participation rate in the incomplete count <br>
             ", placement = "bottom", trigger = "hover")
  
  
  
  output$template.csv <- downloadHandler("template.csv",  function(file){file.copy("./temp.csv", file)} , contentType = "text/csv")
  output$manual.pdf <- downloadHandler("manual.pdf",  function(file){file.copy("./man.pdf", file)} , contentType = "text/csv")
  output$manual.html <- downloadHandler("manual.html",  function(file){file.copy("./man.html", file)} , contentType = "text/csv")
  
  dat <- reactive({
    input$start
    
    isolate({
      
      if(is.null(input$data))
      {
        d = DATA
      }
      
      else
      {
        uploaded = read.csv(input$data$datapath)
        uploaded = data.frame(uploaded[,1], rep(NA, length(uploaded[,1]) ),  uploaded[,2], uploaded[,3], uploaded[,4], uploaded[,5], uploaded[,6], uploaded[,7] )
        names(uploaded) =  names(DATA[1:8])
        
        d = uploaded
      }
      
      
      validate(
        need(names(d)[1:8] == names(DATA)[1:8], message =  "Improper column labels in uploaded data"),
        need(length(d$District) > 1, message = "Too few districts in uploaded data." ),
        need(is.integer(d$PopSize), message = "Each 'PopSize' entry should be a positive integer in the uploaded data "),
        need(is.integer(d$NASROB) | is.na(d$NASROB), message = "Each 'Incomplete Count' entry should be a positive integer in the uploaded data "),
        need(is.integer(d$BSS) | is.na(d$BSS), message = "Each 'Listing 1' entry should be a positive integer in the uploaded data "),
        need(is.integer(d$NEP) | is.na(d$NEP), message = "Each 'Listing 2' entry should be a positive integer in the uploaded data "),
        need(is.integer(d$Both) | is.na(d$Both), message = "Each 'Overlap' entry should be a positive integer in the uploaded data "),
        need(is.integer(d$RSA) | is.na(d$RSA), message = "Each 'Estimate' entry should be a positive integer in the uploaded data ")
      )
      
      return(d)
    })
    
  })
  
  keyParameters <- reactive({
    
    a = posterior()[[1]]
    b = posterior()[[2]]
    mu = posterior()[[3]]
    
    n.total = rowSums(posterior()[[4]], na.rm = TRUE)
    ephi = a[,1] / (a[,1] + b[,1])
    ep1 =  a[,2] / (a[,2] + b[,2])
    ep2 =  a[,3] / (a[,3] + b[,3])
    etheta =  a[,4] / (a[,4] + b[,4])
    
    mat = cbind(n.total, ephi, mu, ep1, ep2, etheta)
    colnames(mat) = c("n.total", "E(phi)", "mu", "E(p1)", "E(p2)", "E(theta)" )
    
    return(mat)
  })
  
  posterior <- reactive({
    withProgress(message = 'Sampling From Posterior',
                 detail = 'This may take a while...', value = 0,{
                   
                   input$start
                   
                   isolate({
                     
                     validate(
                       need(is.integer(input$iterations), message = "'Number of MCMC Iterations' should be an integer"),
                       need(is.integer(input$burnin), message = "'Number of Burn-In Iterations' should be an interger"),
                       need(is.integer(input$thin), message = "'Keep Every _th Sample' should be an integer")
                     )
                     
                     validate(
                       need(input$iterations > 0, message = "'Number of MCMC Iterations' must be positive"),
                       need(input$burnin > 0, message = "'Number of Burn-In Iterations' should be positive"),
                       need(input$thin > 0, message = "'Keep Every _th Sample' should be positive")
                     )
                     
                     post = sizeestima(dat(), input$iterations*input$thin + input$burnin, input$burnin, input$thin)
                     colnames(post[[1]])
                     
                     return(post)
                   })
                   
                 })})
  
  output$posterior.csv <- downloadHandler("posterior.csv", function(file){write.csv(posterior(), file)}  , contentType = "text/csv")
  
  output$distPlot <- renderPlot({
    
    if(input$type == "Total Population Size")
    {
      populationSize = rowSums(posterior()[[4]], na.rm = TRUE)
      par(mfrow = c(1, 2))
      hist(populationSize)
      boxplot(populationSize)
    }
    if(input$type == "Prevalence by District")
    {
      phi = posterior()[[5]]
      
      boxplot(phi, outline = FALSE, las = 2)
      
    }
    if(input$type == "Participation Rate For Incomplete Count")
    {
      rate = posterior()[[6]]
      
      keep = which(rowSums(is.na(dat()[,4:8]) ) < 5)
      
      colnames(rate) = dat()$District[keep]
      
      keep = which(!is.na(rate[1,]))
      
      boxplot(rate[,keep], outline = FALSE, las = 2)
      
    }
    if(input$type == "Participation Rate for First Listing")
    {
      rate = posterior()[[7]]
      
      keep = which(rowSums(is.na(dat()[,4:8]) ) < 5)
      
      colnames(rate) = dat()$District[keep]
      
      keep = which(!is.na(rate[1,]))
      
      boxplot(rate[,keep], outline = FALSE, las = 2)
      
    }
    if(input$type == "Participation Rate for Second Listing")
    {
      rate = posterior()[[8]]
      
      keep = which(rowSums(is.na(dat()[,4:8]) ) < 5)
      
      colnames(rate) = dat()$District[keep]
      
      keep = which(!is.na(rate[1,]))
      
      boxplot(rate[,keep], outline = FALSE, las = 2)
      
    }
    if(input$type == "Matrix Plot of Key Parameters")
    {
      
      mat = keyParameters()
      
      ggpairs(data.frame(mat), diag = list(continuous = "barDiag", discrete = "barDiag"),
              upper = list(continuous = customPlot), lower = list(continuous = "cor"),  ) +
        theme(axis.text.x = element_text(angle = 270))
      
      
    }
    
  })
  
  output$table <- renderTable({
    
    if(input$tabletype == "Key Parameter Estimates")
    {
      
      mat = keyParameters()
      
      Mean = colMeans(mat)
      Quants = apply(mat, MARGIN = 2, quantile, probs = c(.025, .5, .975))
      
      Estimates = data.frame(cbind(Mean, t(Quants) ))
      
      colnames(Estimates) = c("Mean", ".025 Quantile", "Median" , ".975 Quantile")
      rownames(Estimates) = colnames(mat)
      
      
      return(Estimates)
      
      
    }
    
    if(input$tabletype == "Prevalence by District")
    {
      phi = posterior()[[5]]
  
      Mean = colMeans(phi, na.rm = TRUE)
      Quants = apply(phi, MARGIN = 2, quantile, probs = c(.025, .5, .975), na.rm = TRUE)
      
      Estimates = data.frame(cbind(Mean, t(Quants) ))
      
      colnames(Estimates) = c("Mean", ".025 Quantile", "Median" , ".975 Quantile")
      rownames(Estimates) = colnames(phi)
      
      return(Estimates)
      
    }
    if(input$tabletype == "Participation Rate For Incomplete Count")
    {
      rate = posterior()[[6]]
      
      Mean = colMeans(rate, na.rm = TRUE)
      Quants = apply(rate, MARGIN = 2, quantile, probs = c(.025, .5, .975), na.rm = TRUE)
      
      Estimates = data.frame(cbind(Mean, t(Quants) ))
      
      colnames(Estimates) = c("Mean", ".025 Quantile", "Median" , ".975 Quantile")
      rownames(Estimates) = colnames(rate)
      
      keep = which(!is.na(rate[1,]))
      
      return(Estimates[keep,])

    }
    if(input$tabletype == "Participation Rate for First Listing")
    {
      rate = posterior()[[7]]
      
      Mean = colMeans(rate, na.rm = TRUE)
      Quants = apply(rate, MARGIN = 2, quantile, probs = c(.025, .5, .975), na.rm = TRUE)
      
      Estimates = data.frame(cbind(Mean, t(Quants) ))
      
      colnames(Estimates) = c("Mean", ".025 Quantile", "Median" , ".975 Quantile")
      rownames(Estimates) = colnames(rate)
      
      keep = which(!is.na(rate[1,]))
      
      return(Estimates[keep,])

    }
    if(input$tabletype == "Participation Rate for Second Listing")
    {
      rate = posterior()[[8]]
      
      Mean = colMeans(rate, na.rm = TRUE)
      Quants = apply(rate, MARGIN = 2, quantile, probs = c(.025, .5, .975), na.rm = TRUE)
      
      Estimates = data.frame(cbind(Mean, t(Quants) ))
      
      colnames(Estimates) = c("Mean", ".025 Quantile", "Median" , ".975 Quantile")
      rownames(Estimates) = colnames(rate)
      
      keep = which(!is.na(rate[1,]))
      
      return(Estimates[keep,])

    }
    
  }, rownames = TRUE, digit = 4)
  
}



# Run the application 
shinyApp(ui = ui, server = server)
