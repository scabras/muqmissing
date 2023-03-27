#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(mice)
library(BayesVarSel)
library(DT)
library(LaplacesDemon)
library(Matrix)

MC.imputation<- function(Y, X, nMC=500, seed=runif(1,0,1000)){
  #(Works for continuous covariates)
  
  #results:
  rX.imput<- array(0, dim=c(dim(X), nMC))
  rSigma<- array(0, dim=c(dim(X)[2], dim(X)[2], nMC))
  rmu<- matrix(0, nr=dim(X)[2], nc=nMC)
  #####
  
  p<- dim(X)[2]; n<-dim(X)[1]
  if (p<2) stop("It does not work with p<2")
  
  O<- 1*(!is.na(X))
  #the following units have at least one NA
  these<- which(rowSums(O)<p)
  
  #Imputamos el valor inicial usando mice
  imputed<- mice(X)
  X.full<- complete(imputed)
  X.full<-as.matrix(X.full)
  ### Gibbs sampler
  set.seed(seed)
  
  
  #ss0<- sum((Y-mean(Y))^2)
  #BF.our<- rep(NA, nMC)
  cont<- 0
  for(s in 1:nMC)
  {
    
    #posterior dist. with Jeffreys independent prior
    Sigma<- rinvwishart(nu=n-1, S=(n-1)*var(X.full))
    mu<- rmvn(1, colMeans(X.full), Sigma/n)
    
    ###update missing data
    for(i in these)
    { 
      b <- ( O[i,]==0 )
      a <- ( O[i,]==1 )
      solve.Sigma.a.a<- solve(Sigma[a,a])
      thetab.mid.a<- mu[b]+Sigma[b,a]%*%solve.Sigma.a.a%*%(X.full[i,a]-mu[a])
      Sigmab.mid.a<- Sigma[b,b] - as.matrix(forceSymmetric(Sigma[b,a]%*%solve.Sigma.a.a%*%Sigma[a,b]))
      X.full[i,b]<- rmvn(1, as.vector(thetab.mid.a), Sigmab.mid.a)
    }
    
    rX.imput[,,s]<- X.full
    rSigma[,,s]<- Sigma
    rmu[,s]<- mu
    
  }
  
  return(list(rX.imput=rX.imput, rSigma=rSigma, rmu=rmu))
  
}

Bvs.missing<- function(Y, X, imputation.list){
  
  p<- dim(X)[2]
  modelsprob<- matrix(0, nc=p+2, nr=2^p)
  colnames(modelsprob)<- c(colnames(X), "postprob", "BF")
  
  #Null model
  modelsprob[2^p, p+1]<- 1
  for (i in 1:(2^p-1)){
    M1.bin<- BayesVarSel:::integer.base.b_C(i, p)
    modelsprob[i,1:p]<- M1.bin
    BF<- mean(gBF.missing.var.sel(Y=Y, X=X, M1=which(M1.bin==1), imputation.list=imputation.list))
    modelsprob[i,p+2]<- BF
    modelsprob[i,p+1]<- exp(log(BF)-lchoose(p, sum(M1.bin)))
  }	
  
  #renormalize
  modelsprob[,p+1]<- modelsprob[,p+1]/sum(modelsprob[,p+1])
  
  return(modelsprob)		
}

gBF.missing.var.sel<- function(Y, X, M1, imputation.list){
  #X is nxp with missing values coded as NA
  #Y is a n vector without missing
  #M1 is a vector of indexes in {1,2,...,p}
  ss0<- sum((Y-mean(Y))^2)
  n<- length(Y)
  nMC<- dim(imputation.list$rSigma)[3]
  BF.our<- rep(0, nMC)
  
  if (length(M1)>1){
    for(s in 1:dim(imputation.list$rSigma)[3])
    {
      
      #posterior dist. with Jeffreys independent prior
      Sigma<- imputation.list$rSigma[,,s]
      mu<- imputation.list$rmu[,s]
      X.full<- imputation.list$rX.imput[,,s]
      
      #X.center<- t(apply(X.full[,M1], 1, function(x) {x-mu[M1]}))
      X.center<- apply(X.full[,M1], 2, function(x) {x-mean(x)})
      
      tX.center.X.center<- t(X.center)%*%X.center
      
      Sigma11<- Sigma[M1,M1]
      
      #In BF.our we save the BFs obtained with our proposal, where m_1 depends through pi_1(beta | nu), 
      #in this case it depends on Sigma:  
      BF.our[s]<- exp(-.5*(n-1)*log(1-t(Y)%*%X.center%*%solve((tX.center.X.center+Sigma11))%*%t(X.center)%*%Y/ss0)-
                        .5*determinant(tX.center.X.center%*%solve(Sigma11)+diag(rep(1,length(M1))), log=T)$modulus[1])
    }
  }
  else
    for(s in 1:dim(imputation.list$rSigma)[3])
    {
      
      #posterior dist. with Jeffreys independent prior
      Sigma<- imputation.list$rSigma[,,s]
      mu<- imputation.list$rmu[,s]
      X.full<- imputation.list$rX.imput[,,s]
      
      #X.center<- X.full[,M1]-mu[M1]
      X.center<- X.full[,M1]-mean(X.full[,M1])
      
      tX.center.X.center<- sum(X.center^2)
      Sigma11<- Sigma[M1,M1]
      
      #In BF.our we save the BFs obtained with our proposal, where m_1 depends through pi_1(beta | nu), 
      #in this case it depends on Sigma:  
      BF.our[s]<- exp(-.5*(n-1)*log(1-t(Y)%*%X.center%*%solve((tX.center.X.center+Sigma11))%*%t(X.center)%*%Y/ss0)-
                        .5*determinant(tX.center.X.center%*%solve(Sigma11)+diag(rep(1,length(M1))), log=T)$modulus[1])
    }
  
  
  return(BF.our)
  
}	


# Define UI for application that draws a histogram
ui <- (fluidPage(
    
    tags$head(tags$style(HTML("
    .shiny-text-output {
      background-color:#fff;
    }
  "))),
    h1("Bayesian variable selection with missing data", 
       style = "font-family: 'Source Sans Pro';
        color: #fff; text-align: center;
        background-image: url('texturebg.png');
        padding: 20px"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        p("Bayesian variable selection with missing data."),
        h4("Upload your data."),
        fileInput("dat",label= "",
                  accept =c("txt/csv", "text/comma-separated-values,text/plain", ".csv")),
        p("Upload a CSV/Text file with 2+p columns with the following labels:"),
        tags$li("id: subject label;"), 
        tags$li("Y: response continuous variables;"), 
        tags$li("others: numerical covariates (no factors allowed)."),
        tags$br(),
        a(href="example.csv","You can try our CSV example."),
        tags$br(),  tags$br(),
        numericInput("nMC", "Number of imputations:", 500, min = 10, max = 1000),
        actionButton("run", "run")
      ),
      
      mainPanel(
        h2("RESULTS"),
        br(),
        br(),
        h3("All inclusion probabilities are in this table."),
        p("(Those above 50% are in the Median Probability Model)"),
        fluidRow(DT::dataTableOutput("tab.incprob")),
        downloadButton('dwn.incprob', 'Download inclusion probabilities'),
      )
    )
    
    
    
  ))

  

# Define server logic required to draw a histogram
server <- (function(input, output) {
  
  dat <- reactive({
    dat=read.table(file=input$dat$datapath,sep=";",dec=".",
                   header = TRUE,row.names = "id")
    return(dat)})
  
  res=eventReactive(eventExpr = input$run, {
    dd=dat()
    Y=dd$Y
    X.miss=dd[,colnames(dd)!="Y"]
    p=ncol(X.miss)
    simul<- MC.imputation(Y=Y, X=X.miss, nMC=input$nMC)
    modelsprob<- Bvs.missing(Y=Y, X=X.miss, imputation.list=simul)
    ip.imp=colSums(modelsprob[,1:p]*modelsprob[,p+1])
    names(ip.imp)=colnames(X.miss)
    return(ip.imp)
  })
  
  inclusion.probability <- reactive({
    return(sort(c(round(res(),3)),decreasing = TRUE))
  })
  
  output$tab.incprob <- DT::renderDataTable({
    ii=inclusion.probability()
    dd=data.frame(vv=names(ii),ii=round(ii,4))
    return(datatable(dd,rownames = FALSE,colnames = c("Variable","Inclusion Probability")) %>% 
             formatStyle('ii',target = 'row',backgroundColor = styleInterval(c(0.5,1),c("white","lightgreen","lightgreen"))) 
    )
  })
  
  
  output$dwn.incprob <- downloadHandler(
    filename = "inclusion-probabilities.csv",
    content = function(file) {
      ii=inclusion.probability()
      dd=data.frame(var=names(ii),inc.prob=ii)
      write.table(dd, file,sep=";",dec=".",row.names = FALSE)
    }
  )
  
})


# Run the application 
shinyApp(ui = ui, server = server)
