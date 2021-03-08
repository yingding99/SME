# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# jch 2014.Sep.23
# hml 2014.Oct.06 prevalence of marker negative group will be updated
# hml 2014.Oct.08 remove error messages (Have to reopen the app if you want to run another dataset)
# hml 2014.Oct.27 put in a Button upfront for the user to choose, either Difference, or Ratio
# hml 2014.Oct.28 allow for categorical "Treatment" and "Marker" variables
# hml 2014.Nov.13 continuous Biomarker
# jch 2014.Nov.14 changed prompts appropriate for cintinuout biomarker, changed PFS to time
# jch 2014.Nov.18 changed text "Surivival Time" to "Time-to-Event"
# hml 2014.Dec.04 speparate C0 and C1 into two slider bars
# jch 2014.Dec.07 minor text edits
# hml 2015.Feb.16 added download button
# hml 2015.Apr.15 added condition "length(levels(factor(analysisData()[, input$markerVarSelect])))<2" to avoid unnecessary error message
# jch 2015.Aug.23 added na.rm=TRUE and na.last=NA in slider codes to handle missing values, changed length to sum in updating slider

library(shiny)

list.of.packages <- c("eha", "survival", "rootSolve", "mvtnorm", "plotrix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(eha)
library(survival)
library(rootSolve)
library(mvtnorm)
library(plotrix)

source('Weibull_Median_RatDiff_2grp.R', local=TRUE)
source('App_RatDiffPlot.R', local=TRUE)

shinyServer(function(input, output, clientData, session) {

  ## ---- Input data set as reactive object ----
  inputData <- reactive({
    
    validate(
      need(input$file1 != "", "Please upload an input data set")
    )
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)

    input_data <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
    input_data
  })
  
  ## ---- variable names as a reactive object ----
  varNames <- reactive({
    names(inputData())
  })
  
  ## ---- Patient ID UI ----
  output$patVar <- renderUI({
    varlist <- c("", varNames())
    selectInput("patVarSelect", "Patient ID Variable", varlist)
  })
  
  ## ---- Treatment variable UI ----
  output$trtVar <- renderUI({
    varlist <- c("", varNames())
    selectInput("trtVarSelect", "Treatment Variable", varlist)
  })
  
  ## ---- Treatment levels reactive object ----
  trt_levels <- reactive({
    my_levels <- sort(unique(inputData()[, input$trtVarSelect]))
  })
  
  ## ---- Control Value UI ----
  output$controlValue <- renderUI({
    validate(
      need(input$trtVarSelect != "", "Please select Treatment Variable")
    )
    trt_levels <- c("", c("", sort(unique(inputData()[, input$trtVarSelect]))))
    selectInput("controlValueSelect", "Control Value", trt_levels)
  })

  ## ---- Biomarker Variable UI ----
  output$markerVar <- renderUI({
    varlist <- c("", varNames())
    selectInput("markerVarSelect", "Biomarker Variable", varlist)
  })
  
  ##-------------------------------##
  ## ---- make dynamic slider ---- ##
  ##-------------------------------##
  output$sliderC0 <- renderUI({
    validate(
      need(input$markerVarSelect != "", "Please select a Biomarker Variable")
    )
       
    sliderInput("inSliderC0", "Data filter cut-point c0 (subjects w/ value < C0 are excluded)", 
                min=min(inputData()[, input$markerVarSelect],na.rm=TRUE), max=as.numeric(quantile(inputData()[, input$markerVarSelect], 0.50,na.rm=TRUE)),
                value=min(inputData()[, input$markerVarSelect],na.rm=TRUE)) # na.rm=TRUE for missing values
  
  })
  
  output$sliderC1 <- renderUI({
    validate(
      need(input$markerVarSelect != "", "Please select a Biomarker Variable")
    )
    
    trt<-inputData()[inputData()[, input$markerVarSelect]>=input$inSliderC0, input$trtVarSelect]
    Marker<-inputData()[inputData()[, input$markerVarSelect]>=input$inSliderC0, input$markerVarSelect]
    CMarker<-Marker[trt==input$controlValueSelect]
    RxMarker<-Marker[trt!=input$controlValueSelect]
    ind.C.max<-order(-CMarker,na.last=NA)[2] # na.last=NA for missing values
    ind.Rx.max<-order(-RxMarker,na.last=NA)[2]
    ind.C.min<-order(CMarker,na.last=NA)[2]
    ind.Rx.min<-order(RxMarker,na.last=NA)[2]   
    
    
    if(!is.na(max(CMarker[ind.C.min], RxMarker[ind.Rx.min]))){
    
      if(CMarker[order(CMarker,na.last=NA)[2]]==CMarker[order(CMarker,na.last=NA)[3]]){
        min.inSliderC1<-max(CMarker[ind.C.min], RxMarker[ind.Rx.min],na.rm=TRUE)+0.01
        max.inSliderC1<-min(CMarker[ind.C.max], RxMarker[ind.Rx.max],na.rm=TRUE)
      }else{
        min.inSliderC1<-max(CMarker[ind.C.min], RxMarker[ind.Rx.min],na.rm=TRUE)
        max.inSliderC1<-min(CMarker[ind.C.max], RxMarker[ind.Rx.max],na.rm=TRUE)
      }
  
    sliderInput("inSliderC1", "Marker +/- cut-point c1 (subjects w/ value < C1 are marker -)", 
                min=min.inSliderC1, max=max.inSliderC1,
                value=median(c(CMarker, RxMarker)))
    
    }
    
  })
  
  # ---- Time to event Variable UI ----
  output$timeVar <- renderUI({
    varlist <- c("", varNames())
    selectInput("timeVarSelect", "Time to Event Outcome Variable", varlist)
  })
  
  ## ---- censoring variable UI ----
  output$censorVar <- renderUI({
    varlist <- c("", varNames())
    selectInput("censorVarSelect", "Censor Indicator Variable", varlist)
  })
  
  ## ---- censoring value UI ----
  output$censorVal <- renderUI({
    validate(
      need(input$censorVarSelect != "", "Please select Censor indicator Variable")
    )
    censor_levels <- c("", sort(unique(inputData()[, input$censorVarSelect])))
    selectInput("censorValueSelect", "Censor Value", censor_levels)
  })

  ## ---- analysis data reactive object ----
  analysisData <- reactive({
    if (is.null(input$patVarSelect)|  
          is.null(input$trtVarSelect)|  
          is.null(input$controlValueSelect)|
          is.null(input$markerVarSelect)|
          #is.null(input$M_ValueSelect)|
          is.null(input$timeVarSelect)|
          is.null(input$censorVarSelect)|
          is.null(input$censorValueSelect)) { return() }
    tmp <- inputData()
    # assemble all variable names selected from various inputs
    cols_select <- c(input$patVarSelect, 
                     input$trtVarSelect,
                     input$markerVarSelect,
                     input$timeVarSelect,
                     input$censorVarSelect)
    
    analysis_set <- tmp[, cols_select]
    if(class(analysis_set[, input$trtVarSelect]) != "character") {
      analysis_set[, input$trtVarSelect] <- as.character(analysis_set[, input$trtVarSelect])
    }
    
    ## --------------------------------------------##
    ## -- remove data with Marker value below C0 --##
    ## --------------------------------------------##
    
    if(sum(inputData()[, input$markerVarSelect]<input$inSliderC0,na.rm=TRUE)>0){
      analysis_set<- analysis_set[inputData()[, input$markerVarSelect]>=input$inSliderC0,]
    }
    
    analysis_set[, input$markerVarSelect]<-factor(ifelse(analysis_set[, input$markerVarSelect] < input$inSliderC1, "Neg", "Pos"))
    
    analysis_set
  })


  ## ---- update slider bar by SAMPLE prevalence of marker negative subgroup
  observe({
    if (is.null(input$patVarSelect)|  
          is.null(input$trtVarSelect)|  
          is.null(input$controlValueSelect)|
          is.null(input$markerVarSelect)|
          #is.null(input$M_ValueSelect)|
          is.null(input$timeVarSelect)|
          is.null(input$censorVarSelect)|
          is.null(input$censorValueSelect)) { return() }
    
    prev.neg<-sum(analysisData()[, input$markerVarSelect]=="Neg",na.rm=TRUE)/sum(!is.na(analysisData()[, input$markerVarSelect]=="Neg"),na.rm=TRUE)*100
    updateSliderInput(session, "p.neg",  value = prev.neg) # length replaced by sum(!is.na)
  })
  

  ## ---- model fit reactive object ----
  WeibullFit <- reactive({
        if (length(unique(analysisData()[, input$trtVarSelect]))>10 |
              is.null(input$patVarSelect)|  
              is.null(input$trtVarSelect)|  
              is.null(input$controlValueSelect)|
              is.null(input$markerVarSelect)|
              length(levels(factor(analysisData()[, input$markerVarSelect])))<2|
              #is.null(input$M_ValueSelect)|
              is.null(input$timeVarSelect)|
              is.null(input$censorVarSelect)|
              is.null(input$censorValueSelect)) { return() }
        
        # censor indicator (0=alive, 1=dead) ( TRUE/FALSE w/ TRUE = death) (1/2 w/ 2=death)
        
        outcomeEventInd <- analysisData()[, input$censorVarSelect]
        outcomeEventInd[outcomeEventInd == input$censorValueSelect] <- 0
        outcomeEventInd[outcomeEventInd != input$censorValueSelect] <- 1    
        
        # define survival object
        response_vec <- Surv(time = as.numeric(analysisData()[, input$timeVarSelect]), 
                             event = as.numeric(outcomeEventInd))
    
        Trt <- analysisData()[, input$trtVarSelect]
    
        M <- analysisData()[, input$markerVarSelect]
        
    
        # fit Weibull model that reacts to event in user interface 
        fit <- weib(formula = response_vec ~ 
                    relevel(as.factor(Trt),input$controlValueSelect)+
                    relevel(as.factor(M),"Neg")+
                    relevel(as.factor(Trt),input$controlValueSelect)*relevel(as.factor(M),"Neg"), 
                  data=analysisData()) 

        # log.param are the log scaled parameter estimates from weibull fitting
        log.param.est <- fit$log.coef
        log.param.var <- diag(fit$log.var)
        # param are the original scale (exponential of log.param) parameter estimates
        param.est <- fit$coef
        param.var <- fit$var
        param.var.diag <- diag(param.var)
        
        # get prevalence of marker negative group from slider 
        p.neg = as.double(input$p.neg)/100
        
        #-------------------------------------------------------------------------------#
        # compute median survival times and their associated variance covariance matrix #
        #-------------------------------------------------------------------------------#
        comp.median <- median.surv(lambda=param.est[1],k=param.est[2],theta=param.est[3:5],p.neg=p.neg)
        (median.surv.est <- unlist(comp.median))
    
        median.surv.C <- median.surv.est[c(1,2,5)]
        median.surv.Rx <- median.surv.est[c(3,4,6)]
    
        median.surv.estiMates <- matrix(ncol=3, rbind(median.surv.Rx, median.surv.C))
        Negative <- median.surv.estiMates[,1]
        Positive <- median.surv.estiMates[,2]
        Mixture <- median.surv.estiMates[,3]
        Estimates <- data.frame(Negative, Positive, Mixture)
        row.names(Estimates) <- c("Rx", "C")
    
        median.var.est <-var.median.surv.dt.dim6(med.C.neg=median.surv.est[1], med.C.pos=median.surv.est[2], 
                                             med.Rx.neg=median.surv.est[3], med.Rx.pos=median.surv.est[4],
                                             med.C=median.surv.est[5],med.Rx=median.surv.est[6],p.neg=p.neg,
                                             lambda=param.est[1],k=param.est[2],theta=param.est[3:5],
                                             param.var=fit$var)$median.var

  
        
        #----------------------------------------#
        # Calculate ratio of median and their CI #
        #----------------------------------------#  
        #----------------------------------------------------------#
        # begin log scale delta method
        comp.ratio <- log.ratio.median.surv.dt.dim3(median.surv=median.surv.est,median.var=median.var.est)
        log.ratio <- c(comp.ratio$log.ratio.neg,comp.ratio$log.ratio.pos,comp.ratio$log.ratio.mix)
        var.log.ratio <- comp.ratio$var.log.ratio
        
        original.ratio <- c(comp.ratio$ratio.neg,comp.ratio$ratio.pos,comp.ratio$ratio.mix)
        ratio.negative <- original.ratio[1]
        ratio.positive <- original.ratio[2]
        ratio.mixture  <- original.ratio[3]
        Ratios <- data.frame(ratio.negative, ratio.positive, ratio.mixture)
        row.names(Ratios) <- c("Rx/C")
        colnames(Ratios) <- c("Negative", "Positive", "Mixture")
        
        comp.ratio.CI <- log.ratio.simu.CI(seed=12345,alpha=as.double(input$alpha)/100,log.ratio=log.ratio,var.log.ratio=var.log.ratio)
        ratio.CI.lbd <- comp.ratio.CI$ratio.CI.lbd
        ratio.CI.ubd <- comp.ratio.CI$ratio.CI.ubd
        
        # end of log delta method
        #----------------------------------------------------------#
        # CIs from lg scale delta method
        ratio.CI <- matrix(ncol=3, rbind(ratio.CI.ubd, ratio.CI.lbd))
        
        ratio.CI.negative <- ratio.CI[,1]
        ratio.CI.positive <- ratio.CI[,2]
        ratio.CI.mixture <- ratio.CI[,3]
        Ratio.CI <- data.frame(ratio.CI.negative, ratio.CI.positive, ratio.CI.mixture)
        row.names(Ratio.CI) <- c("Upper","Lower")
        colnames(Ratio.CI) <- c("Negative", "Positive", "Mixture")
        
    
    
        #---------------------------------------------#
        # Calculate difference of median and their CI #
        #---------------------------------------------#
        #----------------------------------------------------------#
        # begin of delta method
        comp.dif <- dif.median.surv.dt.dim3(median.surv=median.surv.est,median.var=median.var.est)
        dif.median <- c(comp.dif$dif.neg, comp.dif$dif.pos,comp.dif$dif.mix)
        var.dif.median <- comp.dif$var.dif
        
        dif.median.negative <- dif.median[1]
        dif.median.positive <- dif.median[2]
        dif.median.mixture <- dif.median[3]
        Diff <- data.frame(dif.median.negative, dif.median.positive, dif.median.mixture)
        row.names(Diff) <- c("Rx-C")
        colnames(Diff) <- c("Negative", "Positive", "Mixture")
        
        comp.dif.CI <- dif.median.simu.CI(seed=12345,alpha=as.double(input$alpha)/100, dif.median=dif.median, var.dif.median=var.dif.median)
        median.dif.CI.lbd <- comp.dif.CI$dif.median.CI.lbd
        median.dif.CI.ubd <- comp.dif.CI$dif.median.CI.ubd
    
        # end of delta method
        #----------------------------------------------------------#
        # CIs from delta method
        median.dif.CI <- matrix(ncol=3, rbind(median.dif.CI.ubd, median.dif.CI.lbd))
        
        median.dif.CI.negative <- median.dif.CI[,1]
        median.dif.CI.positive <- median.dif.CI[,2]
        median.dif.CI.mixture <- median.dif.CI[,3]
        Diff.CI <- data.frame(median.dif.CI.negative, median.dif.CI.positive, median.dif.CI.mixture)
        row.names(Diff.CI) <- c("Upper","Lower")
        colnames(Diff.CI) <- c("Negative", "Positive", "Mixture")
        
        #---------------#
        # return retuls #
        #---------------#    
        WeibullResult <- list(Estimates = Estimates,
                              Ratios = Ratios,
                              Ratio.CI = Ratio.CI,
                              Diff = Diff,
                              Diff.CI = Diff.CI,
                              median.surv.estiMates = median.surv.estiMates)
        return(WeibullResult)
  })


  ##---------------------------##
  ## ---- display results ---- ##
  ##---------------------------##
  output$Estimates <- renderPrint({
      if (length(unique(analysisData()[, input$trtVarSelect]))>10 |
            is.null(input$patVarSelect)|  
            is.null(input$trtVarSelect)|  
            is.null(input$controlValueSelect)|
            is.null(input$markerVarSelect)|
            length(levels(factor(analysisData()[, input$markerVarSelect])))<2|
            is.null(input$timeVarSelect)|
            is.null(input$censorVarSelect)|
            is.null(input$censorValueSelect)) { return() }
  
      round(WeibullFit()$Estimates,3)
  
  })
  
  output$text1 <- renderText({ 
    paste("Estimated", input$method, "Time-to-Event")
  })
  
  output$Ratios <- renderPrint({
      if (length(unique(analysisData()[, input$trtVarSelect]))>10 |
            is.null(input$patVarSelect)|  
            is.null(input$trtVarSelect)|  
            is.null(input$controlValueSelect)|
            is.null(input$markerVarSelect)|
            length(levels(factor(analysisData()[, input$markerVarSelect])))<2|
            is.null(input$timeVarSelect)|
            is.null(input$censorVarSelect)|
            is.null(input$censorValueSelect)) { return() }
      if(input$method=="Difference of Median"){
        round(WeibullFit()$Diff, 3)
      } else {
        round(WeibullFit()$Ratios,3)
      }
  })
  
  output$text2 <- renderText({ 
    paste("Confidence Interval for", input$method, "Time-to-Event")
  })
  
  output$logCI <- renderPrint({
      if (length(unique(analysisData()[, input$trtVarSelect]))>10 |
            is.null(input$patVarSelect)|  
            is.null(input$trtVarSelect)|  
            is.null(input$controlValueSelect)|
            is.null(input$markerVarSelect)|
            length(levels(factor(analysisData()[, input$markerVarSelect])))<2|
            is.null(input$timeVarSelect)|
            is.null(input$censorVarSelect)|
            is.null(input$censorValueSelect)) { return() }
      if(input$method=="Difference of Median"){
        round(WeibullFit()$Diff.CI,3)
      } else {
        round(WeibullFit()$Ratio.CI,3)
      }
    
  })
  
  Results <- reactive({
    if (is.null(input$patVarSelect)|  
          is.null(input$trtVarSelect)|  
          is.null(input$controlValueSelect)|
          is.null(input$markerVarSelect)|
          length(levels(factor(analysisData()[, input$markerVarSelect])))<2|
          #is.null(input$M_ValueSelect)|
          is.null(input$timeVarSelect)|
          is.null(input$censorVarSelect)|
          is.null(input$censorValueSelect)) { return() }
    
    if(input$method=="Difference of Median"){
      Est.Eff<-round(WeibullFit()$Diff, 3)
      CI.Eff<-round(WeibullFit()$Diff.CI,3)
    } else {
      Est.Eff<-round(WeibullFit()$Ratios,3)
      CI.Eff<-round(WeibullFit()$Ratio.CI,3)
    }
    
    tmp <- list(Estimated.Median.Survival=round(WeibullFit()$Estimates,3), 
                Estimated.Efficacy=Est.Eff,
                Confidence.Interval.for.Efficacy=CI.Eff)
  })
  
  
  output$downloadTables <- downloadHandler(
    filename = function() {
      paste('Results-', input$method, '.csv', sep="")
    },
    content = function(file) {
      write.table(t("Estimated Median Time-to-Event"), file, col.names=FALSE, row.names=FALSE, sep=",")
      write.table(t(c("", "Negative","Positive", "Mixture")), file, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
      write.table(Results()$Estimated.Median.Survival, file, col.names=FALSE, sep=",", append=TRUE)
      write.table("", file, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
      write.table(paste("Estimated", input$method, "Time-to-Event"), file, col.names=FALSE, row.names=FALSE,sep=",", append=TRUE)
      write.table(t(c("", "Negative","Positive", "Mixture")), file, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
      write.table(Results()$Estimated.Efficacy, file, col.names=FALSE, sep=",", append=TRUE)
      write.table("", file, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
      write.table(paste("Confidence Interval for", input$method, "Time-to-Event"), file, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
      write.table(t(c("", "Negative","Positive", "Mixture")), file, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
      write.table(Results()$Confidence.Interval.for.Efficacy, file, col.names=FALSE, sep=",", append=TRUE)
    }
  )
  
    
  output$logRatDiff <- renderPlot({
      if (length(unique(analysisData()[, input$trtVarSelect]))>10 |
            is.null(input$patVarSelect)|  
            is.null(input$trtVarSelect)|  
            is.null(input$controlValueSelect)|
            is.null(input$markerVarSelect)|
            length(levels(factor(analysisData()[, input$markerVarSelect])))<2|
            is.null(input$timeVarSelect)|
            is.null(input$censorVarSelect)|
            is.null(input$censorValueSelect)) { return() }

    # Generate the M&M plots
    par(pty="s") 
    # plot the ratio of median M&M plot (
    plotLimit <- max(WeibullFit()$median.surv.estiMates)*1.3
    
    if(input$method=="Ratio of Median"){
      RatPlot(xlim=c(0,plotLimit),ylim=c(0,plotLimit), # range of the data
              xlab=expression(paste(nu[C]," = median time for C")),
              ylab=expression(paste(nu[Rx]," = median time for Rx")), # label of two axes
              med.surv.Rx=WeibullFit()$median.surv.estiMates[1,], # median survival time for the Rx (sub)groups
              med.surv.C=WeibullFit()$median.surv.estiMates[2,], # median survival time for the C (sub)groups
              # radius=c(sqrt(20),sqrt(25),sqrt(30)), # the length of the arcs
              angle.1=as.numeric(WeibullFit()$Ratio.CI[2,]), # angle.1 (a vector) = the upper bounds of simultaneous CIs for the ratios
              angle.2=as.numeric(WeibullFit()$Ratio.CI[1,]), # angle.2 (a vector) = the lower bounds of simultaneous CIs for the ratios
              color=c(2,3,4), # the color of the lines/arcs
              pch=c(21,22,23), # the shape of points (to represent the median survival times)
              legend=c(paste(":    ", "g-"), 
                       paste(":    ", "g+"),
                       paste(": {", "g-", "," , "g+", "}")), # the legend of the plot
              axis.1=c(expression(paste(nu[C],"-")),expression(paste(nu[C],"+")), expression(paste(nu[C],"*"))),
              axis.2=c(expression(paste(nu[Rx],"-")),expression(paste(nu[Rx],"+")), expression(paste(nu[Rx],"*")))# labels of the two axes
      )
    }
    
    if(input$method=="Difference of Median"){
      DifPlot(xlim=c(0,plotLimit),ylim=c(0,plotLimit), # range of the data
              xlab=expression(paste(nu[C]," = median time for C")),
              ylab=expression(paste(nu[Rx]," = median time for Rx")), # label of two axes
              med.surv.Rx=WeibullFit()$median.surv.estiMates[1,], # median survival time for the Rx (sub)groups
              med.surv.C=WeibullFit()$median.surv.estiMates[2,], # median survival time for the C (sub)groups
              UCI=WeibullFit()$Diff.CI[2,], # UCI (a vector) = the upper bounds of simultaneous CIs for the differences
              LCI=WeibullFit()$Diff.CI[1,], # LCI (a vector) = the lower bounds of simultaneous CIs for the differences
              color=c(2,3,4), # the color of the lines/arcs
              pch=c(21,22,23), # the shape of points (to represent the median survival times)
              legend=c(paste(":    ", "g-"), 
                       paste(":    ", "g+"),
                       paste(": {", "g-", "," , "g+", "}")), # the legend of the plot
              axis.1=c(expression(paste(nu[C],"-")),expression(paste(nu[C],"+")), expression(paste(nu[C],"*"))),
              axis.2=c(expression(paste(nu[Rx],"-")),expression(paste(nu[Rx],"+")), expression(paste(nu[Rx],"*")))# labels of the two axes
      )
    }
  })
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() { 
      paste('MMplot-', input$method, '.pdf', sep="")
    },
    content = function(file) {
      pdf(file)
      # Generate the M&M plots
      par(pty="s") 
      # plot the ratio of median M&M plot (
      plotLimit <- max(WeibullFit()$median.surv.estiMates)*1.3
      
      if(input$method=="Ratio of Median"){
        RatPlot(xlim=c(0,plotLimit),ylim=c(0,plotLimit), # range of the data
                xlab=expression(paste(nu[C]," = median time for C")),
                ylab=expression(paste(nu[Rx]," = median time for Rx")), # label of two axes
                med.surv.Rx=WeibullFit()$median.surv.estiMates[1,], # median survival time for the Rx (sub)groups
                med.surv.C=WeibullFit()$median.surv.estiMates[2,], # median survival time for the C (sub)groups
                # radius=c(sqrt(20),sqrt(25),sqrt(30)), # the length of the arcs
                angle.1=as.numeric(WeibullFit()$Ratio.CI[2,]), # angle.1 (a vector) = the upper bounds of simultaneous CIs for the ratios
                angle.2=as.numeric(WeibullFit()$Ratio.CI[1,]), # angle.2 (a vector) = the lower bounds of simultaneous CIs for the ratios
                color=c(2,3,4), # the color of the lines/arcs
                pch=c(21,22,23), # the shape of points (to represent the median survival times)
                legend=c(paste(":    ", "g-"), 
                         paste(":    ", "g+"),
                         paste(": {", "g-", "," , "g+", "}")), # the legend of the plot
                axis.1=c(expression(paste(nu[C],"-")),expression(paste(nu[C],"+")), expression(paste(nu[C],"*"))),
                axis.2=c(expression(paste(nu[Rx],"-")),expression(paste(nu[Rx],"+")), expression(paste(nu[Rx],"*")))# labels of the two axes
        )
      }
      
      if(input$method=="Difference of Median"){
        DifPlot(xlim=c(0,plotLimit),ylim=c(0,plotLimit), # range of the data
                xlab=expression(paste(nu[C]," = median time for C")),
                ylab=expression(paste(nu[Rx]," = median time for Rx")), # label of two axes
                med.surv.Rx=WeibullFit()$median.surv.estiMates[1,], # median survival time for the Rx (sub)groups
                med.surv.C=WeibullFit()$median.surv.estiMates[2,], # median survival time for the C (sub)groups
                UCI=WeibullFit()$Diff.CI[2,], # UCI (a vector) = the upper bounds of simultaneous CIs for the differences
                LCI=WeibullFit()$Diff.CI[1,], # LCI (a vector) = the lower bounds of simultaneous CIs for the differences
                color=c(2,3,4), # the color of the lines/arcs
                pch=c(21,22,23), # the shape of points (to represent the median survival times)
                legend=c(paste(":    ", "g-"), 
                         paste(":    ", "g+"),
                         paste(": {", "g-", "," , "g+", "}")), # the legend of the plot
                axis.1=c(expression(paste(nu[C],"-")),expression(paste(nu[C],"+")), expression(paste(nu[C],"*"))),
                axis.2=c(expression(paste(nu[Rx],"-")),expression(paste(nu[Rx],"+")), expression(paste(nu[Rx],"*")))# labels of the two axes
        )
      }
      dev.off()
    }
  )

})



