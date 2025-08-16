 {

      read.cb <- function(header=TRUE,...) read.table("clipboard", header = header,
                                                      sep ="\t",...)
      ##########################
      # Write x to clipbord
      write.cb =
        function(x, row.names=TRUE, col.names=TRUE, comment=FALSE, text=NULL, ...){
          datafile <- file(paste0("clipboard-", object.size(x)), open='wt')
          on.exit(close(datafile))
          if(comment == TRUE)   {
            if(is.null(comment(x))) warning("There is no comment for x! first add one by comment(x) = '...'") else
              writeLines(comment(x), con=datafile)}
          write.table(x, file = datafile, sep = "\t", row.names = row.names,
                      col.names = col.names, ...)
          if(!is.null(text))   {writeLines(text , con=datafile)}
        }

plot.UPG    = function(x         = NULL,
                       ...,
                       sort      = FALSE,           # sort coefficients by average effect size
                       names     = NULL,            # provide names for variables, alternatively
                       groups    = NULL,            # provide names for groups except baseline
                       xlab      = NULL,            # provide x axis label
                       ylab      = NULL,            # provide y axis label
                       q         = c(0.025, 0.975), # credible intervals
                       include   = NULL,              # which variables to include? default:all (numeric vector)
                       OR        = TRUE    ){


  if(is.null(include)) include = 2:ncol(x$inputs$X)
  if(is.null(names))   names = colnames(x$inputs$X[,include,drop=F])
  if(is.null(names))   names = paste0("Variable", 1:ncol(x$inputs$X[,include,drop=F]))
  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")


  #create some global variables such that r cmd check is happy
  variable = iter = value = NULL

  c.point  = apply(x$posterior$beta[,include,,drop=F], c(2,3), mean)
  c.upper  = apply(x$posterior$beta[,include,,drop=F], c(2,3), quantile, q[2])
  c.lower  = apply(x$posterior$beta[,include,,drop=F], c(2,3), quantile, q[1])


  if(OR == TRUE){
    c.point  = apply(exp(x$posterior$beta[,include,,drop=F]), c(2,3), mean)
    c.upper  = apply(exp(x$posterior$beta[,include,,drop=F]), c(2,3), quantile, q[2])
    c.lower  = apply(exp(x$posterior$beta[,include,,drop=F]), c(2,3), quantile, q[1])
  }

  if(nrow(c.point) == 1){
    c.upper = matrix(c.upper,nrow=1)
    c.lower = matrix(c.lower,nrow=1)}

  if(is.null(groups)) groups = paste(x$posterior$groups[-ncol(c.point)])
  if(length(groups) != ncol(c.point)-1) stop("Wrong number of group names supplied. Need K-1 names where K is the number of choices.")


  #kick baseline
  c.upper = c.upper[,-ncol(c.upper),drop=F]
  c.lower = c.lower[,-ncol(c.lower),drop=F]
  c.point = c.point[,-ncol(c.point),drop=F]

  #add variable names
  c.upper = data.frame(c.upper, names)
  c.lower = data.frame(c.lower, names)
  c.point = data.frame(c.point, names)

  #add group names
  colnames(c.upper) = colnames(c.lower) = colnames(c.point) = c(groups, "names")

  #add measurement
  c.upper$measure   = "c.upper"
  c.lower$measure   = "c.lower"
  c.point$measure   = "c.point"

  plot.df = rbind(c.upper,c.lower,c.point)
  plot.df = reshape2::melt(plot.df, id.vars = c("names","measure"))
  plot.df = reshape2::dcast(plot.df, "names + variable ~ measure")

  #sorting (bit more complicated in MNL, I use average point estimate over all groups)
  if(length(names)>1){
    if(sort){

      average       = aggregate(plot.df$c.point, by=list(plot.df$names), FUN=mean)
      lvls          = unique(plot.df$names)
      plot.df$names = factor(plot.df$names, levels = lvls[order(average$x)])
    }
    if(!sort){
      # plot in order of appearance in X
      plot.df$names = factor(plot.df$names, levels = rev(names))
    }
  }
  #axis labeling, take care b/c of coord_flip
  if(is.null(ylab)) ylab = ""
  if(is.null(xlab)) xlab = "Posterior Estimate"

  final =  ggplot(plot.df, aes(x=names, y=c.point, shape = variable, group = variable)) +
    geom_errorbar(aes(x=names, ymin=c.lower, ymax=c.upper, group = variable),
                  position = position_dodge(width=.5),
                  width=0, col="grey60") +
    geom_point(position = position_dodge(width=.5)) +
    theme_minimal()   +
    coord_flip() +
    xlab(ylab)     +
    ylab(xlab) +
    theme(legend.position = "bottom",
          legend.title = element_blank())

  if(!isTRUE(OR))
  {final = final + geom_hline(yintercept = 0, col="red",lty=2)}


  if(isTRUE(OR))
  {final = final +  geom_hline(yintercept = 1, col="red",lty=2) +
    geom_text(mapping = aes (x = 0.6, y = 1),
              label = "OR = 1", color = "red")+
    labs(y = "Posterior Estimate of OR \n(95% Credible interval)")}

  return(final)
}

multinomial = function(data, formula2,
                       ref = NULL,
                       newdata = NULL,
                       ggplot.mapping = NULL ,
                       bayesian = FALSE,
                       OR = TRUE){

  # if(!(exists("%f%", envir = globalenv()) & exists("%+%", envir = globalenv()))){
  # source("https://raw.githubusercontent.com/ahadalizadeh/utility_fun/master/utility_fun.R")
  # }
  require(nnet)
  require(ggplot2)
  require(reshape2)
  formula2 <<- formula2
  data<-data
  summerBayesian <- NULL
  ################# change ref ---------
  if(!is.null(ref)){
    DV =   all.vars(update.formula(formula2, .~1)  )
    # Data$DV = Data[[DV]]
    data[[DV]] <- relevel(as.factor(data[[DV]]), ref = as.character(ref))
    # formula = update.formula(formula, DV2~.)
  }

  ################## run model   ---------
  "%+%" <- function(a, b) {
    if(!(is.atomic(a) | is.atomic(b)| is.matrix(a) | is.matrix(b)))
      stop("a and b must be matrix or atomic")
    if(is.atomic(a) & is.atomic(b)){
      z = paste0(a, b)
    }
    if(is.matrix(a) & is.matrix(b)){
      if( all (dim(a) == dim(b))){
        z = paste0(a, b)
        dim(z) = dim(a)}
    }
    if((is.matrix(a) & is.atomic(b))){
      if(length(b) == 1){
        z = paste0(a, b)
        dim(z) = dim(a)}
    }
    if((is.matrix(b) & is.atomic(a))){
      if(length(a) == 1){
        z = paste0(a, b)
        dim(z) = dim(b)}
    }
    z}
  # formula2 = formula
  test <- multinom(formula = formula2, data = data)
  summ <-nnet:::summary.multinom(test)
  co=summ$coefficients
  se=summ$standard.errors
  z <- co/se
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  res<- (2 %f%    (coef(test)))%+% " (" %+%(3 %f% p) %+% ")"
  Effect = "Beta"
  ci <- round(confint(test,  level=0.95),2)

  summerclassic = list(type = "classic",
                       Effect = Effect,
                       est = coef(test),
                       ci = confint(test,  level=0.95),
                       p = p)
  if(OR)  {
    Effect = "OR"
    res<- (2 %f% exp(coef(test)))%+% " (" %+%(3%f% p) %+% ")"
    res1<- (2 %f% exp(coef(test)))
    res2<- (3 %f% p)
    ci <- confint(test,  level=0.95)
    ci = round(exp(ci),2)


    summerClassic = list(type = "classic",
                         Effect = Effect,
                         est = exp(coef(test)),
                         ci = exp(confint(test,  level=0.95)),
                         p = p)
  }

  Effect =paste0(Effect, " Level ")

  row.names(res) = row.names(p)
  row.names(res1) = row.names(p)
  row.names(res2) = row.names(p)
  colnames(res) = colnames(p)
  colnames(res1) = colnames(p)
  colnames(res2) = colnames(p)
  ress = as.data.frame (t(res))
  ress1 = as.data.frame (t(res1))
  ress2 = as.data.frame (t(res2))
  na = names(ress1)
  for(j in 1:length(na)){
    if(j ==1) { temp = data.frame(paste0(ress1[,na[j]]," (",ci[,1,j],", ",ci[,2,j], ")") , ress2[,na[j]])
    names(temp) = c(paste0(Effect ,na[j], "(95% CI)"), paste0("P-value ",na[j]) )
    }else{
      temp2= data.frame(paste0(ress1[,na[j]]," (",ci[,1,j],", ",ci[,2,j], ")") , ress2[,na[j]])
      names(temp2) =  c(paste0(Effect ,na[j], "(95% CI)"), paste0("P-value ",na[j]) )
      temp = cbind(temp,temp2)
      rm(temp2)
    }
    res = temp
    row.names(res) = row.names(ress2)
  }

  ################# new data  ------------
  plot = NULL
  bayesian.plot  = NULL
  pp.write = NULL
  lpp = NULL
  if(!is.null(newdata)){
    pp.write <- cbind(newdata, predict(test, newdata = newdata, type = "probs", se = TRUE))
    lpp <- melt(pp.write, id.vars =names(newdata), value.name = "probability")

  }

  if(!is.null(ggplot.mapping)){
    if(is.null(newdata)){
      stop("ggplot needs to use `newdata` arqument.")}

    lpp <- melt(pp.write, id.vars =names(newdata), value.name = "probability")

    plot= ggplot(lpp, mapping = ggplot.mapping) +
      geom_line() +
      facet_grid(variable ~ ., scales = "free")
    print(plot)
  }


  bayesian.plot = NULL
  bayesian.result = NULL

  co = NULL
  results.mnl = NULL
  summary.mnl = NULL
  if(bayesian == TRUE){
    library(UPG)
    Date.temp =  model.matrix( formula2, data)
    y =  model.frame( formula2, as.data.frame(data))[[1]]
    y =  as.factor(y)
    y <- relevel(y, ref = as.character(ref))

    # if(is.factor(y)) stop("dependent variable (y) must be factor.")
    X = cbind(1, Date.temp[,-1])
    # if (is.null(ref)) {
    #   tt = table(y)
    #   ref = names(tt[which.max(tt)])
    # }
    # pos.bl = which(levels(y) == ref)
    # new.lvls = c(levels(y)[-pos.bl], ref)
    # y.mnl = factor(y, levels = new.lvls)
    groups = levels(y)
    if (is.null(ref)) {
      ref = groups[1]
    }

    results.mnl <- UPG(y = y,
                       X = X,
                       model = 'mnl',
                       # A0 = 1,  # Intercept variance decreased from 4 to 1
                       # B0 = 1,  # Coefficient variance decreased from 4 to 1
                       # G0 = 10,
                       verbose = TRUE,
                       draws       = 11000,
                       burnin      = 4000,
                       baseline = ref)
    summary.mnl =  ( summary(results.mnl))
    colQuantile_2.5 = function(x) apply(x,c(2,3), quantile, p = 0.025)
    colQuantile_97.5 = function(x) apply(x,c(2,3), quantile, p = 0.975)
    p_value = function(x, null_value = 0){
      min(  2 * min(
        mean(x <= null_value),
        mean(x >= null_value)
      ),1)}

    colPvalue = function(x) apply(x,c(2,3), p_value )

    coMean = colMeans(  (results.mnl$posterior$beta))
    co2.5 = colQuantile_2.5( (results.mnl$posterior$beta ))
    co97.5 = colQuantile_97.5( (results.mnl$posterior$beta ))
    PPP = colPvalue( (results.mnl$posterior$beta ))
    Effect = "Beta"
    if(OR == TRUE){
      Effect = "OR"
      coMean = colMeans( exp(results.mnl$posterior$beta))
      co2.5  = colQuantile_2.5(exp(results.mnl$posterior$beta ))
      co97.5 = colQuantile_97.5(exp(results.mnl$posterior$beta ))
    }

    summerBayesian = list(type = "bayesian",
                          Effect = Effect,
                          est = coMean,
                          lower = co2.5,
                          upper =  co97.5,
                          p = PPP,
                          groups = results.mnl$posterior$groups)
    Effect =paste0(Effect, " Level ")


    co= (2 %f% coMean)%+% " (" %+%(2 %f% co2.5) %+% ", " %+% (2 %f% co97.5)%+% ")"
    rownames(co) = colnames(X)
    rownames(PPP) = colnames(X)
    colnames(co) = results.mnl$posterior$groups
    colnames(PPP) = results.mnl$posterior$groups
    TE = data.frame(co[,1],PPP[,1])
    names(TE) = c(paste0(Effect,results.mnl$posterior$groups[1] ),
                  paste0("P-value ",results.mnl$posterior$groups[1] ))

    for(k in 2:dim(PPP)[2]){
      temP = data.frame(co[,k],PPP[,k])
      names(temP) = c(paste0(Effect,results.mnl$posterior$groups[k] ),
                      paste0("P-value ",results.mnl$posterior$groups[k] ))
      TE = cbind.data.frame(TE, temP)
      rm(temP)
    }
    bayesian.result = TE


    bayesian.plot  = plot.UPG(x=results.mnl, OR = OR)
  }
  tt =   list(results = as.data.frame( (res)),meltData = lpp, ggplot = plot ,

              bayesian.plot = bayesian.plot, bayesian.result = bayesian.result,
              summerClassic = summerClassic,summerBayesian = summerBayesian,
              formula=formula2, data = data,
              results.mnl=results.mnl

  )
  class(tt) = union("Multinomial",class(tt))
  tt
}





###### Example

MultinomialPlot <-  R6::R6Class("Multinomial",
                                lock_object = FALSE,

                                public = list(
                                  legend_value = NULL,
                                  legend_name = NULL,
                                  group.n = 0,
                                  group = NULL,
                                  dat = NULL,
                                  initialize = function(object = NULL,
                                                        legend_name = "Group",
                                                        legend_value =  NULL) {

                                    if(!("Multinomial" %in% class(object)) & !is.null(object))  stop("The class of object is not Multinomial.")
                                    self$legend_name = legend_name
                                    self$legend_value = legend_value
                                    if(!is.null(object)){
                                      self$dat = data.frame(name = row.names(object$bayesian.result) )
                                      private$bayesEstExtracer(object,legend_value = legend_value)
                                      private$bayesLowerExtracer(object)
                                      private$bayesUpperExtracer(object)
                                      # self$plotterbayes()
                                    }
                                  },
                                  add = function(object,
                                                 legend_name = NULL,
                                                 legend_value =  NULL) {
                                    if(is.null(self$dat) ) self$dat = data.frame(name = row.names(object$bayesian.result))

                                    if(!("Multinomial" %in% class(object)))  stop("The class of object is not Multinomial.")
                                    if(is.null(self$legend_name) & !is.null(legend_name)) self$legend_name = legend_name
                                    if(is.null(self$legend_value) & !is.null(legend_value)) self$legend_value = legend_value

                                    private$bayesEstExtracer(object,legend_value = legend_value)
                                    private$bayesLowerExtracer(object)
                                    private$bayesUpperExtracer(object)
                                    # self$plotterbayes()
                                    invisible(self)
                                  },
                                  getRowNames = function(){
                                    self$dat$name
                                  },
                                  setRowNames = function(rownames){
                                    self$dat$name  =  rownames
                                  }
                                  ,

                                  addNewRow = function(new.names ){
                                    if(TRUE){
                                      dat.name = self$dat[,"name"]
                                      dat <-   self$dat

                                      lowerB = private$lowerB
                                      estB = private$estB
                                      upperB = private$upperB
                                      group.n = length(estB)
                                      n.estB = dim(estB[[1]])[2]
                                      kk = 0
                                      addRow = function(data , .before  , reName= NA  ){
                                        r.d = dim(data)[1]
                                        c.d = dim(data)[2]
                                        flag = FALSE
                                        if(c.d == 1){
                                          data[[c.d + 1]] = NA
                                          c.d = c.d + 1
                                          flag = TRUE
                                        }
                                        if(.before == 1){
                                          res =   rbind.data.frame(rep(reName,c.d),data)
                                        } else {
                                          temp1 = data[1:(.before-1), ]
                                          temp3 = data[(.before:r.d), ]
                                          res =   rbind.data.frame(temp1,rep(reName,c.d),temp3)
                                        }
                                        if(flag) res[[c.d]] = NULL
                                        res
                                      }
                                      while (TRUE) {
                                        kk = kk+1
                                        onRow = dat.name[kk]
                                        re =  readline(paste0(kk," - ", "Do you want to add new column before '",onRow ,
                                                              "' (yes, no, cancel)?    "))
                                        if(!(re  %in% c("","y", "n", "c","no", "yes", "cancel", 0 , 1))) {
                                          cat("\nWrong answer, back to previous question.")
                                          kk <- kk-1
                                        }
                                        if(re %in% c("y", "yes", 1)) {
                                          reName =  readline(paste0(kk," - ", "What was the name?"))
                                          n.dat = dim(dat)[1]
                                          i.onRow <- which(dat$name == onRow)
                                          dat =   addRow(dat, .before = i.onRow, reName= reName  )

                                          for (j in 1:group.n) {
                                            estB[[j]]     =   addRow (estB[[j]],   .before = i.onRow, reName= NA)
                                            upperB[[j]]   =   addRow (upperB[[j]], .before = i.onRow, reName= NA)
                                            lowerB[[j]]   =   addRow (lowerB[[j]], .before = i.onRow, reName= NA)
                                          }

                                        }
                                        if(re %in% c("c", "cancel")) {
                                          break}

                                        if(length(dat.name) == kk) break
                                      }
                                      self$dat       =  dat
                                      private$estB = estB
                                      private$lowerB = lowerB
                                      private$upperB = upperB


                                      # cat("New Rows:",
                                      #     paste0(self$dat$name   , collapse = "\n  ")
                                      # )
                                    }


                                  },


                                  rowRemover = function()  {
                                    if(TRUE){
                                      dat.name = self$dat[,"name"]
                                      dat.name.include = rep(TRUE, length(dat.name))
                                      kk = 0
                                      while (TRUE) {
                                        kk = kk+1
                                        re =  readline(paste0(kk," - ", "Do you want to include '",dat.name[kk] ,
                                                              "' in the figure (yes, no, cancel)?    "))
                                        if(!(re  %in% c("","y", "n", "c","no", "yes", "cancel", 0 , 1))) {
                                          cat("\nWrong answer, back to previous question.")
                                          kk <- kk-1
                                        }
                                        if(re %in% c("n", "no","No", "NO", 0)) {dat.name.include[kk] = FALSE}
                                        if(re %in% c("c", "cancel")) {
                                          dat.name.include = rep(TRUE, length(dat.name))
                                          break}

                                        if(length(dat.name) == kk) break
                                      }
                                      self$dat       =  data.frame(name = self$dat[dat.name.include,])
                                      for (i in 1: length(private$estB)) {
                                        private$estB[[i]]   = private$estB[[i]][dat.name.include,]
                                        private$lowerB[[i]] = private$lowerB[[i]][dat.name.include,]
                                        private$upperB[[i]] = private$upperB[[i]][dat.name.include,]
                                      }

                                      cat("Selected:",
                                          paste0(self$dat$name   , collapse = "\n  ")
                                      )
                                    }
                                  }
                                  ,
                                  plotterbayes = function(show = FALSE,
                                                          ref_line = NULL,
                                                          xlim = c(0, 2),
                                                          ticks_at = c(0, 0.5, 1, 1.5, 2),
                                                          nudge_y = 0.4 ,
                                                          plotCellSize = 20,
                                                          padding = 2  ,
                                                          name  = NULL,
                                                          tm = NULL
                                  ){
                                    if(is.null(ref_line)) {
                                      if(private$effect == "OR") ref_line = 1
                                      if(private$effect == "Beta") ref_line = 0
                                    }
                                    est = list()
                                    k = 0
                                    for (i in 1:(self$group.n-1)) {
                                      for (j in 1:length(private$formula)) {
                                        k = k+1
                                        est[[k]] = private$estB[[i]][,j]
                                      }
                                    }
                                    ########
                                    lower = list()
                                    k = 0
                                    for (i in 1:(self$group.n-1)) {
                                      for (j in 1:length(private$formula)) {
                                        k = k+1
                                        lower[[k]] = private$lowerB[[i]][,j]
                                      }
                                    }
                                    ##########
                                    upper = list()
                                    k = 0
                                    for (i in 1:(self$group.n-1)) {
                                      for (j in 1:length(private$formula)) {
                                        k = k+1
                                        upper[[k]] = private$upperB[[i]][,j]
                                      }
                                    }
                                    ################
                                    library(forestploter)
                                    library(grid)



                                    f= private$formula
                                    f.n = length(f)
                                    y.names = c()

                                    for(i in 1:(2*f.n)){
                                      if(i%%2==1)     y.names[i] = paste(" ", collapse = "")
                                      if(i%%2==0)     y.names[i] = all.vars(update.formula(f[[i/2]] , .~1))
                                    }
                                    # y.names[4] = "ff"
                                    # y.names[2] = "ff2"
                                    dat = as.data.frame(self$dat)
                                    for (i in 1:(2*f.n)) {
                                      if(i%%2==0)  dat[[i+1]] <- paste(rep(" ", plotCellSize), collapse = " ")
                                      if(i%%2==1)  dat[[i+1]] <- paste(rep(" ", padding), collapse = " ")
                                    }
                                    names(dat)[-1]  = y.names
                                    # names(dat)
                                    dat[[2]] = NULL
                                    if(!is.null(name)){
                                      if(length(name) != f.n) stop("The length of name is not correct.")
                                      names(dat)[(1:f.n)*2] =  name
                                    }
                                    if(f.n ==1)  ci_c = 2
                                    if(f.n !=1)  ci_c =  (1:f.n)*2
                                    names(dat)[1]="Factors"
                                    private$bayesianPlot <- forest( dat,
                                                                    est = est,
                                                                    lower = lower,
                                                                    upper = upper,
                                                                    ci_column = ci_c,
                                                                    ref_line = ref_line,
                                                                    xlim = xlim,
                                                                    xlab = rep(private$effect, f.n),
                                                                    # vert_line = c(0.5, 2),
                                                                    ticks_at = ticks_at,
                                                                    nudge_y = nudge_y ,
                                                                    theme = tm
                                    )
                                    if(show)    (plot(private$bayesianPlot))
                                  }
                                  ,
                                  saveBayesPlot = function(file = NULL, width= NULL,
                                                           height= NULL, dpi= NULL){
                                    if(is.null(private$bayesianPlot)){
                                      cat("\nNo Bayesian plot available.")
                                    } else {
                                      if(is.null(dpi))     dpi    = 1000
                                      if(is.null(width))   width  = forestploter::get_wh(private$bayesianPlot)[1]
                                      if(is.null(height))  height = forestploter::get_wh(private$bayesianPlot)[2]
                                      if(is.null(file))    file   = paste0(as.numeric(Sys.time()),"-OR.jpeg")
                                      ggsave(file, private$bayesianPlot, width = width,
                                             height = height, dpi = dpi)
                                    }
                                  }
                                ),
                                private = list(
                                  estB = NULL,
                                  lowerB = NULL,
                                  upperB = NULL,
                                  estC = NULL,
                                  lowerC = NULL,
                                  upperC = NULL,
                                  formula = list(),
                                  data = list(),
                                  effect = NULL,
                                  bayesianPlot = NULL,
                                  calsicalPlot = NULL,
                                  bayesEstExtracer =   function(object,legend_value){
                                    if(!is.null(object$bayesian.result)){
                                      summerBayesian = object$summerBayesian
                                      if(self$group.n == 0){
                                        self$group.n = dim(summerBayesian$est)[2]
                                        self$group = summerBayesian$group
                                      } else {
                                        if(self$group.n != dim(summerBayesian$est)[2])
                                          stop("number of group levels is not the same as previouse object.")
                                      }
                                      if(is.null(legend_value))   self$legend_value = summerBayesian$groups
                                      private$formula = append(private$formula, object$formula)
                                      private$data = append(private$data, object$data)
                                      private$effect =  summerBayesian$Effect
                                      estB =  summerBayesian$est
                                      if(is.null( private$estB)){
                                        for(i in 1:dim(estB)[2]){
                                          private$estB[[i]] = data.frame()
                                        }}
                                      for(i in 1:dim(estB)[2]){
                                        if( dim(private$estB[[i]])[2] ==0){
                                          private$estB[[i]] =   data.frame(estB[,i])
                                        }else {
                                          private$estB[[i]] = cbind.data.frame(private$estB[[i]],estB[,i])
                                        }
                                      }

                                    }},

                                  bayesLowerExtracer =   function(object){
                                    if(!is.null(object$bayesian.result)){
                                      summerBayesian = object$summerBayesian

                                      lowerB =  summerBayesian$lower
                                      if(is.null( private$lowerB)){
                                        for(i in 1:dim(lowerB)[2]){
                                          private$lowerB[[i]] = data.frame()
                                        }}
                                      for(i in 1:dim(lowerB)[2]){
                                        if( dim(private$lowerB[[i]])[2] ==0){
                                          private$lowerB[[i]] =   data.frame(lowerB[,i])
                                        }else {
                                          private$lowerB[[i]] = cbind.data.frame(private$lowerB[[i]],lowerB[,i])
                                        }
                                      }

                                    }},

                                  bayesUpperExtracer =   function(object){
                                    if(!is.null(object$bayesian.result)){
                                      summerBayesian = object$summerBayesian

                                      upperB =  summerBayesian$upper
                                      if(is.null( private$upperB)){
                                        for(i in 1:dim(upperB)[2]){
                                          private$upperB[[i]] = data.frame()
                                        }}
                                      for(i in 1:dim(upperB)[2]){
                                        if( dim(private$upperB[[i]])[2] ==0){
                                          private$upperB[[i]] =   data.frame(upperB[,i])
                                        }else {
                                          private$upperB[[i]] = cbind.data.frame(private$upperB[[i]],upperB[,i])
                                        }
                                      }

                                    }}


                                )
)
}
# joined_data =  readRDS(choose.files())
###############################################bayes model
# joined_data$IPAQ_TotalScoreCat3 = as.factor(joined_data$IPAQ_TotalScoreCat3 )

# var = c(
#   "N_Subjective_sleep_quality",
#   "N_Sleep_latency",
#   "N_Sleep_duration",
#   "N_Sleep_efficiency",
#   "N_Sleep_disturbance",
#   "N_Daytime_dysfunction",
#   "N_Use_of_sleep_medication")
#
# Haschronic = c("N_HasDiabet",
#                "N_HasFattyLiver",
#                "N_HasChronicLungDisease",
#                "N_HasHypertension",
#                "N_HasThyroid",
#                "N_HasDepression",
#                "N_HasCardiacIschemic" )
#
# selectedData=  joined_data
# for (i in 1:length(Haschronic)) {
#  selectedData[[Haschronic[i]]] = as.character(selectedData[[Haschronic[i]]])
#  selectedData[[Haschronic[i]]] [which(selectedData[[Haschronic[i]]]%in% c("DontKnow","Unknown"))] = "No"
# }
#
#
# RES = list()
# for (i in 1:7 ) {
#   # selectedData=na.omit(joined_data[which(!(joined_data[[Haschronic[i]]] %in% c("DontKnow","Unknown")) ),allVar])
#
#
#   cat(var[i],"\n")
#   ######## TODO : ReOrder
#   formula2 =as.formula(paste0(
#     var[i]," ~ Age+
#     BMI+
#    GenderID +
#    N_MaritalStatus+
#    N_Education +
#    SES3 +
#    N_MobileUseCat +
#    IPAQ_TotalScoreCat3+
#    SmokingNoCat2+
#    ERIQ_cat +
#    N_HasCardiacIschemic +
#    N_HasChronicLungDisease +
#    N_HasHypertension +
#    N_HasThyroid +
#    N_HasDepression +
#    N_HasHypertension +
#    N_HasDiabet +
#    N_HasFattyLiver+
#    N_CancersValue"))
#   set.seed(12345)
#   data = selectedData
#   RES[[i]]=multinomial(data = selectedData  ,
#                        formula2 = formula2 ,
#                        ref = "0" ,
#                        bayesian = TRUE)
#
#
#
# }
 # return(RES)

    # , args = list(joined_data=joined_data))

# MUL$is_alive()
# MUL$get_result()

######################################### Plot
# row.names(RES[[1]]$bayesian.result)[1] = "intercept"
# as.data.frame(RES[[1]]$bayesian.result) %>%
#   write.cb(all.vars(update.formula(RES[[1]]$formula, .~1)) )




# self = MultinomialPlot$new(RES[[1]])
# self$add(RES[[2]])
# self$add(RES[[3]])
# self$add(RES[[4]])
# self$add(RES[[5]])
# self$add(RES[[6]])
# self$add(RES[[7]])
