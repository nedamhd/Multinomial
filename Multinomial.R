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
    source("https://raw.githubusercontent.com/ahadalizadeh/utility_fun/master/utility_fun.R")
  # }
  require(nnet)
  require(ggplot2)
  require(reshape2)

  ################# change ref ---------
if(!is.null(ref)){
  DV =   all.vars(update.formula(formula2, .~1)  )
  # Data$DV = Data[[DV]]
  data[[DV]] <- relevel(as.factor(data[[DV]]), ref = as.character(ref))
  # formula = update.formula(formula, DV2~.)
}

  ################## run model   ---------
# formula2 = formula
  test <- multinom(formula = formula2, data = data)
  summ <-nnet:::summary.multinom(test)
co=summ$coefficients
se=summ$standard.errors
  z <- co/se
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  res<- (2 %f%    (coef(test)))%+% " (" %+%("3.P" %f% p) %+% ")"
  Effect = "Beta Level "
  ci <- round(confint(test,  level=0.95),2)

  if(OR)  {
    Effect = "OR Level "
    res<- (2 %f% exp(coef(test)))%+% " (" %+%("3.P" %f% p) %+% ")"
    res1<- (2 %f% exp(coef(test)))
    res2<- ("3.P" %f% p)
    ci <- confint(test,  level=0.95)
    ci = round(exp(ci),2)
    }
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
    y =  model.frame( formula2, data)[[1]]
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
                       verbose = TRUE,
                       draws       = 5000,
                       burnin      = 2000,
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
    Effect = "Beta Level "
    if(OR == TRUE){
      Effect = "OR Level "
      coMean = colMeans( exp(results.mnl$posterior$beta))
      co2.5 = colQuantile_2.5(exp(results.mnl$posterior$beta ))
      co97.5 = colQuantile_97.5(exp(results.mnl$posterior$beta ))
    }
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
  list(results = as.data.frame( (res)),meltData = lpp, ggplot = plot ,

       bayesian.plot = bayesian.plot, bayesian.result = bayesian.result,
       formula=formula2
  )}





###### Example 


N_Subjective_sleep_quality=multinomial(data =
                                        joined_data  ,
                                       formula2 =N_Subjective_sleep_quality ~ Age+
                                    GenderID +BMI+N_MaritalStatus+
                                    N_Education + SES3 +
                                    N_MobileUseCat +
                                    ERIQ_cat + IPAQ_TotalScoreCat3+
                                    SmokingNoCat2+
                                    N_HasHypertension +
                                    N_HasDiabet +
                                  N_HasHypertension +
                                    N_HasDiabet +
                                    N_HasFattyLiver+
                                    N_HasCardiacIschemic +
                                    N_HasChronicLungDisease +
                                    N_HasThyroid +
                                    N_HasDepression +
                                    N_CancersValue,


                                  ref = "0" ,
                                  # newdata =  newdata1 ,
                                  bayesian = TRUE
                                  # ggplot.mapping =
                                  #   aes(x = Age,
                                  #       y = probability ,
                                  #       colour = as.factor(GenderID),
                                  #       group = as.factor(GenderID)
                                  # )
)

# N_Subjective_sleep_quality $results %>% write.cb()
N_Subjective_sleep_quality $bayesian.result %>% write.cb()
