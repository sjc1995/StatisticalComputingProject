# This function deals with two variables, when number of group elements is 2,if not normal, use wilcox-test
# and var of two group is same, use var equal t test, and var different use var different t test, if group 
# elements is more than two, if normal use anova test, if not use krustal test

single.variable.analysis <- function(x,group,dataclean = 'Remove'){
  if(!all(!is.na(x))){
    warning('There is some NA values in data')
  }
  if(!all(!is.na(group))){
    stop('There is NA value in group')
  }
  if(tolower(dataclean) == 'remove'){
    data <- data.frame(group=group,x=x)
    data <- na.omit(data)
    x <- data$x
    group <- data$group
  }
  if(length(group) != length(x)) stop("The length of data and group are different")
  if(!is.numeric(x)) stop("data is not only containing numbers")
  d <- length(unique(group))
  for(i in 1:d){
    if(tolower(dataclean) == 'mean'){
      x[group==unique(group)[i]&is.na(x)] <- mean(na.omit(x[group==unique(group)[i]]))
    }else if(tolower(dataclean) == 'median'){
      x[group==unique(group)[i]&is.na(x)] <- median(na.omit(x[group==unique(group)[i]]))
    }
  }
  if(d < 2){
    stop("The group of data is less than 2")
  } else if(d == 2){
    if(length(x[group==unique(group)[1]])<30){
      w1 <- shapiro.test(x[group==unique(group)[1]])$p.value
    }else{
      w1 <- 1
    }
    if(length(x[group==unique(group)[2]])<30){
      w2 <- shapiro.test(x[group==unique(group)[2]])$p.value
    }else{
      w2 <- 1
    }
    if(w1 > 0.05 & w2 >0.05){
      method <- 't-test'
      var <- var.test(x[group==unique(group)[1]],x[group == unique(group)[2]])
      if(var$p.value > 0.05){
        tmp <- t.test(x[group==unique(group)[1]],x[group == unique(group)[2]],var.equal = T)
      }else{
        tmp <- t.test(x[group==unique(group)[1]],x[group == unique(group)[2]],var.equal = F)
      }
      statistic <- tmp$statistic
      p.value <- tmp$p.value
    }else{
      method <- 'wilcox-test'
      tmp <- wilcox.test(x[group==unique(group)[1]],x[group == unique(group)[2]])
      statistic <- tmp$statistic
      p.value <- tmp$p.value
    }
  }else{
    w <- rep(1,d)
    for(i in 1:d){
      if(length(x[group==unique(group)[i]])<30){
        w[i] <- shapiro.test(x[group==unique(group)[i]])$p.value
      }
    }
    if(all(w>0.05)){
      method <- 'ANOVA-test'
      var <- bartlett.test(x~factor(group))$p.value
      if(var < 0.05){
        warning('The variance between group in ANOVA test may not same')
      }
      tmp <- anova(lm(x~factor(group)))
      statistic <- tmp$'F value'[1]
      p.value <- tmp$'Pr(>F)'[1]
    }else{
      method <- 'Kruskal-test'
      tmp <- kruskal.test(x~factor(group))
      statistic <- tmp$statistic
      p.value <- tmp$p.value
    }
  }
  result <- list(method = method,statistic = statistic, p.value = p.value, data =data.frame(group = group, data = x), model = tmp)
  class(result) <- 'singlevariable'
  return(result)
}

print.singlevariable <- function(x){
  print(x$model)
}


oneframe.analysis <- function(x,dataclean = 'Remove'){
  if(!is.data.frame(x)){
    if(!is.matrix(x)){
      stop('Input is not a dataframe')
    }else{
      x <- data.frame(x)
    }
  }
  if('group' %in% tolower(colnames(x))){
    x <- data.frame('group'=x[,tolower(colnames(x))=='group'],x[,tolower(colnames(x))!='group'])
  }
  group <- list(group = x[,1])
  data <- as.list(x[,-1])
  result <- lapply(FUN=function(x) single.variable.analysis(x,group,dataclean),data)
  tmp <- function(x){
    x$p.value
  }
  p.value <- sapply(FUN=tmp,result)
  tmp <- function(x){
    x$statistic
  }
  statistic <- sapply(FUN=tmp,result)
  tmp <- function(x){
    x$method
  }
  method <- sapply(FUN=tmp,result)
  result <- list(method = method, statistic = statistic, p.value = p.value, data = x)
  class(result) <- 'oneframe'
  return(result)
}



print.oneframe <- function(x){
  cat('\n         One frame analysis\n')
  sig <- NA
  sig[x$p.value<0.001] <- '***'
  sig[x$p.value>0.001&x$p.value<0.01] <- '**'
  sig[x$p.value>0.01&x$p.value<0.05] <- '*'
  sig[x$p.value>0.05&x$p.value<0.1] <- '.'
  sig[x$p.value>0.1] <- ' '
  res <- data.frame(method = x$method, statistic = x$statistic, p.value=x$p.value, sig)
  colnames(res)[4] <- ' '
  print(res)
  cat('\n---')
  cat('\nSignif. codes:  0 \'***\' 0.001 \'**\' 0.01 \'*\' 0.05 \'.\' 0.1 \' \' 1\n')
}



meta.analysis <- function(data1 = NULL,data2 = NULL,data3 = NULL,data4 = NULL,data5 = NULL,method = "Fisher",dataclean = 'Remove'){
  data.name <- list(as.character(substitute(data1)),as.character(substitute(data2)),as.character(substitute(data3)),as.character(substitute(data4)),as.character(substitute(data5)))
  data <- list(data1, data2, data3, data4, data5)
  for( i in c(5:1)){
    if(!is.data.frame(data[[i]])){
      if(is.character(data[[i]])){
        if( tolower(data[[i]]) %in% c('fisher','stouffer' ,'min','minp','minpvalue','min pvalue','max','maxp','maxpvalue','max pvalue')){
          method <- data[[i]]
          data[[i]] <- NULL
        }else if(tolower(data[[i]]) %in% c('remove','mean','median')){
          cleandata <- data[[i]]
          data[[i]] <- NULL
        }else{
          warning('unexpected value for method and cleandata,Fiser method is used')
          data[[i]] <- NULL
        }
      }else if(is.null(data[[i]])){
        data[[i]] <- NULL
      }else if(is.matrix(data[[i]])){
        if('group' %in% tolower(colnames(data[[i]]))){
          data[[i]] <- data.frame('group'=data[[i]][,tolower(colnames(data[[i]]))=='group'],data[[i]][,tolower(colnames(data[[i]]))!='group'])
          names(data)[i] <- data.name[i]
        }else{
          data[[i]] <- data.frame(data[[i]]) 
          names(data)[i] <- data.name[i]
        }
      }else{
        stop("Input data is not a dataframe")
      }
    }else{
      if('group' %in% tolower(colnames(data[[i]]))){
        data[[i]] <- data.frame('group'=data[[i]][,tolower(colnames(data[[i]]))=='group'],data[[i]][,tolower(colnames(data[[i]]))!='group'])
        names(data)[i] <- data.name[i]
      }
      names(data)[i] <- data.name[i]
    }
  }
  k <- length(data)
  if(k < 1){
    stop('There is no data input')
  }
  if(k==1){
    result <- oneframe.analysis(data[[1]],dataclean)
    warning('There is only one dataframe')
    return(result)
  }
  d <- NA
  for(i in c(1:k)){
    d[i] <- dim(data[[i]])[2]
  }
  if(!all(d == rep(d[1]))) stop('The dim is wrong')
  pv <- lapply(data, function(x){oneframe.analysis(x,dataclean)})
  tmp <- function(x){
    x$p.value
  }
  pv <- sapply(FUN=tmp,pv)
  pv <- data.frame(pv)
  statistic <- NA
  p.value <- NA
  if(tolower(method) == 'fisher'){
    method <- 'Fisher'
    for(i in 1:(d[1]-1)){
      statistic[i] <- -2*sum(log(pv[i,]))
      p.value[i] <- pchisq(statistic[i],(2*k),lower.tail = F)
    }
  }else if(tolower(method) == 'stouffer'){
    method <- 'Stouffer'
    for(i in 1:(d[1]-1)){
      statistic[i] <- sum(sapply(pv[i,],qnorm))/sqrt(k)
      p.value[i] <- pnorm(statistic[i])
    }
  }else if(tolower(method) %in% c('min','minp','minpvalue','min pvalue')){
    method <- 'MinP'
    for(i in 1:(d[1]-1)){
      statistic[i] <- min(pv[i,])
      p.value[i] <- pbeta(statistic[i],1,k)
    }
  }else if(tolower(method) %in% c('max','maxp','maxpvalue','max pvalue')){
    method <- 'MaxP'
    for(i in 1:(d[1]-1)){
      statistic[i] <- max(pv[i,])
      p.value[i] <- pbeta(statistic[i],k,1)
    }
  }else{
    method <- 'Fisher'
    for(i in 1:(d[1]-1)){
      statistic[i] <- -2*sum(log(pv[i,]))
      p.value[i] <- pchisq(statistic[i],(2*k),lower.tail = F)
    }
    warning('unexpected value for method, Fisher method is used')
  }
  result <- list(statistics = statistic, p.value = p.value , data = data, method = method)
  class(result) <- 'meta'
  return(result)
}



print.meta <- function(x){
  cat('\n         Meta-analysis for',length(x$data),'dataframes\n')
  cat('Data:  ',names(x$data),'\n')
  cat('Method:',x$method)
  result <- data.frame(stat = x$statistic, pv = x$p.value)
  rownames(result) <- names(x$data[[1]])[-1]
  cat('\nStatistics and p-values:\n')
  print(t(result)[,1:length(x$statistic)])
  cat('\nalternative hypothesis: true difference in means in different groups and dataframe is not equal to 0\n')
}



summary.meta <- function(x){
  result <- data.frame(statistic = x$statistic,p.value = x$p.value)
  rownames(result) <- names(x$data[[1]])[-1]
  samplesize <- sapply(FUN = dim,x$data)[1,]
  sig1 <- result[x$p.value<0.01,]
  sig2 <- result[x$p.value<0.05,]
  sig3 <- result[x$p.value<0.1,]
  res <- list(result= result,samplesize=samplesize,sig1=sig1,sig2=sig2,sig3=sig3,data=x$data,method=x$method)
  class(res) <- 'summary.meta'
  return(res)
}



print.summary.meta <- function(x){
  cat('\n         Meta-analysis for',length(x$data),'dataframes\n')
  cat('Data:  ',names(x$data),'\n')
  cat('Method:',x$method,'\n')
  printCoefmat(x$result,P.values = T, has.Pvalue = T)
  cat('\nSample size:\n')
  print(x$samplesize)
  if(dim(x$sig1)[1]>0){
    cat('\nSignificant at 0.01:\n')
    print(rownames(x$sig1))
  }
  if(dim(x$sig2)[1]>0){
    cat('\nSignificant at 0.05:\n')
    print(rownames(x$sig2))
  }
  if(dim(x$sig3)[1]>0){
    cat('\nSignificant at 0.10:\n')
    print(rownames(x$sig3))
  }
}

