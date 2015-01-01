add.zscore <- function(data, set) {
  if (set %in% c('F210Ihom', 'M323Khom')) cond <- 'conditionHom'
  if (set %in% c('Zfp106')) cond <- 'conditionWT'
  if (set %in% c('F210Ihet', 'M323Khet')) cond <- 'conditionHet'
  if (set %in% c('Mbnl2')) cond <- 'conditionTreat'
  if (set %in% c('FUS', 'TDP43')) cond <- 'conditionKD'
  
  data$basic.pval <- ifelse (  data$basic.pval < 10^-15, 10^-15, data$basic.pval )
  data$logP <- -log(data$basic.pval)
  data$zscore <- qnorm(p = 1 - data$basic.pval/2)
  data$zscore<- ifelse ( data[, cond] > 0, abs(data$zscore), - abs(data$zscore))
  return(data)
}

get.slope.estimate <- function(pos.reads.kb = NULL, start.intron, end.intron, strand, plot = FALSE,
                               warren.junctions.pos = c(), warren.junctions.IDs = c(), 
                               extra.exons.pos = c(),
                               my.title = '',
                               add.lines = TRUE,
                               interaction.term = FALSE,
                               min.intron.size = 50,
                               max.count = 10^6
                               ) {

  n.warren.junctions <- length(warren.junctions.pos)
  if (length(warren.junctions.IDs) != n.warren.junctions) stop('IDs must match positions')
  
  n.extra.exons <- length(extra.exons.pos)
  res <- list(n.reads = length(pos.reads.kb), length.kb = end.intron - start.intron, strand = strand,
              n.warren.junctions = n.warren.junctions,
              n.extra.exons = n.extra.exons)

  my.hist <- hist(x = pos.reads.kb,
                  breaks = seq(from = start.intron, to = end.intron+5, by = 5),
                  plot = FALSE)
  
  ##### Main table of counts
  for.regression <- data.frame(x = my.hist$mids,
                               y = my.hist$counts)
  
  smooth <- TRUE
  if (smooth) {
    sd.loc <- sd(for.regression$y)
    for.regression$outlier <- for.regression$y > median(for.regression$y) + 7*sd.loc
    sd.loc <- sd(subset(for.regression$y, ! for.regression$outlier))
    for.regression$outlier <- for.regression$y > median(for.regression$y) + 7*sd.loc
      
    max.normal <- max(subset(for.regression$y, ! for.regression$outlier))
    for.regression$y <- ifelse (  for.regression$outlier, max.normal, for.regression$y)

    my.hist$counts <- for.regression$y
  }

    
  ### mod0 will contain the basic regression analysis, no Warren junction added at this stage
  for.regression$xclean <- my.hist$mids - min(my.hist$mids)
  base.formula <- 'y ~ xclean'
  
  if (n.extra.exons > 0) {
     if (strand == 1)  extra.exons.pos <- sort(extra.exons.pos, decreasing = FALSE)
     if (strand == -1) extra.exons.pos <- sort(extra.exons.pos, decreasing = TRUE)  ##start with proper first intron
    
     extra.exons.IDs <- paste('E', round(1000*extra.exons.pos), sep = '')
     

     if ((n.extra.exons >= 2) && (interaction.term)) {
       intron.length <- abs(diff(extra.exons.pos))
       index <- 1 + min(which(intron.length > min.intron.size))
       base.formula <- paste('y ~ xclean*', extra.exons.IDs[ index ], sep = '')
       print(base.formula)
       coeff.inter <- paste('xclean:', extra.exons.IDs[ index ], 'TRUE', sep = '')
     }
    
    for (i in 1:n.extra.exons) {
      for.regression[, extra.exons.IDs[ i ] ]  <- for.regression$x > extra.exons.pos[ i ]
      base.formula <- paste(base.formula, '+', extra.exons.IDs[ i ], sep = ' ')
    }
     res$extra$formula <- base.formula
  }


  mod0 <- lm (data = for.regression, formula = base.formula)
  modbase <- mod0

  if ((n.extra.exons >= 2) && (interaction.term)) {
    print( summary(modbase) )
    res$slope.basic <- summary(modbase)[[4]]['xclean','Estimate']
    res$slope.inter <- summary(modbase)[[4]]['xclean','Estimate']  + summary(modbase)[[4]][coeff.inter, 'Estimate']
    res$evidence.that.slopes.differ.pval <-  summary(modbase)[[4]][coeff.inter, 'Pr(>|t|)']
    res$junction.for.interaction <- extra.exons.IDs[ index ]
  }
  
  ###### Now actually run the stepwise regression
  junctions.selected <- c()
  res$my.overall.P <- 1

  message('Number of junctions provided by Warren: ', n.warren.junctions)
  if (n.warren.junctions > 0) {

    for (i in 1:n.warren.junctions) {
      for.regression[, warren.junctions.IDs[ i ] ]  <- as.numeric(for.regression$x > warren.junctions.pos[ i ])
    }
    
    first <- TRUE
    while (TRUE)  {
      min.junc <- -1
      min.P <- 1
      
      for (j in 1:n.warren.junctions) {
        new.formula <- paste(base.formula, ' + ', warren.junctions.IDs[ j ])
        mod1 <- lm (data = for.regression, formula = as.formula(new.formula ))
        print(summary(mod1))
        my.P <- anova(mod0, mod1, test = 'F')[2,6]

        if (my.P < min.P && !is.na(my.P)) {
          min.P <- my.P
          min.junc <- j
          updated.formula <- new.formula
        }
      }

      if ( (min.junc != -1) && (first | min.P < 10^-3) ) {
        message('Selecting ', warren.junctions.IDs[ min.junc ], ' with P = ', min.P)
        junctions.selected <- c(junctions.selected, warren.junctions.IDs[ min.junc ])
        base.formula <- updated.formula
        mod0 <- lm (data = for.regression, formula = base.formula)
      }
      
      first <- FALSE
      if (min.P > 10^-3) break
    }
  
    res$my.overall.P <- anova(modbase, mod0, test = 'F')[2,6]
    if (is.na(res$my.overall.P)) res$my.overall.P <- 1
  } else message('No stepwise algorithm') #### end of the stepwise regression step


  if (res$my.overall.P < 0.001) {fitted.ordered <- as.numeric(fitted(mod0))} else {fitted.ordered<- as.numeric(fitted(modbase))}
  

  if (plot) {
    message('Plot')
    my.hist$counts <- pmin( my.hist$counts, max.count)
    plot(my.hist,
         xlab = 'Position (kb)',
         ylab = 'Read count',
         xlim = c(start.intron, end.intron),
         col = 'red',
         main = '')
    points ( x = for.regression$x, y = fitted.ordered, type = 'l', lty = 1, col = 'black')
    
    if (length(warren.junctions.IDs) > 0) {
      selected <- warren.junctions.IDs %in% junctions.selected
      abline(v = warren.junctions.pos, col = ifelse( selected, 'blue', 'black'))
      abline(v = subset(warren.junctions.pos, selected), col = 'blue', lwd = 2)  ##selected junctions in blue
      title(sub = paste('Strand:', res$strand, ', junction(s): ',  junctions.selected, 'regression P: ', signif(res$my.overall.P, 3)))
    } else {
      title(sub = paste('Strand:', res$strand))
    }

    if ((add.lines) && (length(extra.exons.pos) > 0)) {abline(v = extra.exons.pos, lwd = 2, col = 'black')}



    title(main = my.title)
  }
  message('Done with plot')

  res$extra$regression.table <- for.regression
  res$extra$histogram <- my.hist
  res$extra$modbase <- modbase
  
  res$junctions.selected <- paste(junctions.selected, collapse = ',')
  
  res$r2.before <- summary(modbase)$r.squared
  res$r2.after <- summary(mod0)$r.squared
  res$slope.before <-  coef(modbase)[['xclean']]
  res$slope.after <-  coef(mod0)[['xclean']]
  #res$formula <- base.formula 
  
  if ( res$strand == 1 ) {fitted.ordered <- rev(fitted.ordered)}

  nbins <- length(fitted.ordered)
  #if (nbins > 10) {res$norm.sum.reads.50kb <- sum(fitted.ordered[1:10])}
  #if (nbins > 20) {res$norm.sum.reads.100kb <- sum(fitted.ordered[1:20])}
  
  #if (nbins > 30) {
  #  res$value.150kb <- fitted.ordered[30]
  #  res$norm.sum.reads.150kb <- sum(fitted.ordered[1:30])
  #}
  
  #if (nbins > 60) {res$value.300kb <- fitted.ordered[60]}

  res$value.end.intron <- as.numeric(fitted.ordered[1])
    
  return (res)
}


copy.list <- function (list, data.frame, row) {
  for (n in subset(names(list), names(list) != 'extra')) {
    if (! n %in% names(data.frame)) data.frame[, n ] <- NA
    data.frame[row, n ] <- list[[ n ]] 
  }
  return (data.frame)
}
