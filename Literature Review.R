# -----------------------------------------------------------------------------------
# ---- Set Up -----------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# set the path of where the input files are
mywd0 = "C:/Users/Nick Morris/Downloads/ABP/Budget-Uncertainty-Papers/Literature-Review/Vaccine Budget Uncertainty Abstracts"
mywd1 = "C:/Users/Nick Morris/Downloads/ABP/Budget-Uncertainty-Papers/Literature-Review/Vaccine Budget Uncertainty Abstracts/vaccine budget uncertainty ScienceDirect"
mywd2 = "C:/Users/Nick Morris/Downloads/ABP/Budget-Uncertainty-Papers/Literature-Review/Vaccine Budget Uncertainty Abstracts/vaccine budget uncertainty JSTOR"
mywd3 = "C:/Users/Nick Morris/Downloads/ABP/Budget-Uncertainty-Papers/Literature-Review/Vaccine Budget Uncertainty Abstracts/vaccine budget uncertainty IEEE"
mywd4 = "C:/Users/Nick Morris/Downloads/ABP/Budget-Uncertainty-Papers/Literature-Review/Vaccine Budget Uncertainty Abstracts/vaccine budget uncertainty WebOfScience"
mywd5 = "C:/Users/Nick Morris/Downloads/ABP/Budget-Uncertainty-Papers/Literature-Review/Vaccine Budget Uncertainty Abstracts/vaccine budget uncertainty EngineeringVillage"

# open up a graphics window
windows()

# -----------------------------------------------------------------------------------
# ---- Packages ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# data handling
require(data.table)
require(tm)
require(gtools)
require(stringdist)
require(stringr)
require(zoo)
require(missRanger)
require(xml2)
require(slam)

# plotting
require(ggplot2)
require(gridExtra)
require(GGally)
require(scales)
require(corrplot)
require(plot3D)
require(wordcloud)

# modeling
require(fpc)
require(caret)
require(cluster)
require(h2o)
require(MLmetrics)
require(NLP)
require(factoextra)

# parallel computing
require(foreach)
require(doSNOW)

}

# -----------------------------------------------------------------------------------
# ---- Functions --------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- prints the data types of each column in a data frame -------------------------

types = function(dat)
{
  require(data.table)
  
  # make dat into a data.table
  dat = data.table(dat)
  
  # get the column names
  column = names(dat)
  
  # get the class of the columns
  dataType = sapply(1:ncol(dat), function(i) class(unlist(dat[, i, with = FALSE])))
  
  # compute the number of levels for each column
  levels = sapply(1:ncol(dat), function(i) ifelse(dataType[i] == "factor", length(levels(droplevels(unlist(dat[, i, with = FALSE])))), 0))
  
  # compute the number of unique values for each column
  uniqueValues = sapply(1:ncol(dat), function(i) length(unique(unname(unlist(dat[, i, with = FALSE])))))
  
  # compute the portion of missing data
  missing = sapply(1:ncol(dat), function(i) nrow(na.omit(dat[, i, with = FALSE], invert = TRUE)) / nrow(dat))
  
  # build the output table 
  output = data.table(column, id = 1:length(column), dataType, levels, uniqueValues, missing)
  
  # order output by dataType
  output = output[order(dataType)]
  
  return(output)
}

# ---- converts all columns to a character data type --------------------------------

tochar = function(dat)
{
  require(data.table)
  
  # make dat into a data.frame
  dat = data.table(dat)
  
  # get the column names
  column = names(dat)
  
  # get the values in the columns and convert them to character data types
  values = lapply(1:ncol(dat), function(i) as.character(unname(unlist(dat[, i, with = FALSE]))))
  
  # combine the values back into a data.frame
  dat = data.table(do.call("cbind", values), stringsAsFactors = FALSE)
  
  # give dat its column names
  setnames(dat, column)
  
  return(dat)
}

# ---- a qualitative color scheme ---------------------------------------------------

qcolor = function(n, a = 1)
{
  require(grDevices)
  require(scales)
  return(alpha(colorRampPalette(c("#e41a1c", "#0099ff", "#4daf4a", "#984ea3", "#ff7f00", "#ff96ca", "#a65628"))(n), 
               a))
}

# ---- the ggplot2 color scheme -----------------------------------------------------

ggcolor = function(n, a = 1)
{
  require(grDevices)
  require(scales)
  return(alpha(hcl(h = seq(15, 375, length = n + 1), 
                   l = 65, c = 100)[1:n], 
               a))
}

# ---- prints out a dat file object in ampl syntax ----------------------------------

ampl = function(dat, object = "param", name = "c")
{
  tochar = function(dat)
  {
    require(data.table)
    
    # make dat into a data.frame
    dat = data.table(dat)
    
    # get the column names
    column = names(dat)
    
    # get the values in the columns and convert them to character data types
    values = lapply(1:ncol(dat), function(i) as.character(unname(unlist(dat[, i, with = FALSE]))))
    
    # combine the values back into a data.frame
    dat = data.table(do.call("cbind", values), stringsAsFactors = FALSE)
    
    # give dat its column names
    setnames(dat, column)
    
    return(dat)
  }
  
  # make sure the data is a data frame object
  dat = tochar(dat)
  
  # every parameter/set object in an ampl dat file must end with a semicolon
  # so set up 1 semicolon to give to dat
  semicolon = c(";", rep(" ", ncol(dat) - 1))
  
  # add this semicolon as the last row of the data frame
  result = data.frame(rbind(dat, semicolon))
  
  # every parameter/set object in an ample dat file must begin with the name of the object and what it equals
  # for example: param c := 
  # so set up a header to give to dat
  header = c(paste(object, name, ":="), rep(" ", ncol(dat) - 1))
  
  # update the column names of dat to be the header we created
  colnames(result) = header
  
  # print out the result without any row names
  # print out the result left adjusted
  # print(result, right = FALSE, row.names = FALSE)
  
  return(result)	
}

# ---- compares the quantiles of emprical data against the quantiles of any statistical distribution 

ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), alpha = 0.33, basefont = 20, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
{
  require(ggplot2)
  
  # compute the sample quantiles and theoretical quantiles
  q.function = eval(parse(text = paste0("q", distribution)))
  d.function = eval(parse(text = paste0("d", distribution)))
  x = na.omit(x)
  ord = order(x)
  n = length(x)
  P = ppoints(length(x))
  df = data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  # compute the quantile line
  Q.x = quantile(df$ord.x, c(probs[1], probs[2]))
  Q.z = q.function(c(probs[1], probs[2]), ...)
  b = diff(Q.x) / diff(Q.z)
  coef = c(Q.x[1] - (b * Q.z[1]), b)
  
  # compute the confidence interval band
  zz = qnorm(1 - (1 - conf) / 2)
  SE = (coef[2] / d.function(df$z, ...)) * sqrt(P * (1 - P) / n)
  fit.value = coef[1] + (coef[2] * df$z)
  df$upper = fit.value + (zz * SE)
  df$lower = fit.value - (zz * SE)
  
  # plot the qqplot
  p = ggplot(df, aes(x = z, y = ord.x)) + 
    geom_point(color = "blue", alpha = alpha) +
    geom_abline(intercept = coef[1], slope = coef[2], size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    coord_cartesian(ylim = c(min(df$ord.x), max(df$ord.x))) + 
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # conditional additions
  if(main != "")(p = p + ggtitle(main))
  
  return(p)
}

# ---- plots 4 residual plots -------------------------------------------------------

residplots = function(actual, fit, binwidth = NULL, from = NULL, to = NULL, by = NULL, histlabel.y = -10, basefont = 20)
{
  require(ggplot2)
  
  residual = actual - fit 
  DF = data.frame("actual" = actual, "fit" = fit, "residual" = residual)
  
  rvfPlot = ggplot(DF, aes(x = fit, y = residual)) + 
    geom_point(na.rm = TRUE) +
    stat_smooth(method = "loess", se = FALSE, na.rm = TRUE, color = "blue") +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Fitted values") +
    ylab("Residuals") +
    ggtitle("Residual vs Fitted Plot") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), note = TRUE, alpha = 0.33, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
  {
    # compute the sample quantiles and theoretical quantiles
    q.function = eval(parse(text = paste0("q", distribution)))
    d.function = eval(parse(text = paste0("d", distribution)))
    x = na.omit(x)
    ord = order(x)
    n = length(x)
    P = ppoints(length(x))
    DF = data.frame(ord.x = x[ord], z = q.function(P, ...))
    
    # compute the quantile line
    Q.x = quantile(DF$ord.x, c(probs[1], probs[2]))
    Q.z = q.function(c(probs[1], probs[2]), ...)
    b = diff(Q.x) / diff(Q.z)
    coef = c(Q.x[1] - (b * Q.z[1]), b)
    
    # compute the confidence interval band
    zz = qnorm(1 - (1 - conf) / 2)
    SE = (coef[2] / d.function(DF$z, ...)) * sqrt(P * (1 - P) / n)
    fit.value = coef[1] + (coef[2] * DF$z)
    DF$upper = fit.value + (zz * SE)
    DF$lower = fit.value - (zz * SE)
    
    # plot the qqplot
    p = ggplot(DF, aes(x = z, y = ord.x)) + 
      geom_point(color = "black", alpha = alpha) +
      geom_abline(intercept = coef[1], slope = coef[2], size = 1, color = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
      coord_cartesian(ylim = c(min(DF$ord.x), max(DF$ord.x))) + 
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # conditional additions
    if(main != "")(p = p + ggtitle(main))
    
    return(p)
  }
  
  qqPlot = ggqq(residual, 
                alpha = 1,				  
                main = "Normal Q-Q Plot", 
                xlab = "Theoretical Quantiles", 
                ylab = "Residuals")
  
  rvtPlot = ggplot(data.frame("x" = 1:length(DF$residual), "y" = DF$residual), aes(x = x, y = y)) + 
    geom_line(na.rm = TRUE) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Obs. Number") +
    ylab("Residuals") +
    ggtitle("Residual Time Series") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  test = t.test(DF$residual)
  
  CI = data.frame("x" = test$estimate, 
                  "LCB" = test$conf.int[1], 
                  "UCB" = test$conf.int[2], 
                  row.names = 1)
  
  histPlot = ggplot(DF, aes(x = residual)) +
    geom_histogram(color = "white", fill = "black", binwidth = binwidth) +
    geom_segment(data = CI, aes(x = LCB, xend = LCB, y = 0, yend = Inf), color = "blue") +
    geom_segment(data = CI, aes(x = UCB, xend = UCB, y = 0, yend = Inf), color = "blue") +
    annotate("text", x = CI$x, y = histlabel.y, 
             label = "T-Test C.I.", size = 5, 
             color = "blue", fontface = 2) + 
    ggtitle("Residual Histogram") +
    labs(x = "Residuals", y = "Frequency") +
    theme_bw(base_size = basefont) +
    theme(legend.key.size = unit(.25, "in"),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(class(from) != "NULL" & class(to) != "NULL" & class(by) != "NULL") (histPlot = histPlot + scale_x_continuous(breaks = seq(from = from, to = to, by = by)))
  
  return(list("rvfPlot" = rvfPlot, 
              "qqPlot" = qqPlot, 
              "rvtPlot" = rvtPlot,  
              "histPlot" = histPlot))
}

# ---- builds a square confusion matrix ---------------------------------------------

confusion = function(ytrue, ypred)
{
  require(gtools)
  
  # make predicted and actual vectors into factors, if they aren't already
  if(class(ytrue) != "factor") ytrue = factor(ytrue)
  if(class(ypred) != "factor") ypred = factor(ypred)
  
  # combine their levels into one unique set of levels
  common.levels = mixedsort(unique(c(levels(ytrue), levels(ypred))))
  
  # give each vector the same levels
  ytrue = factor(ytrue, levels = common.levels)
  ypred = factor(ypred, levels = common.levels)
  
  # build the confusion matrix
  output = table(ytrue, ypred)
  
  # add a per class accuracy column
  output = cbind(output, "Accuracy" = unname(100 * (diag(output) / rowSums(output))))
  
  # name the dimensions (a name for the rows, and a name for the columns)
  names(dimnames(output)) = c("Actual", "Predicted")
  
  # return the confusion matrix
  return(output)
}

# ---- generates a logarithmically spaced sequence ----------------------------------

lseq = function(from, to, length.out)
{
  return(exp(seq(log(from), log(to), length.out = length.out)))
}

# ---- compute weight sum of squares for an hclust model ----------------------------

WSS = function(k, hc, dat) 
{
  # compute wss for model hc grouping data d into k clusters 
  wss = function(d) 
  {
    return(sum(scale(d, scale = FALSE)^2))
  }
  
  # determine the clusters
  cl = cutree(hc, k)
  
  # split the data accroding to the clusters
  spl = split(dat, cl)
  
  # compute the value of wss 
  value = sum(sapply(spl, wss))
  
  return(value)
}

}

# -----------------------------------------------------------------------------------
# ---- Prepare the Data -------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# choose how much of your CPU to use (cpu = 1/2 -> 50% use)
# cpu = 2/3
cpu = 0
workers = max(c(1, round(cpu * getDTthreads(), 0)))

# do we need to join the article data or has it been done already?
join.data = FALSE

if(join.data)
{
  # ---- Science Direct ----
  
  # set the work directory
  setwd(mywd1)
  
  # get the names of the xml files
  filenames = list.files()
  filenames = filenames[which(grepl(".xml", filenames))]
  
  # determine the number of files to join
  tasks = length(filenames)
  
  # create a log file to keep track of parallel computing progress
  myfile = "log.txt"
  file.create(myfile)
  
  # set up a cluster if workers > 1, otherwise don't set up a cluster
  if(workers > 1)
  {
    # setup parallel processing
    cl = makeCluster(workers, type = "SOCK", outfile = "")
    registerDoSNOW(cl)
    
    # define %dopar%
    `%fun%` = `%dopar%`
    
    # write out start time to log file
    sink(myfile, append = TRUE)
    cat("\n------------------------------------------------\n")
    cat("joining articles from Science Direct\n")
    cat(paste(workers, "workers started at", Sys.time(), "\n"))
    sink()
    
  } else
  {
    # define %do%
    `%fun%` = `%do%`
    
    # write out start time to log file
    sink(myfile, append = TRUE)
    cat("\n------------------------------------------------\n")
    cat("joining articles from Science Direct\n")
    cat(paste("task 1 started at", Sys.time(), "\n"))
    sink()
  }
  
  # loop through each xml file and collect article information
  SCID = foreach(i = 1:tasks) %fun%
  {
    # load packages we need for our tasks
    require(data.table)
    require(xml2)
    
    # read in xml data
    XML = read_xml(file(filenames[i]))
    
    # convert xml data into a list
    XML = as_list(XML)
    
    # get the total number of articles in XML
    n = length(XML$Sources)
    
    # get information on publication year, SOURCE, title, and abstract for each article in XML
    YEAR = unlist(lapply(1:n, function(k) 
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$Year[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, as.numeric(value))
      
      return(value)
    }))
    
    SOURCE = unlist(lapply(1:n, function(k) 
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$JournalName[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, value)
      
      return(value)
    }))
    
    TITLE = unlist(lapply(1:n, function(k)
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$Title[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, value)
      
      return(value)
    }))
    
    ABSTRACT = unlist(lapply(1:n, function(k) 
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$BIBTEX_Abstract[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, value)
      
      return(value)
    }))
    
    # combine the information into a table
    DT = data.table(DATABASE = "Science Direct", YEAR = YEAR, SOURCE = SOURCE, TITLE = TITLE, ABSTRACT = ABSTRACT)
    
    # export progress information
    sink(myfile, append = TRUE)
    cat(paste("task", i, "of", tasks, "finished at", Sys.time(), "\n"))
    sink()
    
    return(DT)
  }
  
  # write out end time to log file
  sink(myfile, append = TRUE)
  cat(paste(tasks, "tasks finished at", Sys.time(), "\n"))
  sink()
  
  # end the cluster if it was set up
  if(workers > 1)
  {
    stopCluster(cl)
  }
  
  # combine the list of tables into one table
  SCID = rbindlist(SCID)
  
  # ---- JSTOR ---- 
  
  # set the work directory
  setwd(mywd2)
  
  # get the names of the xml files
  filenames = list.files()
  filenames = filenames[which(grepl(".xml", filenames))]
  
  # determine the number of files to join
  tasks = length(filenames)
  
  # create a log file to keep track of parallel computing progress
  myfile = "log.txt"
  file.create(myfile)
  
  # set up a cluster if workers > 1, otherwise don't set up a cluster
  if(workers > 1)
  {
    # setup parallel processing
    cl = makeCluster(workers, type = "SOCK", outfile = "")
    registerDoSNOW(cl)
    
    # define %dopar%
    `%fun%` = `%dopar%`
    
    # write out start time to log file
    sink(myfile, append = TRUE)
    cat("\n------------------------------------------------\n")
    cat("joining articles from JSTOR\n")
    cat(paste(workers, "workers started at", Sys.time(), "\n"))
    sink()
    
  } else
  {
    # define %do%
    `%fun%` = `%do%`
    
    # write out start time to log file
    sink(myfile, append = TRUE)
    cat("\n------------------------------------------------\n")
    cat("joining articles from JSTOR\n")
    cat(paste("task 1 started at", Sys.time(), "\n"))
    sink()
  }
  
  # loop through each xml file and collect article information
  JSR = foreach(i = 1:tasks) %fun%
  {
    # load packages we need for our tasks
    require(data.table)
    require(xml2)
    
    # read in xml data
    XML = read_xml(file(filenames[i]))
    
    # convert xml data into a list
    XML = as_list(XML)
    
    # get information on publication year, SOURCE, title, and abstract for each article in XML
    YEAR = tryCatch(XML$article$front$`article-meta`$`pub-date`$year[[1]], error = function(...) return(NULL))
    SOURCE = tryCatch(XML$article$front$`journal-meta`$`journal-title-group`$`journal-title`[[1]], error = function(...) return(NULL))
    TITLE = tryCatch(XML$article$front$`article-meta`$`title-group`$`article-title`[[1]], error = function(...) return(NULL))
    ABSTRACT = tryCatch(XML$article$front$`article-meta`$abstract$p[[1]], error = function(...) return(NULL))
    
    # replace NULL's with NA's
    YEAR = ifelse(is.null(YEAR), NA, as.numeric(YEAR))
    SOURCE = ifelse(is.null(SOURCE), NA, SOURCE)
    TITLE = ifelse(is.null(TITLE), NA, TITLE)
    ABSTRACT = ifelse(is.null(ABSTRACT), NA, ABSTRACT)
    
    # combine the information into a table
    DT = data.table(DATABASE = "JSTOR", YEAR = YEAR, SOURCE = SOURCE, TITLE = TITLE, ABSTRACT = ABSTRACT)
    
    # export progress information
    sink(myfile, append = TRUE)
    cat(paste("task", i, "of", tasks, "finished at", Sys.time(), "\n"))
    sink()
    
    return(DT)
  }
  
  # write out end time to log file
  sink(myfile, append = TRUE)
  cat(paste(tasks, "tasks finished at", Sys.time(), "\n"))
  sink()
  
  # end the cluster if it was set up
  if(workers > 1)
  {
    stopCluster(cl)
  }
  
  # combine the list of tables into one table
  JSR = rbindlist(JSR)
  
  # ---- IEEE ---- 
  
  # set the work directory
  setwd(mywd3)
  
  # get the names of the xml files
  filenames = list.files()
  filenames = filenames[which(grepl(".xml", filenames))]
  
  # determine the number of files to join
  tasks = length(filenames)
  
  # create a log file to keep track of parallel computing progress
  myfile = "log.txt"
  file.create(myfile)
  
  # set up a cluster if workers > 1, otherwise don't set up a cluster
  if(workers > 1)
  {
    # setup parallel processing
    cl = makeCluster(workers, type = "SOCK", outfile = "")
    registerDoSNOW(cl)
    
    # define %dopar%
    `%fun%` = `%dopar%`
    
    # write out start time to log file
    sink(myfile, append = TRUE)
    cat("\n------------------------------------------------\n")
    cat("joining articles from IEEE\n")
    cat(paste(workers, "workers started at", Sys.time(), "\n"))
    sink()
    
  } else
  {
    # define %do%
    `%fun%` = `%do%`
    
    # write out start time to log file
    sink(myfile, append = TRUE)
    cat("\n------------------------------------------------\n")
    cat("joining articles from IEEE\n")
    cat(paste("task 1 started at", Sys.time(), "\n"))
    sink()
  }
  
  # loop through each xml file and collect article information
  IEX = foreach(i = 1:tasks) %fun%
  {
    # load packages we need for our tasks
    require(data.table)
    require(xml2)
    
    # read in xml data
    XML = read_xml(file(filenames[i]))
    
    # convert xml data into a list
    XML = as_list(XML)
    
    # get the total number of articles in XML
    n = length(XML$Sources)
    
    # get information on publication year, SOURCE, title, and abstract for each article in XML
    YEAR = unlist(lapply(1:n, function(k) 
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$Year[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, as.numeric(value))
      
      return(value)
    }))
    
    SOURCE = unlist(lapply(1:n, function(k) 
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$ConferenceName[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, value)
      
      return(value)
    }))
    
    TITLE = unlist(lapply(1:n, function(k)
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$Title[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, value)
      
      return(value)
    }))
    
    ABSTRACT = unlist(lapply(1:n, function(k) 
    {
      # get the value for article k
      value = tryCatch(XML$Sources[[k]]$BIBTEX_Abstract[[1]], error = function(...) return(NULL))
      
      # if the value is a NULL, make it a NA
      value = ifelse(is.null(value), NA, value)
      
      return(value)
    }))
    
    # combine the information into a table
    DT = data.table(DATABASE = "IEEE", YEAR = YEAR, SOURCE = SOURCE, TITLE = TITLE, ABSTRACT = ABSTRACT)
    
    # export progress information
    sink(myfile, append = TRUE)
    cat(paste("task", i, "of", tasks, "finished at", Sys.time(), "\n"))
    sink()
    
    return(DT)
  }
  
  # write out end time to log file
  sink(myfile, append = TRUE)
  cat(paste(tasks, "tasks finished at", Sys.time(), "\n"))
  sink()
  
  # end the cluster if it was set up
  if(workers > 1)
  {
    stopCluster(cl)
  }
  
  # combine the list of tables into one table
  IEX = rbindlist(IEX)
  
  # ---- Web of Science ---- 
  
  # set the work directory
  setwd(mywd4)
  
  # get the names of the csv files
  filenames = list.files()
  filenames = filenames[which(grepl(".csv", filenames))]
  
  # determine the number of files to join
  tasks = length(filenames)
  
  # read in csv data
  WOS = lapply(1:tasks, function(i)
  {
    # read in csv data
    DT = fread(filenames[i])
    
    # get information on publication year, SOURCE, title, and abstract for each article in XML
    DT = DT[,.(DATABASE = "Web of Science", YEAR = PY, SOURCE = SO, TITLE = TI, ABSTRACT = AB)]
    
    return(DT)
  })
  
  # combine the list of tables into one table
  WOS = rbindlist(WOS)
  
  # ---- Engineering Village ---- 
  
  # set the work directory
  setwd(mywd5)
  
  # get the names of the csv files
  filenames = list.files()
  filenames = filenames[which(grepl(".csv", filenames))]
  
  # determine the number of files to join
  tasks = length(filenames)
  
  # read in csv data
  EV = lapply(1:tasks, function(i)
  {
    # read in csv data
    DT = fread(filenames[i])
    
    # get information on publication year, SOURCE, title, and abstract for each article in XML
    DT = DT[,.(DATABASE = "Engineering Village", YEAR = `Publication year`, SOURCE = Source, TITLE = Title, ABSTRACT = Abstract)]
    
    return(DT)
  })
  
  # combine the list of tables into one table
  EV = rbindlist(EV)
  
  # combine all database tables into one single table
  dat = rbind(SCID, WOS, IEX, EV, JSR)
  
  # remove any duplicates
  dat = dat[!duplicated(dat[,!"DATABASE"]),]
  
  # remove any papers without an abstract
  dat = dat[!is.na(ABSTRACT)]
  dat = dat[!(ABSTRACT == "")]
  
  # write out dat
  setwd(mywd0)
  fwrite(dat, "vaccine budget uncertainty abstracts.csv")
  
  # remove objects we no longer need
  rm(`%fun%`, ABSTRACT, DT, EV, filenames, i, IEX, 
     JSR, myfile, n, SCID, SOURCE, tasks, TITLE, WOS, XML, YEAR)
  
  # free up RAM
  gc()
  
} else
{
  # read in the article data
  setwd(mywd0)
  dat = fread("vaccine budget uncertainty abstracts.csv")
}

}

# -----------------------------------------------------------------------------------
# ---- Model the Data ---------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# do we need to model the data or has it already been done?
model.data = FALSE

if(model.data)
{
  # ---- plot the original data ----
  
  # keep an original version of dat
  dat.original = data.table(dat)
  
  # get counts on the number of articles by DATABASE and YEAR
  counts = data.table(dat[, .(COUNT = .N), by = .(YEAR, DATABASE)])
  
  # order counts by YEAR and COUNT
  counts = counts[order(YEAR, -COUNT)]
  
  # plot counts v. DATABASE and YEAR
  counts.plot = ggplot(data = counts, aes(x = YEAR, y = COUNT, color = DATABASE, group = DATABASE)) + 
    geom_line(size = 1.25) + 
    geom_point(size = 2) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 1) + 
    scale_y_continuous(label = comma) + 
    ggtitle("Query Results\n'vaccine budget uncertainty'") + 
    labs(x = "Years", y = "No. of Papers", color = "Database") + 
    facet_wrap(~DATABASE, scales = "free_y") + 
    theme_bw(base_size = 25) +
    theme(legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  counts.plot
  
  # ---- build a term-document matrix ----
  
  # build a corpus of ABSTRACT, but first remove special characters
  corp = gsub("[^0-9A-Za-z///' ]", " ", dat$ABSTRACT, ignore.case = TRUE)
  corp = Corpus(VectorSource(corp))
  
  # remove punctuation, extra white space, and stop words
  # convert words to lowercase letters
  # stem words to their root letters
  corp = tm_map(corp, removePunctuation)
  corp = tm_map(corp, stripWhitespace)
  corp = tm_map(corp, content_transformer(tolower))
  corp = tm_map(corp, removeWords, stopwords(kind = "english"))
  corp = tm_map(corp, removeWords, stopwords(kind = "SMART"))
  corp = tm_map(corp, stemDocument, language = "english")
  
  # build a term document matrix (tdm) on the corpus
  tdm = TermDocumentMatrix(corp)
  inspect(tdm)
  
  # get the matrix form of tdm
  mat = as.matrix(tdm)
  
  # get the count for each word
  counts = rowSums(mat)
  
  # determine the 1st and 99th percentiles of counts
  cuts = round(as.numeric(quantile(x = counts, probs = c(0.01, 0.99))), 0)
  
  # update the minimum cut to be at least 5
  cuts[1] = max(c(cuts[1], 5))
  
  # determine which words to keep
  keep = as.numeric(which(counts >= cuts[1] & counts <= cuts[2]))
  keep = as.numeric(which(counts >= cuts[1]))
  
  # determine which words to remove
  rmwords = names(counts)[-keep]
  
  # set up the indexing for a k set of words within rmwords
  kmin = 1
  kinc = 50
  kmax = kmin + kinc - 1
  total = length(rmwords)
  
  # remove subsets of words until we have iterated through all words in rmwords
  while(kmin < total)
  {
    # get k words
    k = rmwords[kmin:kmax]
    
    # remove k words from corp
    corp = tm_map(corp, removeWords, k)
    
    # index for the next k words
    kmin = kmin + kinc
    kmax = min(c(total, kmax + kinc))
  }
  
  # update dat with WORDS
  dat[, WORDS := sapply(corp, as.character)]
  
  # build a term document matrix (tdm) on the corpus
  tdm = TermDocumentMatrix(corp)
  inspect(tdm)
  
  # ---- cluster documents by cosine distance ----
  
  # compute the cosine distance between all documents
  cdm = 1 - crossprod_simple_triplet_matrix(tdm)/(sqrt(col_sums(tdm^2) %*% t(col_sums(tdm^2))))
  
  # https://stackoverflow.com/questions/25959385/how-to-draw-the-plot-of-within-cluster-sum-of-squares-for-a-cluster
  # https://cran.r-project.org/web/packages/textmineR/vignettes/b_document_clustering.html
  
  # convert cdm to a distance matrix type
  dmat = as.dist(cdm)
  
  # plot dmat
  distance.matrix.plot = fviz_dist(dmat, gradient = list(low = "blue", mid = "white", high = "red")) + 
    ggtitle("\nAbstract Cosine Distance Matrix") + 
    theme_void(25) +
    labs(fill = "Cosine Distance") + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top", 
          legend.key.size = unit(.25, "in"),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  # distance.matrix.plot
  
  # determine which hierarchical clustering method best summarizes dmat
  # choice of hclust methods
  methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  
  # build cluster methods
  hc = lapply(methods, function(j) hclust(dmat, method = j))
  
  # compute correlations
  cors = rbindlist(lapply(1:length(methods), function(j) 
    data.table(method = methods[j], value = cor(dmat, cophenetic(hc[[j]])))))
  
  # order by value
  cors = cors[order(value, decreasing = TRUE)]
  
  # make method a factor for plotting purposes
  cors[, method := factor(method, levels = unique(method))]
  
  # plot correlation v. method
  cophenetic.plot = ggplot(data = cors, aes(x = method, y = value, color = value, group = 1)) + 
    geom_smooth(size = 1.5, method = 'loess', color = "black", linetype = "dashed", fill = NA) + 
    geom_point(size = 7) + 
    scale_color_continuous(low = "royalblue", high = "orangered") + 
    scale_y_continuous(label = comma) + 
    ggtitle("Hierarchical Cluster Analysis") + 
    labs(x = "Link Method", y = "Cophenetic Correlations", color = "Correlation") + 
    theme_bw(base_size = 25) +
    theme(legend.position = "none", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  cophenetic.plot
  
  # lets go with the best method
  hc.method = as.character(cors$method[1])
  
  # update the hc model
  hc = hc[[which(methods == hc.method)]]
  
  # determine a set of clusters to try
  min.k = 4  # lower bound of the set
  max.k = 100  # upper bound of the set
  tot.k = 20  # size of the set
  set.k = unique(round(lseq(from = min.k, to = max.k, length.out = tot.k), 0))  # a log-spaced sequence
  
  # compute WSS for each cluster in set.k
  hclust.wss = sapply(set.k, function(i) WSS(k = i, hc = hc, dat = cdm))
  
  # make hclust.wss into a table
  hclust.wss = data.table(k = set.k, tot_withinss = hclust.wss)
  
  # plot tot_withinss v. k
  hc.wss.plot = ggplot(data = hclust.wss, aes(x = k, y = tot_withinss, color = tot_withinss)) + 
    geom_smooth(size = 1.5, method = 'loess', color = "black", linetype = "dashed", fill = NA) + 
    geom_point(size = 7) + 
    scale_color_continuous(low = "royalblue", high = "orangered") + 
    scale_y_continuous(label = comma) + 
    ggtitle("Hierarchical Cluster Analysis") + 
    labs(x = "No. Clusters", y = "Total Within Sum of Squares", color = "TWSS") + 
    theme_bw(base_size = 25) +
    theme(legend.position = "none", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  hc.wss.plot
  
  # choose a number of clusters to move forward with
  num.clusters = 84
  clustering = cutree(hc, num.clusters)
  
  # plot the Dendrogram
  hcd = as.dendrogram(hc)
  plot(hcd, main = "Hierarchical Clustering of Abstracts",
       ylab = "", xlab = "", yaxt = "n")
  rect.hclust(hc, num.clusters, border = "red")
  
  # compute the presence of each word, proportional to eachother, across all documents
  tdm.mat = t(as.matrix(tdm))
  p_words = colSums(tdm.mat) / sum(tdm.mat)
  
  # determine the presence of each word in each cluster
  cluster_words = lapply(sort(unique(clustering)), function(x)
  {
    # get the documents which appear in cluster x
    rows = as.matrix(tdm.mat[unique(which(clustering == x)),], ncol = ncol(tdm.mat))
    
    # get the presence of each word in cluster x
    value = (tryCatch(colSums(rows), error = function(...) return(0)) / sum(rows)) - p_words[colnames(rows)]
    
    return(value)
  })
  
  # build a summary table for the most frequent words in each cluster
  cluster_summary = data.table(cluster = sort(unique(clustering)),
                               size = as.numeric(table(clustering)),
                               top_words = sapply(cluster_words, function(d){
                                 paste(
                                   names(d)[ order(d, decreasing = TRUE) ][1:5], 
                                   collapse = ", ")
                               }),
                               stringsAsFactors = FALSE)
  
  cluster_summary[size >= 10]
  
  # build a word cloud for words that appear frequently in each cluster
  cluster.pick = which.max(cluster_summary$size)[1]
  cluster.pick = 1
  hc.wc = wordcloud(words = names(cluster_words[[cluster.pick]]), scale = c(12, 1.5),
                    freq = cluster_words[[cluster.pick]], 
                    max.words = 50, 
                    random.order = FALSE, 
                    colors = c("royalblue", "goldenrod2", "orangered"),
                    main = paste("Top words in Cluster", cluster.pick))
  
  hc.wc
  
  # convert tdm to have values based on term frequency - inverse document frequency (tfidf)
  tdm = weightTfIdf(tdm)
  inspect(tdm)
  
  # convert tdm to a matrix
  tdm = as.matrix(tdm)
  
  # transpose the matrix so documents are rows and words are columns
  tdm = t(tdm)
  
  # make tdm into a data table
  tdm = data.table(tdm)
  
  # add the clusters to tdm
  tdm[, Cluster := as.numeric(clustering)]
  
  # determine which clusters to keep based on their size
  min.size = 10
  keep.clusters = unname(unlist(cluster_summary[size >= min.size, .(cluster)]))
  
  # remove clusters that fail to satisfy min.size
  tdm[!(Cluster %in% keep.clusters), Cluster := NA]
  
  # split up tdm into clean and missing versions
  tdm.clean = data.table(na.omit(tdm))
  tdm.missing = data.table(na.omit(tdm, invert = TRUE))
  
  # ---- rank words with random forests ----
  
  # choose the number of workers for parallel processing
  workers = getDTthreads() - 2
  
  # initialize the h2o instance
  h2o.init(nthreads = workers, max_mem_size = "9g")
  
  # remove any objects in the h2o instance
  h2o.removeAll()
  
  # remove the progress bar when model building
  h2o.no_progress()
  
  # convert Cluster to a factor in tdm.clean
  tdm.clean[, Cluster := factor(Cluster, levels = unique(Cluster))]
  
  # identify predictors (x) and response (y)
  y = "Cluster"
  x = names(tdm.clean[, !"Cluster", with = FALSE])
  
  # make an h2o version of tdm.clean
  tdm.clean.h2o = as.h2o(tdm.clean)
  
  # compute the max class weight
  max.class.weight = table(unname(unlist(tdm.clean[, y, with = FALSE])))
  max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  
  # set up hyperparameters of interest
  rf.hyper.params = list(ntrees = 250,
                         # min_rows = c(1, 11, 25),
                         # max_depth = c(20, 40, 60),
                         min_rows = 11,
                         max_depth = 40,
                         stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 30
  rf.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # lets run a grid search for a good model
  h2o.rm("rf.random.grid")
  rf.random.grid = h2o.grid(algorithm = "randomForest",
                            grid_id = "rf.random.grid",
                            y = y,
                            x = x,
                            training_frame = tdm.clean.h2o,
                            # stopping_rounds = 3,
                            # histogram_type = "RoundRobin",
                            nfolds = 3,
                            fold_assignment = "Stratified",
                            seed = 3,
                            # score_each_iteration = TRUE,
                            balance_classes = TRUE,
                            max_after_balance_size = max.class.weight,
                            hyper_params = rf.hyper.params,
                            search_criteria = rf.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # pick the top model from all grid searches
  imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
  
  # extract variable importance
  imp = data.table(imp.rf@model$variable_importances)
  
  # make variable into a factor
  imp[, variable := factor(variable, levels = unique(variable))]
  
  # pick a cut off value
  # this value should show where importance drops the most (look for the center of the "knee")
  cutoff = 0.08
  
  # plot a barplot of variable importance
  imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
    ggtitle("GINI Importance\nRandom Forests") +
    labs(x = "Variable", y = "Scaled Importance") +
    scale_y_continuous(labels = percent) +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_color_gradient(low = "yellow", high = "red") +
    theme_dark(25) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.grid.major.x = element_blank())
  
  # plot a density plot of variable importance
  imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
    geom_density(fill = "cornflowerblue", alpha = 2/3) +
    geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
    ggtitle("GINI Importance\nRandom Forests") +
    labs(x = "Scaled Importance", y = "Density") +
    scale_x_continuous(labels = percent) +
    coord_flip() + 
    theme_bw(25) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # show the importance plots
  grid.arrange(imp.plot1, imp.plot2, nrow = 1)
  
  # check out imp
  imp
  
  # find which indicators meet the cutoff value
  keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))
  
  # reverse the order of the levels of variable
  imp[, variable := factor(variable, levels = rev(levels(variable)))]
  
  # plot imp
  imp.plot3 = ggplot(imp[variable %in% keep.indicators], aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
    ggtitle("Random Forests: GINI Importance") +
    labs(x = "Variable", y = "Scaled Importance") +
    scale_y_continuous(labels = percent) +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_color_gradient(low = "yellow", high = "red") +
    coord_flip() + 
    theme_dark(20) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  imp.plot3
  
  # update tdm.clean and tdm.missing to only contain indicators of interest
  tdm.clean = tdm.clean[, c(keep.indicators, y), with = FALSE]
  tdm.missing = tdm.missing[, c(keep.indicators), with = FALSE]
  
  # ---- learn clusters with random forests ----
  
  # make an h2o version of tdm.clean and tdm.missing
  tdm.clean.h2o = as.h2o(tdm.clean)
  tdm.missing.h2o = as.h2o(tdm.missing)
  
  # update features
  x = names(tdm.missing)
  
  # set up hyperparameters of interest
  rf.hyper.params = list(ntrees = 250,
                         # min_rows = c(1, 11, 25),
                         # max_depth = c(20, 40, 60),
                         min_rows = 11,
                         max_depth = 40,
                         stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 30
  rf.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("rf.random.grid")
  rf.random.grid = h2o.grid(algorithm = "randomForest",
                            grid_id = "rf.random.grid",
                            y = y,
                            x = x,
                            training_frame = tdm.clean.h2o,
                            # stopping_rounds = 3,
                            # histogram_type = "RoundRobin",
                            nfolds = 3,
                            fold_assignment = "Stratified",
                            seed = 3,
                            # score_each_iteration = TRUE,
                            balance_classes = TRUE,
                            max_after_balance_size = max.class.weight,
                            hyper_params = rf.hyper.params,
                            search_criteria = rf.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # pick the top model from all grid searches
  imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
  
  # predict clusters for tdm.clean
  clean.predictions = as.data.frame(predict(imp.rf, newdata = tdm.clean.h2o[-which(names(tdm.clean.h2o) == y)])[,1])[,1]
  
  # create a confusion matrix
  clean.conf = confusion(ytrue = tdm.clean$Cluster, ypred = clean.predictions)
  clean.conf
  
  # predict clusters for tdm.missing
  predictions = as.data.frame(predict(imp.rf, newdata = tdm.missing.h2o)[,1])[,1]
  tdm.missing[, Cluster := predictions]
  
  # give tdm an ID column
  tdm[, ID := 1:nrow(tdm)]
  
  # add ID to tdm.clean and tdm.missing
  tdm.clean = cbind(tdm.clean, na.omit(tdm)[,.(ID)])
  tdm.missing = cbind(tdm.missing, na.omit(tdm, invert = TRUE)[,.(ID)])
  
  # get tdm back without any missing values
  tdm = data.table(rbind(tdm.clean, tdm.missing))
  
  # order tdm by ID and add Cluster to dat.original
  tdm = tdm[order(ID)]
  dat.original[, CLUSTER := tdm$Cluster]
  
  # remove ID from tdm
  tdm[, ID := NULL]
  
  # ---- summarize clusters ----
  
  # compute the presence of each word, proportional to eachother, across all documents
  clustering = tdm$Cluster
  p_words = colSums(tdm.mat) / sum(tdm.mat)
  
  # determine the presence of each word in each cluster
  cluster_words = lapply(sort(unique(clustering)), function(x)
  {
    # get the documents which appear in cluster x
    rows = as.matrix(tdm.mat[unique(which(clustering == x)),], ncol = ncol(tdm.mat))
    
    # get the presence of each word in cluster x
    value = (tryCatch(colSums(rows), error = function(...) return(0)) / sum(rows)) - p_words[colnames(rows)]
    
    return(value)
  })
  
  # build a summary table for the most frequent words in each cluster
  cluster_summary = data.table(cluster = sort(unique(clustering)),
                               size = as.numeric(table(clustering)),
                               top_words = sapply(cluster_words, function(d){
                                 paste(
                                   names(d)[ order(d, decreasing = TRUE) ][1:5], 
                                   collapse = ", ")
                               }),
                               stringsAsFactors = FALSE)
  
  cluster_summary
  
  # build a word cloud for words that appear frequently in each cluster
  cluster.pick = 11
  hc.wc = wordcloud(words = names(cluster_words[[cluster.pick]]), scale = c(7, 1),
                    freq = cluster_words[[cluster.pick]], 
                    max.words = 50, 
                    random.order = FALSE, 
                    colors = c("royalblue", "goldenrod2", "orangered"),
                    main = paste("Top words in Cluster", cluster.pick))
  
  hc.wc
  
  # ---- sample documents to read ----
  
  # create a column to sample documents by DATABASE and CLUSTER
  dat.original[, CATEGORY := as.factor(paste(DATABASE, CLUSTER, sep = "_"))]
  
  # whats the total number documents you would like to read to interpret the clustering?
  num.docs = 100
  
  # build the fold assignment
  set.seed(42)
  k.folds = round(nrow(dat.original) / num.docs, 0)
  folds = createFolds(y = dat.original$CATEGORY, k = k.folds)
  
  # determine the number of digits that define a fold number
  digits = nchar(k.folds)
  
  # make folds into a vector
  folds = sort(unlist(folds))
  
  # update the names of folds
  names(folds) = substr(names(folds), start = 1, stop = 4 + digits)
  
  # update folds to just be the names
  folds = names(folds)
  
  # add folds to dat.original
  dat.original[, FOLD := folds]
  
  # remove CATEGORY
  dat.original[, CATEGORY := NULL]
  
  # write out dat.original
  setwd(mywd0)
  fwrite(dat.original, "clustered vaccine budget uncertainty abstracts.csv")
  
} else
{
  # read in dat.original
  setwd(mywd0)
  dat.original = fread("clustered vaccine budget uncertainty abstracts.csv")
}
  

#########################################
# DONT RUN THE CODE BELOW
# CODE BELOW IS REFERENCE MATERIAL

# ---- Modeling Options ----

# build a word to vector model on WORDS?
build.w2v = FALSE

if(build.w2v)
{
  # get the words to model
  words = as.character(as.h2o(dat$WORDS))
  
  # build a word2vec (w2v) model on words
  w2v = h2o.word2vec(training_frame = words,
                     vec_size = 30,
                     min_word_freq = 5,
                     window_size = 5,
                     init_learning_rate = 0.025,
                     sent_sample_rate = 0.001,
                     epochs = 10) 
  
  # get values for each word from w2v
  
  
  # get synonyms for each word from w2v
  
  
  # join word values and synonyms onto tdm
  
}

# build a kmeans model on tfidf?
build.kmeans = TRUE

if(build.kmeans)
{
  # get an h2o version of tdm
  DT.h2o = as.h2o(cdm)
  
  # get the column names of DT.h2o
  x = names(DT.h2o)
  
  # determine a set of clusters to try
  min.k = 4  # lower bound of the set
  max.k = 36  # upper bound of the set
  tot.k = 12  # size of the set
  set.k = unique(round(lseq(from = min.k, to = max.k, length.out = tot.k), 0))  # a log-spaced sequence
  
  # set up hyperparameters of interest
  km.hyper.params = list(k = set.k, max_iterations = 30)
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 20
  km.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # lets run a grid search for a good model
  h2o.rm("km.random.grid")
  km.random.grid = h2o.grid(algorithm = "kmeans",
                            grid_id = "km.random.grid",
                            x = x,
                            training_frame = DT.h2o,
                            standardize = FALSE,
                            nfolds = 3,
                            fold_assignment = "Modulo",
                            seed = 21,
                            hyper_params = km.hyper.params,
                            search_criteria = km.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  km.grid = h2o.getGrid("km.random.grid", sort_by  = "tot_withinss", decreasing = FALSE)
  
  # get the summary table of the grid search
  DT.km.grid = data.table(km.grid@summary_table)
  
  # convert k, max_iterations, tot_withinss to numeric data types
  DT.km.grid[, k := as.numeric(k)]
  DT.km.grid[, max_iterations := as.numeric(max_iterations)]
  DT.km.grid[, tot_withinss := as.numeric(tot_withinss)]
  
  # plot tot_withinss v. k
  twss.plot = ggplot(data = DT.km.grid, aes(x = k, y = tot_withinss, color = tot_withinss)) + 
    geom_smooth(size = 1.5, method = 'loess', color = "black", linetype = "dashed", fill = NA) + 
    geom_point(size = 7) + 
    scale_color_continuous(low = "royalblue", high = "orangered") + 
    scale_y_continuous(label = comma) + 
    ggtitle("Cluster Analysis") + 
    labs(x = "No. Clusters", y = "Total Within Sum of Squares", color = "TWSS") + 
    theme_bw(base_size = 25) +
    theme(legend.position = "none", 
          legend.key.size = unit(.25, "in"),
          plot.title = element_text(hjust = 0.5),
          # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
  
  twss.plot
  
  # pick a model from all grid searches
  p = 3
  km.mod = h2o.getModel(km.grid@model_ids[[p]])
  
  # get clusters
  clus = predict(km.mod, newdata = DT.h2o)
  colnames(clus) = "Cluster"
}

# rank features based on cluster predictions?
rank.features = TRUE

if(rank.features)
{
  # identify predictors (x) and response (y)
  y = "Cluster"
  x = names(DT[, !"Cluster", with = FALSE])
  
  # make an h2o version of DT
  DT.h2o = as.h2o(DT)
  
  # compute the max class weight
  max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
  max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
  
  # set up hyperparameters of interest
  rf.hyper.params = list(ntrees = 250,
                         # min_rows = c(1, 11, 25),
                         # max_depth = c(20, 40, 60),
                         min_rows = 11,
                         max_depth = 40,
                         stopping_metric = "mean_per_class_error")
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 30
  rf.search.criteria = list(strategy = "RandomDiscrete", 
                            max_runtime_secs = minutes * 60, 
                            # max_models = 100, 
                            seed = 42)
  
  # lets run a grid search for a good model, without drop out ratios
  h2o.rm("rf.random.grid")
  rf.random.grid = h2o.grid(algorithm = "randomForest",
                            grid_id = "rf.random.grid",
                            y = y,
                            x = x,
                            training_frame = DT.h2o,
                            # stopping_rounds = 3,
                            # histogram_type = "RoundRobin",
                            # nfolds = 3,
                            # fold_assignment = "Stratified",
                            seed = 3,
                            # score_each_iteration = TRUE,
                            balance_classes = TRUE,
                            max_after_balance_size = max.class.weight,
                            hyper_params = rf.hyper.params,
                            search_criteria = rf.search.criteria)
  
  # free up RAM
  gc()
  
  # rank each model in the random grids
  rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
  
  # pick the top model from all grid searches
  imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
  
  # extract variable importance
  imp = data.table(imp.rf@model$variable_importances)
  
  # make variable into a factor
  imp[, variable := factor(variable, levels = unique(variable))]
  
  # pick a cut off value
  # this value should show where importance drops the most (look for the center of the "knee")
  cutoff = 0.25
  
  # plot a barplot of variable importance
  imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
    ggtitle("GINI Importance\nRandom Forests") +
    labs(x = "Variable", y = "Scaled Importance") +
    scale_y_continuous(labels = percent) +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_color_gradient(low = "yellow", high = "red") +
    theme_dark(25) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.grid.major.x = element_blank())
  
  # plot a density plot of variable importance
  imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
    geom_density(fill = "cornflowerblue", alpha = 2/3) +
    geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
    ggtitle("GINI Importance\nRandom Forests") +
    labs(x = "Scaled Importance", y = "Density") +
    scale_x_continuous(labels = percent) +
    coord_flip() + 
    theme_bw(25) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # show the importance plots
  grid.arrange(imp.plot1, imp.plot2, nrow = 1)
  
  # check out imp
  imp
  
  # reverse the order of the levels of variable
  imp[, variable := factor(variable, levels = rev(levels(variable)))]
  
  # plot imp
  cutoff = 0.25
  imp.plot3 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
    ggtitle("Random Forests: GINI Importance") +
    labs(x = "Variable", y = "Scaled Importance") +
    scale_y_continuous(labels = percent) +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_color_gradient(low = "yellow", high = "red") +
    coord_flip() + 
    theme_dark(20) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  imp.plot3
  
  # find which indicators meet the cutoff value
  keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))
  
  # find the top m indicators
  m = 15
  keep.indicators = as.character(imp$variable)[1:m]
  
  # update DT to only contain indicators of interest
  DT = DT[, c(y, keep.indicators), with = FALSE]
}

# run PCA on the data to create more features?
build.pca = TRUE

if(build.pca)
{
  # pick the number of PC to retain
  r = 5
  
  # build PC
  pca.mod = h2o.prcomp(training_frame = DT.h2o, x = x, seed = 21, 
                       k = r, pca_method = "GramSVD",
                       impute_missing = TRUE, transform = "NORMALIZE")
  
  
  # get eigenvalues of pca.mod
  eigens = pca.mod@model$model_summary[1,]^2
  row.names(eigens) = "Eigenvalue"
  
  # show importance of PC
  pca.imp = rbind(pca.mod@model$model_summary, eigens)
  pca.imp
  
  # get PC
  pca = predict(pca.mod, newdata = DT.h2o)
  
  # add clus to pca
  pca = as.data.table(h2o.cbind(pca, clus))
  
  # make Cluster into a factor
  pca[, Cluster := factor(Cluster, levels = sort(unique(pca$Cluster)))]
  
  # build a simple plot to grab a legend from
  legend.plot = ggplot(pca, aes(x = PC1, y = PC2, color = Cluster)) + 
    geom_point() + 
    theme(legend.position = "top")
  
  # plot pca v. clusters
  pca.plot = ggpairs(pca, columns = names(pca)[-ncol(pca)], 
                     mapping = aes(color = Cluster), # axisLabels = "internal",
                     legend = grab_legend(legend.plot) + theme_bw(25),
                     upper = list(continuous = wrap("points", alpha = 1/3)),
                     # upper = list(continuous = wrap(ggally_cor, size = 10)),
                     lower = list(continuous = wrap("points", alpha = 1/3)),
                     # diag = list(continuous = wrap("diagAxis", labelSize = 15)),
                     diag = list(continuous = wrap("densityDiag", fill = "black", alpha = 1/3)),
                     title = "Cluster Analysis with PCA") + 
    theme_bw(25) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # pca.plot
  
  # pick x, y, and z variables
  x = pca$PC1
  y = pca$PC2
  z = pca$PC3
  
  # make Cluster into a numeric
  pca[, Cluster := as.numeric(as.character(Cluster))]
  
  # create Cluster.factor and Cluster.number columns
  pca[, Cluster.factor := factor(Cluster, levels = sort(unique(pca$Cluster)))]
  pca[, Cluster.number := as.numeric(Cluster.factor)]
  
  # plot PC v. clusters
  # https://rpubs.com/yoshio/95844
  scatter3D(x = x, y = -y, z = -z, # xlim = c(0, 11), zlim = c(0, 90),
            main = "Cluster Analysis with PCA",
            colvar = pca$Cluster.number, 
            col = jet.col(n = length(unique(pca$Cluster.number)), alpha = 2/3),
            # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
            colkey = list(at = sort(unique(pca$Cluster.number)), labels = levels(pca$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
            # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
            xlab = "\nPC1\n22.4%", ylab = "\n\nPC2\n13.5%", zlab = "\n\nPC3\n7.3%", clab = "Cluster", 
            # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
            theta = 20, phi = 40,
            pch = 16, ticktype = "detailed", type = "p", bty = "b2",
            # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
            cex = 1, cex.main = 2.5, cex.axis = 2, cex.lab = 2.5)
  
  x = pca$PC1
  y = pca$PC2
  z = pca$PC4
  
  scatter3D(x = x, y = -y, z = z, # xlim = c(0, 11), zlim = c(0, 90),
            main = "Cluster Analysis with PCA",
            colvar = pca$Cluster.number, 
            col = jet.col(n = length(unique(pca$Cluster.number)), alpha = 2/3),
            # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
            colkey = list(at = sort(unique(pca$Cluster.number)), labels = levels(pca$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
            # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
            xlab = "\nPC1\n22.4%", ylab = "\n\nPC2\n13.5%", zlab = "\n\nPC4\n5.5%", clab = "Cluster", 
            # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
            theta = 20, phi = 40,
            pch = 16, ticktype = "detailed", type = "p", bty = "b2",
            # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
            cex = 1, cex.main = 2.5, cex.axis = 2, cex.lab = 2.5)
  
}

# run autoencoders on the data to create more features?
build.encoder = TRUE

if(build.encoder)
{
  # convert Cluster into a factor
  DT[, Cluster := factor(Cluster, levels = sort(unique(Cluster)))]
  
  # identify predictors (x) and response (y)
  y = "Cluster"
  x = names(DT[, !"Cluster", with = FALSE])
  
  # get the predictors
  X.h2o = as.h2o(DT[, x, with = FALSE])
  
  # set up hyperparameters of interest for the autoencoder
  auto.hyper.params = list(hidden = list(c(round((2/3) * length(x), 0), round((2/3)^2 * length(x), 0), round((2/3) * length(x), 0))),
                           epochs = c(30),
                           activation = "Tanh",
                           l1 = 1e-5,
                           l2 = 0,
                           rho = 0.95,
                           epsilon = 1e-8,
                           adaptive_rate = TRUE)
  
  # lets use a random grid search and specify a time limit and/or model limit
  minutes = 30
  auto.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
  
  # run a random grid search for a good model
  h2o.rm("auto.random.grid")
  auto.random.grid = h2o.grid(algorithm = "deeplearning",
                              grid_id = "auto.random.grid",
                              autoencoder = TRUE,
                              x = x,
                              training_frame = X.h2o,
                              # score_each_iteration = TRUE,
                              seed = 3,
                              hyper_params = auto.hyper.params,
                              search_criteria = auto.search.criteria)
  
  # rank each model in the random grid
  auto.grid = h2o.getGrid("auto.random.grid", sort_by = "rmse", decreasing = FALSE)
  
  # check out model ranking
  auto.grid
  
  # get the best model from our grid search
  auto.mod = h2o.getModel(auto.grid@model_ids[[1]])
  
  # get features from auto.mod
  features = lapply(1:3, function(i) 
    h2o.deepfeatures(object = auto.mod, data = X.h2o, layer = i))
  
  # combine features into a single table
  features = do.call("h2o.cbind", features)
  
  # add features to DT
  DT = cbind(DT, as.data.table(features))
  
  #########################
  
  # create Cluster.factor and Cluster.number columns
  DT.numeric[, Cluster.factor := factor(Cluster, levels = sort(unique(Cluster)))]
  DT.numeric[, Cluster.number := as.numeric(Cluster.factor)]
  
  # plot some features colored by cluster
  x = DT.numeric$Birth_Cohort
  y = DT.numeric$GNIpc
  z = DT.numeric$DF.L3.C2
  
  # plot features for the clusters
  scatter3D(x = -x / 1e6, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
            main = "Features for Clusters",
            colvar = DT.numeric$Cluster.number, 
            col = ggcolor(n = length(unique(DT.numeric$Cluster.number)), a = 1/3),
            colkey = list(at = sort(unique(DT.numeric$Cluster.number)), labels = levels(DT.numeric$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
            xlab = "\n\nAnnual Births\n(Millions)", ylab = "\n\nGNIpc\n(USD)", zlab = "\n\nDeep Feature\nL3C2", clab = "Cluster", 
            theta = 50, phi = 20,
            pch = 16, ticktype = "detailed", type = "p", bty = "g",
            cex = 1.5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
  
  # plot some features colored by cluster
  x = DT.numeric$Reservation_Price
  y = DT.numeric$Country_Risk
  z = DT.numeric$DF.L1.C6
  
  # plot features for the clusters
  scatter3D(x = -x, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
            main = "Features for Clusters",
            colvar = DT.numeric$Cluster.number, 
            col = ggcolor(n = length(unique(DT.numeric$Cluster.number)), a = 1/3),
            colkey = list(at = sort(unique(DT.numeric$Cluster.number)), labels = levels(DT.numeric$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
            xlab = "\n\n\nReservation\nPrice\n(USD)", ylab = "\nRisk", zlab = "\n\nDeep Feature\nL1C6", clab = "Cluster", 
            theta = 50, phi = 20,
            pch = 16, ticktype = "detailed", type = "p", bty = "g",
            cex = 1.5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
  
  
}

# ---- feature selection ----

{

# identify predictors (x) and response (y)
y = "Cluster"
x = names(DT[, !"Cluster", with = FALSE])

# make an h2o version of DT
DT.h2o = as.h2o(DT)

# compute the max class weight
max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))

# set up hyperparameters of interest
rf.hyper.params = list(ntrees = 250,
                       # min_rows = c(1, 11, 25),
                       # max_depth = c(20, 40, 60),
                       min_rows = 11,
                       max_depth = 40,
                       stopping_metric = "mean_per_class_error")

# lets use a random grid search and specify a time limit and/or model limit
minutes = 30
rf.search.criteria = list(strategy = "RandomDiscrete", 
                          max_runtime_secs = minutes * 60, 
                          # max_models = 100, 
                          seed = 42)

# lets run a grid search for a good model, without drop out ratios
h2o.rm("rf.random.grid")
rf.random.grid = h2o.grid(algorithm = "randomForest",
                          grid_id = "rf.random.grid",
                          y = y,
                          x = x,
                          training_frame = DT.h2o,
                          # stopping_rounds = 3,
                          # histogram_type = "RoundRobin",
                          # nfolds = 3,
                          # fold_assignment = "Stratified",
                          seed = 3,
                          # score_each_iteration = TRUE,
                          balance_classes = TRUE,
                          max_after_balance_size = max.class.weight,
                          hyper_params = rf.hyper.params,
                          search_criteria = rf.search.criteria)

# free up RAM
gc()

# rank each model in the random grids
rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)

# pick the top model from all grid searches
imp.rf = h2o.getModel(rf.grid@model_ids[[1]])

# extract variable importance
imp = data.table(imp.rf@model$variable_importances)

# make variable into a factor
imp[, variable := factor(variable, levels = unique(variable))]

# pick a cut off value
# this value should show where importance drops the most (look for the center of the "knee")
cutoff = 0.25

# plot a barplot of variable importance
imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
  ggtitle("GINI Importance\nRandom Forests") +
  labs(x = "Variable", y = "Scaled Importance") +
  scale_y_continuous(labels = percent) +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_dark(25) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank())

# plot a density plot of variable importance
imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
  geom_density(fill = "cornflowerblue", alpha = 2/3) +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
  ggtitle("GINI Importance\nRandom Forests") +
  labs(x = "Scaled Importance", y = "Density") +
  scale_x_continuous(labels = percent) +
  coord_flip() + 
  theme_bw(25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# show the importance plots
grid.arrange(imp.plot1, imp.plot2, nrow = 1)

# check out imp
imp

# reverse the order of the levels of variable
imp[, variable := factor(variable, levels = rev(levels(variable)))]

# plot imp
cutoff = 0.25
imp.plot3 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
  ggtitle("Random Forests: GINI Importance") +
  labs(x = "Variable", y = "Scaled Importance") +
  scale_y_continuous(labels = percent) +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_color_gradient(low = "yellow", high = "red") +
  coord_flip() + 
  theme_dark(20) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

imp.plot3

# find which indicators meet the cutoff value
keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))

# find the top m indicators
m = 15
keep.indicators = as.character(imp$variable)[1:m]

# update DT to only contain indicators of interest
DT = DT[, c(y, keep.indicators), with = FALSE]

# determine which features to categorize
exclude.columns = c(y)
features = data.table(DT[, !exclude.columns, with = FALSE])
feature.names = names(features)

# go through each column of features and categorize it based on 10 quantiles
features = lapply(1:ncol(features), function(i)
{
  # get column i from features
  v = unname(unlist(features[, i, with = FALSE]))
  
  # categorize v
  v = cut(x = v, 
          breaks = unique(c(min(v) - 1e-6, 
                            as.numeric(quantile(v, probs = seq(0.1, 0.9, 0.1))),
                            max(v) + 1e-6)), 
          ordered_result = TRUE)
  
  # make v into a column again
  v = data.table(v)
  
  # give v its name back
  setnames(v, feature.names[i])
  
  return(v)
})

# combine features into a single table
features = do.call("cbind", features)

# convert features into binary values
features = data.table(model.matrix(~ ., data = features, 
                                   contrasts.arg = lapply(features, contrasts, contrasts = FALSE))[,-1])

# keep a numeric version of DT
DT.numeric = data.table(DT)

# plot some features colored by cluster
DT.plot = ggpairs(DT.numeric, columns = names(DT.numeric)[2:6], 
                  mapping = aes(color = Cluster),
                  upper = list(continuous = wrap("density", alpha = 1)),
                  lower = list(continuous = wrap("points", alpha = 1/3, size = 1.5)),
                  diag = list(continuous = wrap("densityDiag", alpha = 1)),
                  title = "Features for Clusters") + 
  theme_bw(20) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# DT.plot

# create Cluster.factor and Cluster.number columns
DT.numeric[, Cluster.factor := factor(Cluster, levels = sort(unique(Cluster)))]
DT.numeric[, Cluster.number := as.numeric(Cluster.factor)]

# plot some features colored by cluster
x = DT.numeric$Birth_Cohort
y = DT.numeric$GNIpc
z = DT.numeric$DF.L3.C2

# plot features for the clusters
scatter3D(x = -x / 1e6, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
          main = "Features for Clusters",
          colvar = DT.numeric$Cluster.number, 
          col = ggcolor(n = length(unique(DT.numeric$Cluster.number)), a = 1/3),
          colkey = list(at = sort(unique(DT.numeric$Cluster.number)), labels = levels(DT.numeric$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
          xlab = "\n\nAnnual Births\n(Millions)", ylab = "\n\nGNIpc\n(USD)", zlab = "\n\nDeep Feature\nL3C2", clab = "Cluster", 
          theta = 50, phi = 20,
          pch = 16, ticktype = "detailed", type = "p", bty = "g",
          cex = 1.5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)

# plot some features colored by cluster
x = DT.numeric$Reservation_Price
y = DT.numeric$Country_Risk
z = DT.numeric$DF.L1.C6

# plot features for the clusters
scatter3D(x = -x, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
          main = "Features for Clusters",
          colvar = DT.numeric$Cluster.number, 
          col = ggcolor(n = length(unique(DT.numeric$Cluster.number)), a = 1/3),
          colkey = list(at = sort(unique(DT.numeric$Cluster.number)), labels = levels(DT.numeric$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
          xlab = "\n\n\nReservation\nPrice\n(USD)", ylab = "\nRisk", zlab = "\n\nDeep Feature\nL1C6", clab = "Cluster", 
          theta = 50, phi = 20,
          pch = 16, ticktype = "detailed", type = "p", bty = "g",
          cex = 1.5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)

# update DT to only have binary features
DT = cbind(DT[, exclude.columns, with = FALSE], features)

# identify predictors (x) and response (y)
y = "Cluster"
x = names(DT[, !"Cluster", with = FALSE])

# make an h2o version of DT
DT.h2o = as.h2o(DT)

# compute the max class weight
max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))

# set up hyperparameters of interest
rf.hyper.params = list(ntrees = 250,
                       # min_rows = c(1, 11, 25),
                       # max_depth = c(20, 40, 60),
                       min_rows = 11,
                       max_depth = 40,
                       stopping_metric = "mean_per_class_error")

# lets use a random grid search and specify a time limit and/or model limit
minutes = 30
rf.search.criteria = list(strategy = "RandomDiscrete", 
                          max_runtime_secs = minutes * 60, 
                          # max_models = 100, 
                          seed = 42)

# lets run a grid search for a good model, without drop out ratios
h2o.rm("rf.random.grid")
rf.random.grid = h2o.grid(algorithm = "randomForest",
                          grid_id = "rf.random.grid",
                          y = y,
                          x = x,
                          training_frame = DT.h2o,
                          # stopping_rounds = 3,
                          # nfolds = 3,
                          # fold_assignment = "Stratified",
                          seed = 3,
                          # score_each_iteration = TRUE,
                          
                          # control how many bins are created for splitting on a feature
                          nbins_cats = 3, 
                          nbins = 3,
                          nbins_top_level = 3,
                          
                          balance_classes = TRUE,
                          max_after_balance_size = max.class.weight,
                          hyper_params = rf.hyper.params,
                          search_criteria = rf.search.criteria)

# free up RAM
gc()

# rank each model in the random grids
rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)

# pick the top model from all grid searches
imp.rf = h2o.getModel(rf.grid@model_ids[[1]])

# extract variable importance
imp = data.table(imp.rf@model$variable_importances)

# make variable into a factor
imp[, variable := factor(variable, levels = unique(variable))]

# pick a cut off value
# this value should show where importance drops the most (look for the center of the "knee")
cutoff = 0.23

# plot a barplot of variable importance
imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
  ggtitle("GINI Importance\nRandom Forests") +
  labs(x = "Variable", y = "Scaled Importance") +
  scale_y_continuous(labels = percent) +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_dark(25) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank())

# plot a density plot of variable importance
imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
  geom_density(fill = "cornflowerblue", alpha = 2/3) +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
  ggtitle("GINI Importance\nRandom Forests") +
  labs(x = "Scaled Importance", y = "Density") +
  scale_x_continuous(labels = percent) +
  coord_flip() + 
  theme_bw(25) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# show the importance plots
grid.arrange(imp.plot1, imp.plot2, nrow = 1)

# check out imp
imp

# reverse the order of the levels of variable
imp[, variable := factor(variable, levels = rev(levels(variable)))]

# plot imp
cutoff = 0.17
imp.plot3 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
  ggtitle("Random Forests: GINI Importance") +
  labs(x = "Variable", y = "Scaled Importance") +
  scale_y_continuous(labels = percent) +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_color_gradient(low = "yellow", high = "red") +
  coord_flip() + 
  theme_dark(20) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

imp.plot3

# find which indicators meet the cutoff value
keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))

# find the top m indicators
m = 50
keep.indicators = as.character(imp$variable)[1:m]

# update DT to only contain indicators of interest
DT = DT[, c(y, keep.indicators), with = FALSE]

}

# ---- machine learning ----

{

# identify predictors (x) and response (y)
y = "Cluster"
x = names(DT[, !"Cluster", with = FALSE])

# build the fold assignment
set.seed(42)
k.folds = 5
folds = createFolds(y = unname(unlist(DT[, y, with = FALSE])), k = k.folds)

# split up DT into train, valid, and test
train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
train = data.table(DT[train.rows])

valid.rows = unname(unlist(folds[[k.folds - 1]]))
valid = data.table(DT[valid.rows])

test.rows = unname(unlist(folds[[k.folds]]))
test = data.table(DT[test.rows])

# split up YX.h2o into train, valid, and test
train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])

# write out train, valid, test, and DT.numeric
fwrite(train, "train.csv")
fwrite(valid, "valid.csv")
fwrite(test, "test.csv")
fwrite(DT.numeric, "numeric.csv")

# ---- GLM ----


# ---- RF ----


# ---- GB ----


# ---- NNET ----


# ---- STACK ----


}





}



