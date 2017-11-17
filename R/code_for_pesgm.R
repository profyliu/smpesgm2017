
#' Compute the past 24-hour and 6-hour mean and sd for each point in range.
#'
#' @param df the input data frame.
#' @param range the datetime range for calculation.
#'
calc.T2 <- function(df, range){
  df$H24.sd = NA
  df$H6.sd = NA
  beg_row = which(grepl(range[1], df$datetime))
  end_row = which(grepl(range[2], df$datetime))
  a = vector(mode="numeric",length=48)
  b = vector(mode="numeric",length=12)
  for(t in beg_row:end_row){
    for(h in 1:48){
      a[h] = df[t-h, 'Total']
      if(h <= 12){
        b[h] = df[t-h, 'Total']
      }
    }
    df[t,'H24.sd'] = stats::sd(a)
    df[t,'H6.sd'] = stats::sd(b)
  }
  return(df)
}

#' Compute the D-day same-time avearge for the T-hour sd for each point in range
#' and do single-point anomaly detection.
#'
#' @param df the input data frame.
#' @param range the datetime range for calculation.
#' @param alpha confidence level in the Grubb's test.
#' @param w weight vector for the three features.
#'
calc.DOT <- function(df, range, alpha=0.9, w=c(0.6,0.3,0.1)){
  D = 30
  c = 48
  df$p = NA
  df$y = NA
  df$thresh = (D-1)/sqrt(D)*sqrt(stats::qt(alpha/(2*D),D-2)^2/(D-2+stats::qt(alpha/(2*D),D-2)^2))
  df$anomaly = 0
  beg_row = which(grepl(range[1], df$datetime))
  end_row = which(grepl(range[2], df$datetime))
  m = data.frame(matrix(ncol=3,nrow=D))
  for(t in beg_row:end_row){
    for(d in 1:D){
      m[d,] = df[t-c*d, c('H24.sd','H6.sd','Total')]
    }
    mean.m = colSums(m)/D
    cov.m = stats::cov(m*w)
    val = df[t,c('H24.sd','H6.sd','Total')]
    p = mvtnorm::dmvnorm(val, mean.m, cov.m, log=T)
    dm = as.matrix((val-mean.m)*w, nrow=1)
    y = sqrt(dm %*% solve(cov.m) %*% t(dm))
    df[t,'p'] = p
    df[t,'y'] = y
    df[t,'anomaly'] = ifelse(df[t,'y'] > df[t,'thresh'], 1, 0)
  }
  return(df)
}

#' Count number of anomalies within a time range.
#'
#' @param df the input data frame.
#' @param range the datetime range, a vector of two elements.
#' @param K the minimum number of contiguous single-point anomalies needed to conclude an anomalous sequence.
#'
count_anomaly <- function(df, range, K=6){
  beg_row = which(grepl(range[1], df$datetime))
  end_row = which(grepl(range[2], df$datetime))
  n_anomaly = 0
  count = 0
  for(t in beg_row:end_row){
    if(df[t,'anomaly'] == 1){
      count = count + 1
    } else{
      if(count >= K){
        n_anomaly = n_anomaly + 1
      }
      count = 0
    }
  }
  return(n_anomaly)
}

#' Find the first anomaly after a given time point
#' and return the number of intervals past until detection.
#'
#' @param df the input data frame.
#' @param start.time the starting datetime.
#' @param K the minimum number of contiguous single-point anomalies needed to conclude an anomalous sequence.
#' @param timeout maximum acceptable period (number of intervals) for a useful detection.
#'
find_anomaly <- function(df, start.time, K=6, timeout=48){
  beg_row = which(grepl(start.time, df$datetime))
  count = 0
  detection.time = beg_row + timeout
  for (t in beg_row:(beg_row + timeout)){
    if(df[t,'anomaly'] == 1 & t < (beg_row + timeout)){
      if(df[t-1,'anomaly'] == 0){detection.time = t}
      count = count + 1
    } else {
      if(count >=K){
        #print(paste("detection.time =", detection.time))
        return(detection.time + K - beg_row)
      }
      count = 0
      #print(paste("reset count at t=",t))
    }
  }
  return(timeout)
}


#' Sample N occurrence time within range.
#'
#' @param df the input data frame.
#' @param range the datetime range within which samples are drawn. It is a
#'   vector of two elements, which mark the start and end datetime of the range.
#' @param N number of samples to draw.
#' @param seed random seed, an integer.
#'
sample.occurrence <- function(df, range, N=10, seed=2017){
  set.seed(seed)
  beg_row = which(grepl(range[1], df$datetime))
  end_row = which(grepl(range[2], df$datetime))
  row_nums = sample(beg_row:end_row, N)
  ocvec = df[row_nums, "datetime"]
  return(ocvec)
}

#' Given DOT begin datetime, calculate the T2 begin datetime.
#'
#' @param df the input data frame.
#' @param cur.datetime the reference datetime.
#' @param c number of time intervals per day.
#' @param D number of days to look back.
calc.T2.datetime <- function(df, cur.datetime, c=48, D=30){
  cur_index = which(grepl(cur.datetime, df$datetime))
  beg_index = cur_index - (D)*c
  beg_time = df[beg_index, 'datetime']
  return(beg_time)
}

#' Given DOT begin datetime, calculate the overall begin datetime.
#'
#' @param df the input data frame.
#' @param cur.datetime DOT begin datetime.
#' @param c number of time intervals per day.
#' @param D number of days to look back.
#'
calc.begin.datetime <- function(df, cur.datetime, c=48, D=30){
  cur_index = which(grepl(cur.datetime, df$datetime))
  beg_index = cur_index - (D+1)*c
  beg_time = df[beg_index, 'datetime']
  return(beg_time)
}

#' Given an anomaly occurrence time, calculate the detection period end time.
#'
#' @param df the input data frame.
#' @param oc.datetime anomaly occurrence time.
#' @param timeout maximum acceptable period (number of intervals) for a useful detection.
#'
calc.detect.end <- function(df, oc.datetime, timeout=48){
  oc_index = which(grepl(oc.datetime, df$datetime))
  end_index = oc_index + timeout
  end_time = df[end_index, 'datetime']
  return(end_time)
}

#' Make an anomaly data case given occurrence time.
#'
#' @param df.in the input data frame.
#' @param occurrence.time time at which the anomaly is to occur.
#'
makeAbCase <- function(df.in, occurrence.time){
  begin.datetime = calc.begin.datetime(df.in, occurrence.time)
  end.datetime = calc.detect.end(df.in, occurrence.time)
  beg_row = which(grepl(begin.datetime, df.in$datetime))
  end_row = which(grepl(end.datetime, df.in$datetime))
  ab = df.in[beg_row:end_row,c("datetime","date","FridgeRange","KitchenLights","BedroomLights","ElectricRange")]
  ab[which(ab$datetime >= occurrence.time), c("KitchenLights","BedroomLights","ElectricRange")] = ab[which(ab$datetime == occurrence.time),c("KitchenLights","BedroomLights","ElectricRange")]
  ab$Total = rowSums(ab[colnames(ab)[3:6]])
  return(ab)
}

#' Run an experiment.
#'
#' Generate N random cases in range and report detection results.
#' @param df the input data frame.
#' @param range time range in which simulated anomalies can occur.
#' @param N number of detection experiments to run.
#'
detect.simu <- function(df, range, N){
  ocvec = sample.occurrence(df,range,N)
  result = data.frame(occurred = rep("when",N), detect.time = vector(mode="integer",length=N), detected = vector(length=N), stringsAsFactors=FALSE)
  for(i in 1:N){
    oc = strptime(ocvec[i], format="%Y-%m-%d %H:%M:%S")
    print(paste("Evaluating oc = ", oc, "i =", i))
    ab1 = makeAbCase(df,oc)
    end.datetime = calc.detect.end(ab1,oc)
    range.T2 = c(calc.T2.datetime(ab1,oc), end.datetime)
    range.DOT = c(oc, end.datetime)
    ab1 = calc.T2(ab1, range.T2)
    ab1 = calc.DOT(ab1, range.DOT)
    dt = find_anomaly(ab1,oc)
    result[i, 'occurred'] = paste(oc)
    result[i, 'detect.time'] = dt
    result[i, 'detected'] = ifelse(dt < 48, 1, 0)
  }
  return(result)
}


