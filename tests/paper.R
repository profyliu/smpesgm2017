library(smpesgm2017)

# Load and preprocessing raw data
data(ew) # Load the data
ew$Total = rowSums(ew[colnames(ew)[2:5]])
ew[,2:6] = ew[,2:6]*1000 # change from kw to watt
ew$datetime = strptime(ew$Timestamp, format="%m/%d/%Y %H:%M")
ew$datetime = ew$datetime + 1 # Add one second here to bypass the strptime issue caused by 00:00:00 time instances
ew$date = as.POSIXct(strftime(ew$datetime, format="%Y-%m-%d"),format="%Y-%m-%d")


# Run a detection pass on the normal data
beg.datetime = "2014-02-02 00:00:01"
end.datetime = "2014-11-29 23:30:01"
range.T2 = c("2014-01-02 00:00:01", end.datetime)
range.DOT = c(beg.datetime, end.datetime)
ew = calc.T2(ew, range.T2)
ew = calc.DOT(ew, range.DOT)
# False positive rate on single time point
mean(ew[which(ew$date >= '2014-02-02' & ew$date <= '2014-11-29'),'anomaly'])
# False positive rate on k-sequence
range.count = c(beg.datetime, end.datetime)
count_anomaly(ew, range.count, K = 6)
# Save the results to file
#write.csv(file="normal_case_full.csv", ew[which(ew$date >= '2014-02-02' & ew$date <= '2014-11-30'),])


# Run simulated detection experiments
simu.range = c("2014-02-02 00:00:01", "2014-11-29 23:30:01")
start.time = Sys.time()
res = detect.simu(ew, simu.range, N=10) # N = 1000 used in the paper
print(paste("Run time is",round(Sys.time() - start.time,digits=1),"minutes."))
# Save the results to file
#write.csv(file="experimentN1000.csv", res)


