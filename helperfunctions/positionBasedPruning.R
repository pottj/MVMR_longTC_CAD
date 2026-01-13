positionBasedPruning = function(data,distance,col_chr,col_pos,col_prio,smallerIsBetter){
  # data = copy(set1_fem)
  # distance = 1000000
  # col_chr = "CHR"
  # col_pos = "POS"
  # col_prio = "tval_abs"
  # smallerIsBetter = F
  
  # check column names
  stopifnot(col_chr %in% names(data))
  stopifnot(col_pos %in% names(data))
  stopifnot(col_prio %in% names(data))
  
  # rename columns if necessary
  if(col_chr != "CHR") data[,CHR := get(col_chr)]
  if(col_pos != "POS") data[,POS := get(col_pos)]
  if(col_prio != "PRIO") data[,PRIO := get(col_prio)]
  
  # order data set and get all possible chromosomes
  setorder(data,CHR,POS)
  myCHRs = unique(data$CHR)
  
  # run loop per chromosome for filtering
  result.22 = foreach(s2 = myCHRs) %do% {
    # s2 = myCHRs[1]
    subdata2 = copy(data)
    subdata2 = subdata2[CHR == s2, ]
    setkey(subdata2, POS)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, NR_SNPs := 0]
    }else{
      subdata2[, keep := NA]
      subdata2[, NR_SNPs := as.numeric(NA)]
      
      smallestDist = getSmallestDist(subdata2[, POS])
      while(smallestDist < distance) {
        #minP = min(subdata2[is.na(keep), p])
        if(smallerIsBetter == T){
          myPrio = min(subdata2[is.na(keep), PRIO])
        }else{
          myPrio = max(subdata2[is.na(keep), PRIO])
        }
        myPOS = subdata2[myPrio == PRIO & is.na(keep), POS]
        if(length(myPOS)>1){
          myPOS = myPOS[1]
        }
        subdata2[POS == myPOS, keep := T]
        
        #filter for SNPs that can stay within the set (outside the +- 1MB range or keep==T)
        myFilt = (subdata2[, POS] < (myPOS - distance)) | 
          (subdata2[, POS] > (myPOS + distance)) | 
          subdata2[, keep] 
        myFilt[is.na(myFilt)] = FALSE
        subdata2 = subdata2[myFilt == TRUE, ]
        
        subdata2[POS == myPOS, NR_SNPs := sum(myFilt==F)]
        smallestDist = getSmallestDist(subdata2[, POS])
      }
      
      #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
      subdata2[is.na(keep), NR_SNPs := 0]
      subdata2[is.na(keep), keep := TRUE]
    }
    
    subdata2
  }
  data2 = rbindlist(result.22)
  
  return(data2)
  
}