library(doParallel)
library(foreach)


parallelFunc = function(i){
  print("heyBabe")
  if(i ==1){
    print("hey")
    a = 1:10000000
    print("finish a")
  }
  if(i ==3){
    print("hey1")
    b = 1:1000000000
    print("finish b")
  }
  if(i ==2){
    print("hey2")
    c= 1:100000000
    print("finish c")
  }
  return(0)
}

foreach(i=1:3) %dopar% parallelFunc(i)
registerDoParallel(cores = 4)
stime = system.time( foreach(i=1:3, .combine = cbind) %dopar% parallelFunc(i))
stime
rm(stime)
print("good")
