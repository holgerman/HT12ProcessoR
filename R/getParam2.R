getParam2 <- function(x, myparam = param) {
    
    res = myparam[myparam$variable == x, "value"]
    if (res == "auto" & is.na(res) == F) 
        res = myparam[myparam$variable == x, "givenvalue"]
    res
}
