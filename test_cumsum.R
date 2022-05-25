s = c(1:3)

test_func = function(s){
  x = s
  cumsum(x)
}
test_func(s) #returns cumsum(x)

test_func2 = function(s){
  x = s
}
test_func2(s) #returns nothing