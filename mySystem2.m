function F = mySystem2(vars,params)

    x1 = tickMap(vars,params);
    x2 = tickMap(x1,params);
  
    F = x2 - vars;

end