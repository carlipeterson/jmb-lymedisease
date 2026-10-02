function F = mySystem2Infection(vars,params)

    x1 = tickMapInfection(vars,params);
    x2 = tickMapInfection(x1,params);
  
    F = x2 - vars;

end