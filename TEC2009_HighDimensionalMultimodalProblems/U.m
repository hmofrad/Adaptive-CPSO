function out = U(x)
    a = 10;
    k = 100;
    m = 4;
    if (x>a)
        out = k*(x - a)^m;
    elseif (x>=-a && x<=a)
        out = 0;
    elseif  (x<-a)
        out = k*(-x - a)^m;
    end