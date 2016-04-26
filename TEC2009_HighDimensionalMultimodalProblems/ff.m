function fit = ff(c,x)
% c as the fitness function number
    fit_func = c{1};
    fit = 0; xtmp = 0; xtmp1 = 0;i=0;
    switch (fit_func)
        case {'sphere',1}
            fit = sum(x.^2);
        case {'Schwefels1.2',2}
            for i=1:length(x)
                xtmp = sum(x(1:i))^2;
                fit  = fit+xtmp;
            end
        case {'Schwefels2.21',3}
            fit = max(abs(x));
        case {'Rosenbrock',4}
            for i=1:length(x)-1
                xtmp = 100*((x(i+1)-x(i)^2)^2) + (x(i) - 1)^2;
                fit = fit+xtmp;
            end
        case {'Ackley',5}
        	fit = -20*exp(-0.2*sqrt(sum(x.^2)/length(x))) - exp(sum(cos(2*pi*x))/length(x)) + 20 + exp(1);
        case {'Griewank',6}  
            sum1 = 0;
            prod1 = 1;
            for i = 1:1:length(x)
                sum1 = sum1 + (x(i)^2/4000);
                prod1 = prod1*cos(x(i)/sqrt(i));
            end
            fit = sum1 - prod1 + 1;
        case {'Penalized1',7}
            y = 1 + ((x+1)/4);
            for i=1:length(x)-1
                xtmp = xtmp + (y(i)-1)^2 * (1+10*sin(pi*y(i+1))^2);
                xtmp1 = xtmp1 + U(x(i));
            end
            xtmp1 = xtmp1 + U(x(end));
            fit = (pi/length(x))*(10*sin(pi*y(1))^2 + xtmp + (y(end)-1)^2) + xtmp1;   
       case {'Penalized2',8}
            y = 1 + ((x+1)/4);
            for i=1:length(x)-1
                xtmp = xtmp + (y(i)-1)^2 * (1+10*sin(pi*y(i+1))^2);
                xtmp1 = xtmp1 + U(x(i));
            end
            xtmp1 = xtmp1 + U(x(end));
            fit = 0.1*(sin(pi*y(1))^2 + xtmp + (y(end)-1)^2) + xtmp1;
    end
end