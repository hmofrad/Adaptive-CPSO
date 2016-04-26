% fitness function
function f = fit_func(c,x)
    % c as the fitness function number
    func = c{1}; % benchmark number
    [nop D]=size(x);
    switch (func)
        case {'Schewfel',8}
            x=x-(-12569.5);
            f = sum(-x.*sin(sqrt(abs(x))),2);
            f = abs(f);
        case {'Rastrigin',9}
            f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
        case {'Ackley',10}
            f=sum(x.^2,2);
            f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
        case {'Griewank',11}
            f=1;
            for i=1:D
                f=f.*cos((x(:,i)-100)./sqrt(i));
            end
            f=sum((x-100).^2,2)./4000-f+1;
        case {'Penalized1',12}
            y = 1 + ((x+1)./4);     
            f = (pi/D)*(10.*sin(pi.*y(:,1)).^2 + ...
            sum(((y(:,1:D-1)-1).^2).*(1+10*sin(pi.*y(:,2:D)).^2),2) + ...
            (y(:,D)-1).^2);
            for i=1:nop
                for j=1:D
                    u(j) = U(x(i,j));
                end
                f(i) = f(i) + sum(u);
            end
       case {'Penalized2',13}     
            f = 0.1*(10.*sin(3*pi.*x(:,1)).^2 + ...
            sum(((x(:,1:D-1)-1).^2).*(1+sin(3*pi.*x(:,2:D)).^2),2) + ...
            (x(:,D)-1).^2);
            for i=1:nop
                for j=1:D
                    u(j) = U(x(i,j));
                end
                f(i) = f(i) + sum(u);
            end
    end
    function u = U(x)
        a = 10; k = 100; m = 4;
        if (x>a)
            u = k*(x - a)^m;
        elseif (x>=-a && x<=a)
            u = 0;
        elseif  (x<-a)
            u = k*(-x - a)^m;
        end
    end
end