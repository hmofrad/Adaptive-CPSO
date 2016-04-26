% fitness function
function f = fit_func(c,x)
    global orthm1 orthm2 swarm1 swarm2 K1 K2 K1len K2len s
    swarmNum = s;
    % c as the fitness function number
    func = c{1}; % benchmark number
    [nop D]=size(x);
    greal=[0 1 0 0 0 0 4.209687462275036e+002 0 0 0 0 0 0 0];
    x=x-greal(func);
    if K1 == 0
        for i=1:swarmNum
            x(:,swarm2(i,:))=x(:,swarm2(i,:))*orthm2(:,swarm2(i,:));
        end
    else
        for i=1:K1
             x(:,swarm1(i,:))=x(:,swarm1(i,:))*orthm1(:,swarm1(i,:));
        end
        if K2len ~= 1
            for i=1:K2
                x(:,swarm2(i,:)-K1*K1len)=x(:,swarm2(i,:)-K1*K1len)*orthm2(:,swarm2(i,:)-K1*K1len);
            end
        end
    end
    
    x=x+greal(func);

    switch (func)
        case {'Sphere',1}
            f=sum(x.^2,2);   
            
        case {'Rosenbrock',2}
            f=sum(100.*(x(:,1:D-1).^2-x(:,2:D)).^2+(x(:,1:D-1)-1).^2,2);
        case {'Ackley',3,9}
            f=sum(x.^2,2);
            f=20-20.*exp(-0.2.*sqrt(f./D))-exp(sum(cos(2.*pi.*x),2)./D)+exp(1);
        case {'Griewank',4,10}
            f=1;
            for i=1:D
                f=f.*cos(x(:,i)./sqrt(i));
            end
            f=sum(x.^2,2)./4000-f+1;
        case {'Weierstrass',5,11}
            x=x+0.5;
            a = 0.5;
            b = 3;
            kmax = 20;
            c1(1:kmax+1) = a.^(0:kmax);
            c2(1:kmax+1) = 2*pi*b.^(0:kmax);
            f=0;
            c=-w(0.5,c1,c2);
            for i=1:D
            f=f+w(x(:,i)',c1,c2);
            end
            f=f+c*D;
        case {'rastrigin', 6,12}
            f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
        case {'Rastrigin_noncont',7,13}
            x=(abs(x)<0.5).*x+(abs(x)>=0.5).*(round(x.*2)./2);
            f=sum(x.^2-10.*cos(2.*pi.*x)+10,2);
        case {'Schewfel',8,14}
            f=0;
            for i=1:D
                f=f-(abs(x(:,i))<=500).*(x(:,i).*sin(sqrt(abs(x(:,i)))))+(abs(x(:,i))>500).*0.001.*(500-abs(x(:,i))).^2;
            end
            f=4.189828872724338e+002*D+f;
    end
    
    function y = w(x,c1,c2)
        y = zeros(length(x),1);
        for k = 1:length(x)
            y(k) = sum(c1 .* cos(c2.*x(:,k)));
        end
    end
end