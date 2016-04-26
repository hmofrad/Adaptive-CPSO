% function [x_std  gbest_fit gbest_hist] = pso(f,bnd,dim,nop,endgen)
clear,clc
    f= 8;nop = 50;dim=30;
    endgen = 10^4;    % maximum Generation
    max_fe = 10^6; % current fitness evaluation
    if f == 8;  Ub = 500;     end % Schwefel
    if f == 9;  Ub = 32.768;  end % Rastrigin
    if f == 10; Ub = 600;     end % Ackley
    if f == 11; Ub = 0.5;     end % Griewank
    if f == 12; Ub = 5.12;    end % Penalized1
    if f == 13; Ub = 5.12;    end % Penalized1
    w = 0.72;
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity    nop = 10; dim = 9;
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    x(:,end+1)=0;
    c1 = 1.49445;    c2 = 1.49445;
    x(:,end) = fit_func({f},x(:,1:end-1));
    pbest = x; %initialize Best Particle Position
    [mn ind] = min(x(:,end));
    gbest = x(ind,:);
    gbest_hist = [];
    fe = 0;
    for i=1:endgen
        for j=1:nop
            if (x(j,end)<pbest(j,end))
                pbest(j,:) = x(j,:);
            end
            if (pbest(j,end)<gbest(end))
                gbest = pbest(j,:);
            end 
            spd(j,:) = w.* spd(j,:)+c1.*rand(1,dim).*(pbest(j,1:end-1)-x(j,1:end-1))...
                                   +c2.*rand(1,dim).*(gbest(1:end-1) - x(j,1:end-1));
%             ind = find(abs(spd>Vmax));
%             spd(ind) = Vmax.*sign(spd(ind));
            spd=(spd>Vmax).*Vmax+(spd<=Vmax).*spd; 
            spd=(spd<(Vmin)).*(Vmin)+(spd>=(Vmin)).*spd;

            x(j,1:end-1) = x(j,1:end-1)+spd(j,:);
        end
        fe = fe+nop;
%         x(:,end) = abs(benchmark_func(x(:,1:end-1),f)-f_bias);
        x(:,end) = fit_func({f},x(:,1:end-1));
%         if mod(i,200) == 0,fprintf('iteration=%u,gbest=%e\n',i,gbest(end));end
        gbest_hist = [gbest_hist gbest(end)];
        fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end))
    end
% end