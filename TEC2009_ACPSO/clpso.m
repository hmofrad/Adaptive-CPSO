function [gbest_hist fe_hist] = clpso(f,nop,dim)
%     clear,clc
%     f= 11;nop = 50;dim=300;
    endgen = 10^5;    % maximum Generation
    max_fe =  7*10^6; % 3750000 rcurrent fitness evaluation

    if f == 8;  Ub = 500;    end % Schwefel
    if f == 9;  Ub = 5.12;   end % Rastrigin
    if f == 10; Ub = 32.768; end % Ackley
    if f == 11; Ub = 600;    end % Griewank
    if f == 12; Ub = 50;     end % Penalized1
    if f == 13; Ub = 50;     end % Penalized2
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    x(:,end+1)=0; % calculate each particle fittness
    c = 1.49445;
    w_max = 0.9;  w_min = 0.4;
    for i=1:nop
        x(i,end) = fit_func({f},x(i,1:end-1));
    end

    pbest = x; %initialize Best Particle Position
    [bst ind] = min(x(:,end));
    gbest = x(ind,:); % initialize global best position
    gbest_hist = [];
    % initialise pbestf
    pbestf = []; pbestftmp = [];
    for j = 1:nop
        for i = 1:dim
            Pind = randi([1 nop]);  % Particle Index
            Dind = randi([1 dim]);  % Dimension Index
            pbestftmp = [pbestftmp  pbest(Pind,Dind)];
        end
        pbestf = [pbestf; pbestftmp];
        pbestftmp = [];
    end
    fe_hist = []; gbest_hist = [];fe=0;m=7;
    while i<=endgen && fe<=max_fe
        w = w_max-(w_max-w_min)*i/endgen;
        for j=1:nop
            w = w_max*((w_max-w_min)*i/endgen);
            if mod(i, m) == 0
                pbestf(j,:) = exemplar(x,pbest,j);
            end
            if (x(j,end) < pbest(j,end))
                pbest(j,:) = x(j,:);
            end
            if (pbest(j,end) < gbest(end))
                gbest = pbest(j,:);
            end
        end
        spd = w.*spd+c.*rand(nop,dim).*(pbestf-x(:,1:end-1));
        spd=(spd>Vmax).*Vmax+(spd<=Vmax).*spd; 
        spd=(spd<(Vmin)).*(Vmin)+(spd>=(Vmin)).*spd;
        x(:,1:end-1) = x(:,1:end-1)+spd;
        x(:,end) = fit_func({f},x(:,1:end-1));
        fe = fe+nop;
        gbest_hist = [gbest_hist gbest(end)]; fe_hist = [fe_hist fe];
        if mod(i,1000)==0,fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end)),end
        if ((i>1300 && gbest_hist(length(gbest_hist)-500) == gbest(end)) || gbest(end) ==0),break,end
        i = i + 1;
    end
end