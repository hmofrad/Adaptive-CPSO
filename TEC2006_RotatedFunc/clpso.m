function [gbest_hist fe_hist] = clpso(f,nop,dim,swarmNum)
%     clear,clc
%     f= 14;nop = 10;dim=10;swarmNum=6;
    endgen = 10^4;    % maximum Generation
%     if dim == 10, max_fe = 2*10^4; end
%     if dim == 30, max_fe = 8*10^5; end
    if dim == 10, max_fe = 3*10^4; end
    if dim == 30, max_fe = 3*10^5; end


    global  orthm1 orthm2 swarm1 swarm2 K1 K2 K1len K2len s
    if f == 9;  Ub = 32.768; end % Rotated Ackley
    if f == 10; Ub = 600;    end % Rotated Griewank
    if f == 11; Ub = 0.5;    end % Rotated Weierstrass
    if f == 12; Ub = 5.12;   end % Rotated Rastrigin
    if f == 13; Ub = 5.12;   end % Rotated Rastrigin_noncont
    if f == 14; Ub = 500;    end % Rotated Schewfel
    % =================  Rotation =======================
    s = swarmNum;
    K1 = mod(dim,swarmNum); K1len = ceil(dim/swarmNum);    orthm1 = [];
    K2 = swarmNum - K1;     K2len = floor((dim/swarmNum)); orthm2 = [];
    rc=1:dim; swarm1=[]; swarm2=[];
    if K1 == 0
        for i=1:K2
           orthm2 = [orthm2 orthm_generator(K2len)];
        end
    else
        for i=1:K1
            orthm1 = [orthm1 orthm_generator(K1len)];
        end
        for i=1:K1
            swarm1(i,:)=rc((i*K1len-K1len+1):(i*K1len));
        end
        if K2len == 1
            orthm2 = [];
        else
            for i=1:K2
               orthm2 = [orthm2 orthm_generator(K2len)];
            end
        end
    end
    for i=1:K2
        swarm2(i,:)=rc((K1*K1len+i*K2len-K2len+1):(K1*K1len+i*K2len));
    end
    
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