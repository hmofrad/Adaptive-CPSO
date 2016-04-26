%   Adaptive Cooperative Particle Swarm Optimization
% function [gbest_hist fe_hist]=ACPSO(f,nop,dim,swarmNum,alpha,beta)
clear,clc
f= 16;nop = 30;dim=30;swarmNum=6;alpha=0.1;beta=0.1;
%     [gh fe]=ACPSO(15,30,30,3,0.01,0.01);
    endgen = 10^4;    % maximum Generation
    max_fe = 3*10^5; % current fitness evaluation
    global initial_flag action 
    initial_flag = 0;
    if f == 1; Bounds=[-100,100]; f_bias=-450;end
    if f == 2; Bounds=[-100,100]; f_bias=-450;end
    if f == 3; Bounds=[-100,100]; f_bias=-450;end
    if f == 4; Bounds=[-100,100]; f_bias=-450;end
    if f == 5; Bounds=[-100,100]; f_bias=-310;end
    if f == 6; Bounds=[-100,100]; f_bias=+390;end
    if f == 7; Bounds=[0,600];    f_bias=-180;end
    if f == 8; Bounds=[-32, 32];  f_bias=-140;end
    if f == 9; Bounds=[-5,5];     f_bias=-330;end
    if f == 10;Bounds=[-5,5];     f_bias=-330;end
    if f == 11;Bounds=[-0.5,0.5]; f_bias=+90;end
    if f == 12;Bounds=[-100,100]; f_bias=-460;end
    if f == 13,Bounds=[-3,1];     f_bias=-130;end
    if f == 14;Bounds=[-100,100]; f_bias=-300;end
    if f == 15;Bounds=[-5,5];     f_bias= 120;end
    if f == 16;Bounds=[-5,5];     f_bias= 120;end
    if f == 17;Bounds=[-5,5];	  f_bias= 120;end
    if f == 18;Bounds=[-5,5];	  f_bias=10;end
    if f == 19;Bounds=[-5,5];	  f_bias=10;end
    if f == 20;Bounds=[-5,5];     f_bias=10;end
    if f == 21;Bounds=[-5,5];     f_bias=360;end
    if f == 22;Bounds=[-5,5];	  f_bias=360;end
    if f == 23;Bounds=[-5,5]; 	  f_bias=360;end
    if f == 24;Bounds=[-5,5];	  f_bias=260;end
    if f == 25;Bounds=[-2,5]; 	  f_bias=260;end
    Ub=Bounds(2);
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    tr_num = 500; % train epochs
    c1 = 1.49445;    c2 = 1.49445;
    w_max = 0.9;  w_min = 0.4;
    % ======================================================
    pbest = x; %initialize Best Particle Position
    imp_tag = 0; % improvement TAG
    gbest_hist = []; % swarm best history
    fe_hist = [];
    % ====== Learning Automata initialization =======
    action = 1:swarmNum; r = length(action);
    p = repmat(1/r,swarmNum,dim); % action probability vector
    swarmTable = automataActSel(p);
    % ===== calculate sbest of each swarm =====
    B=[];
    for i=1:nop
        B = [B benchmark_func(x(i,:),f)-f_bias];
    end
	[mn ind] = min(B);
    gbest = x(ind,:); gbest(end+1) = mn;
    swarm_ind= []; dim_ind = [];
    i = 1; fe = 0; cpy=10;
        %   Train Learning Automata Phase
    for i=1:tr_num
        w = w_max-(w_max-w_min)*i/endgen;
        for j = 1:swarmNum
            swarmDim = find(swarmTable(j,:) == 1);
            swarmLength = length(swarmDim);
            if ~isempty(swarmLength)
                for k = 1:nop
                    if benchmark_func(b(gbest(1:end-1),swarmDim,x(k,swarmDim)),f)-f_bias<...
                            benchmark_func(b(gbest(1:end-1),swarmDim,pbest(k,swarmDim)),f)-f_bias
                        pbest(k,swarmDim) = x(k,swarmDim);
%                         imp_tag = imp_tag+1;
                    end
                    if benchmark_func(b(gbest(1:end-1),swarmDim,pbest(k,swarmDim)),f)-f_bias < gbest(end)
                        gbest(swarmDim) = pbest(k,swarmDim);
                        gbest(end) = benchmark_func(gbest(1:end-1),f)-f_bias;
                        imp_tag = 1;
                    end
                end
                spd(:,swarmDim) = w.* spd(:,swarmDim)...
                                +c1.*rand(nop,length(swarmDim)).*(pbest(:,swarmDim)-x(:,swarmDim))...
                                +c2.*rand(nop,length(swarmDim)).*(repmat(gbest(swarmDim),nop,1)-x(:,swarmDim));
                spd(:,swarmDim)=(spd(:,swarmDim)>Vmax).*Vmax+(spd(:,swarmDim)<=Vmax).*spd(:,swarmDim); 
                spd(:,swarmDim)=(spd(:,swarmDim)<(Vmin)).*(Vmin)+(spd(:,swarmDim)>=(Vmin)).*spd(:,swarmDim);
                x(:,swarmDim) = x(:,swarmDim)+spd(:,swarmDim);
                if i>1 && j == swarm_ind
                    p = automataProbUp(p,imp_tag,j,dim_ind,swarmTable,alpha,beta);
                end
                imp_tag = 0;
            end
        end 
        fe = fe + swarmNum*nop;
        dim_ind = mod(i,dim);
        if dim_ind == 0
            dim_ind = dim;
        end
        rc=randperm(swarmNum);   swarm_ind=rc(1);
        if swarmTable(swarm_ind,dim_ind) == 1
            swarm_ind=rc(2);
        end
        swarmTable(:,dim_ind) = 0;
        swarmTable(swarm_ind,dim_ind) = 1;
        if mod(i,cpy)==0,fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end)),end
    end
%             % ========= Play Learning Automata Phase ==========
    while i<=endgen && fe<=max_fe
        w = w_max-(w_max-w_min)*i/endgen;
        swarmTable = automataActSel(p);
        for j = 1:swarmNum
            swarmDim = find(swarmTable(j,:) == 1);
            swarmLength = length(swarmDim);
            if ~isempty(swarmLength)
                for k = 1:nop
                    if benchmark_func(b(gbest(1:end-1),swarmDim,x(k,swarmDim)),f)-f_bias<...
                            benchmark_func(b(gbest(1:end-1),swarmDim,pbest(k,swarmDim)),f)-f_bias
%                         imp_tag = imp_tag+1;
                    end
                    if benchmark_func(b(gbest(1:end-1),swarmDim,pbest(k,swarmDim)),f)-f_bias < gbest(end)
                        gbest(swarmDim) = pbest(k,swarmDim);
                        gbest(end) = benchmark_func(gbest(1:end-1),f)-f_bias;
                        imp_tag = 1;
                    end
                end
                spd(:,swarmDim) = w.* spd(:,swarmDim)...
                                +c1.*rand(nop,length(swarmDim)).*(pbest(:,swarmDim)-x(:,swarmDim))...
                                +c2.*rand(nop,length(swarmDim)).*(repmat(gbest(swarmDim),nop,1)-x(:,swarmDim));
                spd(:,swarmDim)=(spd(:,swarmDim)>Vmax).*Vmax+(spd(:,swarmDim)<=Vmax).*spd(:,swarmDim); 
                spd(:,swarmDim)=(spd(:,swarmDim)<(Vmin)).*(Vmin)+(spd(:,swarmDim)>=(Vmin)).*spd(:,swarmDim);
                x(:,swarmDim) = x(:,swarmDim)+spd(:,swarmDim);
                p = automataProbUp(p,imp_tag,j,swarmDim,swarmTable,alpha,beta);
                imp_tag = 0;
            end
        end
        fe = fe + swarmNum*nop;
       gbest_hist = [gbest_hist gbest(end)]; fe_hist = [fe_hist fe];
        if mod(i,cpy)==0,fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end)),end
        if ((i>1300 && gbest_hist(length(gbest_hist)-500) == gbest(end)) || gbest(end) ==0)
%            gbest_hist(end+1) =  gbest(end); fe_hist(end+1) = max_fe;
            break
        end
        i = i + 1;
    end
 
%     gbest(end)
%     plot(fe_hist,log(gbest_hist))
% end
