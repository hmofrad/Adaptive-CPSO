% function [gbest_hist fe_hist]=ACPSO(f,nop,dim,swarmNum,alpha,beta)
    clear,clc
    f= 11;nop = 40;dim=30;swarmNum=3;alpha=0.1;beta=0.0;
    %[gh fe]=ACPSO(13,30,30,3,0.01,0.01);
    endgen = 10^4;    % maximum Generation
    max_fe = 10^6; % current fitness evaluation
    global  orthm1 orthm2 swarm1 swarm2 K1 K2 K1len K2len action
    if f == 9;  Ub = 32.768; end % Rotated Ackley
    if f == 10; Ub = 600;    end % Rotated Griewank
    if f == 11; Ub = 0.5;    end % Rotated Weierstrass
    if f == 12; Ub = 5.12;   end % Rotated Rastrigin
    if f == 13; Ub = 5.12;   end % Rotated Rastrigin_noncont
    if f == 14; Ub = 500;    end % Rotated Schewfel
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    tr_num = 500; % train epochs
    c1 = 1.49445;    c2 = 1.49445;
    w_max = 0.9;  w_min = 0.4;
    
    % =================  Rotation =======================
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
	[mn ind] = min(fit_func({f},x));
    gbest = x(ind,:); gbest(end+1) = mn;
    swarm_ind= []; dim_ind = [];
    i = 1; fe = 0;
    while i<=endgen && fe<=max_fe
        w = w_max-(w_max-w_min)*i/endgen;
        %   Train Learning Automata Phase
        if i < tr_num
            for j = 1:swarmNum
                swarmDim = find(swarmTable(j,:) == 1);
                swarmLength = length(swarmDim);
                if ~isempty(swarmLength)
                    for k = 1:nop
                        if fit_func({f},b(gbest(1:end-1),swarmDim,x(k,swarmDim)))<fit_func({f},b(gbest(1:end-1),swarmDim,pbest(k,swarmDim)))
                            pbest(k,swarmDim) = x(k,swarmDim);
                        end
                        if fit_func({f},b(gbest(1:end-1),swarmDim,pbest(k,swarmDim))) < gbest(end)
                            gbest(swarmDim) = pbest(k,swarmDim);
                            gbest(end) = fit_func({f},gbest(1:end-1));
                            imp_tag = 1;
                        end
                    end
                    spd(:,swarmDim) = w.* spd(:,swarmDim)...
                                    +c1.*rand(nop,length(swarmDim)).*(pbest(:,swarmDim)-x(:,swarmDim))...
                                    +c2.*rand(nop,length(swarmDim)).*(repmat(gbest(swarmDim),nop,1)-x(:,swarmDim));
                    spd(:,swarmDim)=(spd(:,swarmDim)>Vmax).*Vmax+(spd(:,swarmDim)<=Vmax).*spd(:,swarmDim); 
                    spd(:,swarmDim)=(spd(:,swarmDim)<(Vmin)).*(Vmin)+(spd(:,swarmDim)>=(Vmin)).*spd(:,swarmDim);
                    x(:,swarmDim) = x(:,swarmDim)+spd(:,swarmDim);
                    fe = fe + nop;
                    if i>1 && j == swarm_ind
                        p = automataProbUp(p,imp_tag,j,dim_ind,swarmTable,alpha,beta);
                    end
                    imp_tag = 0;
                end
            end 
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
        else
            % ========= Play Learning Automata Phase ==========
            swarmTable = automataActSel(p);
            for j = 1:swarmNum
                swarmDim = find(swarmTable(j,:) == 1);
                swarmLength = length(swarmDim);
                if ~isempty(swarmLength)
                    for k = 1:nop
                        if fit_func({f},b(gbest(1:end-1),swarmDim,x(k,swarmDim)))<fit_func({f},b(gbest(1:end-1),swarmDim,pbest(k,swarmDim)))
                            pbest(k,swarmDim) = x(k,swarmDim);
                        end
                        if fit_func({f},b(gbest(1:end-1),swarmDim,pbest(k,swarmDim))) < gbest(end)
                            gbest(swarmDim) = pbest(k,swarmDim);
                            gbest(end) = fit_func({f},gbest(1:end-1));
                            imp_tag = 1;
                        end
                    end
                    spd(:,swarmDim) = w.* spd(:,swarmDim)...
                                    +c1.*rand(nop,length(swarmDim)).*(pbest(:,swarmDim)-x(:,swarmDim))...
                                    +c2.*rand(nop,length(swarmDim)).*(repmat(gbest(swarmDim),nop,1)-x(:,swarmDim));
                    spd(:,swarmDim)=(spd(:,swarmDim)>Vmax).*Vmax+(spd(:,swarmDim)<=Vmax).*spd(:,swarmDim); 
                    spd(:,swarmDim)=(spd(:,swarmDim)<(Vmin)).*(Vmin)+(spd(:,swarmDim)>=(Vmin)).*spd(:,swarmDim);
                    x(:,swarmDim) = x(:,swarmDim)+spd(:,swarmDim);
                    fe = fe + nop;
                    p = automataProbUp(p,imp_tag,j,swarmDim,swarmTable,alpha,beta);
                    imp_tag = 0;
                end
            end    
        end
        gbest_hist = [gbest_hist gbest(end)]; fe_hist = [fe_hist fe];
        if mod(i,500)==0,fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end)),end
        if ((i>1300 && gbest_hist(length(gbest_hist)-500) == gbest(end)) || gbest(end) ==0)
%            gbest_hist(end+1) =  gbest(end); fe_hist(end+1) = max_fe;
            break
        end
        i = i + 1;
    end
    gbest(end)
%     plot(fe_hist,log(gbest_hist))
% end