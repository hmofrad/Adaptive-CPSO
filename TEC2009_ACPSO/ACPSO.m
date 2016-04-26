function [gbest_hist fe_hist]=ACPSO(f,nop,dim,swarmNum,alpha,beta)
%     clear,clc
%     f= 8;nop = 50;dim=300;swarmNum=60;alpha=0.1;beta=0.1;
    endgen = 10^5;    % maximum Generation
    max_fe =  5*10^6; % 3750000 rcurrent fitness evaluation
    global action 
    if f == 8;  Ub = 500;    end % Schwefel
    if f == 9;  Ub = 5.12;   end % Rastrigin
    if f == 10; Ub = 32.768; end % Ackley
    if f == 11; Ub = 600;    end % Griewank
    if f == 12; Ub = 50;     end % Penalized1
    if f == 13; Ub = 50;     end % Penalized2
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    tr_num = 300; % train epochs
    c1 = 1.49445;    c2 = 1.49445;
    w_max = 0.9;  w_min = 0.4;
    % =================  Rotation =======================
%     orthm = orthm_generator(dim);
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
    i = 0; fe = 0;
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
                    if i>0 && j == swarm_ind
                        p = automataProbUp(p,imp_tag,j,dim_ind,swarmTable,alpha,beta);
                    end
                    imp_tag = 0;
                end
            end 
            dim_ind = mod(i,dim)+1;
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
    if f == 8, gbest_hist(end) = gbest_hist(end) + (-12569.5); end
    gbest(end)
