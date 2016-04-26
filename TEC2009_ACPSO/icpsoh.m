%   Idealized CPSO-H
function [gbest_hist fe_hist] = icpsoh(f,nop,dim,swarmNum)
%     clear,clc
%     f= 11;nop = 50;dim=300;swarmNum=60;alpha=0.1;beta=0.1;
    endgen = 10^5;    % maximum Generation
    max_fe =  3750000; % 3750000 rcurrent fitness evaluation
    if f == 8;  Ub = 500;    end % Schwefel
    if f == 9;  Ub = 5.12;   end % Rastrigin
    if f == 10; Ub = 32.768; end % Ackley
    if f == 11; Ub = 600;    end % Griewank
    if f == 12; Ub = 50;     end % Penalized1
    if f == 13; Ub = 50;     end % Penalized2

    % =================  Rotation =======================
    s = swarmNum;
    K1 = mod(dim,swarmNum); K1len = ceil(dim/swarmNum);
    K2 = swarmNum - K1;     K2len = floor((dim/swarmNum));
    rc=1:dim;swarm=[];
    for i=1:K2
        swarm(i,:)=rc((K1*K1len+i*K2len-K2len+1):(K1*K1len+i*K2len));
    end
    
    Lb = -Ub;    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    c1 = 1.49445; c2 = 1.49445;
    w_max = 0.9;  w_min = 0.4;
    [bst ind]=min(fit_func({f},x));
    gbest = x(ind,:); gbest(end+1) = bst;
    pbest = x; %initialize Best Particle Position
    gbest_hist = []; fe_hist=[];
    % ========================
    % initilise Q Swarm parameters
    % ========================
    x1 = Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    spd1 = Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    x1(:,end+1)=0;
    for i=1:nop
        x1(i,end) = fit_func({f},x1(i,1:end-1));
    end
    pbest1 = x1; %initialise Best Particle Position
    [bst ind] = min(x1(:,end));
    gbest1 = x1(ind,:); % initialise global best position
    fe =0;i=1;
    while i<=endgen && fe<=max_fe
        w = w_max-(w_max-w_min)*i/endgen;
        if mod(i,2) == 1
            for ii=1:K2
                for k=1:nop
                    if fit_func({f},b(gbest(1:end-1),swarm(ii,:),x(k,swarm(ii,:)))) <...
                            fit_func({f},b(gbest(1:end-1),swarm(ii,:),pbest(k,swarm(ii,:))))
                        pbest(k,swarm(ii,:)) = x(k,swarm(ii,:));
                    end
                    if fit_func({f},b(gbest(1:end-1),swarm(ii,:),pbest(k,swarm(ii,:)))) < gbest(end)
                        gbest(swarm(ii,:)) = pbest(k,swarm(ii,:));
                        gbest(end)=fit_func({f},gbest(1:end-1));
                    end
                end
                spd(:,swarm(ii,:)) = w.* spd(:,swarm(ii,:))...
                    +c1.*rand(nop,K2len).*(pbest(:,swarm(ii,:))-x(:,swarm(ii,:)))...
                    +c2.*rand(nop,K2len).*(repmat(gbest(swarm(ii,:)),nop,1)-x(:,swarm(ii,:)));
                spd(:,swarm(ii,:))=(spd(:,swarm(ii,:))>Vmax).*Vmax+(spd(:,swarm(ii,:))<=Vmax).*spd(:,swarm(ii,:)); 
                spd(:,swarm(ii,:))=(spd(:,swarm(ii,:))<(Vmin)).*(Vmin)+(spd(:,swarm(ii,:))>=(Vmin)).*spd(:,swarm(ii,:));
                x(:,swarm(ii,:)) = x(:,swarm(ii,:))+spd(:,swarm(ii,:)); 
            end
            fe = fe + swarmNum*nop;
            rc=randperm(round(nop/2));   ind=rc(1);
            if x1(ind,end) == gbest1(end)
                ind=rc(2);
            end
            x1(ind,:) =  gbest;
        else
            % Q SWARM
            for j=1:nop
                if (x1(j,end) < pbest1(j,end))
                    pbest1(j,:) = x1(j,:);
                end
                if (pbest1(j,end) < gbest1(end))
                    gbest1 = pbest1(j,:);
                end
            end
            spd1 = w.* spd1+c1.*rand(nop,dim).*(pbest1(:,1:end-1)-x1(:,1:end-1))...
                               +c2.*rand(nop,dim).*(repmat(gbest1(1:end-1),nop,1) - x1(:,1:end-1));
            spd1=(spd1>Vmax).*Vmax+(spd1<=Vmax).*spd1; 
            spd1=(spd1<(Vmin)).*(Vmin)+(spd1>=(Vmin)).*spd1;
            x1(:,1:end-1) = x1(:,1:end-1)+spd1;
            x1(:,end) = fit_func({f},x1(:,1:end-1));
            fe = fe+nop;
        % ==== Interval THREE ====
            for ii=1:K2
                rc=randperm(round(nop/2));   ind=rc(1);
                if x(ind,swarm(ii,:)) == gbest(swarm(ii,:))
                    ind=rc(2);
                end
                x(ind,swarm(ii,:)) = gbest1(swarm(ii,:));
            end
        end
        gbest_hist = [gbest_hist gbest(end)]; fe_hist = [fe_hist fe];
        if mod(i,100)==0,fprintf('fun=%u,Gene=%u,Fit_Eval=%u,Gbest=%e\n',f,i,fe,gbest(end)),end
        if ((i>1300 && gbest_hist(length(gbest_hist)-500) == gbest(end)) || gbest(end) ==0),break,end
        i=i+1;
    end
end