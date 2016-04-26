cd 'C:\Users\M0RBiD\Documents\MATLAB\PSO implementation\CPSOLA_Final\CEC2005'
clear
clc
global initial_flag orthm swarmNum action alpha beta 
% f = 1; Bounds=[-100,100]; f_bias=-450;
% f = 2; Bounds=[-100,100]; f_bias=-450;
% f = 3; Bounds=[-100,100]; f_bias=-450;
% f = 4; Bounds=[-100,100]; f_bias=-450;
% f = 5; Bounds=[-100,100]; f_bias=-310;
% f = 6; Bounds=[-100,100]; f_bias=+390;
% f = 7; Bounds=[0,600];    f_bias=-180;
% f = 8; Bounds=[-32, 32];  f_bias=-140;
% f = 9; Bounds=[-5,5];     f_bias=-330;
% f = 10;Bounds=[-5,5];     f_bias=-330;
% f = 11;Bounds=[-0.5,0.5]; f_bias=90;
% f = 12;Bounds=[-100,100]; f_bias=-460; 
% f = 13;Bounds=[-3,1];     f_bias=-130;
% f = 14;Bounds=[-100,100]; f_bias=-300;
% f = 15;Bounds=[-5,5];     f_bias= 120;
% f = 16;Bounds=[-5,5];     f_bias= 120;
% f = 17;Bounds=[-5,5];	  f_bias= 120;
% f = 18;Bounds=[-5,5];	  f_bias=10;
% f = 19;Bounds=[-5,5];	  f_bias=10;
% f = 20;Bounds=[-5,5];     f_bias=10;
f = 21;Bounds=[-5,5];     f_bias=360;
% f = 22;Bounds=[-5,5];	  f_bias=360;
% f = 23;Bounds=[-5,5]; 	  f_bias=360;
% f = 24;Bounds=[-5,5];	  f_bias=260;
% f = 25;Bounds=[-2,5]; 	  f_bias=260;

nop = 30;dim = 25;swarmNum = 3;
    Ub = Bounds(2);
    Lb = -Ub;	bnd = [Lb Ub];     vm = Ub;
    Vmin = -0.2*(Ub-Lb);	Vmax = -Vmin;
    spd=Vmin+2.*Vmax.*rand(nop,dim);  % initialize velocity
    endgen = 10^3;    % maximum Generation
    max_fe = 10^6; % current fitness evaluation
    fe = 0;
    w=0.9-(1:endgen)*(0.7/endgen);
    s = 1*10^5; % Maximum Stagnation epochs
    tr_num = 1*10^5; % train epochs
    alpha = 0.01; beta = 0.01;
    c1 = 1.49445;    c2 = 1.49445;
    w_max = 0.9;  w_min = 0.4;

    x=Lb+(Ub-Lb).*rand(nop,dim);  % initialize position
    
    pbest = x; %initialize Best Particle Position

    sbest_fit =  []; %  swarm best fitness
    sbest_hist = []; % swarm best history

    % ===== calculate sbest of each swarm =====
    [mn ind] = min(abs(benchmark_func(x,f)-f_bias));
    sbest = x(ind,:);

    % ====== Learning Automata initialization =======
    action = 1:swarmNum; r = length(action);
    p = repmat(1/r,swarmNum,dim); % action probability vector
%     p(1,1:10) = 0.5;p(1,11:20) = 0.25;p(1,21:30) = 0.25;
%     p(2,1:10) = 0.25;p(2,11:20) = 0.5;p(2,21:30) = 0.25;
%     p(3,1:10) = 0.25;p(3,11:20) = 0.25;p(3,21:30) = 0.5;
    i = 1;
    while i<=endgen && fe<=max_fe
        w = w_max-(w_max-w_min)*i/endgen;
        %   Train Learning Automata Phase
%         if i < tr_num
            
%         else
        %   Play Learning Automata Phase
            swarmTable = automataActSel(p);
            for j = 1:swarmNum
                swarmDim = find(swarmTable(j,:) == 1);
                swarmLength = length(swarmDim);
                if ~isempty(swarmDim)
                    sbest_tmp = benchmark_func(sbest,f)+f_bias;
                    for jj =1:swarmLength
                        for k = 1:nop
                            
                            if benchmark_func(b(sbest,swarmDim(jj),x(k,swarmDim(jj))),f)+f_bias...
                                    < benchmark_func(b(sbest,swarmDim(jj),x(k,swarmDim(jj))),f)+f_bias
                                pbest(k,swarmDim(jj)) = x(k,swarmDim(jj));
                            end
                            
                            if benchmark_func(b(sbest,swarmDim(jj),pbest(k,swarmDim(jj))),f)+f_bias < benchmark_func(sbest,f)+f_bias
                                sbest(swarmDim(jj)) = pbest(k,swarmDim(jj));
                            end
                            spd(k,swarmDim(jj)) = w.* spd(k,swarmDim(jj))...
                                                +c1.*rand(1,length(swarmDim(jj))).*(pbest(k,swarmDim(jj))-x(k,swarmDim(jj)))...
                                                +c2.*rand(1,length(swarmDim(jj))).*(sbest(swarmDim(jj))  -x(k,swarmDim(jj)));
%                             spd(k,swarmDim(jj))=(spd(k,swarmDim(jj))>Vmax).*Vmax+(spd(k,swarmDim(jj))<=Vmax).*spd(k,swarmDim(jj)); 
%                             spd(k,swarmDim(jj))=(spd(k,swarmDim(jj))<(Vmin)).*(Vmin)+(spd(k,swarmDim(jj))>=(Vmin)).*spd(k,swarmDim(jj));
                            spdtmp = spd(k,swarmDim(jj));               ind = find(abs(spdtmp>vm));
                            spdtmp(ind) = vm.*sign(spdtmp(ind));	spd(k,swarmDim(jj)) = spdtmp;
                            x(k,swarmDim(jj)) = x(k,swarmDim(jj))+spd(k,swarmDim(jj));
                            sbest_fit = benchmark_func(sbest,f)+f_bias;
                            sbest_hist = [sbest_hist sbest_fit];
                        end
                        fe = fe + nop;
                    end
                    p = automataProbUp(p,swarmDim,swarmTable,sbest_fit,sbest_tmp);
                end
            end    
%         end

    %  Ackley   Griewank Rastrigin_N   Schewfel
%     sigma = sigma * e^ tau .* normrnd(0,1,[nop dim]);
%     x = x + sigma .* normrnd(0,1,[nop dim]);
    
        fprintf('Fitness Evaluation=%u,SBEST=%e\n',fe,sbest_fit)
        if (mod(fe,s) == 0 && (sbest_hist(fe-s+1) == sbest_fit)) || sbest_fit ==0, break, end
        i = i + 1;
    end
    sbest_hist(end:max_fe) = sbest_hist(end);
    plot(log(sbest_hist))
