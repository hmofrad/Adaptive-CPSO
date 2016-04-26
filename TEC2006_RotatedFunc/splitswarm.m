
    % =================  Rotation =======================
    K1 = mod(dim,swarmNum); K1len = ceil(dim/swarmNum);    orthm1 = [];
    K2 = swarmNum - K1;     K2len = floor((dim/swarmNum)); orthm2 = [];
    rc=1:dim;
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
    