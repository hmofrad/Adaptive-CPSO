function [out] = b(sbest,swarmind,pos)
    sbest(swarmind) = pos;
    out = sbest;
end