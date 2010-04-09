%G_BLOCK_IF Compute the reconstruction matrix for multiple IAF time decoder.
%   G = G_BLOCK_IF(TI,TJ) computes the reconstruction matrix
%   G[i,j] = <phi_k^i,psi_l^j> used to decode a signal in L2 space
%   that was encoded by a population of ideal IAF neurons. TI and
%   TJ contain the spike times from neurons i and j, respectively.
%
%   The calculation is described in further detail in Equation 29 of the
%   Consistent Recovery paper mentioned in the toolbox references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function Gb = G_block_IF(ti,tj)

li=length(ti)-1;
lj=length(tj)-1;

Gb = zeros(li,lj);

if isequal(ti,tj)
    for i=1:li
        for j=1:li
            tmz=ti(min(i,j));
            tmp=ti(min(i,j)+1);
            tpz=ti(max(i,j));
            tpp=ti(max(i,j)+1);
            Gb(i,j)=((tpp-tmz)^5-(tpp-tmp)^5+(tpz-tmp)^5-(tpz-tmz)^5)/20;
        end
    end
else
    for k=1:li
        for l=1:lj
            tip=ti(k+1);
            tim=ti(k);
            tjp=tj(l+1);
            tjm=tj(l);
            if tjp<=tim
                Gb(k,l)=((tip-tjm)^5+(tim-tjp)^5-(tim-tjm)^5-(tip-tjp)^5)/20;
            elseif (tjm<=tim)+(tim<=tjp)+(tjp<=tip)==3
                Gb(k,l)=((tip-tjm)^5-(tip-tjp)^5+(tjp-tim)^5-(tim-tjm)^5)/20;
            elseif (tjm<=tim)+(tip<tjp)==2
                Gb(k,l)=((tip-tjp)^5+(tip-tjm)^5-(tim-tjp)^5-(tim-tjm)^5)/20;
            elseif (tim<=tjm)+(tjm<=tip)==2
                Gb(k,l)=((tip-tjm)^5+(tip-tjm)^5-(tip-tjp)^5-(tim-tjp)^5)/20;
            elseif (tim<=tjm)+(tjm<=tip)+(tip<=tjp)==3
                Gb(k,l)=((tip-tjp)^5+(tip-tjm)^5-(tim-tjp)^5+(tim-tjm)^5)/20;
            else
                Gb(k,l)=((tip-tjp)^5+(tim-tjm)^5-(tim-tjp)^5-(tip-tjm)^5)/20;
            end
        end
    end
end
