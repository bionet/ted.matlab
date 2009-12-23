function Gb = Gblock_IF(ti,tj)

li=length(ti)-1;
lj=length(tj)-1;

Gb = zeros(li,lj);

if isequal(ti,tj) %sum(ti(1:min(li,lj)+1)==tj(1:min(li,lj)+1))==length(ti)
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