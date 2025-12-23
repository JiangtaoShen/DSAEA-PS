function [x] = IfDominate(a,b,m)

t=0;
q=0;
p=0;
for i=1:m
    if a(1,i)<=b(1,i)
        t=t+1;
    elseif a(1,i)>= b(1,i)
        q=q+1;
    else a(1,i)== b(1,i)
        p=p+1;
    end
end

if t==m&p~=m     
    x=1;
elseif q==m&p~=m 
    x=2;
else
    x=3;
end
end

