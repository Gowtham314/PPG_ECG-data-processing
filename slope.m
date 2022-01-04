function y = slope(x,n)
for i=1:length(x)-n;
    if x(i+n)> x(i)
       y(i)=x(i+n)-x(i);
    else
       y(i)=0;
    end
    
end
    end

