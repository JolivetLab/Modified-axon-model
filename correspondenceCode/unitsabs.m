function out=unitsabs(in)
% units
units={'s','m','A','V','S','F','O'};
scale={'m','m','u','m','m','u','k'};
cal1={'G','M','k',' ','d','c','m','u','n','p','f'};
cal2={1e9,1e6,1e3,1,1e-1,1e-2,1e-3,1e-6,1e-9,1e-12,1e-15};

out=zeros(1,in{1});
for i=1:in{1}
    a=cal2(strcmp(cal1,in{i+1}(1)));
    b=scale(strcmp(units,in{i+1}(2)));
    c=cal2(strcmp(cal1,b));
    out(i)=a{1}/c{1};
end
out=prod(out.^in{end});