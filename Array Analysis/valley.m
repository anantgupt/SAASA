function out=valley(in)
temp=in(end);
for i=length(in):-1:1
    if in(i)>temp
            temp=in(i);
    end
out(i)=temp;
end
end