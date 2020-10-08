%没用，该程序是补零后进行fft操作
%插值程序里用不到
disp(dft([1 2 8]));
function f=deepdft(a,c)
    if c==0
        f=a;
    else
       f1=deepdft(a(1:2:end),c-1);
       f2=deepdft(a(2:2:end),c-1).*exp((0:2^(c-1)-1)*-pi*1i/2^(c-1));
       f=[f1+f2,f1-f2]; 
    end
end
function f=dft(a)
    arguments
        a (1,:) double
    end
    b=size(a,2);
    c=ceil(log2(b));
    a=[a,zeros(2^c-b)];
    f=deepdft(a,c);
end