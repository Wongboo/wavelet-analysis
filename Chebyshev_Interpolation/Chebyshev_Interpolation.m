global x;
global z;
syms x;
z=[1 x zeros(1,78)];
for i=3:80
    z(i)=expand(2*x*z(i-1)-z(i-2));
end
%将z和x设置为全局变量，防止多次计算切比雪夫多项式
chebyshevplot(@f1);
chebyshevplot(@f2);
y=log(2:80);
[y1,w1]=loglerror(@f1);
[y2,w2]=loglerror(@f2);
figure;
plot(y,y1-y*3)
%对于f1，只能验证定理2，不能验证定理3，为3次收敛
figure;
plot(y,y2-y*10)
% 对于f2，定理2说明其收敛速度快于任意有限次收敛
% 显然这不好验证，这里我们只验证了快于十次收敛
% 幸运的是，对f2，验证定理3就意味这验证了定理2
figure;
plot((2:80),w2/log(1.2)-(2:80))
function y=f1(s)
    y=abs(sin(6*s)).^3-cos(5*exp(s));
end
function y=f2(s)
    y=(1+25*s.^2).^-1-sin(20*s);
end
function  y=dft(x)
    arguments
        x (1,:) double 
    end
    a=size(x,2); y=zeros(1,a);
    for j=1:a
        y(j)=x*exp(pi*2i*(j-1)*(0:a-1)/a)';
    end
end   
function f=interpolation(b,c)
    arguments
        b
        c (1,1) {mustBePositive,mustBeInteger}
    end
    global z; global x;
    d=b(cos(linspace(0,pi,c+1)));
    d=real(dft([d,flip(d(:,2:c),2)]));
    a=d(:,1:c+1)/c;
    a(1)=a(1)/2; a(c+1)=a(c+1)/2;
    e(x)=dot(a,z(:,1:c+1));
    f=e((-1:0.01:1));
end
function chebyshevplot(b)
    t=-1:0.01:1;
    b0=b(t);
    b1=interpolation(b,24);
    b2=interpolation(b,39);
    b3=interpolation(b,64);
    b4=interpolation(b,79);
    figure;
    plot(t,b0,t,b1,t,b2,t,b3,t,b4);
    legend('原曲线','插值25点曲线','插值40点曲线','插值65点曲线','插值80点曲线')
end
function [y,w]=loglerror(b)
    t=b(-1:0.01:1); y=zeros(1,79); w=zeros(1,79);
    for i=1:79
        ts=interpolation(b,i);
        ts=((t-ts).^2);
        y(i)=-log(mean(ts))/2;
        w(i)=-log(max(ts))/2;
        %不用内置的integral函数，利用mean函数实现
        %更精确点，y(i)=-log(mean(ts)*2)/2;只差一个常数，可以忽略
    end
end

    
    

