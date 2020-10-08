%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
syms x;
f={@(x)exp(x),@(x)1./(1+x.^2)};
I0={exp(1)-exp(-1),pi/2};
%I0记录标准积分值
for n = 2:7
    I1={0 0}; I2={0 0};
    S1=linspace(-1,1,n);
    %该步骤生成了一个以-1为首，1为尾的n项等差数列,由于太简单直接用内置函数了
    S2=sort(double(solve(diff(((x^2-1)^n),n),x)));
    %算正交多项式的根，得到Gauss点
    for i=1:floor((1+n)/2)
        g1=[1]; g2=[1];
        for j = 1:n
            if i~=j
                g1=conv(g1,[1/(S1(i)-S1(j)),-S1(j)/(S1(i)-S1(j))]);
                g2=conv(g2,[1/(S2(i)-S2(j)),-S2(j)/(S2(i)-S2(j))]);
                %conv是多项式相乘，以得到Lagrange多项式，也可以不调用内置函数写做
                %g1=([g1,0]-[0,g1]*S1(j))/(S1(i)-S1(j));
                %或者不用循环，用内置的deconv函数
                %[g1,s1]=deconv(w1,[1,-S1(i)]);
                %g1=g1/polyval（g1,S1(i));
                %其中w1为ppt中的w函数
            end
        end
        k1=0; k2=0;
        for j=0:ceil(n/2)-1
            k1=k1+g1(n-2*j)*2/(2*j+1);
            k2=k2+g2(n-2*j)*2/(2*j+1);
        end    
        %多项式积分可以直接用相应单项式系数乘以相应单项式积分值再相加即可
        %数值积分这里用到了S1(i),S1(n+1-i)的Gauss(或Newton)系数相同（对称性），故只算了一半    
        for m=1:2
            if i*2==n+1
                I1{m}=I1{m}+f{m}(S1(i))*k1;
                I2{m}=I2{m}+f{m}(S2(i))*k2;
            else
                I1{m}=I1{m}+(f{m}(S1(i))+f{m}(S1(n+1-i)))*k1;
                I2{m}=I2{m}+(f{m}(S2(i))+f{m}(S2(n+1-i)))*k2;
            end
        end
    end
    for m=1:2
        fprintf('f%d在%d点的Newton-Cotes积分下的值为%f，误差为%f\n',m,n,I1{m},abs(I1{m}-I0{m}));
        fprintf('f%d在%d点的Gauss积分下的值为%f，误差为%f\n',m,n,I2{m},abs(I2{m}-I0{m}));    
    end
end