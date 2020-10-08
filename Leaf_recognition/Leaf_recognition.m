A=zeros(8,8);
A(1,1:1:7)=Feature_extraction('1-1.jpg');
A(2,1:1:7)=Feature_extraction('1-2.jpg');
A(3,1:1:7)=Feature_extraction('1-3.jpg');
A(4,1:1:7)=Feature_extraction('1-4.jpg');
A(5,1:1:7)=Feature_extraction('2-1.jpg');
A(6,1:1:7)=Feature_extraction('2-2.jpg');
A(7,1:1:7)=Feature_extraction('2-3.jpg');
A(8,1:1:7)=Feature_extraction('2-4.jpg');
%前四个为银杏的向量值，后四个为枫叶的向量值
%几张图的分辨率都调到了1024*1024.据论文要求
%实际只需要四的倍数就行了，但分辨率最好一致，不然还要添加处理
B=zeros(8);
disp('向量值');
disp(A);
for i=1:8
    for j=1:i
        B(i,j)=myint(A(i,1:1:7),A(j,1:1:7));
        B(j,i)=B(i,j);
    end
end
disp('差异值');
disp(B);

function a=Feature_extraction(x)
    y=perprocessing(x);
    [d,h,v,cd,ch,cv,cc]=wavelet_decompose(y);
    a1=fractal_dimension(d);
    a2=fractal_dimension(h);
    a3=fractal_dimension(v);
    a4=fractal_dimension(cd);
    a5=fractal_dimension(ch);
    a6=fractal_dimension(cv);
    a7=fractal_dimension(cc);
    a=[a1 a2 a3 a4 a5 a6 a7];
end

function w=perprocessing(x) 
    w=double(imread(x));
    w=rgbtogray(w);
    %可以把w再转化为uint8，以显示，但没必要
end

function y=rgbtogray(x)
    a=size(x,1);
    b=size(x,2);
    y=zeros(a,b);
    for i=1:a
        for j=1:b
            y(i,j)=x(i,j,1)*0.2989+x(i,j,2)*0.5870+x(i,j,3)*0.1141;
            %有些地方用的是0.1140
        end
    end
end 

function [d,h,v,cd,ch,cv,cc]=wavelet_decompose(y)
    arguments
        y (1024,1024)
    end
    [d,h,v,c]=lifting_wavelet(y);
    [cd,ch,cv,cc]=lifting_wavelet(c);
    %按照原论文的论述，做两次小波后，识别正确率就已足够
end

function [d,h,v,c]=lifting_wavelet(y)
    a=size(y,1);
    b=zeros(a);
    for i=1:a
        b(i,1:1:end)=[y(i,1:2:a),y(i,2:2:a)];%懒变换
        b(i,a/2+1:1:a)=b(i,a/2+1:1:a)-(b(i,1:1:a/2)+[b(i,2:1:a/2),b(i,a/2)])/2;%预测
        b(i,1:1:a/2)=b(i,1:1:a/2)+(b(i,a/2+1:1:a)+[b(i,a/2+2:1:a),b(i,a)])/4;%更新
    end
    for i=1:a
        b(1:1:end,i)=[b(1:2:a,i);b(2:2:a,i)];%懒变换
        b(a/2+1:1:a,i)=b(a/2+1:1:a,i)-(b(1:1:a/2,i)+[b(2:1:a/2,i);b(a/2,i)])/2;%预测
        b(1:1:a/2,i)=b(1:1:a/2,i)+(b(a/2+1:1:a,i)+[b(a/2+2:1:a,i);b(a,i)])/4;%更新
    end
    c=b(1:1:a/2,1:1:a/2);
    d=b(a/2+1:1:a,a/2+1:1:a);
    v=b(a/2+1:1:a,1:1:a/2);
    h=b(1:1:a/2,a/2+1:1:a);
end

function D=fractal_dimension(I)
    J=zeros([size(I),5]);
    a=size(I,1);
    b=size(I,2);
    J(1:1:a,1:1:b,1)=I;
    K=J;
    L=zeros(1,5);
    A=zeros(1,4);
    for i=2:5
        for j=1:a
            for k=1:b
                if j==1
                    c0=J(j+1,k,i-1);
                elseif j==a
                    c0=J(j-1,k,i-1);
                else
                    c0=[J(j+1,k,i-1),J(j-1,k,i-1)];
                end
                if k==1
                    d0=J(j,k+1,i-1);
                elseif k==b
                    d0=J(j,k-1,i-1);
                else
                    d0=[J(j,k+1,i-1),J(j,k-1,i-1)];
                end
                %这两段用intersect函数写更简洁，但不知为何，运行速度太慢
                J(j,k,i)=max([J(j,k,i-1)+1,c0,d0]);
                K(j,k,i)=min([K(j,k,i-1)-1,c0,d0]);
                L(i)=L(i)+J(j,k,i)-K(j,k,i);
            end
        end
     A(i-1)=(L(i)-L(i-1))/2;
    end
    %按照原论文提供的blanket method计算分形维数
    D=polyfit(log(1:1:4),log(A),1);
    D=-D(1)+2;
end

function c=myint(a,b)
    arguments
        a (1,7)
        b (1,7)
    end
    d=zeros(1,7);
    for i=1:7
        d(i)=abs(a(i)-b(i));
        %无法定义在这个向量上定义完全合定义的内积
        %分形维数的加法没有意义，内积无线性性
    end
    c=-d*[4 4 4 1 1 1 1]'/16;
    %赋权时候假设horizontal detail images(h), vertical detail images(v), 
    %diagonal detail images(d) and coarse images(c)各占1/4的细节
    %因此对h,v,d,ch,cv,cd,cc各赋权1/4,1/4,1/4,1/16,1/16,1/16,1/16
end
    