for i=1:8
    A(i,8)=sqrt(myint(A(i,1:1:7),A(i,1:1:7)));
    for j=1:i
        B(i,j)=myint(A(i,1:1:7),A(j,1:1:7))/(A(i,8)*A(j,8));
        B(j,i)=B(i,j);
    end
end
disp('相似度');
disp(B);

function c=myint(a,b)
    arguments
        a (1,7)
        b (1,7)
    end
    d=zeros(1,7);
    for i=1:7
        d(i)=(a(i)-2.5)*(b(i)-2.5);
        %可改成d(i)=(a(i))*(b(i));
        %或者自己定义的“内积”
    end
    c=d*[4 4 4 1 1 1 1]'/16;
end