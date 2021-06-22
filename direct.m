%MIT License
%Copyright (c) [2021] [Shuang Qin]
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.
function x = direct(X,U,t,bias)
% The first row of X is the real position of target
% The rest are sensor coordinates vector, we suppose 6 sensors here
%X=[10,-20,20,20,-20,0,20;8,-20,-20,20,20,-20,0]
% t=U*(measurements of each sensor)
% U=[-1,1,0,0,0,0;0,1,-1,0,0,0;0,1,0,-1,0,0;0,1,0,0,-1,0;0,1,0,0,0,-1]
% 1 corresponds to the reference sensor and
% is also the sensor with the largest measured value
% bias is pre-set
di = length(X(1,:)); %dimension of position
xr = [];
txr = [];
x0 = mean(X(2:end,:),1)';%expasion potin
for i =1:size(U, 1)+1
    a(i) = norm(x0-X(i+1,:)');%distance from sensor to expansion point
end
for i =1:size(U, 1)
    l = find(U(i,:) ==  1);
    m = find(U(i,:) == -1);
    if t(i) >= 0
        xr = [xr,[m;l]];
        txr = [txr,t(i)];
    else
        xr = [xr,[l;m]]; 
        txr = [txr,-t(i)];
    end
end
xc = unique(xr);%Subscripts of all participating sensors
%% construct the augmented matrix
m = length(xr);
B = zeros(di+5*m,1);
B(1:m) = 1;
for i =1:m
    j = find(xc==xr(1,i));
    r = txr(i);
    si = X(xr(2,i)+1,:)';
    sj = X(xr(1,i)+1,:)';
    c1 = bias^2+2*bias*r+r^2+norm(sj)^2-norm(si)^2;
    c2 = -bias^2+2*bias*r-r^2-norm(sj)^2+norm(si)^2;
    B(di+2*m+i)=-c1;
    B(di+3*m+i)=-c2;
    B(di+4*m+i)=(x0-sj)'*x0 - a(j)^2;
end
BB = B;
A = zeros(di+5*m,di+5*m);
for i=1:m
    A(i,di+2*m+i) = 1;
    A(i,di+3*m+i) = 1;
end
for i=1:m
    j = find(xc==xr(1,i));
    r = txr(i);
    si = X(xr(2,i)+1,:)';
    sj = X(xr(1,i)+1,:)';
    A(m+1:m+di,di+2*m+i)=2*(si-sj);
    A(m+1:m+di,di+3*m+i)=-2*(si-sj);
    A(m+1:m+di,di+4*m+i)=x0-sj;
end
for i=1:m
    j = find(xc==xr(1,i));
    r = txr(i);
    si = X(xr(2,i)+1,:)';
    sj = X(xr(1,i)+1,:)';
    A(di+m+i,di+2*m+i)=2*(bias+r);
    A(di+m+i,di+3*m+i)=2*(bias-r);
    A(di+m+i,di+4*m+i)=-a(j);
end

for i=1:m
    j = find(xc==xr(1,i));
    r = txr(i);
    si = X(xr(2,i)+1,:)';
    sj = X(xr(1,i)+1,:)';
    A(di+2*m+i,1:di)=2*(si-sj)';
    A(di+2*m+i,di+i)=2*(bias+r);
    A(di+2*m+i,di+m+i)=-1;
end
for i=1:m
    j = find(xc==xr(1,i));
    r = txr(i);
    si = X(xr(2,i)+1,:)';
    sj = X(xr(1,i)+1,:)';
    A(di+3*m+i,1:di)=-2*(si-sj)';
    A(di+3*m+i,di+i)=2*(bias-r);
    A(di+3*m+i,di+m+i)=-1;
end
for i=1:m
    j = find(xc==xr(1,i));
    r = txr(i);
    si = X(xr(2,i)+1,:)';
    sj = X(xr(1,i)+1,:)';
    A(di+4*m+i,1:di)=(x0-sj)';
    A(di+4*m+i,di+i)=-a(j);
end
AA = A;
%% solve
for i1=1:4 %6 sensors
    for i2=1:4
        for i3=1:4
            for i4=1:4
                for i5=1:4
                    res=[i1,i2,i3,i4,i5];
                    B1 = ones(di+5*m,1);
                    B2 = ones(di+5*m,1);
                    for ii = 1:length(res)
                        if res(ii)==4
                            continue;
                        elseif res(ii)==1
                            B1(di+2*m+ii) = 0;
                            B2(di+2*m+ii) = 0;
                        elseif res(ii)==2
                            B1(di+3*m+ii) = 0;
                            B2(di+3*m+ii) = 0;
                        elseif res(ii)==3
                            B1(di+4*m+ii) = 0;
                            B2(di+4*m+ii) = 0;
                        else
                            disp('error');
                            break;
                        end
                    end
                    B1=find(B1==0);
                    B2=find(B2==0);
                    A = AA;
                    B = BB;
                    A(B1,:)=[];
                    A(:,B2)=[];
                    B(B2)=[];
                    if rank(A) ~= length(A) %singular matrix, enter next cycle
                        continue
                    else
                        F=inv(A)*B;
                    end
                    if ~all(F(di+1:end)>=0) %Variable is less than zero, enter the next cycle
                        continue
                    end
                    y = F(1:di+m);%feasible region verification
                    y = [y(1:di+find(xc==xr(2,1))-1);0;y(di+find(xc==xr(2,1)):end)];
                    lam = F(di+m+1:di+2*m);
                    flag = 0;
                    for ii =1:length(xr)
                        j = find(xc==xr(1,ii));
                        r = txr(ii);
                        si = X(xr(2,ii)+1,:)';
                        sj = X(xr(1,ii)+1,:)';
                        c1 = bias^2+2*bias*r+r^2+norm(sj)^2-norm(si)^2;
                        c2 = -bias^2+2*bias*r-r^2-norm(sj)^2+norm(si)^2;
                        u1 = [2*(si-sj);zeros(length(xc),1)];
                        u1(di+j) = 2*(bias+r);
                        u2 = [-2*(si-sj);zeros(length(xc),1)];
                        u2(di+j) = 2*(bias-r);
                        if c1 + u1'*y - lam(ii) > 10^-2
                            flag = 1;
                            break
                        elseif c2 + u2'*y - lam(ii) > 10^-2
                            flag = 1;
                            break
                        elseif a(j)^2-(x0-sj)'*x0 + (x0-sj)'*y(1:di) -...
                                a(j)*y(di+j) > 10^-2
                            flag = 1;
                            break
                        end
                    end
                    if flag == 0 %the result satisfies the feasible region
                       x = F(1:di);
                       return;
                    end
                end
            end
        end
    end
end
end