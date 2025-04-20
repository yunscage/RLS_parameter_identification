function [Uoc,R0,R1,Tau,alpha] = RLS_Para(T,Ut,Io,NoFF)
%Parameter 识别一阶电池模型
%  Uo Output voatage,Io Output current
% H=[1 Ut(i) -Io(i+1) -Io(i)];
% thata=[(1-a)*Uoc ; a ; R0 ; R1(1-a)-a*R0];
% NoFF=1; % 是否使用遗忘因子 1: 使用;   0: 不使用
if nargin<4
    NoFF=0;
end

r=0.993;%遗忘因子
w=1/r;
L1=length(Ut);
L2=length(Io);
if L1 ~= L2 || L1<10
    R0=-1;
    R1=-1;
    Tau=-1;
    disp('error!');
    return;
end
Ut=Ut(:);
Io=Io(:);



if NoFF==1
    HK=min(50,L2-1);
else
    HK=L2-1;
end
H=[ones(HK,1) Ut(1:HK) -Io(2:HK+1) -Io(1:HK)];
Z=Ut(2:HK+1);
% 使用 SVD 来计算等价于 (H'*H)\(h'*w*e)
M=H'*H;
[U, S, V] = svd(M);% 计算 M 的 SVD 分解
S_pinv = diag(1./diag(S));% 计算 S 的伪逆
S_pinv(S == 0) = 0; % 确保奇异值不为零
Theta=V * S_pinv * U'*H'*Z;
w=1;
if  NoFF==1
estimate=zeros(length(Theta),L1);
Estore=zeros(1,L1);
M=1000*eye(length(Theta));
estimate(:,1)=Theta;
for i=1:L1-1
    h=[1 Ut(i) -Io(i+1) -Io(i)];
    e=Ut(i+1)-h*Theta;
    M=M+h'*w*h;
    [U, S, V] = svd(M);% 计算 M 的 SVD 分解
    S_pinv = diag(1./diag(S));% 计算 S 的伪逆
    S_pinv(S == 0) = 0; % 确保奇异值不为零
    % 使用 SVD 来计算等价于 M \ (h' * w * e)
    Theta=Theta + V * S_pinv * U' * (h' * w * e);
    % Theta=Theta+M\(h'*w*e);
    estimate(:,i+1)=Theta;
    Estore(:,i)=e;
    % 自适应调整w
    if i>HK
        e_temp =Estore(i-9:i);
        variance=mean(e_temp.*e_temp)/20;
        r=0.95+0.04*(exp(-variance));
        w=1/r;
    end
end
Estore(:,i+1)=e;
else
    estimate=Theta;
    Estore=Z-H*Theta;
end
alpha=estimate(2,:);

%% 方案求解1
% Uoc = estimate(1,:) ./ (1 - alpha);
% R0=estimate(3,:);
% R1 = (estimate(4,:) + alpha .* R0) ./ (1 - alpha);
% Tau=-T./log(alpha);



%% 双线性变换方案求解
Uoc = estimate(1,:) ./ (1 - estimate(2,:));
R0=(estimate(3,:)-estimate(4,:))./(1+estimate(2,:));
Tau=T/2*(1+estimate(2,:))./(1-estimate(2,:));
R1 =2*(estimate(3,:)-R0)./(1-estimate(2,:));

alpha = max(0.001, min(alpha, 0.999));
Uoc=max(0.8*min(Ut),min(Uoc,1.2*max(Ut)));
R0 = max(1e-5, min(R0, 0.5));
R1 = max(3e-5, min(R1, 1));
Tau = max(0.12, min(Tau, 500));
end
