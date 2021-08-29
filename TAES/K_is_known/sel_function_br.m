%SDR2
%基于距离法
clc;clear;close all;
warning off
num_sensor=20;
sensor_sel=5;
dim=2;
Source_filed_po=unifrnd(0,num_sensor*2.5,dim,1);
Sensor_posi=unifrnd(0,num_sensor*2.5,dim,num_sensor);
Noise_cvm=eye(num_sensor);
dim=length(Source_filed_po(:,1));
lambda=0.01;
%--------------------------------------------------------------------------
%发送时间未知的TOA模型（TDOA转化）
    for j=1:num_sensor
        Range(j,:)=norm((Source_filed_po-Sensor_posi(:,j)),2);
    end
        RANGE=Range*ones(1,dim);
        H=[((Source_filed_po*ones(1,num_sensor))'-(Sensor_posi)')./RANGE,ones(num_sensor,1)];
Range_a=Range./(sum(Range)/num_sensor);
%凸优化
cvx_clear
cvx_begin sdp
cvx_solver sdpt3
cvx_precision best
cvx_quiet(1)
variable w(num_sensor,1);
variable Z(dim+1,dim+1);
variable W(num_sensor,num_sensor);
minimize trace(Z(1:dim,1:dim))+lambda*(Range_a'*w)
subject to
[H'*diag(w)*inv(Noise_cvm)*H,eye(dim+1);eye(dim+1),Z]>=0;%线性优于矩阵
[W,w;w',1]>=0;
trace(W)==sensor_sel;
diag(W)==w;
cvx_end
%传感器选择向量恢复
w_sort=sort(w,'descend');
w_max=w_sort(sensor_sel); 
w(w<w_max)=0;
w(w>0)=1;
%计算CRLB
Fisher_information_matrix_inverse=inv(H'*diag(w)/(Noise_cvm)*diag(w)*H);
crlb=trace(Fisher_information_matrix_inverse(1:dim,1:dim));



