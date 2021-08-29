%SDR2
%基于距离法
%多目标优化
%BR方法修改：
%1.约束中加入trace(W)>=3;
%2.传感器恢复中加入选择的传感器数量为trace(W)（向下取整）
function ret=sel_function_br_MOP_Revision(Noise_var,Source_filed_po,Sensor_posi,gamma)
warning off
num_sensor=length(Sensor_posi);
lambda=0.01;
Noise_cvm=diag(Noise_var);
dim=length(Source_filed_po(:,1));
%--------------------------------------------------------------------------
%发送时间未知的TOA模型（TDOA转化）
    for j=1:num_sensor
        Range(j,:)=norm((Source_filed_po-Sensor_posi(:,j)),2);
    end
        RANGE=Range*ones(1,dim);
        H=[((Source_filed_po*ones(1,num_sensor))'-(Sensor_posi)')./RANGE,ones(num_sensor,1)];
Range_a=Range./(sum(Range)/num_sensor);%lambda的N乘在这里了
%凸优化
cvx_clear
cvx_begin sdp
cvx_solver sdpt3
cvx_precision best
cvx_quiet(1)
variable w(num_sensor,1);
variable Z(dim+1,dim+1);
variable W(num_sensor,num_sensor);
minimize trace(Z(1:dim,1:dim))+gamma*sum(w)+lambda*(Range_a'*w)
subject to
[H'*diag(w)*inv(Noise_cvm)*H,eye(dim+1);eye(dim+1),Z]>=0;%线性优于矩阵
[W,w;w',1]>=0;
trace(W)>=dim+1;
diag(W)==w;
cvx_end
w1=w;
w1_sort=sort(w1,'descend');
w1_max=w1_sort(floor(sum(w1+1e-4)));
w1(w1<w1_max)=0;
w1(w1>0)=1;

w2=w;
w2=round(w2);

%计算CRLB
Fisher_information_matrix_inverse=inv(H'*diag(w1)/(Noise_cvm)*diag(w1)*H);
crlb1=trace(Fisher_information_matrix_inverse(1:dim,1:dim));

Fisher_information_matrix_inverse=inv(H'*diag(w2)/(Noise_cvm)*diag(w2)*H);
crlb2=trace(Fisher_information_matrix_inverse(1:dim,1:dim));

if (abs(crlb1)+gamma*sum(w1))<=(abs(crlb2)+gamma*sum(w2))
ret=[w1;crlb1];
else
ret=[w2;crlb2];
end

