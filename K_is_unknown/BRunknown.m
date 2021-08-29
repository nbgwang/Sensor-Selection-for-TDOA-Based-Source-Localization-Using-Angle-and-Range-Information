%传感器选择与目标定位
clc;clear;close all;
warning off
sensor_number=20;
dim=2;
Source_filed_po=unifrnd(0,sensor_number*2.5,dim,1);
SENSOR_coordinate=unifrnd(0,sensor_number*2.5,dim,sensor_number);
Noise_variance=ones(sensor_number,1);
dim=length(Source_filed_po(:,1));
gamma=1;
%SDR-TOA+BR方法
%--------------------------------------------------------------------------
%传感器选择
%--------------------------------------------------------------------------
%根据上一时刻位置选出的传感器
S=sel_function_br_MOP_Revision(Noise_variance,Source_filed_po,SENSOR_coordinate,gamma);
CRB=S(sensor_number+1);
NUM_sensor=sum(S(1:sensor_number));

