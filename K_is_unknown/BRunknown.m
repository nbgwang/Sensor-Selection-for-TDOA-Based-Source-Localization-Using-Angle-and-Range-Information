%������ѡ����Ŀ�궨λ
clc;clear;close all;
warning off
sensor_number=20;
dim=2;
Source_filed_po=unifrnd(0,sensor_number*2.5,dim,1);
SENSOR_coordinate=unifrnd(0,sensor_number*2.5,dim,sensor_number);
Noise_variance=ones(sensor_number,1);
dim=length(Source_filed_po(:,1));
gamma=1;
%SDR-TOA+BR����
%--------------------------------------------------------------------------
%������ѡ��
%--------------------------------------------------------------------------
%������һʱ��λ��ѡ���Ĵ�����
S=sel_function_br_MOP_Revision(Noise_variance,Source_filed_po,SENSOR_coordinate,gamma);
CRB=S(sensor_number+1);
NUM_sensor=sum(S(1:sensor_number));

