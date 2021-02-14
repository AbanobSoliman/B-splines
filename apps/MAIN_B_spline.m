%% B-spline Trajectory (Start date: 09/02/2021)
% Algorithm Description 
% The main idea is to optimize RGB-D-T cameras configuration on a
% B-spline continuous time trajectory based on:
% - Rigid body motion kinematics based-on Lie group SE(3)
% Developed by Abanob SOLIMAN, PhD Student, IBISC Laboratory, France
% Email: abanob.soliman@univ-evry.fr
% Under the supervision of:
% Prof. "Samia Bouchafa Bruneau" <samia.bouchafabruneau@univ-evry.fr>
% Prof. "Dro Désiré Sidibie" <drodesire.sidibie@univ-evry.fr>
% Dr. "fabien bonardi" <fabien.bonardi@univ-evry.fr>
%% Initializing Environment
clc;close all;clear all;
addpath '../data/b-spline';
addpath '../src';
%% Loading the Poses T=[R|t] --> Testing data
tic
% Loading EKF data
load('dP.mat','pqpos','pqorient')
% Reading OKVIS Data 
Tbl_OKVIS = readtable('OKVIS_state_estimation.csv');
Tbl_cam = readtable('data_cam.csv');
camTimeSteps = Tbl_cam.x_timestamp_ns_;
camTimeSteps  = (camTimeSteps - camTimeSteps(1))*1e-9;
OKVISPositions    = [Tbl_OKVIS.p_x,Tbl_OKVIS.p_y,Tbl_OKVIS.p_z];
OKVISQuaternions  = quaternion([Tbl_OKVIS.q_w,Tbl_OKVIS.q_x,Tbl_OKVIS.q_y,Tbl_OKVIS.q_z]);
OKVISOrientations = quat2eul(OKVISQuaternions,'ZYX')*180/pi; % In Degrees
OKVISSteps = size(Tbl_OKVIS,1);
% Reading the GroundTruth Data 
Tbl_gt = readtable('data_gt.csv');
gtTimeSteps = Tbl_gt.x_timestamp_ns_;
gtTimeSteps  = (gtTimeSteps - gtTimeSteps(1))*1e-9;
gtPositions    = [Tbl_gt.p_RS_R_x_m_,Tbl_gt.p_RS_R_y_m_,Tbl_gt.p_RS_R_z_m_];
gtQuaternions  = quaternion([Tbl_gt.q_RS_w__,Tbl_gt.q_RS_x__,Tbl_gt.q_RS_y__,Tbl_gt.q_RS_z__]);
gtOrientations = quat2eul(gtQuaternions,'ZYX')*180/pi; % In Degrees
gtSteps = size(Tbl_gt,1);
toc
%% B-spline simulation for R(3) Pose representation -X,Y,Z Translations ONLY-
tic
Poses = OKVISPositions';
CP = size(Poses,2); % number of control points = poses
n = CP - 1;
k = [2,3]; % Degree of b-spline (Quadratic/Cubic)
Pr = 50; %Precision of curve
t = linspace(0,1,Pr);
S2_R3 = b_splineR3(Poses,t,k(1));
T2_R3 = linspace(0,camTimeSteps(CP),size(S2_R3,2));
S3_R3 = b_splineR3(Poses,t,k(2));
T3_R3 = linspace(0,camTimeSteps(CP),size(S3_R3,2));
toc
%% Plotting Quadratic/Cubic B-spline in R(3)
tic
figure
plot3(Poses(1,:),Poses(2,:),Poses(3,:),'DisplayName','OKVIS');
hold all
plot3(S2_R3(1,:),S2_R3(2,:),S2_R3(3,:),'DisplayName','Quad-B-spline');
hold all
plot3(S3_R3(1,:),S3_R3(2,:),S3_R3(3,:),'DisplayName','Cubic-B-spline');
hold all
xlabel('X cm')
ylabel('Y cm')
zlabel('Z cm')
legend
title('X-Y-Z Path B-spline')
grid on

figure
plot(Poses(1,:),Poses(2,:),'DisplayName','OKVIS');
hold all
plot(S2_R3(1,:),S2_R3(2,:),'DisplayName','Quad-B-spline');
hold all
plot(S3_R3(1,:),S3_R3(2,:),'DisplayName','Cubic-B-spline');
hold all
xlabel('X cm')
ylabel('Y cm')
legend
title('X-Y Path B-spline')
grid on
toc
%% Comulative B-spline on R(3), SO(3) & SE(3) Lie group 
tic
Pr = 50; % spline precision
u = linspace(0,1,Pr);

% R(3) Representation -Translations (3x1 vector) ONLY-
P  = OKVISPositions';
n = 3; % spline order
S3_CR3 = comul_b_splineR3(P,u,n);
T3_CR3 = linspace(0,camTimeSteps(size(P,2)),size(S3_CR3,2));
n = 2; % spline order
S2_CR3 = comul_b_splineR3(P,u,n);
T2_CR3 = linspace(0,camTimeSteps(size(P,2)),size(S2_CR3,2));

% SO(3) Representation -Rotations (3x3 matrix) ONLY-
O = OKVISOrientations;
Q = OKVISQuaternions;
n = 3; % spline order
S3_SO3 = comul_b_splineSO3(Q,u,n);
T3_SO3 = linspace(0,camTimeSteps(size(Q,1)),size(S3_SO3,2));
QS3_SO3 = S3_SO3';
ES3_SO3 = quat2eul(QS3_SO3,'ZYX')*180/pi; % deg
n = 2; % spline order
S2_SO3 = comul_b_splineSO3(Q,u,n);
T2_SO3 = linspace(0,camTimeSteps(size(Q,1)),size(S2_SO3,2));
QS2_SO3 = S2_SO3';
ES2_SO3 = quat2eul(QS2_SO3,'ZYX')*180/pi; % deg

% SE(3) Representation -Rotations & Translations (4x4 matrix)-
T = [P;compact(Q)'];
n = 3; % spline order
S3_SE3 = comul_b_splineSE3(T,u,n);
T3_SE3 = linspace(0,camTimeSteps(size(T,2)),size(S3_SE3,2));
QS3_SE3 = S3_SE3(4:end,:)';
ES3_SE3 = quat2eul(QS3_SE3,'ZYX')*180/pi; % deg
n = 2; % spline order
S2_SE3 = comul_b_splineSE3(T,u,n);
T2_SE3 = linspace(0,camTimeSteps(size(T,2)),size(S2_SE3,2));
QS2_SE3 = S2_SE3(4:end,:)';
ES2_SE3 = quat2eul(QS2_SE3,'ZYX')*180/pi; % deg
toc
%% Plotting Comulative Cubic/Quadratic B-spline in R(3)/SO(3)/SE(3)
tic
figure
plot(P(1,:),P(2,:),'DisplayName','OKVIS');
hold all
plot(S2_R3(1,:),S2_R3(2,:),'DisplayName','Quad. B-spline');
hold all
plot(S2_CR3(1,:),S2_CR3(2,:),'DisplayName','Quad. comul R3');
hold all
plot(S2_SE3(1,:),S2_SE3(2,:),'DisplayName','Quad. comul SE3');
hold all
plot(S3_R3(1,:),S3_R3(2,:),'DisplayName','Cubic B-spline');
hold all
plot(S3_CR3(1,:),S3_CR3(2,:),'DisplayName','Cubic comul R3');
hold all
plot(S3_SE3(1,:),S3_SE3(2,:),'DisplayName','Cubic comul SE3');
hold all
xlabel('X cm')
ylabel('Y cm')
legend
title('X-Y Path Comulative B-spline in R3,SE3')
grid on

figure
subplot(3,1,1)
plot(camTimeSteps(1:size(P,2)),P(1,:),'DisplayName','OKVIS');
hold all
plot(T2_R3,S2_R3(1,:),'DisplayName','Quad. B-spline');
hold all
plot(T3_R3,S3_R3(1,:),'DisplayName','Cubic B-spline');
hold all
plot(T2_CR3,S2_CR3(1,:),'DisplayName','Quad. comul-R3');
hold all
plot(T3_CR3,S3_CR3(1,:),'DisplayName','Cubic comul-R3');
hold all
plot(T2_SE3,S2_SE3(1,:),'DisplayName','Quad. comul-SE3');
hold all
plot(T3_SE3,S3_SE3(1,:),'DisplayName','Cubic comul-SE3');
hold all
ylabel('X cm')
legend
title('X-Y-Z-t Path Comulative B-spline in R3,SE3')
grid on
subplot(3,1,2)
plot(camTimeSteps(1:size(P,2)),P(2,:),'DisplayName','OKVIS');
hold all
plot(T2_R3,S2_R3(2,:),'DisplayName','Quad. B-spline');
hold all
plot(T3_R3,S3_R3(2,:),'DisplayName','Cubic B-spline');
hold all
plot(T2_CR3,S2_CR3(2,:),'DisplayName','Quad. comul-R3');
hold all
plot(T3_CR3,S3_CR3(2,:),'DisplayName','Cubic comul-R3');
hold all
plot(T2_SE3,S2_SE3(2,:),'DisplayName','Quad. comul-SE3');
hold all
plot(T3_SE3,S3_SE3(2,:),'DisplayName','Cubic comul-SE3');
hold all
ylabel('Y cm')
legend
grid on
subplot(3,1,3)
plot(camTimeSteps(1:size(P,2)),P(3,:),'DisplayName','OKVIS');
hold all
plot(T2_R3,S2_R3(3,:),'DisplayName','Quad. B-spline');
hold all
plot(T3_R3,S3_R3(3,:),'DisplayName','Cubic B-spline');
hold all
plot(T2_CR3,S2_CR3(3,:),'DisplayName','Quad. comul-R3');
hold all
plot(T3_CR3,S3_CR3(3,:),'DisplayName','Cubic comul-R3');
hold all
plot(T2_SE3,S2_SE3(3,:),'DisplayName','Quad. comul-SE3');
hold all
plot(T3_SE3,S3_SE3(3,:),'DisplayName','Cubic comul-SE3');
hold all
ylabel('Z cm')
xlabel('t sec')
legend
grid on

figure
subplot(3,1,1)
plot(camTimeSteps(1:size(O,1)),O(:,1),'DisplayName','OKVIS');
hold all
plot(T2_SO3,ES2_SO3(:,1),'DisplayName','Quad. comul SO3');
hold all
plot(T3_SO3,ES3_SO3(:,1),'DisplayName','Cubic comul SO3');
hold all
plot(T2_SE3,ES2_SE3(:,1),'DisplayName','Quad. comul SE3');
hold all
plot(T3_SE3,ES3_SE3(:,1),'DisplayName','Cubic comul SE3');
hold all
ylabel('\phi°')
legend
title('\phi-\theta-\psi Rotations Comulative B-spline in SO3,SE3')
grid on
subplot(3,1,2)
plot(camTimeSteps(1:size(O,1)),O(:,2),'DisplayName','OKVIS');
hold all
plot(T2_SO3,ES2_SO3(:,2),'DisplayName','Quad. comul SO3');
hold all
plot(T3_SO3,ES3_SO3(:,2),'DisplayName','Cubic comul SO3');
hold all
plot(T2_SE3,ES2_SE3(:,2),'DisplayName','Quad. comul SE3');
hold all
plot(T3_SE3,ES3_SE3(:,2),'DisplayName','Cubic comul SE3');
hold all
ylabel('\theta°')
legend
grid on
subplot(3,1,3)
plot(camTimeSteps(1:size(O,1)),O(:,3),'DisplayName','OKVIS');
hold all
plot(T2_SO3,ES2_SO3(:,3),'DisplayName','Quad. comul SO3');
hold all
plot(T3_SO3,ES3_SO3(:,3),'DisplayName','Cubic comul SO3');
hold all
plot(T2_SE3,ES2_SE3(:,3),'DisplayName','Quad. comul SE3');
hold all
plot(T3_SE3,ES3_SE3(:,3),'DisplayName','Cubic comul SE3');
hold all
ylabel('\psi°')
xlabel('t sec')
legend
grid on
toc
%% end script (Last update date: 14/02/2021)