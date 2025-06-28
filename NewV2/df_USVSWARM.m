function dz=df_USVSWARM(t,Data)
gain_gamma=50;
gain_alpha=50;
gain_beta=130;
gain_kappa=190;

global N1 N2 L_1 L_2 L_w O;

% Parameters
X_min=-1000; X_max=1000;
Y_min=-1000; Y_max=1000;
phi_min=0;  phi_max=2*pi;
eta_max=[X_max;Y_max;phi_max];
b1=kron(ones(6,1),eta_max);
eta_min=[X_min;Y_min;phi_min];
delta1=kron(ones(6,1),eta_min);

b2=b1;
delta2=delta1;

h_1=[500 500 300 300 400 500]'+ones(N1,1)*X_min;
h_2=[500 500 500 300 300 400]'-ones(N2,1)*X_max;

global R_1 R_2 G_1 G_2;

%% data introduce
eta1=Data(1:3*N1);
vartheta1=Data(3*N1+1:6*N1);
omega1=Data(6*N1+1:9*N1);
varpi1=Data(9*N1+1:12*N1);

lambda1=Data(12*N1+1:13*N1);
rho1=Data(13*N1+1:14*N1);

xi1=Data(14*N1+1:14*N1+N1^2*3);
zeta1=Data(14*N1+N1^2*3+1:14*N1+2*N1^2*3);

Data_Num=14*N1+2*N1^2*3;

eta2=Data(Data_Num+1:Data_Num+3*N2);
vartheta2=Data(Data_Num+3*N2+1:Data_Num+6*N2);
omega2=Data(Data_Num+6*N2+1:Data_Num+9*N2);
varpi2=Data(Data_Num+9*N2+1:Data_Num+12*N2);

lambda2=Data(Data_Num+12*N2+1:Data_Num+13*N2);
rho2=Data(Data_Num+13*N2+1:Data_Num+14*N2);

xi2=Data(Data_Num+14*N2+1:Data_Num+14*N2+N2^2*3);
zeta2=Data(Data_Num+14*N2+N2^2*3+1:Data_Num+14*N2+2*N2^2*3);

S=Data(Data_Num+14*N2+2*N2^2*3+1:Data_Num+14*N2+2*N2^2*3+432);

% Gradient
S1=S(1:216);
Ext_grd=Grad_J1(S1(1:36));    dJ_11_w=Ext_grd(1:18,1); clear Ext_grd;
Ext_grd=Grad_J1(S1(37:72));   dJ_12_w=Ext_grd(19:36,1); clear Ext_grd;
Ext_grd=Grad_J1(S1(73:108));   dJ_13_w=Ext_grd(37:54,1); clear Ext_grd;
Ext_grd=Grad_J1(S1(109:144));   dJ_14_w=Ext_grd(55:72,1); clear Ext_grd;
Ext_grd=Grad_J1(S1(145:180));   dJ_15_w=Ext_grd(73:90,1); clear Ext_grd;
Ext_grd=Grad_J1(S1(181:216));   dJ_16_w=Ext_grd(91:108,1); clear Ext_grd;
dJ1=[dJ_11_w;dJ_12_w;dJ_13_w;dJ_14_w;dJ_15_w;dJ_16_w];

S2=S(217:432);
Ext_grd=Grad_J2(S2(1:36));    dJ_21_w=Ext_grd(1:18,1); clear Ext_grd;
Ext_grd=Grad_J2(S2(37:72));   dJ_22_w=Ext_grd(19:36,1); clear Ext_grd;
Ext_grd=Grad_J2(S2(73:108));   dJ_23_w=Ext_grd(37:54,1); clear Ext_grd;
Ext_grd=Grad_J2(S2(109:144));   dJ_24_w=Ext_grd(55:72,1); clear Ext_grd;
Ext_grd=Grad_J2(S2(145:180));   dJ_25_w=Ext_grd(73:90,1); clear Ext_grd;
Ext_grd=Grad_J2(S2(181:216));   dJ_26_w=Ext_grd(91:108,1); clear Ext_grd;
dJ2=[dJ_21_w;dJ_22_w;dJ_23_w;dJ_24_w;dJ_25_w;dJ_26_w];

% Updation
df_eta1=vartheta1;
df_vartheta1=-gain_alpha*vartheta1-(R_1*xi1+G_1'*proj(lambda1)+proj(omega1+eta1-b1)-proj(varpi1-eta1+delta1));
df_omega1=-omega1+proj(omega1+eta1-b1)+vartheta1;
df_varpi1=-varpi1+proj(varpi1-eta1+delta1)+vartheta1;
df_lambda1=-lambda1+proj(lambda1)+G_1*eta1-h_1-kron(L_1,eye(1))*proj(lambda1)-kron(L_1,eye(1))*rho1;
df_rho1=kron(L_1,eye(1))*proj(lambda1);
d_xi1=-gain_beta*( xi1+kron(L_1,eye(18))*xi1+kron(L_1,eye(18))*zeta1-dJ1 );
d_zeta1=gain_beta*( kron(L_1,eye(18))*xi1 );

df_eta2=vartheta2;
df_vartheta2=-gain_alpha*vartheta2-(R_2*xi2+G_2'*proj(lambda2)+proj(omega2+eta2-b2)-proj(varpi2-eta2+delta2));
df_omega2=-omega2+proj(omega2+eta2-b2)+vartheta2;
df_varpi2=-varpi2+proj(varpi2-eta2+delta2)+vartheta2;
df_lambda2=-lambda2+proj(lambda2)+G_2*eta2-h_2-kron(L_2,eye(1))*proj(lambda2)-kron(L_2,eye(1))*rho2;
df_rho2=kron(L_2,eye(1))*proj(lambda2);
d_xi2=-gain_beta*( xi2+kron(L_2,eye(18))*xi2+kron(L_2,eye(18))*zeta2-dJ2 );
d_zeta2=gain_beta*( kron(L_2,eye(18))*xi2 );

eta=[eta1;eta2];
d_S=-gain_kappa*( (kron(L_w,eye(3*12))+O)*S -O*kron(ones(12,1),eta) );

%% USV
vartheta=[varpi1;varpi2];
df_vartheta=[df_vartheta1;df_vartheta2];
df_eta=[df_eta1;df_eta2];
m_usv=23.8;
xg_usv=0.046;
iz_usv=1.76;
Xu_usv=-0.723;
Yv_usv=-0.861;
Yr_usv=0.108;
Nv_usv=0.105;
Nr_usv=-1.9;
Xdu_usv=-2;
Ydv_usv=-10;
Ydr_usv=0;
Ndv_usv=0;
Ndr_usv=-1;

NUM_USV=Data_Num+14*N2+2*N2^2*3+432;
X_usv=Data(NUM_USV+1:NUM_USV+36);
PsiV=Data(NUM_USV+37:NUM_USV+72);
Hat_mu=Data(NUM_USV+73:NUM_USV+192);
hat_d=Data(NUM_USV+193:NUM_USV+204);
%dz2=[df_X_usv,dfdf_X_usv,df_Hat_mu,df_hat_x];
tau=zeros(36,1);
dfdf_X_usv=zeros(36,1);
df_hat_d=zeros(12,1);

df_X_usv=PsiV;
df_hat_x=vartheta-(X_usv-eta);
difdif_hat_x=df_vartheta-(PsiV-df_eta);

e=df_X_usv-df_hat_x;
N=12;
for i=1:N
    dot_x_i=PsiV((i-1)*3+1:i*3);
    psi_usv=X_usv(i*3);
    Psi=[cos(psi_usv)    -sin(psi_usv)    0;
        sin(psi_usv)     cos(psi_usv)    0;
        0             0         1];
    vi=Psi'*dot_x_i;
    tau_varrho=vi(1);
    tau_delta=vi(2);
    tau_theta=vi(3);

    mathcal_H=[0  -tau_theta   0;
        tau_theta   0   0;
        0   0   0];

    M=[m_usv-Xdu_usv          0         0;
        0            m_usv-Ydv_usv      m_usv*xg_usv-Ydr_usv;
        0            m_usv*xg_usv-Ndv_usv   iz_usv-Ndr_usv];

    c13_usv=-(m_usv-Ydv_usv)*tau_delta-(m_usv*xg_usv-Ydr_usv)*tau_theta;
    c23_usv=(m_usv-Xdu_usv)*tau_varrho;
    Pi_usv=[0     0   c13_usv;
        0     0   c23_usv;
        -c13_usv  -c23_usv  0];

    Phi_usv=[-Xu_usv   0    0;
        0    -Yv_usv   -Yr_usv;
        0    -Nv_usv   -Nr_usv];

    E_usv(:,:,i)=Psi*M*Psi';
    C_usv(:,:,i)=Psi*(Pi_usv-M*mathcal_H)*Psi';
    D_usv(:,:,i)=Psi*Phi_usv*Psi';
    R_usv(:,:,i)=Psi;


    y=difdif_hat_x((i-1)*3+1:i*3);
    z=df_hat_x((i-1)*3+1:i*3);
    Upsilon_usv(:,:,i)=[-cos(psi_usv)^2*z(1) - cos(psi_usv)*sin(psi_usv)*z(2), cos(psi_usv)*sin(psi_usv)*z(2) - sin(psi_usv)^2*z(1) ...
        sin(psi_usv)*z(3), 0 ...
        0, -cos(psi_usv)^2*tau_theta*z(2) + cos(psi_usv)*tau_theta*sin(psi_usv)*z(1) - cos(psi_usv)^2*y(1) - cos(psi_usv)*sin(psi_usv)*y(2) + sin(psi_usv)*tau_varrho*z(3)...
        -cos(psi_usv)*tau_theta*sin(psi_usv)*z(1) - sin(psi_usv)^2*tau_theta*z(2) + cos(psi_usv)*sin(psi_usv)*y(2) + cos(psi_usv)*tau_delta*z(3) - sin(psi_usv)^2*y(1), cos(psi_usv)*tau_theta*z(3) + sin(psi_usv)*y(3) ...
        0, 0;
        -cos(psi_usv)*sin(psi_usv)*z(1) - sin(psi_usv)^2*z(2), -cos(psi_usv)^2*z(2) + cos(psi_usv)*sin(psi_usv)*z(1) ...
        -cos(psi_usv)*z(3), 0 ...
        0, -cos(psi_usv)*tau_theta*sin(psi_usv)*z(2) + sin(psi_usv)^2*tau_theta*z(1) - cos(psi_usv)*sin(psi_usv)*y(1) - cos(psi_usv)*tau_varrho*z(3) - sin(psi_usv)^2*y(2) ...
        cos(psi_usv)^2*tau_theta*z(1) + cos(psi_usv)*tau_theta*sin(psi_usv)*z(2) - cos(psi_usv)^2*y(2) + cos(psi_usv)*sin(psi_usv)*y(1) + sin(psi_usv)*tau_delta*z(3), sin(psi_usv)*tau_theta*z(3) - cos(psi_usv)*y(3) ...
        0, 0;
        0, 0 ...
        0, -cos(psi_usv)*z(2) + sin(psi_usv)*z(1), -z(3), cos(psi_usv)*tau_varrho*z(2) - sin(psi_usv)*tau_varrho*z(1) ...
        -cos(psi_usv)*tau_delta*z(1) - sin(psi_usv)*tau_delta*z(2), -cos(psi_usv)*tau_theta*z(1) - sin(psi_usv)*tau_theta*z(2) ...
        cos(psi_usv)*tau_theta*z(1) + sin(psi_usv)*tau_theta*z(2) - cos(psi_usv)*y(2) + sin(psi_usv)*y(1), -y(3)];

    distb(1)=sin(2*t);
    distb(2)=cos(3*t);
    distb(3)=sin(0.1*t);

    tau((i-1)*3+1:i*3,end)=-gain_gamma*e((i-1)*3+1:i*3,end)+Upsilon_usv(:,:,i)*Hat_mu((i-1)*10+1:i*10,end);
    dfdf_X_usv((i-1)*3+1:i*3,end)=(E_usv(:,:,i))\(-(C_usv(:,:,i)+D_usv(:,:,i))*dot_x_i+tau((i-1)*3+1:i*3,end) +hat_d(i)*sign_vector(e((i-1)*3+1:i*3,end)) +Psi*distb');

    df_hat_d(i,end)=e((i-1)*3+1:i*3,end)'*sign_vector(e((i-1)*3+1:i*3,end));

end
Upsilon=blkdiag(Upsilon_usv(:,:,1),Upsilon_usv(:,:,2),Upsilon_usv(:,:,3),Upsilon_usv(:,:,4),Upsilon_usv(:,:,5),Upsilon_usv(:,:,6), ...
    Upsilon_usv(:,:,7),Upsilon_usv(:,:,8),Upsilon_usv(:,:,9),Upsilon_usv(:,:,10),Upsilon_usv(:,:,11),Upsilon_usv(:,:,12));
% MM=blkdiag(E_usv(:,:,1),E_usv(:,:,2),E_usv(:,:,3),E_usv(:,:,4),E_usv(:,:,5),E_usv(:,:,6),...
%     E_usv(:,:,7),E_usv(:,:,8),E_usv(:,:,9),E_usv(:,:,10),E_usv(:,:,11),E_usv(:,:,12));
% CC=blkdiag(C_usv(:,:,1),C_usv(:,:,2),C_usv(:,:,3),C_usv(:,:,4),C_usv(:,:,5),C_usv(:,:,6), ...
%     C_usv(:,:,7),C_usv(:,:,8),C_usv(:,:,9),C_usv(:,:,10),C_usv(:,:,11),C_usv(:,:,12));
% DD=blkdiag(D_usv(:,:,1),D_usv(:,:,2),D_usv(:,:,3),D_usv(:,:,4),D_usv(:,:,5),D_usv(:,:,6), ...
%     D_usv(:,:,7),D_usv(:,:,8),D_usv(:,:,9),D_usv(:,:,10),D_usv(:,:,11),D_usv(:,:,12));
% RR=blkdiag(R_usv(:,:,1),R_usv(:,:,2),R_usv(:,:,3),R_usv(:,:,4),R_usv(:,:,5),R_usv(:,:,6), ...
%     R_usv(:,:,7),R_usv(:,:,8),R_usv(:,:,9),R_usv(:,:,10),R_usv(:,:,11),R_usv(:,:,12));

% randong=zeros(N*3,1);
% for i=1:N
%     randong(3*i-2)=30*sin(i*t);
%     randong(3*i-1)=30*cos(i*t);
%     randong(3*i)=30*sin(i*t);
% end

df_Hat_mu=-Upsilon'*e;
df_hat_x=vartheta-(X_usv-eta);

%dz2=[df_X_usv,dfdf_X_usv,df_Hat_mu,df_hat_x];

%% different
dz=[df_eta1;df_vartheta1;df_omega1;df_varpi1;df_lambda1;df_rho1;d_xi1;d_zeta1; ...
    df_eta2;df_vartheta2;df_omega2;df_varpi2;df_lambda2;df_rho2;d_xi2;d_zeta2;d_S;df_X_usv; dfdf_X_usv; df_Hat_mu; df_hat_x; df_hat_d];

end