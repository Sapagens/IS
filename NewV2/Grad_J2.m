function Grad=Grad_J2(Z)
Z1=Z(1:18);
Z2=Z(19:36);
for i=1:6
    X_1(i)=Z1(3*i-2);X_2(i)=Z2(3*i-2);
    Y_1(i)=Z1(3*i-1);Y_2(i)=Z2(3*i-1);
    phi_1(i)=Z1(3*i);phi_2(i)=Z2(3*i);
end
% parameter
X_2TRG=1000;Y_2TRG=0;
X_1TRG=-1000;Y_1TRG=0;
iota2_11=0.01; iota2_12=5;   iota2_13=1;   iota2_14=10;
iota2_21=0.01; iota2_22=5;   iota2_23=1;   iota2_24=10;
iota2_31=0.01; iota2_32=5;   iota2_33=1;   iota2_34=10;   
iota2_41=2; iota2_42=1;   iota2_43=30;   iota2_44=1;   iota2_45=3;
iota2_51=2; iota2_52=1;   iota2_53=30;   iota2_54=1;   iota2_55=3;
iota2_61=2; iota2_62=1;   iota2_63=30;   iota2_64=1;   iota2_65=3;
% Formation
DeltaX_21=0; DeltaX_22=0; DeltaX_23=0; DeltaX_24=0; DeltaX_25=0; DeltaX_26=0;
DeltaY_21=100; DeltaY_22=0; DeltaY_23=-100; DeltaY_24=100; DeltaY_25=0; DeltaY_26=-100;
% DeltaY_21=0; DeltaY_22=0; DeltaY_23=0; DeltaY_24=0; DeltaY_25=0; DeltaY_26=0;
% DeltaX_21=100; DeltaX_22=0; DeltaX_23=-100; DeltaX_24=100; DeltaX_25=0; DeltaX_26=-100;
%% Gradinet_X
X=[X_1',X_2']'; %横着排列的，第i行表示第i个联盟

X_1sum=sum(X_1);
X_2ivsum=X_2(1)+X_2(2)+X_2(3);
X_2dfsum=X_2(4)+X_2(5)+X_2(6);

dJ21_X21=-iota2_11*(6*X(2,1)-X_1sum) +iota2_12*(X(2,1)-1/3*X_2ivsum-DeltaX_21)*(1-1/3) +iota2_14*(X(2,1)-X_1TRG);
dJ21_X22=iota2_12*(X(2,1)-1/3*X_2ivsum-DeltaX_21)*(-1/3);
dJ21_X23=iota2_12*(X(2,1)-1/3*X_2ivsum-DeltaX_21)*(-1/3);
dJ21_X24=0; dJ21_X25=0; dJ21_X26=0;

dJ22_X21=iota2_22*(X(2,2)-1/3*X_2ivsum-DeltaX_22)*(-1/3);
dJ22_X22=-iota2_21*(6*X(2,2)-X_1sum) +iota2_22*(X(2,2)-1/3*X_2ivsum-DeltaX_22)*(1-1/3) +iota2_24*(X(2,2)-X_1TRG);
dJ22_X23=iota2_22*(X(2,2)-1/3*X_2ivsum-DeltaX_22)*(-1/3);
dJ22_X24=0; dJ22_X25=0; dJ22_X26=0;

dJ23_X21=iota2_32*(X(2,3)-1/3*X_2ivsum-DeltaX_23)*(-1/3);
dJ23_X22=iota2_32*(X(2,3)-1/3*X_2ivsum-DeltaX_23)*(-1/3);
dJ23_X23=-iota2_31*(6*X(2,3)-X_1sum) +iota2_32*(X(2,3)-1/3*X_2ivsum-DeltaX_23)*(1-1/3) +iota2_34*(X(2,3)-X_1TRG);
dJ23_X24=0; dJ23_X25=0; dJ23_X26=0;

dJ24_X21=0; dJ24_X22=0; dJ24_X23=0;
dJ24_X24=iota2_41*(X(2,4)-1/3*X_2dfsum-DeltaX_24)*(1-1/3) +iota2_42*(1/3*X_2dfsum-X_2TRG)*(1/3)    +iota2_43*(X(2,4)-1/2*(X(1,1)+X_2TRG))   +iota2_44*(X(2,4)-X(1,1));
dJ24_X25=iota2_41*(X(2,4)-1/3*X_2dfsum-DeltaX_24)*(-1/3)  +iota2_42*(1/3*X_2dfsum-X_2TRG)*(1/3);
dJ24_X26=iota2_42*(X(2,4)-1/3*X_2dfsum-DeltaX_24)*(-1/3)  +iota2_42*(1/3*X_2dfsum-X_2TRG)*(1/3);

dJ25_X21=0; dJ25_X22=0; dJ25_X23=0;
dJ25_X24=iota2_51*(X(2,5)-1/3*X_2dfsum-DeltaX_25)*(-1/3)  +iota2_52*(1/3*X_2dfsum-X_2TRG)*(1/3);
dJ25_X25=iota2_51*(X(2,5)-1/3*X_2dfsum-DeltaX_25)*(1-1/3) +iota2_52*(1/3*X_2dfsum-X_2TRG)*(1/3)    +iota2_53*(X(2,5)-1/2*(X(1,1)+X_2TRG))   +iota2_54*(X(2,5)-X(1,1));
dJ25_X26=iota2_51*(X(2,5)-1/3*X_2dfsum-DeltaX_25)*(-1/3)  +iota2_52*(1/3*X_2dfsum-X_2TRG)*(1/3);

dJ26_X21=0; dJ26_X22=0; dJ26_X23=0;
dJ26_X24=iota2_61*(X(2,6)-1/3*X_2dfsum-DeltaX_26)*(-1/3)  +iota2_62*(1/3*X_2dfsum-X_2TRG)*(1/3);
dJ26_X25=iota2_61*(X(2,6)-1/3*X_2dfsum-DeltaX_26)*(-1/3)  +iota2_62*(1/3*X_2dfsum-X_2TRG)*(1/3);
dJ26_X26=iota2_61*(X(2,6)-1/3*X_2dfsum-DeltaX_26)*(1-1/3) +iota2_62*(1/3*X_2dfsum-X_2TRG)*(1/3)    +iota2_63*(X(2,6)-1/2*(X(1,2)+X_2TRG))   +iota2_64*(X(2,6)-X(1,2));


%% Gradinet_Y
Y=[Y_1',Y_2']'; %横着排列的，第i行表示第i个联盟

Y_1sum=sum(Y_1);
Y_2ivsum=Y_2(1)+Y_2(2)+Y_2(3);
Y_2dfsum=Y_2(4)+Y_2(5)+Y_2(6);

dJ21_Y21=-iota2_11*(6*Y(2,1)-Y_1sum) +iota2_12*(Y(2,1)-1/3*Y_2ivsum-DeltaY_21)*(1-1/3) +iota2_14*(Y(2,1)-Y_1TRG);
dJ21_Y22=iota2_12*(Y(2,1)-1/3*Y_2ivsum-DeltaY_21)*(-1/3);
dJ21_Y23=iota2_12*(Y(2,1)-1/3*Y_2ivsum-DeltaY_21)*(-1/3);
dJ21_Y24=0; dJ21_Y25=0; dJ21_Y26=0;

dJ22_Y21=iota2_22*(Y(2,2)-1/3*Y_2ivsum-DeltaY_22)*(-1/3);
dJ22_Y22=-iota2_21*(6*Y(2,2)-Y_1sum) +iota2_22*(Y(2,2)-1/3*Y_2ivsum-DeltaY_22)*(1-1/3) +iota2_24*(Y(2,2)-Y_1TRG);
dJ22_Y23=iota2_22*(Y(2,2)-1/3*Y_2ivsum-DeltaY_22)*(-1/3);
dJ22_Y24=0; dJ22_Y25=0; dJ22_Y26=0;

dJ23_Y21=iota2_32*(Y(2,3)-1/3*Y_2ivsum-DeltaY_23)*(-1/3);
dJ23_Y22=iota2_32*(Y(2,3)-1/3*Y_2ivsum-DeltaY_23)*(-1/3);
dJ23_Y23=-iota2_31*(6*Y(2,3)-Y_1sum) +iota2_32*(Y(2,3)-1/3*Y_2ivsum-DeltaY_23)*(1-1/3) +iota2_34*(Y(2,3)-Y_1TRG);
dJ23_Y24=0; dJ23_Y25=0; dJ23_Y26=0;

dJ24_Y21=0; dJ24_Y22=0; dJ24_Y23=0;
dJ24_Y24=iota2_41*(Y(2,4)-1/3*Y_2dfsum-DeltaY_24)*(1-1/3) +iota2_42*(1/3*Y_2dfsum-Y_2TRG)*(1/3)    +iota2_43*(Y(2,4)-1/2*(Y(1,1)+Y_2TRG))   +iota2_44*(Y(2,4)-Y(1,1));
dJ24_Y25=iota2_41*(Y(2,4)-1/3*Y_2dfsum-DeltaY_24)*(-1/3)  +iota2_42*(1/3*Y_2dfsum-Y_2TRG)*(1/3);
dJ24_Y26=iota2_42*(Y(2,4)-1/3*Y_2dfsum-DeltaY_24)*(-1/3)  +iota2_42*(1/3*Y_2dfsum-Y_2TRG)*(1/3);

dJ25_Y21=0; dJ25_Y22=0; dJ25_Y23=0;
dJ25_Y24=iota2_51*(Y(2,5)-1/3*Y_2dfsum-DeltaY_25)*(-1/3)  +iota2_52*(1/3*Y_2dfsum-Y_2TRG)*(1/3);
dJ25_Y25=iota2_51*(Y(2,5)-1/3*Y_2dfsum-DeltaY_25)*(1-1/3) +iota2_52*(1/3*Y_2dfsum-Y_2TRG)*(1/3)    +iota2_53*(Y(2,5)-1/2*(Y(1,1)+Y_2TRG))   +iota2_54*(Y(2,5)-Y(1,1));
dJ25_Y26=iota2_51*(Y(2,5)-1/3*Y_2dfsum-DeltaY_25)*(-1/3)  +iota2_52*(1/3*Y_2dfsum-Y_2TRG)*(1/3);

dJ26_Y21=0; dJ26_Y22=0; dJ26_Y23=0;
dJ26_Y24=iota2_61*(Y(2,6)-1/3*Y_2dfsum-DeltaY_26)*(-1/3)  +iota2_62*(1/3*Y_2dfsum-Y_2TRG)*(1/3);
dJ26_Y25=iota2_61*(Y(2,6)-1/3*Y_2dfsum-DeltaY_26)*(-1/3)  +iota2_62*(1/3*Y_2dfsum-Y_2TRG)*(1/3);
dJ26_Y26=iota2_61*(Y(2,6)-1/3*Y_2dfsum-DeltaY_26)*(1-1/3) +iota2_62*(1/3*Y_2dfsum-Y_2TRG)*(1/3)    +iota2_63*(Y(2,6)-1/2*(Y(1,2)+Y_2TRG))   +iota2_64*(Y(2,6)-Y(1,2));


%% Gradinet_phi
phi=[phi_1',phi_2']'; %横着排列的，第i行表示第i个联盟

dJ21_phi21=iota2_13*(phi(2,1)-5*pi/4);dJ21_phi22=0;dJ21_phi23=0;dJ21_phi24=0;dJ21_phi25=0;dJ21_phi26=0;
dJ22_phi21=0;dJ22_phi22=iota2_23*(phi(2,2)-3*pi/2);dJ22_phi23=0;dJ22_phi24=0;dJ22_phi25=0;dJ22_phi26=0;
dJ23_phi21=0;dJ23_phi22=0;dJ23_phi23=iota2_33*(phi(2,3)-7*pi/4);dJ23_phi24=0;dJ23_phi25=0;dJ23_phi26=0;
dJ24_phi21=0;dJ24_phi22=0;dJ24_phi23=0;dJ24_phi24=iota2_45*(phi(2,4)-phi(1,1)-pi);dJ24_phi25=0;dJ24_phi26=0;
dJ25_phi21=0;dJ25_phi22=0;dJ25_phi23=0;dJ25_phi24=0;dJ25_phi25=iota2_55*(phi(2,5)-phi(1,1)-pi);dJ25_phi26=0;
dJ26_phi21=0;dJ26_phi22=0;dJ26_phi23=0;dJ26_phi24=0;dJ26_phi25=0;dJ26_phi26=iota2_65*(phi(2,6)-phi(1,2)-pi);

%% Total_Gradient
Grad_1=[dJ21_X21,dJ21_Y21,dJ21_phi21, ...
dJ21_X22,dJ21_Y22,dJ21_phi22, ...
dJ21_X23,dJ21_Y23,dJ21_phi23, ...
dJ21_X24,dJ21_Y24,dJ21_phi24, ...
dJ21_X25,dJ21_Y25,dJ21_phi25, ...
dJ21_X26,dJ21_Y26,dJ21_phi26]';
Grad_2=[dJ22_X21,dJ22_Y21,dJ22_phi21, ...
dJ22_X22,dJ22_Y22,dJ22_phi22, ...
dJ22_X23,dJ22_Y23,dJ22_phi23, ...
dJ22_X24,dJ22_Y24,dJ22_phi24, ...
dJ22_X25,dJ22_Y25,dJ22_phi25, ...
dJ22_X26,dJ22_Y26,dJ22_phi26]';
Grad_3=[dJ23_X21,dJ23_Y21,dJ23_phi21, ...
dJ23_X22,dJ23_Y22,dJ23_phi22, ...
dJ23_X23,dJ23_Y23,dJ23_phi23, ...
dJ23_X24,dJ23_Y24,dJ23_phi24, ...
dJ23_X25,dJ23_Y25,dJ23_phi25, ...
dJ23_X26,dJ23_Y26,dJ23_phi26]';
Grad_4=[dJ24_X21,dJ24_Y21,dJ24_phi21, ...
dJ24_X22,dJ24_Y22,dJ24_phi22, ...
dJ24_X23,dJ24_Y23,dJ24_phi23, ...
dJ24_X24,dJ24_Y24,dJ24_phi24, ...
dJ24_X25,dJ24_Y25,dJ24_phi25, ...
dJ24_X26,dJ24_Y26,dJ24_phi26]';
Grad_5=[dJ25_X21,dJ25_Y21,dJ25_phi21, ...
dJ25_X22,dJ25_Y22,dJ25_phi22, ...
dJ25_X23,dJ25_Y23,dJ25_phi23, ...
dJ25_X24,dJ25_Y24,dJ25_phi24, ...
dJ25_X25,dJ25_Y25,dJ25_phi25, ...
dJ25_X26,dJ25_Y26,dJ25_phi26]';
Grad_6=[dJ26_X21,dJ26_Y21,dJ26_phi21, ...
dJ26_X22,dJ26_Y22,dJ26_phi22, ...
dJ26_X23,dJ26_Y23,dJ26_phi23, ...
dJ26_X24,dJ26_Y24,dJ26_phi24, ...
dJ26_X25,dJ26_Y25,dJ26_phi25, ...
dJ26_X26,dJ26_Y26,dJ26_phi26]';

Grad=[Grad_1;Grad_2;Grad_3;Grad_4;Grad_5;Grad_6];
end