function Grad=Grad_J1(Z)
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
iota1_11=0.01; iota1_12=5;   iota1_13=1;   iota1_14=10;
iota1_21=0.01; iota1_22=5;   iota1_23=1;   iota1_24=10;
iota1_31=2; iota1_32=1;   iota1_33=20;   iota1_34=1;   iota1_35=1;
iota1_41=2; iota1_42=1;   iota1_43=20;   iota1_44=1;   iota1_45=1;
iota1_51=2; iota1_52=1;   iota1_53=20;   iota1_54=1;   iota1_55=1;
iota1_61=2; iota1_62=1;   iota1_63=20;   iota1_64=1;   iota1_65=1;
% Formation
DeltaX_11=0; DeltaX_12=0; DeltaX_13=0; DeltaX_14=0; DeltaX_15=0; DeltaX_16=0;
DeltaY_11=100; DeltaY_12=-100; DeltaY_13=150; DeltaY_14=50; DeltaY_15=-50; DeltaY_16=-150;
% DeltaY_11=0; DeltaY_12=0; DeltaY_13=0; DeltaY_14=0; DeltaY_15=0; DeltaY_16=0;
% DeltaX_11=100; DeltaX_12=-100; DeltaX_13=150; DeltaX_14=50; DeltaX_15=-50; DeltaX_16=-150;
%% Gradinet_X
X=[X_1',X_2']'; %横着排列的，第i行表示第i个联盟
X_2sum=sum(X_2);
X_1dfsum=X_1(3)+X_1(4)+X_1(5)+X_1(6);
dJ11_X11=-iota1_11*(6*X(1,1)-X_2sum)    +iota1_12*( X(1,1)-1/2*(X(1,1)+X(1,2))-DeltaX_11 )*(1-1/2)    +iota1_14*(X(1,1)-X_2TRG);
dJ11_X12= iota1_12*(X(1,1)-1/2*(X(1,1)+X(1,2))-DeltaX_11)*(-1/2);
dJ11_X13=0; dJ11_X14=0; dJ11_X15=0; dJ11_X16=0;

dJ12_X11= iota1_22*(X(1,2)-1/2*(X(1,1)+X(1,2))-DeltaX_12)*(-1/2);
dJ12_X12= -iota1_21*(6*X(1,2)-X_2sum)   +iota1_22*(X(1,2)-1/2*(X(1,1)+X(1,2))-DeltaX_12)*(1-1/2)      +iota1_24*(X(1,2)-X_2TRG);
dJ12_X13=0; dJ12_X14=0; dJ12_X15=0; dJ12_X16=0;

dJ13_X11=0;dJ13_X12=0;
dJ13_X13=iota1_31*(X(1,3)-1/4*X_1dfsum-DeltaX_13)*(1-1/4) +iota1_32*(1/4*X_1dfsum-X_1TRG)*(1/4)    +iota1_33*(X(1,3)-1/2*(X(2,1)+X_1TRG))   +iota1_34*(X(1,3)-X(2,1));
dJ13_X14=iota1_31*(X(1,3)-1/4*X_1dfsum-DeltaX_13)*(-1/4)  +iota1_32*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ13_X15=iota1_31*(X(1,3)-1/4*X_1dfsum-DeltaX_13)*(-1/4)  +iota1_32*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ13_X16=iota1_31*(X(1,3)-1/4*X_1dfsum-DeltaX_13)*(-1/4)  +iota1_32*(1/4*X_1dfsum-X_1TRG)*(1/4);

dJ14_X11=0;dJ14_X12=0;
dJ14_X13=iota1_41*(X(1,4)-1/4*X_1dfsum-DeltaX_14)*(-1/4)  +iota1_42*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ14_X14=iota1_41*(X(1,4)-1/4*X_1dfsum-DeltaX_14)*(1-1/4) +iota1_42*(1/4*X_1dfsum-X_1TRG)*(1/4)    +iota1_43*(X(1,4)-1/2*(X(2,2)+X_1TRG))   +iota1_44*(X(1,4)-X(2,2));
dJ14_X15=iota1_41*(X(1,4)-1/4*X_1dfsum-DeltaX_14)*(-1/4)  +iota1_42*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ14_X16=iota1_41*(X(1,4)-1/4*X_1dfsum-DeltaX_14)*(-1/4)  +iota1_42*(1/4*X_1dfsum-X_1TRG)*(1/4);


dJ15_X11=0;dJ15_X12=0;
dJ15_X13=iota1_51*(X(1,5)-1/4*X_1dfsum-DeltaX_15)*(-1/4)  +iota1_52*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ15_X14=iota1_51*(X(1,5)-1/4*X_1dfsum-DeltaX_15)*(-1/4)  +iota1_52*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ15_X15=iota1_51*(X(1,5)-1/4*X_1dfsum-DeltaX_15)*(1-1/4) +iota1_52*(1/4*X_1dfsum-X_1TRG)*(1/4)    +iota1_53*(X(1,5)-1/2*(X(2,3)+X_1TRG))   +iota1_54*(X(1,5)-X(2,3));
dJ15_X16=iota1_51*(X(1,5)-1/4*X_1dfsum-DeltaX_15)*(-1/4)  +iota1_52*(1/4*X_1dfsum-X_1TRG)*(1/4);


dJ16_X11=0;dJ16_X12=0;
dJ16_X13=iota1_61*(X(1,6)-1/4*X_1dfsum-DeltaX_16)*(-1/4)  +iota1_62*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ16_X14=iota1_61*(X(1,6)-1/4*X_1dfsum-DeltaX_16)*(-1/4)  +iota1_62*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ16_X15=iota1_61*(X(1,6)-1/4*X_1dfsum-DeltaX_16)*(-1/4)  +iota1_62*(1/4*X_1dfsum-X_1TRG)*(1/4);
dJ16_X16=iota1_61*(X(1,6)-1/4*X_1dfsum-DeltaX_16)*(1-1/4) +iota1_62*(1/4*X_1dfsum-X_1TRG)*(1/4)    +iota1_63*(X(1,6)-1/2*(X(2,3)+X_1TRG))   +iota1_64*(X(1,6)-X(2,3));

%% Gradinet_Y
Y=[Y_1',Y_2']'; %横着排列的，第i行表示第i个联盟
Y_2sum=sum(Y_2);
Y_1dfsum=Y_1(3)+Y_1(4)+Y_1(5)+Y_1(6);

dJ11_Y11=-iota1_11*(6*Y(1,1)-Y_2sum)    +iota1_12*( Y(1,1)-1/2*(Y(1,1)+Y(1,2))-DeltaY_11 )*(1-1/2)    +iota1_14*(Y(1,1)-Y_2TRG);
dJ11_Y12= iota1_12*(Y(1,1)-1/2*(Y(1,1)+Y(1,2))-DeltaY_11)*(-1/2);
dJ11_Y13=0; dJ11_Y14=0; dJ11_Y15=0; dJ11_Y16=0;


dJ12_Y11= iota1_22*(Y(1,2)-1/2*(Y(1,1)+Y(1,2))-DeltaY_12)*(-1/2);
dJ12_Y12= -iota1_21*(6*Y(1,2)-Y_2sum)   +iota1_22*(Y(1,2)-1/2*(Y(1,1)+Y(1,2))-DeltaY_12)*(1-1/2)      +iota1_24*(Y(1,2)-Y_2TRG);
dJ12_Y13=0; dJ12_Y14=0; dJ12_Y15=0; dJ12_Y16=0;


dJ13_Y11=0;dJ13_Y12=0;
dJ13_Y13=iota1_31*(Y(1,3)-1/4*Y_1dfsum-DeltaY_13)*(1-1/4) +iota1_32*(1/4*Y_1dfsum-Y_1TRG)*(1/4)    +iota1_33*(Y(1,3)-1/2*(Y(2,1)+Y_1TRG))   +iota1_34*(Y(1,3)-Y(2,1));
dJ13_Y14=iota1_31*(Y(1,3)-1/4*Y_1dfsum-DeltaY_13)*(-1/4)  +iota1_32*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ13_Y15=iota1_31*(Y(1,3)-1/4*Y_1dfsum-DeltaY_13)*(-1/4)  +iota1_32*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ13_Y16=iota1_31*(Y(1,3)-1/4*Y_1dfsum-DeltaY_13)*(-1/4)  +iota1_32*(1/4*Y_1dfsum-Y_1TRG)*(1/4);


dJ14_Y11=0;dJ14_Y12=0;
dJ14_Y13=iota1_41*(Y(1,4)-1/4*Y_1dfsum-DeltaY_14)*(-1/4)  +iota1_42*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ14_Y14=iota1_41*(Y(1,4)-1/4*Y_1dfsum-DeltaY_14)*(1-1/4) +iota1_42*(1/4*Y_1dfsum-Y_1TRG)*(1/4)    +iota1_43*(Y(1,4)-1/2*(Y(2,2)+Y_1TRG))   +iota1_44*(Y(1,4)-Y(2,2));
dJ14_Y15=iota1_41*(Y(1,4)-1/4*Y_1dfsum-DeltaY_14)*(-1/4)  +iota1_42*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ14_Y16=iota1_41*(Y(1,4)-1/4*Y_1dfsum-DeltaY_14)*(-1/4)  +iota1_42*(1/4*Y_1dfsum-Y_1TRG)*(1/4);


dJ15_Y11=0;dJ15_Y12=0;
dJ15_Y13=iota1_51*(Y(1,5)-1/4*Y_1dfsum-DeltaY_15)*(-1/4)  +iota1_52*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ15_Y14=iota1_51*(Y(1,5)-1/4*Y_1dfsum-DeltaY_15)*(-1/4)  +iota1_52*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ15_Y15=iota1_51*(Y(1,5)-1/4*Y_1dfsum-DeltaY_15)*(1-1/4) +iota1_52*(1/4*Y_1dfsum-Y_1TRG)*(1/4)    +iota1_53*(Y(1,5)-1/2*(Y(2,3)+Y_1TRG))   +iota1_54*(Y(1,5)-Y(2,3));
dJ15_Y16=iota1_51*(Y(1,5)-1/4*Y_1dfsum-DeltaY_15)*(-1/4)  +iota1_52*(1/4*Y_1dfsum-Y_1TRG)*(1/4);


dJ16_Y11=0;dJ16_Y12=0;
dJ16_Y13=iota1_61*(Y(1,6)-1/4*Y_1dfsum-DeltaY_16)*(-1/4)  +iota1_62*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ16_Y14=iota1_61*(Y(1,6)-1/4*Y_1dfsum-DeltaY_16)*(-1/4)  +iota1_62*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ16_Y15=iota1_61*(Y(1,6)-1/4*Y_1dfsum-DeltaY_16)*(-1/4)  +iota1_62*(1/4*Y_1dfsum-Y_1TRG)*(1/4);
dJ16_Y16=iota1_61*(Y(1,6)-1/4*Y_1dfsum-DeltaY_16)*(1-1/4) +iota1_62*(1/4*Y_1dfsum-Y_1TRG)*(1/4)    +iota1_63*(Y(1,6)-1/2*(Y(2,3)+Y_1TRG))   +iota1_64*(Y(1,6)-Y(2,3));

%% Gradinet_phi
phi=[phi_1',phi_2']'; %横着排列的，第i行表示第i个联盟

dJ11_phi11=iota1_13*(phi(1,1)-2*pi/3);dJ11_phi12=0;dJ11_phi13=0;dJ11_phi14=0;dJ11_phi15=0;dJ11_phi16=0;
dJ12_phi11=0;dJ12_phi12=iota1_23*(phi(1,2)-pi/3);dJ12_phi13=0;dJ12_phi14=0;dJ12_phi15=0;dJ12_phi16=0;
dJ13_phi11=0;dJ13_phi12=0;dJ13_phi13=iota1_35*(phi(1,3)-phi(2,1)+pi);dJ13_phi14=0;dJ13_phi15=0;dJ13_phi16=0;
dJ14_phi11=0;dJ14_phi12=0;dJ14_phi13=0;dJ14_phi14=iota1_45*(phi(1,4)-phi(2,2)+pi);dJ14_phi15=0;dJ14_phi16=0;
dJ15_phi11=0;dJ15_phi12=0;dJ15_phi13=0;dJ15_phi14=0;dJ15_phi15=iota1_55*(phi(1,5)-phi(2,3)+pi);dJ15_phi16=0;
dJ16_phi11=0;dJ16_phi12=0;dJ16_phi13=0;dJ16_phi14=0;dJ16_phi15=0;dJ16_phi16=iota1_65*(phi(1,6)-phi(2,3)+pi);

%% Total_Gradient
Grad_1=[dJ11_X11,dJ11_Y11,dJ11_phi11, ...
    dJ11_X12,dJ11_Y12,dJ11_phi12, ...
    dJ11_X13,dJ11_Y13,dJ11_phi13, ...
    dJ11_X14,dJ11_Y14,dJ11_phi14, ...
    dJ11_X15,dJ11_Y15,dJ11_phi15, ...
    dJ11_X16,dJ11_Y16,dJ11_phi16]';
Grad_2=[dJ12_X11,dJ12_Y11,dJ12_phi11, ...
    dJ12_X12,dJ12_Y12,dJ12_phi12, ...
    dJ12_X13,dJ12_Y13,dJ12_phi13, ...
    dJ12_X14,dJ12_Y14,dJ12_phi14, ...
    dJ12_X15,dJ12_Y15,dJ12_phi15, ...
    dJ12_X16,dJ12_Y16,dJ12_phi16]';
Grad_3=[dJ13_X11,dJ13_Y11,dJ13_phi11, ...
    dJ13_X12,dJ13_Y12,dJ13_phi12, ...
    dJ13_X13,dJ13_Y13,dJ13_phi13, ...
    dJ13_X14,dJ13_Y14,dJ13_phi14, ...
    dJ13_X15,dJ13_Y15,dJ13_phi15, ...
    dJ13_X16,dJ13_Y16,dJ13_phi16]';
Grad_4=[dJ14_X11,dJ14_Y11,dJ14_phi11, ...
    dJ14_X12,dJ14_Y12,dJ14_phi12, ...
    dJ14_X13,dJ14_Y13,dJ14_phi13, ...
    dJ14_X14,dJ14_Y14,dJ14_phi14, ...
    dJ14_X15,dJ14_Y15,dJ14_phi15, ...
    dJ14_X16,dJ14_Y16,dJ14_phi16]';
Grad_5=[dJ15_X11,dJ15_Y11,dJ15_phi11, ...
    dJ15_X12,dJ15_Y12,dJ15_phi12, ...
    dJ15_X13,dJ15_Y13,dJ15_phi13, ...
    dJ15_X14,dJ15_Y14,dJ15_phi14, ...
    dJ15_X15,dJ15_Y15,dJ15_phi15, ...
    dJ15_X16,dJ15_Y16,dJ15_phi16]';
Grad_6=[dJ16_X11,dJ16_Y11,dJ16_phi11, ...
    dJ16_X12,dJ16_Y12,dJ16_phi12, ...
    dJ16_X13,dJ16_Y13,dJ16_phi13, ...
    dJ16_X14,dJ16_Y14,dJ16_phi14, ...
    dJ16_X15,dJ16_Y15,dJ16_phi15, ...
    dJ16_X16,dJ16_Y16,dJ16_phi16]';
Grad=[Grad_1;Grad_2;Grad_3;Grad_4;Grad_5;Grad_6];
end