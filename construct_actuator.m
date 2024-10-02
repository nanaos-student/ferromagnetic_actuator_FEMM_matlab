openfemm()
opendocument('femm_template.FEM')
mi_saveas('actuator.fem')

%Load files
load('coil1p.mat')
load('coil2p.mat')
load('coil3p.mat')
load('corep.mat')
load('corep.mat')
load('moverp.mat')
load('coil4p.mat')

%Add nodes
%Core
mi_addnode(corep(:,1),corep(:,2))
for i=1:size(corep,1)
    x = corep(i,1);
    y = corep(i,2);
    mi_selectnode(x,y);
    mi_setnodeprop('None',1);
    mi_clearselected();
end
%Coil1
mi_addnode(coil1p(:,1),coil1p(:,2))
for i=1:size(coil1p,1)
    x = coil1p(i,1);
    y = coil1p(i,2);
    mi_selectnode(x,y);
    mi_setnodeprop('None',6);
    mi_clearselected();
end
%Coil2
mi_addnode(coil2p(:,1),coil2p(:,2))
for i=1:size(coil2p,1)
    x = coil2p(i,1);
    y = coil2p(i,2);
    mi_selectnode(x,y);
    mi_setnodeprop('None',5);
    mi_clearselected();
end
%Coil3
mi_addnode(coil3p(:,1),coil3p(:,2))
for i=1:size(coil3p,1)
    x = coil3p(i,1);
    y = coil3p(i,2);
    mi_selectnode(x,y);
    mi_setnodeprop('None',4);
    mi_clearselected();
end
%Coil4
mi_addnode(coil4p(:,1),coil4p(:,2))
for i=1:size(coil4p,1)
    x = coil4p(i,1);
    y = coil4p(i,2);
    mi_selectnode(x,y);
    mi_setnodeprop('None',3);
    mi_clearselected();
end
%Mover
mi_addnode(moverp(:,1),moverp(:,2))
for i=1:size(moverp,1)
    x = moverp(i,1);
    y = moverp(i,2);
    mi_selectnode(x,y);
    mi_setnodeprop('None',2);
    mi_clearselected();
end

%Add segments
%Core
for i = 1:size(corep, 1)
    mi_addsegment(corep(i, 1), corep(i, 2), corep(mod(i, size(corep, 1)) + 1, 1), corep(mod(i, size(corep, 1)) + 1, 2));
    mi_selectsegment((corep(i, 1) + corep(mod(i, size(corep, 1)) + 1, 1)) / 2, (corep(i, 2) + corep(mod(i, size(corep, 1)) + 1, 2)) / 2);
    mi_setsegmentprop('Corep', 1, 0, 0, corep(i, 3));
    mi_clearselected();
end
%Mover
for i = 1:size(moverp, 1)
    mi_addsegment(moverp(i, 1), moverp(i, 2), moverp(mod(i, size(moverp, 1)) + 1, 1), moverp(mod(i, size(moverp, 1)) + 1, 2));
    mi_selectsegment((moverp(i, 1) + moverp(mod(i, size(moverp, 1)) + 1, 1)) / 2, (moverp(i, 2) + moverp(mod(i, size(moverp, 1)) + 1, 2)) / 2);
    mi_setsegmentprop('Moverp', 1, 0, 0, moverp(i, 3));
    mi_clearselected();
end
%Coil1
for i = 1:size(coil1p, 1)
    mi_addsegment(coil1p(i, 1), coil1p(i, 2), coil1p(mod(i, size(coil1p, 1)) + 1, 1), coil1p(mod(i, size(coil1p, 1)) + 1, 2));
    mi_selectsegment((coil1p(i, 1) + coil1p(mod(i, size(coil1p, 1)) + 1, 1)) / 2, (coil1p(i, 2) + coil1p(mod(i, size(coil1p, 1)) + 1, 2)) / 2);
    mi_setsegmentprop('Coil1p', 1, 0, 0, coil1p(i, 3));
    mi_clearselected();
end
%Coil2
for i = 1:size(coil2p, 1)
    mi_addsegment(coil2p(i, 1), coil2p(i, 2), coil2p(mod(i, size(coil2p, 1)) + 1, 1), coil2p(mod(i, size(coil2p, 1)) + 1, 2));
    mi_selectsegment((coil2p(i, 1) + coil2p(mod(i, size(coil2p, 1)) + 1, 1)) / 2, (coil2p(i, 2) + coil2p(mod(i, size(coil2p, 1)) + 1, 2)) / 2);
    mi_setsegmentprop('coil2p', 1, 0, 0, coil2p(i, 3));
    mi_clearselected();
end
%Coil3
for i = 1:size(coil3p, 1)
    mi_addsegment(coil3p(i, 1), coil3p(i, 2), coil3p(mod(i, size(coil3p, 1)) + 1, 1), coil3p(mod(i, size(coil3p, 1)) + 1, 2));
    mi_selectsegment((coil3p(i, 1) + coil3p(mod(i, size(coil3p, 1)) + 1, 1)) / 2, (coil3p(i, 2) + coil3p(mod(i, size(coil3p, 1)) + 1, 2)) / 2);
    mi_setsegmentprop('coil3p', 1, 0, 0, coil3p(i, 3));
    mi_clearselected();
end
%Coil4
for i = 1:size(coil4p, 1)
    mi_addsegment(coil4p(i, 1), coil4p(i, 2), coil4p(mod(i, size(coil4p, 1)) + 1, 1), coil4p(mod(i, size(coil4p, 1)) + 1, 2));
    mi_selectsegment((coil4p(i, 1) + coil4p(mod(i, size(coil4p, 1)) + 1, 1)) / 2, (coil4p(i, 2) + coil4p(mod(i, size(coil4p, 1)) + 1, 2)) / 2);
    mi_setsegmentprop('coil4p', 1, 0, 0, coil4p(i, 3));
    mi_clearselected();
end

mi_makeABC();

%Add block properties
%Core
mi_addblocklabel(11.0,7.0);
mi_selectlabel(11.0,7.0);
mi_setblockprop('core_nonlinear',1,1,'None',1,1,'None');
mi_clearselected();
%Mover
mi_addblocklabel(59.0,2.0);
mi_selectlabel(59.0,2.0);
mi_setblockprop('core_nonlinear',1,1,'None',1,2,'None');
mi_clearselected();
%Coil1
mi_addblocklabel(45.0,-64.0);
mi_selectlabel(45.0,-64.0);
mi_setblockprop('copper',1,1,'winding_2',1,6,100);
mi_clearselected();
%Coil2
mi_addblocklabel(48.0,-30.0);
mi_selectlabel(48.0,-30.0);
mi_setblockprop('copper',1,1,'winding_2',1,5,-100);
mi_clearselected();
%Coil3
mi_addblocklabel(46.0,31.0);
mi_selectlabel(46.0,31.0);
mi_setblockprop('copper',1,1,'winding_1',1,4,100);
mi_clearselected();
%Coil4
mi_addblocklabel(52.0,64.0);
mi_selectlabel(52.0,64.0);
mi_setblockprop('copper',1,1,'winding_1',1,3,-100);
mi_clearselected();
%Air
mi_addblocklabel(-37.0,-3.0);
mi_selectlabel(-37.0,-3.0);
mi_setblockprop('air',1,1,'None',1,'None',0);
mi_clearselected();

%Refining mesh density for variable air gap
%Core
mi_selectsegment((corep(10,1)+corep(11,1))/2,(corep(10,2)+corep(11,2))/2);
mi_setsegmentprop('Corep',0.5,0,0,1);
mi_clearselected();
%Mover
mi_selectsegment((moverp(1,1)+moverp(4,1))/2,(moverp(1,2)+moverp(4,2))/2);
mi_setsegmentprop('Moverp',0.5,0,0,2);
mi_clearselected();

%Set winding currents
mi_setcurrent('winding_1',10);
mi_setcurrent('winding_2',10);

mi_smartmesh(1)
mi_showmesh()
mi_createmesh()
mi_analyze()
mi_loadsolution()

%Resistances
%FEMM blockintegral for winding 2
mo_groupselectblock(5);
A_2 = mo_blockintegral(5);  %cross-sectional area
V_2 = mo_blockintegral(10);  %volume
mo_clearblock();
L_2 = V_2/A_2;  %active length
R_2 = (100*L_2)/((58*10^6)*((0.6*A_2)/100));  %resistance 
%FEMM get circuit properties for winding 2
CP = mo_getcircuitproperties('winding_2');
R_F = CP(2)/CP(1);  %resistance
%CAD for winding 2
Cross_section_area = 1.656191*10^-3;  %cross-sectional area
active_length = 0.049;  %active length
R_3D = (100*active_length)/((58*10^6)*((0.6*Cross_section_area)/100));  %resistance

%Total power loss for the whole system
Total_power = [];
Total_power_F = [];
Total_power_3D = [];
for i = 1:10
    Total_power = [Total_power, ((i*i)*(2*R_2))];
    Total_power_F = [Total_power_F, ((i*i)*(2*R_F))];
    Total_power_3D = [Total_power_3D, ((i*i)*(2*R_3D))];
end

%Power loss graph
figure;
plot(1:10,Total_power,'r');
hold on;
plot(1:10,Total_power_F,'b');
hold on;
plot(1:10,Total_power_3D,'g');
hold off;
xlabel('Applied Current (A)');
ylabel('Total Power loss (W)');
title('Total power loss against applied current');
legend('FEMM blockintegral','FEMM CP','CAD');

%Inductances 
Inductance_nonlinear = [];  %non-linear
Inductance_linear = [];  %linear
Inductance_no_fringing = [];  %no fringing
Inductance_fringing = [];  %fringing
FE_co_nonlinear = [];  %non-linear co-energy
Psi_all_nonlinear = [];  %non-linear flux linkages
Psi_all_linear = [];  %linear flux linkages
dx = -0.1;  %rate of armature movement
%Non-linear case
for i = 0.0:0.1:4.9 
    mi_selectgroup(2);
    mi_movetranslate(dx,0.0);
    mi_clearselected();
    mi_analyze()
    mi_loadsolution()
    CP = mo_getcircuitproperties('winding_1');
    Inductance_nonlinear = [Inductance_nonlinear, CP(3)/CP(1)];
    %Flux linkages at various positions
    for k = 0.0:2.0:10
        mi_setcurrent('winding_1',k);
        mi_setcurrent('winding_2',k);
        mi_analyze()
        mi_loadsolution()
        CP_nonlinear = mo_getcircuitproperties('winding_1');
        Psi_all_nonlinear = [Psi_all_nonlinear, CP_nonlinear(3)];
    end
    mi_setcurrent('winding_1',10);
    mi_setcurrent('winding_2',10);
    mo_groupselectblock();
    Ec = mo_blockintegral(17);
    mo_clearblock;
    FE_co_nonlinear = [FE_co_nonlinear, Ec/2];
end

%Armature movement
mi_selectgroup(2);
mi_movetranslate(4.9,0.0);
mi_clearselected();

mi_smartmesh(1)
mi_showmesh()
mi_createmesh()
mi_analyze()
mi_loadsolution()

%Linear case
for i = 0.0:0.1:4.9 
    mi_selectlabel(11.0,7.0);
    mi_setblockprop('core_linear',1,1,'None',1,1,'None');
    mi_clearselected();
    mi_selectlabel(59.0,2.0);
    mi_setblockprop('core_linear',1,1,'None',1,2,'None');
    mi_clearselected();
    mi_selectgroup(2);
    mi_movetranslate(dx,0.0);
    mi_clearselected();
    mi_analyze()
    mi_loadsolution()
    CP = mo_getcircuitproperties('winding_1');
    Inductance_linear = [Inductance_linear, CP(3)/CP(1)];
    %Flux linkages at various positions
    for k = 0.0:2.0:10
        mi_setcurrent('winding_1',k);
        mi_setcurrent('winding_2',k);
        mi_analyze()
        mi_loadsolution()
        CP_linear = mo_getcircuitproperties('winding_1');
        Psi_all_linear = [Psi_all_linear, CP_linear(3)];
    end
    mi_setcurrent('winding_1',10);
    mi_setcurrent('winding_2',10);
end

%Core, Armature and Air gap dimensions
L_c = 0.1445;  %flux path in core
L_gf = 5*10^-4;  %flux path in fixed air gap
Permiability_free_space = 1.256637061*10^-6;  %permeability of free space
Permiability = Permiability_free_space*1000;  %relative permeability
A_a = 4*10^-4;  %flux area for armature
A_c = 7.5*10^-4;  %flux area for core
A_gf = 4*10^-4;  %flux area for fixed air gap
L_a = 0.05;  %flux path in armature
W = 0.02;  %width of smallest area
T = 0.02;  %depth of smallest area

%Analytical inductances
for i = 0.0:0.1:4.9
    real_i = i/1000;
    L_no_f = (100*100)/((L_c/(Permiability*A_c))+(L_a/(Permiability*A_a))+(real_i/(Permiability_free_space*A_a))+((L_gf/(Permiability_free_space*A_gf))));
    Inductance_no_fringing = [Inductance_no_fringing, L_no_f];  
    L_f = (100*100)/((L_c/(Permiability*A_c))+(L_a/(Permiability*A_a))+(real_i/(Permiability_free_space*((W+(2*real_i))*(T+(2*real_i)))))+((L_gf/(Permiability_free_space*A_gf))));
    Inductance_fringing = [Inductance_fringing, L_f];  
end

%Winding inductance versus armature position graphs
figure;
plot(0.1:0.1:5.0, flip(Inductance_nonlinear),'b'); 
hold on;
plot(0.1:0.1:5.0, flip(Inductance_linear),'r');
hold on;
plot(0.1:0.1:5.0, Inductance_no_fringing,'g');
hold on;
plot(0.1:0.1:5.0, Inductance_fringing,'m');
xlabel('Distance between the core and armature (mm)');
ylabel('Inductance (H)');
title('Inductance as a function of armature position');
legend('Numerical nonlinear','Numerical linear','Analytical without fringing','Analytical with fringing');

%Analytical
PF = [];  %Fringing flux linkage
PnoF = [];  %No fringing flux linkage
F_co = [];  %Fringing coenergy
Nf_co = [];  %No fringing coenergy
MF = [];  %Fringing mechanical energy
MnF = [];  %No fringing mechanical energy
displace_f = 0.1:0.1:4.9;  %displacement
Ff = [];  %Fringing force
FnF = [];  %No fringing force
%Flux linkage
for i = 1:1:50
    for j = 0.0:2.0:10
        PF = [PF, Inductance_fringing(1,i)*j];
        PnoF = [PnoF, Inductance_no_fringing(1,i)*j];
    end
end
%Coenergy
for w = 1:6:length(PF)
    if w+5 <= length(PF)
        F_co = [F_co, trapz(0.0:2.0:10, PF(w:w+5))];
        Nf_co = [Nf_co, trapz(0.0:2.0:10, PnoF(w:w+5))];
    end
end
%Mechanical energy
for t = 1:1:50
    MF = [MF, F_co(1,t)-F_co(1,50)];
    MnF = [MnF, Nf_co(1,t)-Nf_co(1,50)];
end
%Force
for q = 1:1:length(displace_f)
    Ff = [Ff, MF(1,q)/(displace_f(1,q)/1000)];
    FnF =[FnF, MnF(1,q)/(displace_f(1,q)/1000)];
end
%PSI-I diagrams
%Fringing
figure;
plot(0.0:2.0:10,PF(295:300),'b');
hold on;
plot(0.0:2.0:10,PF(253:258),'r');
hold on;
plot(0.0:2.0:10,PF(187:192),'g');
hold on;
plot(0.0:2.0:10,PF(127:132),'m');
hold on;
plot(0.0:2.0:10,PF(61:66),'k');
hold on;
plot(0.0:2.0:10,PF(1:6),'c');
hold off;
xlabel('Current (A)');
ylabel('Flux Linkage (Wb)');
title('Analytical Psi-I, Fringing');
legend('Open','Intermediate4','Intermediate3','Intermediate2','Intermediate1','Close');
%No Fringing
figure;
plot(0.0:2.0:10,PnoF(295:300),'b');
hold on;
plot(0.0:2.0:10,PnoF(253:258),'r');
hold on;
plot(0.0:2.0:10,PnoF(187:192),'g');
hold on;
plot(0.0:2.0:10,PnoF(127:132),'m');
hold on;
plot(0.0:2.0:10,PnoF(61:66),'k');
hold on;
plot(0.0:2.0:10,PnoF(1:6),'c');
hold off;
xlabel('Current (A)');
ylabel('Flux Linkage (Wb)');
title('Analytical Psi-I, No Fringing');
legend('Open','Intermediate4','Intermediate3','Intermediate2','Intermediate1','Close');
%Force displacement diagrams
figure;
plot(0.3:0.1:4.9, Ff(3:49),'g');
hold on;
plot(0.3:0.1:4.9, FnF(3:49),'r');
xlabel('Displacement (mm)');
ylabel('Force (N)');
title('Force-Displacement Characteristics');
legend('Fringing','No fringing');


%Numerical 
%PSI-I diagrams
%Linear
figure;
plot(0.0:2.0:10,Psi_all_linear(295:300),'b');
hold on;
plot(0.0:2.0:10,Psi_all_linear(253:258),'r');
hold on;
plot(0.0:2.0:10,Psi_all_linear(187:192),'g');
hold on;
plot(0.0:2.0:10,Psi_all_linear(127:132),'m');
hold on;
plot(0.0:2.0:10,Psi_all_linear(61:66),'k');
hold on;
plot(0.0:2.0:10,Psi_all_linear(1:6),'c');
hold off;
xlabel('Current (A)');
ylabel('Flux Linkage (Wb)');
title('Numerical Psi-I, Linear');
legend('Close','Intermediate1','Intermediate2','Intermediate3','Intermediate4','Open');
%Non-linear
figure;
plot(0.0:2.0:10,Psi_all_nonlinear(295:300),'b');
hold on;
plot(0.0:2.0:10,Psi_all_nonlinear(253:258),'r');
hold on;
plot(0.0:2.0:10,Psi_all_nonlinear(187:192),'g');
hold on;
plot(0.0:2.0:10,Psi_all_nonlinear(127:132),'m');
hold on;
plot(0.0:2.0:10,Psi_all_nonlinear(61:66),'k');
hold on;
plot(0.0:2.0:10,Psi_all_nonlinear(1:6),'c');
hold off;
xlabel('Current (A)');
ylabel('Flux Linkage (Wb)');
title('Numerical Psi-I, Non-Linear');
legend('Close','Intermediate1','Intermediate2','Intermediate3','Intermediate4','Open');
%Linear and Nonlinear
Linear_co = [];  %Linear coenergy
Nonlinear_co = [];  %Nonlinear coenergy
Ml = [];  %Linear mechanical energy
Mnl = [];  %Nonlinear mechanical energy
displace = 4.9:-0.1:0.1;  %displacement
Fl = [];  %Linear force
Fnl = [];  %Nonlinear force
%Coenergy
for w = 1:6:length(Psi_all_linear)
    if w+5 <= length(Psi_all_linear)
        Linear_co = [Linear_co, trapz(0.0:2.0:10, Psi_all_linear(w:w+5))];
        Nonlinear_co = [Nonlinear_co, trapz(0.0:2.0:10, Psi_all_nonlinear(w:w+5))];
    end
end
%Mechanical energy
for t = 1:1:50
    Ml = [Ml, Linear_co(1,t)-Linear_co(1,1)];
    Mnl = [Mnl, Nonlinear_co(1,t)-Nonlinear_co(1,1)];
end
%Force
for q = 1:1:length(displace)
    Fl = [Fl, Ml(1,q)/(displace(1,q)/1000)];
    Fnl =[Fnl, Mnl(1,q)/(displace(1,q)/1000)];
end
flip_Fl = flip(Fl);
flip_Fnl = flip(Fnl);
%Force displacement diagrams
figure;
plot(0.3:0.1:4.9, flip_Fl(3:49), 'g');
hold on;
plot(0.3:0.1:4.9, flip_Fnl(3:49), 'r');
xlabel('Displacement (mm)');
ylabel('Force (N)');
title('Force-Displacement Characteristics');
legend('Linear','Non-linear');
%Co-energy trapz()
figure;
plot(0.1:0.1:5.0, flip(Nonlinear_co)); 
xlabel('Armature positions relative to starting position (mm)');
ylabel('Co-energy (J)');
title('Co-energy against position of armature, trapz()');
legend('0.1 - Closed and 4.9 - Open');
%FEMM Co-energy mo_blockintegral()
figure;
plot(0.1:0.1:5.0, flip(FE_co_nonlinear));
xlabel('Armature positions relative to starting position (mm)');
ylabel('Co-energy (J)');
title('Co-energy against position of armature, blockintegral(17)');
legend('0.1 - Closed and 4.9 - Open');