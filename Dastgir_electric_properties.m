clc;clear;close all

% This code is written by Navid Dastgir
% Micromechanical modeling of the functional role of graphene nanoplatelets
% predict magneto-electro-elastic properties using Mori-Tanaka
% as function of piezelectric volume fraction
% Four phase composite
% (PDMS matrix + GNP) + BaTiO3 + CoFe2O4

%% define (matrix + 0.75% GNP) properties (for sample)
c11_f1=3.4334e9; c12_f1=1.4024e9;
c44_f1 = (c11_f1 - c12_f1)/2; 
k_f1 = 8.854e-12; mu_f1 = 4*pi*1e-7; 
L_F1 = zeros(12,12);
L_F1(1:3,1:3) = [c11_f1, c12_f1, c12_f1; c12_f1, c11_f1, c12_f1; c12_f1, c12_f1, c11_f1];
L_F1(4:6,4:6) = diag([c44_f1, c44_f1, c44_f1]);
L_F1(7:9,7:9) = eye(3) * k_f1;
L_F1(10:12,10:12) = eye(3) * mu_f1;
S_F1 = inv(L_F1);
%% define (BaTiO3 piezoelectric) properties
c11_f2=16.6e10; c12_f2=7.7e10; c13_f2=7.8e10; c33_f2=16.2e10; c44_f2=4.3e10; c66_f2=(c11_f2 - c12_f2)/2;
e31_f2=-4.4; e33_f2=18.6; e15_f2=11.6;
k11_f2=112e-10; k33_f2=126e-10; mu11_f2=5e-6; mu33_f2=10e-6;
L_F2 = zeros(12,12);
L_F2(1:3,1:3) = [c11_f2, c12_f2, c13_f2; c12_f2, c11_f2, c13_f2; c13_f2, c13_f2, c33_f2];
L_F2(4:6,4:6) = diag([c44_f2, c44_f2, c66_f2]);
L_F2(9,1:3) = [e31_f2, e31_f2, e33_f2]; L_F2(1:3,9) = L_F2(9,1:3)';
L_F2(8,4) = e15_f2; L_F2(4,8) = e15_f2; L_F2(7,5) = e15_f2; L_F2(5,7) = e15_f2;
L_F2(7:9,7:9) = diag([k11_f2, k11_f2, k33_f2]);
L_F2(10:12,10:12) = diag([mu11_f2, mu11_f2, mu33_f2]);
%% define (CoFe2O4 piezomagnetic) properties
c11_f3=28.6e10; c12_f3=17.3e10; c13_f3=17.05e10; c33_f3=26.95e10; c44_f3=4.53e10; c66_f3=(c11_f3 - c12_f3)/2; 
q31_f3=580.3; q33_f3=699.7; q15_f3=550;
k11_f3=0.8e-10; k33_f3=0.93e-10; 
mu11_f3 = 590e-6; mu33_f3 = 157e-6; 
L_F3 = zeros(12,12);
L_F3(1:3,1:3) = [c11_f3, c12_f3, c13_f3; c12_f3, c11_f3, c13_f3; c13_f3, c13_f3, c33_f3];
L_F3(4:6,4:6) = diag([c44_f3, c44_f3, c66_f3]);
L_F3(12,1:3) = [q31_f3, q31_f3, q33_f3]; L_F3(1:3,12) = L_F3(12,1:3)';
L_F3(11,4) = q15_f3; L_F3(4,11) = q15_f3; L_F3(10,5) = q15_f3; L_F3(5,10) = q15_f3;
L_F3(7:9,7:9) = diag([k11_f3, k11_f3, k33_f3]);
L_F3(10:12,10:12) = diag([mu11_f3, mu11_f3, mu33_f3]);
%% Eshelby Tensor components
c11_mat = L_F1(1,1); c12_mat = L_F1(1,2); c13_mat = L_F1(1,3);
e31_mat = L_F1(9,1); q31_mat = L_F1(12,1); 
s1111 = (5*c11_mat + c12_mat)/(8*c11_mat);
s1122 = (3*c12_mat - c11_mat)/(8*c11_mat);
s1133 = (c13_mat)/(2*c11_mat);
s2211 = (3*c12_mat - c11_mat)/(8*c11_mat);
s2222 = (5*c11_mat + c12_mat)/(8*c11_mat);
s2233 = (c13_mat)/(2*c11_mat);
s2323 = 0.25; s1313 = 0.25; s1212 = (3*c11_mat - c12_mat)/(8*c11_mat);
s1143 = (e31_mat)/(2*c11_mat); s1153 = (q31_mat)/(2*c11_mat);
s2243 = (e31_mat)/(2*c11_mat); s2253 = (q31_mat)/(2*c11_mat);
s4141 = 0.5; s4242 = 0.5; s5151 = 0.5; s5252 = 0.5;
%% Eshelby Matrix Construction
s = zeros(12,12);
s(1,1) = s1111; s(1,2) = s1122; s(1,3) = s1133;
s(2,1) = s2211; s(2,2) = s2222; s(2,3) = s2233;
s(4,4) = 2 * s2323; s(6,6) = 2 * s1212; s(5,5) = 2 * s1313;
s(1,9) = s1143; s(1,12) = s1153;
s(2,9) = s2243; s(2,12) = s2253;
s(7,7) = 2 * s4141; s(8,8) = 2 * s4242;
s(10,10) = 2 * s5151; s(11,11) = 2 * s5252;
%%  Difine Dilute concentration tensor 
I = eye(12);
A_DIL_2 = inv(I + s * S_F1 * (L_F2 - L_F1));
A_DIL_3 = inv(I + s * S_F1 * (L_F3 - L_F1));
%% Difine volume fraction and Mori-Tanaka concentration tensor 
VF_BaTiO3 = 0:0.02:0.6;
y_C11 = zeros(1, length(VF_BaTiO3));
y_C12 = zeros(1, length(VF_BaTiO3));
y_C44 = zeros(1, length(VF_BaTiO3));
y_e31 = zeros(1, length(VF_BaTiO3));
y_e15 = zeros(1, length(VF_BaTiO3));
y_q31 = zeros(1, length(VF_BaTiO3));
y_q15 = zeros(1, length(VF_BaTiO3));
for i = 1:length(VF_BaTiO3)   
    c2 = VF_BaTiO3(i);   % Piezoelectric volume fraction    
    c3 = 0.6 - c2;       % Piezomagnetic volume fraction  
    c1 = 0.4;             % (Matrix + 0.75 GNP) volume fraction
    Denom = c1*I + c2*A_DIL_2 + c3*A_DIL_3;
    A_MOR_2 = A_DIL_2 / Denom; 
    A_MOR_3 = A_DIL_3 / Denom;
    L_eff = L_F1 + c2*(L_F2-L_F1)*A_MOR_2 + c3*(L_F3-L_F1)*A_MOR_3;
    % Report magneto-electro-elastic properties
    y_C11(i) = L_eff(1,1);
    y_C12(i) = L_eff(1,2);
    y_C44(i) = L_eff(4,4);
    y_e31(i) = L_eff(9,1); 
    y_e15(i) = L_eff(7,5); 
    y_q31(i) = L_eff(12,1); 
    y_q15(i) = L_eff(10,5); 
end 
% Report magneto-electro-elastic properties as function of piezoelectric
% volume fraction
y1 = [VF_BaTiO3(:), y_C11(:)/1e9];
y2 = [VF_BaTiO3(:), y_C12(:)/1e9];
y3 = [VF_BaTiO3(:), y_C44(:)/1e9];
y4 = [VF_BaTiO3(:), y_e15(:)];
y5 = [VF_BaTiO3(:), y_e31(:)];
y6 = [VF_BaTiO3(:), y_q15(:)];
y7 = [VF_BaTiO3(:), y_q31(:)];
%% 6. Plotting (Updated X-Axis)
figure('Name', 'Mechanical Stiffness', 'Color', 'w');
subplot(3,1,1); plot(VF_BaTiO3*100, y_C11/1e9, 'b-o', 'LineWidth', 1.5);
ylabel('C_{11} (GPa)'); grid on; title('Longitudinal Stiffness C_{11}');
subplot(3,1,2); plot(VF_BaTiO3*100, y_C12/1e9, 'r-s', 'LineWidth', 1.5);
ylabel('C_{12} (GPa)'); grid on; title('Transverse Stiffness C_{12}'); 
subplot(3,1,3); plot(VF_BaTiO3*100, y_C44/1e9, 'k-^', 'LineWidth', 1.5);
ylabel('C_{44} (GPa)'); xlabel('BaTiO_3 Volume Fraction (%)'); grid on; title('Shear Stiffness C_{44}');

figure('Name', 'Piezoelectric Coefficients', 'Color', 'w');
subplot(2,1,1); plot(VF_BaTiO3*100, y_e31, 'm-o', 'LineWidth', 1.5);
ylabel('e_{31} (C/m^2)'); grid on; title('Piezoelectric Coefficient e_{31}');
subplot(2,1,2); plot(VF_BaTiO3*100, y_e15, 'c-s', 'LineWidth', 1.5);
ylabel('e_{15} (C/m^2)'); xlabel('BaTiO_3 Volume Fraction (%)'); grid on; title('Piezoelectric Coefficient e_{15}');

figure('Name', 'Piezomagnetic Coefficients', 'Color', 'w');
subplot(2,1,1); plot(VF_BaTiO3*100, y_q31, 'g-o', 'LineWidth', 1.5);
ylabel('q_{31} (N/Am)'); grid on; title('Piezomagnetic Coefficient q_{31}');
subplot(2,1,2); plot(VF_BaTiO3*100, y_q15, 'y-s', 'LineWidth', 1.5);
ylabel('q_{15} (N/Am)'); xlabel('BaTiO_3 Volume Fraction (%)'); grid on; title('Piezomagnetic Coefficient q_{15}');