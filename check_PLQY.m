close all
clear all

%% 
h = 6.62607004e-34 ; %m2 kg / s
c = 3e8;           %m/s
q = 1.6e-19; % C
k_B = 1.38e-23; % J/°K
T = 5800 ; %°K
T_solar =300;
QFLS = [];
QFLS_2 = [];
JSC_all = QFLS;
potential_PCE = JSC_all;
PLQYs = QFLS;
QFLS_rad = QFLS;
names = '';
q = 1.6e-19;
I2 = 50;
bandgaps = [];

%%

FAPI = 0;
TCP = 0;
TCP_163 = 1;
TCP_165 = 0;
TCP_167 = 0;
TCP_17 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerlinesIn = 8;

Folder_name = 'PLQY';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Folder = dir(strcat(Folder_name,'\','*_evaluated.dat'));
%%
for i1 = 1:numel(Folder)
filename = Folder(i1).name;
names = [names;{filename(1:end-14)}];
delimiterIn = '\t';

data = importdata(strcat(Folder_name,'\',filename),delimiterIn,headerlinesIn);
text = data.textdata;
string = text(4,1);
string = strrep(string, 'PLQY = ', '');
string = strrep(string, ' ', '');
PLQY = str2num(char(string));



dat = data.data;
wvl = dat(:,1);
spectrum = dat(:,2);

[~,I1]=max(spectrum);
wvl_max = wvl(I1);
E_max = 1240/wvl_max; 


highE_wvl = fliplr(wvl(I2:I1-I2));
highE_spectrum = fliplr(spectrum(I2:I1-I2));
slope = mean(diff(log(highE_spectrum))./diff(1240./highE_wvl));



energy = 1240./wvl;
bb = q*2.*pi./(4.1356e-15.^3*3e8.^2).*energy.^2.*exp(-1.*energy./0.0259);
energy2 = linspace(min(energy),3,1e5);
bb2 = q*2.*pi./(4.1356e-15.^3*3e8.^2).*energy2.^2.*exp(-1.*energy2./0.0259);


energy4 = 1240./max(highE_wvl):mean(diff(energy2)):max(energy2);
%PL_extrap = log(highE_spectrum(end))+(energy4-1240/max(highE_wvl))*slope;


x_q = energy2(energy2<1240./max(highE_wvl));
test_a=interp1(energy(energy<1240./max(highE_wvl)),spectrum(energy<1240./max(highE_wvl)), x_q, 'linear','extrap');
test_abs = test_a./(bb2(energy2<1240./max(highE_wvl)));
test_abs = [test_abs, max(test_abs).*ones(1,numel(bb2)-numel(x_q))];
absorption = test_abs;
absorption = absorption./max(absorption);
absorption = smooth(absorption,5)';

%  diff_spec = diff(absorption)./diff(energy2);
%  diff_spec = smooth(diff_spec./max(diff_spec),100);
%  diff_spec = diff_spec./(20*max(diff_spec));

 j_0_rad = (trapz(energy2, bb2.*absorption));
 j_0_rad_max = max(j_0_rad);
% 
% solar_data1 = importdata('astmg173.xls');
% solar_data = solar_data1.data;
% solar_data_wvl = solar_data(:,1);
% solar_data_wvl  = solar_data_wvl(solar_data_wvl <max(wvl));
% solar_data = solar_data(1:numel(solar_data_wvl),:);
% solar_data_wvl = solar_data(:,1);
% solar_data_sun = solar_data(:,3); %% W /(m2 nm s)
% 
% solar_data_sun =  solar_data_sun.*solar_data_wvl*1e-9./(h*c);
% 
% abso_re_interp = interp1(1240./energy2, absorption , solar_data_wvl, 'linear','extrap');
% 
% JSC = (cumtrapz(solar_data_wvl, solar_data_sun.*abso_re_interp)).*q/10;
% JSC1 = max(JSC);

if FAPI == 1
    
    j_0_rad_max = 2e-20;
    JSC1 = 25.0; 
    Voc_rad = 1.28;
    
end

if TCP_163 == 1 
    
    j_0_rad_max = 9e-22;
    JSC1 = 22.0;
    Voc_rad = 1.32;

end

if TCP_165 == 1 

    JSC1 = 21.5;
    Voc_rad = 1.34;

end

if TCP_167 == 1 

    JSC1 = 21.0;
    Voc_rad = 1.36;

end


if TCP_17 == 1 

    JSC1 = 20.0;
    Voc_rad = 1.375;

end

figure(1);
semilogy(energy,spectrum, 'DisplayName','Sample number XX','Linewidth',3); 

% find peak wavelength

[maximum,index] = max(spectrum);
peak_wavelength = wvl (index);
bandgap = 1240 / peak_wavelength;
bandgaps = [bandgaps;bandgap];


% Create ylabel
ylabel({'Intensity (a.u)'},'FontSize',16);

% Create xlabel
xlabel({'Energy (eV)'},'FontSize',16);

xlim([1.1 2.1])
ylim([1e-6 5])

hold on

%if i1 ==1 
%    plot(energy, bb./q)
%end

%plot(energy2, absorption, 'o')
%plot(energy2(1:end-1),diff_spec)
%legend('PL','\phi_{BB,300K}','norm. abs from reciprocity','dAbs/dE', 'Location', 'SouthEast')
%legend('PL','Location', 'SouthEast')

    
    QFLS = [QFLS;0.0259*log(PLQY*JSC1/j_0_rad_max)];
    QFLS_2 = [QFLS_2; Voc_rad + 0.0259*log(PLQY)];
    PLQYs = [PLQYs;PLQY];
    QFLS_rad = [QFLS_rad;0.0259*log(JSC1/j_0_rad_max)];
    JSC_all = [JSC_all; JSC1];
    potential_PCE = [potential_PCE; 0.825*0.0259*log(PLQY*JSC1/j_0_rad_max)*JSC1.*1];

end
%%

scatterx = [1:1:numel(Folder)];
x = table(names, QFLS, QFLS_2, PLQYs, QFLS_rad,JSC_all, potential_PCE, bandgaps)
writetable(x, strcat(Folder_name,'\','all_analysed.txt'));

figure (3)
scatter (scatterx,QFLS_2)
 