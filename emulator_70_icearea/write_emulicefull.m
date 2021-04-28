clc
clear 
close all

emul1=load('ant_40_geo.emulicefull');
emul2=load('ant_40_geo2.emulicefull');
emul3=load('ant_40_geo3.emulicefull');
% emul4=load('ant_40_geo4.emulicefull');
% emul5=load('ant_40_geo5.emulicefull');
% emul6=load('ant_40_geo6.emulicefull');
% emul7=load('ant_40_geo7.emulicefull');
% emul8=load('ant_40_geo8.emulicefull');
% emul9=load('ant_40_geo9.emulicefull');
% emul10=load('ant_40_geo10.emulicefull');
% % emul11=load('ant_40_geo11.emulicefull');

% emul=[emul1; emul2; emul3; emul4; emul5; emul6; emul7; emul8; emul9; ; ];
emul=[emul1; emul2; emul3;];
% skip the second column
emul(:,2)=[];

figure;
plot(emul)

%write text file
% fid=fopen('emulice_full_dt2000.txt','w');
fid=fopen('emulice_full_dt1000.txt','w');
% fprintf(fid, '%6.0f %6.4f \n', [emul(:,1)' emul(:,2)']);
fprintf(fid, '%6.0f %6.0f \n', [emul']);
fclose(fid);

