clc, clear, close all
% Tref=25;%reference temperature [C]
celldata.data.d0_5C = renamevars(readtable('nca1discharge0_5C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d0_5C.Temperature(1)=(Tref+273.15);
celldata.data.d1_0C = renamevars(readtable('nca1discharge1_0C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d1_0C.Temperature(1)=(Tref+273.15);
celldata.data.d2_5C = renamevars(readtable('nca1discharge2_5C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d2_5C.Temperature(1)=(Tref+273.15);
celldata.data.d5_0C = renamevars(readtable('nca1discharge5_0C','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d5_0C.Temperature(1)=(Tref+273.15);
celldata.data.d7_5C = renamevars(readtable('nca1discharge7_5C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d7_5C.Temperature(1)=(Tref+273.15);
celldata.data.d10_0C = renamevars(readtable('nca1discharge10_0C.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d10_0C.Temperature(1)=(Tref+273.15);
celldata.data.d1_0C45 = renamevars(readtable('nmcdischarge1_0C45Celcius.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d1_0C45.Temperature(1)=(45+273.15);
celldata.data.d1_0C0 = renamevars(readtable('nmcdischarge1_0C0Celcius.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
 % celldata.data.d1_0C0.Temperature(1)=(0+273.15);
celldata.curvepoints=renamevars(readtable('curvepoints.csv','Delimiter','semi'),["x","y"],["DoD","V"]);
%% Data from datasheet
name={};
unit={};
data=[];
% Capacity
name=[name,'Capacity'];
unit=[unit,"Ah"];
data=[data,7.5];
name=[name,'CapacityCurrent'];
unit=[unit,"A"];
data=[data,1.5];
name=[name,'CurvepointCurrent'];
unit=[unit,"A"];
data=[data,7.5*2.5];
% Voltage data
name=[name,'Nom_Voltage'];
unit=[unit,"V"];
data=[data,3.6];
name=[name,'Max_Voltage'];
unit=[unit,"V"];
data=[data,4.2];
name=[name,'Min_Voltage'];
unit=[unit,"V"];
data=[data,3.0];
% Current data
name=[name,'Max_Charge_Current'];
unit=[unit,"A"];
data=[data,120.0];
name=[name,'Max_disCharge_Current'];
unit=[unit,"A"];
data=[data,300.0];
% Temperature limits
name=[name,'Max_Charge_Temperature'];
unit=[unit,"K"];
data=[data,40+273.15];
name=[name,'Min_Charge_Temperature'];
unit=[unit,"K"];
data=[data,0+273.15];
name=[name,'Max_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,60+273.15];
name=[name,'Min_disCharge_Temperature'];
unit=[unit,"K"];
data=[data,-30+273.15];
name=[name,'Reference_Temperature'];
unit=[unit,"K"];
data=[data,25+273.15];
% weight
name=[name,'Weight'];
unit=[unit,"kg"];
data=[data,0.32];
celldata.datatable=array2table(data);
celldata.datatable.Properties.VariableNames=name;
celldata.datatable.Properties.VariableUnits=unit;
celldata.cellname='NCA1';
celldata.chemistry="NCA";
celldata.realcellname='UHP341440 NCA 7.5Ah';
clearvars unit data name Tref
save NCA1DATA.mat
