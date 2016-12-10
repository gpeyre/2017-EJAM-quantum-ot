function [c,C] = check_neuron_class(name)

V = {}; C = {}; 

%% 1 %%
V{end+1} = 'Motoneuron';
C{1} = {'MI' 'NSML' 'NSMR'};

%% 2 %%
V{end+1} = 'Pharyngeal';
C{end+1} = {'I1L' 'I1R' 'I2L' 'I2R' 'I3' 'I4' 'I6' 'I5' 'M1' 'M2L' 'M2R' 'M3L' 'M3R' 'M4' 'M5' 'MCL' 'MCR'};


%% 3 %%
% Somatic nervous system
V{end+1} = 'Amphid';
C{end+1} = {'AIBL' 'AIBR' 'AIYL' 'AIYR' 'AUAL' 'AUAR'};

%% 4 %%
V{end+1} = 'Integrative ring';
C{end+1} = {'ADAL' 'ADAR'};


%% 5 %%
V{end+1} = 'Motoneuron';
C{end+1} = {'AVL' 'DVB' 'PVNL' 'PVNR' 'RIVL' 'RIVR' 'SMBDL' 'SMBDR' 'SMBVL' 'SMBVR' 'SMDDL' 'SMDDR' 'SMDVL' 'SMDVR' 'AS1' 'AS10' 'AS11' 'AS2' 'AS3' 'AS4' 'AS5' 'AS6' 'AS7' 'AS8' 'AS9' 'DA1' 'DA2' 'DA3' 'DA4' 'DA5' 'DA6' 'DA7' 'DA8' 'DA9' 'DB1' 'DB2' 'DB3' 'DB4' 'DB5' 'DB6' 'DB7' 'DD1' 'DD2' 'DD3' 'DD4' 'DD5' 'DD6' 'HSNL' 'HSNR' 'PDB' 'RIML' 'RIMR' 'RMDDL' 'RMDDR' 'RMDL' 'RMDR' 'RMDVL' 'RMDVR' 'RMED' 'RMEL' 'RMER' 'RMEV' 'RMFL' 'RMFR' 'RMHL' 'RMHR' 'URADL' 'URADR' 'URAVL' 'URAVR' 'VA1' 'VA10' 'VA11' 'VA12' 'VA2' 'VA3' 'VA4' 'VA5' 'VA6' 'VA7' 'VA8' 'VA9' 'VB1' 'VB10' 'VB11' 'VB2' 'VB3' 'VB4' 'VB5' 'VB6' 'VB7' 'VB8' 'VB9' 'VC1' 'VC2' 'VC3' 'VC4' 'VC5' 'VC6' 'VD1' 'VD10' 'VD11' 'VD12' 'VD13' 'VD2' 'VD3' 'VD4' 'VD5' 'VD6' 'VD7' 'VD8' 'VD9'};


%% 6 %%
V{end+1} = 'Interneuron';
C{end+1} = {'ALA' 'AVAL' 'AVAR' 'AVBL' 'AVBR' 'AVDL' 'AVDR' 'AVEL' 'AVER' 'AVFL' 'AVFR' 'AVG' 'AVHL' 'AVHR' 'AVJL' 'AVJR' 'AVKL' 'AVKR' 'BDUL' 'BDUR' 'DVC' 'LUAL' 'LUAR' 'PVCL' 'PVCR' 'PVPL' 'PVPR' 'PVQL' 'PVQR' 'PVR' 'PVT' 'PVWL' 'PVWR' 'RIAL' 'RIAR' 'RIBL' 'RIBR' 'RICL' 'RICR' 'RID' 'RIFL' 'RIFR' 'RIGL' 'RIGR' 'RIH' 'RMGL' 'RMGR' 'SAADL' 'SAADR' 'SAAVL' 'SAAVR' 'SABD' 'SABVL' 'SABVR' 'SDQL' 'SDQR' 'SIADL' 'SIADR' 'SIAVL' 'SIAVR' 'SIBDL' 'SIBDR' 'SIBVL' 'SIBVR' 'URBL' 'URBR'};


%% 7 %%
V{end+1} = 'Principal cell';
C{end+1} = {'CANL' 'CANR' 'PLNL' 'PLNR'};


%% 8 %%
V{end+1} = 'Ring';
C{end+1} = {'AIAL' 'AIAR' 'AIML' 'AIMR' 'AINL' 'AINR' 'AIZL' 'AIZR' 'DVA' 'RIR' 'RIS'};


%% 9 %%
V{end+1} = 'Ring/Pharynx';
C{end+1} = {'RIPL' 'RIPR'};

%% 10 %%
V{end+1} = 'Sensory neuron';
C{end+1} = {'IL1DL' 'IL1DR' 'IL1L' 'IL1R' 'IL1VL' 'IL1VR' 'OLQDL' 'OLQDR' 'OLQVL' 'OLQVR' 'URXL' 'URXR' 'ADEL' 'ADER' 'ADFL' 'ADFR' 'ADLL' 'ADLR' 'AFDL' 'AFDR' 'ALML' 'ALMR' 'ALNL' 'ALNR' 'AQR' 'ASEL' 'ASER' 'ASGL' 'ASGR' 'ASHL' 'ASHR' 'ASIL' 'ASIR' 'ASJL' 'ASJR' 'ASKL' 'ASKR' 'AVM' 'AWAL' 'AWAR' 'AWBL' 'AWBR' 'AWCL' 'AWCR' 'BAGL' 'BAGR' 'CEPDL' 'CEPDR' 'CEPVL' 'CEPVR' 'FLPL' 'FLPR' 'IL2DL' 'IL2DR' 'IL2L' 'IL2R' 'IL2VL' 'IL2VR' 'OLLL' 'OLLR' 'PDA' 'PDEL' 'PDER' 'PHAL' 'PHAR' 'PHBL' 'PHBR' 'PHCL' 'PHCR' 'PLML' 'PLMR' 'PQR' 'PVDL' 'PVDR' 'PVM' 'URYDL' 'URYDR' 'URYVL' 'URYVR'};

for i=1:length(C)
    for k=1:length(C{i})
        if strcmp(name, C{i}{k})==1
            c = i; 
            return;
        end
    end
end

% error;
c = Inf; v = Inf;

end