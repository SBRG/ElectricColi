clear

%Load solver - Gurobi offers free licenses for academic use
changeCobraSolver('gurobi','LP',0,-1);
changeCobraSolver('gurobi','QP',0,-1);
changeCobraSolver('gurobi','MILP',0,-1);
changeCobraSolver('gurobi','MIQP',0,-1);

%Change directory to that of the current script
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename))

load('mappedGenes.mat')
mappedGenes;
mappedBnumbers;

load('iML1515_simple.mat');
modelGenes = iML1515;

load('iML1515.mat')
modelIMLMod = iML1515;

solTest = optimizeCbModel(modelIMLMod);

%Convert the base model to anaerobic
%Anaerobic
modelIMLMod = changeRxnBounds(modelIMLMod,'EX_o2_e',0,'l');
solTest0 = optimizeCbModel(modelIMLMod,'max',10^-6);

%Try hard capping THD since it's being used to convert large amounts of
%NADH
modelIMLMod = changeRxnBounds(modelIMLMod,'THD2pp',2,'u');

%% Setting up the knockouts
%Set 1
%Note I changed ccp to its synonym yhjA
geneSet1 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW';'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA'};

bnumberSet1 = cell(length(geneSet1),1);
for i = 1:length(geneSet1)
    curInd = find(strcmp(geneSet1{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSet1{i} = mappedBnumbers{curInd};
    else
        bnumberSet1{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModel1,geneInd1] = ismember(bnumberSet1,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
downRegFraction = 0;
[modelRed1, hasEffect1, constrRxnNames1, deletedGenes1] = deleteModelGenes(modelIMLMod, bnumberSet1(isInModel1));
solTest1 = optimizeCbModel(modelRed1,'max',10^-6);


%Set 2
%Note I changed yieF to chrR
%Note I changed ygc to ygcBEGNOPQRSUW
%Note I changed ydi to ydiBEFHJKLMNOPQRSTUVYZ
%Note I changed fix to fixABCX
geneSet2 = {'azoR';'nfsA';'nfsB';'qorA';'qorB';'wrbA';'chrR';'kefF';'kefC';...
'kefG';'kefB';'nemA';'aegA';'mdaB';'ygiN';'btuE';'yqhD';...%'metF';'cysJ';'cysI'; %Cannot knock out these genes according to the model
%'cysH';... %Cannot knock out this gene either
%'ygcB';'ygcE';'ygcG'; %These three were not mentioned as knocked out
'ygcN';'ygcO';'ygcP';'ygcQ';'ygcR';'ygcS';'ygcU';'ygcW';...
'ydiB';'ydiE';'ydiF';'ydiH';'ydiJ';'ydiK';'ydiL';'ydiM';'ydiN';'ydiO';'ydiP';...
'ydiQ';'ydiR';'ydiS';'ydiT';...
%'ydiU';'ydiV';'ydiY';'ydiZ'; These were not mentioned as knocked out
'fixA';'fixB';'fixC';'fixX'};

bnumberSet2 = cell(length(geneSet2),1);
for i = 1:length(geneSet2)
    curInd = find(strcmp(geneSet2{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSet2{i} = mappedBnumbers{curInd};
    else
        bnumberSet2{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModel2,geneInd2] = ismember(bnumberSet2,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
[modelRed2, hasEffect2, constrRxnNames2, deletedGenes2] = deleteModelGenes(modelRed1, bnumberSet2(isInModel2));
solTest2 = optimizeCbModel(modelRed2,'max',10^-6);

%Individually test reactions that were carrying flux but knocked out
% rxnsProblematic = {'PAPSR';'SULR';'MTHFR2';'PAPSR2';'QMO3';'ACOAD1fr'};
% for i = 1:length(rxnsProblematic)
%     modelTest = changeRxnBounds(modelRed1,rxnsProblematic{i},0,'b');
%     solTemp = optimizeCbModel(modelTest,'max',10^-6);
%     rxnsProblematic{i};
%     solTemp.f;
% end

%SULR and MTHFR2 are reactions that we cannot delete even individually
%SULR is b2763 (cysI) AND b2764 (cysJ)
%MTHFR2 is b3941 (metF)

%Test reactions together that were carrying flux but knocked out
% rxnsProblematic2 = {'PAPSR';'PAPSR2';'QMO3';'ACOAD1fr'};
% modelTest = modelRed1;
% for i = 1:length(rxnsProblematic2)
%     modelTest = changeRxnBounds(modelTest,rxnsProblematic2{i},0,'b');
%     solTemp = optimizeCbModel(modelTest,'max',10^-6);
%     rxnsProblematic2{i}
%     solTemp.f
% end

%CysH also necessary for PAPSR and PAPSR2

%Set 3
%Assuming fdn this means fdnGHI
%Assuming fdo this means fdoGHI
%Assuming fdh means fdhDEF
%Assuming hyc means hycABCDEFGHI
%Assuming hyp means hypABCDEFT
%Added the alternate ethanol producing genes frmA and adhP
%Adding KOs for hya, hyb, hyf because cell can't produce cofactor anymore
%Adding dld to help knock out lactate secretion, as ldhA KO isn't enough
%Adding mhpF to help knock out ethanol secretion, as adhE isn't enough
%Also knocking out methylglyoxyl reductase to remove that as an NADH
%dumping pathway (b1771 (ydjG)) and 12ppd (b3945 (gldA))
%Adding mgsA after seeing it as another escape route
geneSet3 = {'ldhA';'dld';'adhE';'frmA';'adhP';'mhpF';...
    'fdnG';'fdnH';'fdnI';...
    'fdoG';'fdoH';'fdoI';...
    'fdhD';'fdhE';'fdhF';...
    'hycA';'hycB';'hycC';'hycD';'hycE';'hycF';'hycG';'hycH';'hycI';...
    'hypA';'hypB';'hypC';'hypD';'hypE';'hypF';'hypT';...
    'hyaA';'hyaB';'hyaC';'hyaD';'hyaE';'hyaF';...
    'hybO';'hybA';'hybB';'hybC';'hybD';'hybE';'hybF';'hybG';...
    'hyfA';'hyfB';'hyfC';'hyfD';'hyfE';'hyfF';'hyfG';'hyfH';'hyfI';'hyfJ';'hyfR';...
    'ydjG';'gldA';'mgsA'};

bnumberSet3 = cell(length(geneSet3),1);
for i = 1:length(geneSet3)
    curInd = find(strcmp(geneSet3{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSet3{i} = mappedBnumbers{curInd};
    else
        bnumberSet3{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModel3,geneInd3] = ismember(bnumberSet3,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
[modelRed3, hasEffect3, constrRxnNames3, deletedGenes3] = deleteModelGenes(modelRed2, bnumberSet3(isInModel3));
solTest3 = optimizeCbModel(modelRed3,'max',10^-6);

%Try hard capping THD since it's being used to convert large amounts of
%NADH
modelTHD = modelRed3;
modelTHD = changeRxnBounds(modelTHD,'THD2pp',2,'u');

solTest8 = optimizeCbModel(modelTHD,'max',10^-6);
solTest8SP = optimizeCbModel(modelTHD);

%% Check which reactions are supposed to be flagged for deletion and see if any have significant flux still (indicating escape by alternate GPR)

geneSetAll = [geneSet1;geneSet2;geneSet3];
bnumberSetAll = cell(length(geneSetAll),1);
for i = 1:length(geneSetAll)
    curInd = find(strcmp(geneSetAll{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetAll{i} = mappedBnumbers{curInd};
    else
        bnumberSetAll{i} = "not found";
    end
end

rxnSetAll = {};
for i =1:length(bnumberSetAll)
    curGene = bnumberSetAll{i};
    curGeneInd = find(strcmp(curGene,modelGenes.genes));
    curRxnList = modelGenes.rxns(find(modelGenes.rxnGeneMat(:,curGeneInd)));
    rxnSetAll = [rxnSetAll;curRxnList];
end
rxnSetAll = unique(rxnSetAll);

%% Testing bioreactor conditions

% %Make sure the ATP synthase is not going in reverse %Not necessary
% modelRed3 = changeRxnBounds(modelRed3,'ATPS4rpp',0,'l');

%Anaerobic
modelAnaerobic = changeRxnBounds(modelTHD,'EX_o2_e',0,'l');
%Increase glucose
modelAnaerobic = changeRxnBounds(modelAnaerobic,'EX_glc__D_e',-10,'l');
% %Remove the ATPM requirement
% modelAnaerobic = changeRxnBounds(modelAnaerobic,'ATPM',0,'l');
solTest4 = optimizeCbModel(modelAnaerobic,'max',10^-6);

%dhna_c is 1,4-Dihydroxy-2-naphthoate, which is one of the most similar
%forms to Biki's compound 2-Hydroxy-1,4-naphthoquinone

%Add the reactions and metabolite for exogenously supplied HNQ
% nadh_c + HNQ(ox)_c -> nad_c + HNQ(red)_c
% modelAnaerobic = addReaction(modelAnaerobic,'HNQR_new','metaboliteList',{'nadh[c]','nhqox[c]','h[c]','nad[c]','nhqred[c]'},'stoichCoeffList',[-1 -1 -1 1 1], 'reversible',false);
modelHNQ = addReaction(modelAnaerobic,'HNQR_new','metaboliteList',{'nadh[c]','hnqox[c]','nad[c]','hnqred[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false);
modelHNQ = addReaction(modelHNQ,'HNQanode_new','metaboliteList',{'hnqred[c]','hnqox[c]'},'stoichCoeffList',[-1 1], 'reversible',false);
solTest5 = optimizeCbModel(modelHNQ,'max',10^-6);
solTest5SP = optimizeCbModel(modelHNQ);

%Write model to map to Escher
% filename = 'C:\Users\dczie\Dropbox (Personal)\Projects (Collaborations)\Electric coli\modelMod.json';
% filename = 'E:\Dropbox (Personal)\Projects (Collaborations)\Electric coli\modelMod.json';
% writeCbModel(modelAnaerobic, 'fileName',filename);

%% Testing with the peptide supplementation

%Anaerobic
modelPeptide = changeRxnBounds(modelTHD,'EX_o2_e',0,'l');
%Increase glucose
modelPeptide = changeRxnBounds(modelPeptide,'EX_glc__D_e',-4,'l'); %Equal flux of glucose and total AA uptake (of course carbons count is different)
%Add amount of AA uptake
fluxAAUptake = -0.4;
listAA = {'glu__L';'gln__L';'asp__L';'asn__L';'leu__L';'ile__L';'val__L';...
    'ser__L';'gly';'ala__L';'tyr__L';'phe__L';'trp__L';'pro__L';...
    'arg__L';'his__L';'cys__L';'met__L';'lys__L';'thr__L'};
for i = 1:length(listAA)
    curEx = strcat("EX_",listAA{i},"_e");
    modelPeptide = changeRxnBounds(modelPeptide,curEx,fluxAAUptake,'l');
end

%Add the reactions and metabolite for exogenously supplied HNQ
% nadh_c + HNQ(ox)_c -> nad_c + HNQ(red)_c
% modelAnaerobic = addReaction(modelAnaerobic,'HNQR_new','metaboliteList',{'nadh[c]','nhqox[c]','h[c]','nad[c]','nhqred[c]'},'stoichCoeffList',[-1 -1 -1 1 1], 'reversible',false);
modelPeptide = addReaction(modelPeptide,'HNQR_new','metaboliteList',{'nadh[c]','hnqox[c]','nad[c]','hnqred[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false);
modelPeptide = addReaction(modelPeptide,'HNQanode_new','metaboliteList',{'hnqred[c]','hnqox[c]'},'stoichCoeffList',[-1 1], 'reversible',false);
solTestPep = optimizeCbModel(modelPeptide,'max',10^-6);
solTestPepSP = optimizeCbModel(modelPeptide);

%% Test a proton pumping alternate to the transporter

%Add proton pumping to the HNQ export (but not import) and see if that 
% %can power significant ATP synthase compared to substrate level phosphorylation

nProtonsPumped = 1;
modelHNQAlt = addReaction(modelAnaerobic,'HNQR_new','metaboliteList',{'nadh[c]','hnqox[c]','nad[c]','hnqred[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false);
modelHNQAlt = addReaction(modelHNQAlt,'HNQanode_new','metaboliteList',{'hnqred[c]','hnqox[c]','h[c]','h[p]'},'stoichCoeffList',[-1 1 -nProtonsPumped nProtonsPumped], 'reversible',false);
solTest5Alt = optimizeCbModel(modelHNQAlt,'max',10^-6);
solTest5SPAlt = optimizeCbModel(modelHNQAlt);


%% Nitrate respiring model

%Nitrate respiration with KOs
modelNitrateKO = changeRxnBounds(modelAnaerobic,'EX_no3_e',-20,'l');

solTestNitrateKO = optimizeCbModel(modelNitrateKO,'max',10^-6);
solTestNitrateKOSP = optimizeCbModel(modelNitrateKO);

%Nitrate respiration from WT
modelNitrate = changeRxnBounds(modelIMLMod,'EX_o2_e',0,'l');
modelNitrate = changeRxnBounds(modelNitrate,'EX_no3_e',-20,'l');

solTestNitrate = optimizeCbModel(modelNitrate,'max',10^-6);
solTestNitrateSP = optimizeCbModel(modelNitrate);

%Strain 1 eBK31
geneSetNO1 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW'};
bnumberSetNO1 = cell(length(geneSetNO1),1);
for i = 1:length(geneSetNO1)
    curInd = find(strcmp(geneSetNO1{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetNO1{i} = mappedBnumbers{curInd};
    else
        bnumberSetNO1{i} = "not found";
    end
end
%Check if genes are in the metabolic model
[isInModelNO1,geneIndNO1] = ismember(bnumberSetNO1,regexprep(modelIMLMod.genes,'_deleted',''));
%Check the deletion
[modelRedNO1, hasEffectNO1, constrRxnNamesNO1, deletedGenesNO1] = deleteModelGenes(modelNitrate, bnumberSetNO1(isInModelNO1));
solTestNO1 = optimizeCbModel(modelRedNO1,'max',10^-6);


%Strain 2 eBK6
geneSetNO2 = {'ldhA';'adhE';'mhpF';'frdA';'frdB';'frdC';'mgsA'};
bnumberSetNO2 = cell(length(geneSetNO2),1);
for i = 1:length(geneSetNO2)
    curInd = find(strcmp(geneSetNO2{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetNO2{i} = mappedBnumbers{curInd};
    else
        bnumberSetNO2{i} = "not found";
    end
end
%Check if genes are in the metabolic model
[isInModelNO2,geneIndNO2] = ismember(bnumberSetNO2,regexprep(modelIMLMod.genes,'_deleted',''));
%Check the deletion
[modelRedNO2, hasEffectNO2, constrRxnNamesNO2, deletedGenesNO2] = deleteModelGenes(modelNitrate, bnumberSetNO2(isInModelNO2));
solTestNO2 = optimizeCbModel(modelRedNO2,'max',10^-6);


%Strain 3 eBK6
geneSetNO3 = {'fdhD';'fdhE';'fdhF';'fdnG';'fdnH';'fdnI';...
    'fdoG';'fdoH';'fdoI'};
bnumberSetNO3 = cell(length(geneSetNO3),1);
for i = 1:length(geneSetNO3)
    curInd = find(strcmp(geneSetNO3{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetNO3{i} = mappedBnumbers{curInd};
    else
        bnumberSetNO3{i} = "not found";
    end
end
%Check if genes are in the metabolic model
[isInModelNO3,geneIndNO3] = ismember(bnumberSetNO3,regexprep(modelIMLMod.genes,'_deleted',''));
%Check the deletion
[modelRedNO3, hasEffectNO3, constrRxnNamesNO3, deletedGenesNO3] = deleteModelGenes(modelNitrate, bnumberSetNO3(isInModelNO3));
solTestNO3 = optimizeCbModel(modelRedNO3,'max',10^-6);

%KO everything but nar and nap and cap THD in the same way
geneSetNO4 = {'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA';...
'azoR';'nfsA';'nfsB';'qorA';'qorB';'wrbA';'chrR';'kefF';'kefC';...
'kefG';'kefB';'nemA';'aegA';'mdaB';'ygiN';'btuE';'yqhD';...%'metF';'cysJ';'cysI'; %Cannot knock out these genes according to the model
%'cysH';... %Cannot knock out this gene either
%'ygcB';'ygcE';'ygcG'; %These three were not mentioned as knocked out
'ygcN';'ygcO';'ygcP';'ygcQ';'ygcR';'ygcS';'ygcU';'ygcW';...
'ydiB';'ydiE';'ydiF';'ydiH';'ydiJ';'ydiK';'ydiL';'ydiM';'ydiN';'ydiO';'ydiP';...
'ydiQ';'ydiR';'ydiS';'ydiT';...
%'ydiU';'ydiV';'ydiY';'ydiZ'; These were not mentioned as knocked out
'fixA';'fixB';'fixC';'fixX';...
'ldhA';'dld';'adhE';'frmA';'adhP';'mhpF';...
    'fdnG';'fdnH';'fdnI';...
    'fdoG';'fdoH';'fdoI';...
    'fdhD';'fdhE';'fdhF';...
    'hycA';'hycB';'hycC';'hycD';'hycE';'hycF';'hycG';'hycH';'hycI';...
    'hypA';'hypB';'hypC';'hypD';'hypE';'hypF';'hypT';...
    'hyaA';'hyaB';'hyaC';'hyaD';'hyaE';'hyaF';...
    'hybO';'hybA';'hybB';'hybC';'hybD';'hybE';'hybF';'hybG';...
    'hyfA';'hyfB';'hyfC';'hyfD';'hyfE';'hyfF';'hyfG';'hyfH';'hyfI';'hyfJ';'hyfR';...
    'ydjG';'gldA';'mgsA'};
bnumberSetNO4 = cell(length(geneSetNO4),1);
for i = 1:length(geneSetNO4)
    curInd = find(strcmp(geneSetNO4{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetNO4{i} = mappedBnumbers{curInd};
    else
        bnumberSetNO4{i} = "not found";
    end
end
%Check if genes are in the metabolic model
[isInModelNO4,geneIndNO4] = ismember(bnumberSetNO4,regexprep(modelIMLMod.genes,'_deleted',''));
%Check the deletion
[modelRedNO4, hasEffectNO4, constrRxnNamesNO4, deletedGenesNO4] = deleteModelGenes(modelNitrate, bnumberSetNO4(isInModelNO4));
%Cap THD as well
modelRedNO4 = changeRxnBounds(modelRedNO4,'THD2pp',2,'u');
solTestNO4 = optimizeCbModel(modelRedNO4,'max',10^-6);

%KO everything including nar and nap and cap THD in the same way
geneSetNO5 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW';'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA';...
'azoR';'nfsA';'nfsB';'qorA';'qorB';'wrbA';'chrR';'kefF';'kefC';...
'kefG';'kefB';'nemA';'aegA';'mdaB';'ygiN';'btuE';'yqhD';...%'metF';'cysJ';'cysI'; %Cannot knock out these genes according to the model
%'cysH';... %Cannot knock out this gene either
%'ygcB';'ygcE';'ygcG'; %These three were not mentioned as knocked out
'ygcN';'ygcO';'ygcP';'ygcQ';'ygcR';'ygcS';'ygcU';'ygcW';...
'ydiB';'ydiE';'ydiF';'ydiH';'ydiJ';'ydiK';'ydiL';'ydiM';'ydiN';'ydiO';'ydiP';...
'ydiQ';'ydiR';'ydiS';'ydiT';...
%'ydiU';'ydiV';'ydiY';'ydiZ'; These were not mentioned as knocked out
'fixA';'fixB';'fixC';'fixX';...
'ldhA';'dld';'adhE';'frmA';'adhP';'mhpF';...
    'fdnG';'fdnH';'fdnI';...
    'fdoG';'fdoH';'fdoI';...
    'fdhD';'fdhE';'fdhF';...
    'hycA';'hycB';'hycC';'hycD';'hycE';'hycF';'hycG';'hycH';'hycI';...
    'hypA';'hypB';'hypC';'hypD';'hypE';'hypF';'hypT';...
    'hyaA';'hyaB';'hyaC';'hyaD';'hyaE';'hyaF';...
    'hybO';'hybA';'hybB';'hybC';'hybD';'hybE';'hybF';'hybG';...
    'hyfA';'hyfB';'hyfC';'hyfD';'hyfE';'hyfF';'hyfG';'hyfH';'hyfI';'hyfJ';'hyfR';...
    'ydjG';'gldA';'mgsA'};
bnumberSetNO5 = cell(length(geneSetNO5),1);
for i = 1:length(geneSetNO5)
    curInd = find(strcmp(geneSetNO5{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetNO5{i} = mappedBnumbers{curInd};
    else
        bnumberSetNO5{i} = "not found";
    end
end
%Check if genes are in the metabolic model
[isInModelNO5,geneIndNO5] = ismember(bnumberSetNO5,regexprep(modelIMLMod.genes,'_deleted',''));
%Check the deletion
[modelRedNO5, hasEffectNO5, constrRxnNamesNO5, deletedGenesNO5] = deleteModelGenes(modelNitrate, bnumberSetNO5(isInModelNO5));
%Cap THD as well
modelRedNO5 = changeRxnBounds(modelRedNO5,'THD2pp',2,'u');
solTestNO5 = optimizeCbModel(modelRedNO5,'max',10^-6);

%% Other knockout sets for E. coli

%AnoxicNull -
geneSetANull1 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW';'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA'};

geneSetANullAll = geneSetANull1;
bnumberSetANull = cell(length(geneSetANullAll),1);
for i = 1:length(geneSetANullAll)
    curInd = find(strcmp(geneSetANullAll{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetANull{i} = mappedBnumbers{curInd};
    else
        bnumberSetANull{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModelANull,geneIndANull] = ismember(bnumberSetANull,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
downRegFraction = 0;
[modelRedANull, hasEffectANull, constrRxnNamesANull, deletedGenesANull] = deleteModelGenes(modelIMLMod, bnumberSetANull(isInModelANull));
solTestANull = optimizeCbModel(modelRedANull,'max',10^-6);


%QRedNull -
geneSetQNull1 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW';'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA'};
geneSetQNull2 = {'azoR';'nfsA';'nfsB';'qorA';'qorB';'wrbA';'chrR';'kefF';'kefC';...
'kefG';'kefB';'nemA';'aegA';'mdaB';'ygiN';'btuE';'yqhD';...%'metF';'cysJ';'cysI'; %Cannot knock out these genes according to the model
%'cysH';... %Cannot knock out this gene either
%'ygcB';'ygcE';'ygcG'; %These three were not mentioned as knocked out
'ygcN';'ygcO';'ygcP';'ygcQ';'ygcR';'ygcS';'ygcU';'ygcW';...
'ydiB';'ydiE';'ydiF';'ydiH';'ydiJ';'ydiK';'ydiL';'ydiM';'ydiN';'ydiO';'ydiP';...
'ydiQ';'ydiR';'ydiS';'ydiT';...
%'ydiU';'ydiV';'ydiY';'ydiZ'; These were not mentioned as knocked out
'fixA';'fixB';'fixC';'fixX'};

geneSetQNullAll = [geneSetQNull1;geneSetQNull2];
bnumberSetQNull = cell(length(geneSetQNullAll),1);
for i = 1:length(geneSetQNullAll)
    curInd = find(strcmp(geneSetQNullAll{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetQNull{i} = mappedBnumbers{curInd};
    else
        bnumberSetQNull{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModelQNull,geneIndQNull] = ismember(bnumberSetQNull,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
downRegFraction = 0;
[modelRedQNull, hasEffectQNull, constrRxnNamesQNull, deletedGenesQNull] = deleteModelGenes(modelIMLMod, bnumberSetQNull(isInModelQNull));
solTestQNull = optimizeCbModel(modelRedQNull,'max',10^-6);


%EsinkNull -
geneSetENull1 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW';'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA'};
geneSetENull2 = {'azoR';'nfsA';'nfsB';'qorA';'qorB';'wrbA';'chrR';'kefF';'kefC';...
'kefG';'kefB';'nemA';'aegA';'mdaB';'ygiN';'btuE';'yqhD';...%'metF';'cysJ';'cysI'; %Cannot knock out these genes according to the model
%'cysH';... %Cannot knock out this gene either
%'ygcB';'ygcE';'ygcG'; %These three were not mentioned as knocked out
'ygcN';'ygcO';'ygcP';'ygcQ';'ygcR';'ygcS';'ygcU';'ygcW';...
'ydiB';'ydiE';'ydiF';'ydiH';'ydiJ';'ydiK';'ydiL';'ydiM';'ydiN';'ydiO';'ydiP';...
'ydiQ';'ydiR';'ydiS';'ydiT';...
%'ydiU';'ydiV';'ydiY';'ydiZ'; These were not mentioned as knocked out
'fixA';'fixB';'fixC';'fixX'};
geneSetENull3 = {'ldhA';'dld';'adhE';'frmA';'adhP';'mhpF';...
    'fdnG';'fdnH';'fdnI';...
    'fdoG';'fdoH';'fdoI';...
    'fdhD';'fdhE';'fdhF';...
    'hycA';'hycB';'hycC';'hycD';'hycE';'hycF';'hycG';'hycH';'hycI';...
    'hypA';'hypB';'hypC';'hypD';'hypE';'hypF';'hypT';...
    'hyaA';'hyaB';'hyaC';'hyaD';'hyaE';'hyaF';...
    'hybO';'hybA';'hybB';'hybC';'hybD';'hybE';'hybF';'hybG';...
    'hyfA';'hyfB';'hyfC';'hyfD';'hyfE';'hyfF';'hyfG';'hyfH';'hyfI';'hyfJ';'hyfR';...
    'ydjG';'gldA';'mgsA'};

geneSetENullAll = [geneSetENull1;geneSetENull2;geneSetENull3];
bnumberSetENull = cell(length(geneSetENullAll),1);
for i = 1:length(geneSetENullAll)
    curInd = find(strcmp(geneSetENullAll{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetENull{i} = mappedBnumbers{curInd};
    else
        bnumberSetENull{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModelENull,geneIndENull] = ismember(bnumberSetENull,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
downRegFraction = 0;
[modelRedENull, hasEffectENull, constrRxnNamesENull, deletedGenesENull] = deleteModelGenes(modelIMLMod, bnumberSetENull(isInModelENull));
solTestENull = optimizeCbModel(modelRedENull,'max',10^-6);





%FermNull -
geneSetFNull1 = {'napF';'napD';'napA';'napG';'napH';'napB';'napC';'narG';'narH';...
'narI';'narZ';'narY';'narW';'dmsA';'dmsB';'dmsC';'ynfE';'ynfF';'ynfG';'ynfH';...
'torC';'torA';'torD';'torY';'torZ';'frdA';'frdB';'frdC';'nirB';'nirD';'nirC';...
'nrfA';'nrfB';'nrfC';'nrfD';'nrfE';'nrfF';'nrfG';'norV';'norW';'hmp';'yhjA'};
geneSetFNull2 = {'ldhA';'dld';'adhE';'frmA';'adhP';'mhpF';...
    'ydjG';'gldA';'mgsA'};
%The following set is from finding flux through these reactions
geneSetFNull3 = {'fadA';'fadB';'fadI';'fadJ';'fre'};

geneSetFNullAll = [geneSetFNull1;geneSetFNull2;geneSetFNull3];
bnumberSetFNull = cell(length(geneSetFNullAll),1);
for i = 1:length(geneSetFNullAll)
    curInd = find(strcmp(geneSetFNullAll{i},mappedGenes));
    if ~isempty(curInd)
        bnumberSetFNull{i} = mappedBnumbers{curInd};
    else
        bnumberSetFNull{i} = "not found";
    end
end

%Check if genes are in the metabolic model
[isInModelFNull,geneIndFNull] = ismember(bnumberSetFNull,regexprep(modelIMLMod.genes,'_deleted',''));

%Check the deletion
downRegFraction = 0;
[modelRedFNull, hasEffectFNull, constrRxnNamesFNull, deletedGenesFNull] = deleteModelGenes(modelIMLMod, bnumberSetFNull(isInModelFNull));
solTestFNull = optimizeCbModel(modelRedFNull,'max',10^-6);

modelFNullHNQ = addReaction(modelRedFNull,'HNQR_new','metaboliteList',{'nadh[c]','hnqox[c]','nad[c]','hnqred[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false);
modelFNullHNQ = addReaction(modelFNullHNQ,'HNQanode_new','metaboliteList',{'hnqred[c]','hnqox[c]'},'stoichCoeffList',[-1 1], 'reversible',false);

solTestFNullHNQ = optimizeCbModel(modelFNullHNQ,'max',10^-6);

%% Counting redox reactions

metsRedox = {'nad[c]','nadh[c]', 'nadp[c]', 'nadph[c]','q8[c]','q8h2[c]','mql8[c]','mqn8[c]','2dmmql8[c]','2dmmq8[c]','trdox[c]','trdrd[c]','gthrd[c]', 'gthox[c]'};
metIndsRedox = zeros(length(metsRedox),1);
for i = 1:length(metsRedox)
    curInd = find(strcmp(metsRedox{i},modelIMLMod.mets));
    metIndsRedox(i) = curInd;
end

SBin = full(logical(modelIMLMod.S));

rxnsRedox = {};
for i = 1:length(metIndsRedox)
    curMetInd = metIndsRedox(i);
    curRxnInds = find(SBin(curMetInd,:));
    rxnsRedox = union(rxnsRedox,modelIMLMod.rxns(curRxnInds));
end











