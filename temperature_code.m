function Tlist = temperature_code(file)
% Returns the temperature coding table used in experiment file
%  The temperature around the protein is often recorded as the number given
%  by digits 2 and 3 in the Status variable. Tlist has two columns.  Column
%  1 hold the temperature codes from the Status variable, Column 2 holds
%  the corresponding temperatures.

  % Lists that have been used
  lists{1} = [0:2:12; 4,6,10,14,16,20,23]';
  % For SBS-tester 850-808
  lists{2} = [0 2 3 4 6 8 16; 4.4,9.88,12.12,14.89,18.89,21.4,31.7]';
  % New Mexico files:
  lists{3} = [[0:2:8,12]; 6.2,11.1,14.9,18.0, 21.3, 26.6]';
  % Mail from Christian Jun 13, 2023:
  % For Tim's gift':
  lists{4} = [0 2 4 6 8 10 14 12 16 20 24 31; ...
    3.2 6.27 10.16 13.98 16.95 20.12 23.34 26.12 28.30 31.21 33.97 38.03]';
  

  single_T = [
    "20230821/eA_10°.txt",10; ...
    "20230821/eA_16°.txt",16; ...
    "20230821/eA_23°.txt",23; ...
    "20230821/eA_3°.txt",3];

  T = str2double(single_T(single_T(:,1)==file,2));
  if ~isempty(T)  % Read constant T for file from single_T
    Tlist = T;
    return
  end

  files{1} = [
    "02162022/jA.txt"
    "02162022/lA.txt"
    "02162022/mA.txt"
    "02162022/nA.txt"
    "02172022/oA.txt"
    "02182022/qA.txt"
    "02182022/rA.txt"
    "04082022/YA.txt"
    "04082022/jA.txt"
    "07082022/cA.txt"
    "07082022/eA.txt"];

  files{2} = [
    "20220718/bA.txt"   
    "20220718/cA.txt"
    "20220718/cB.txt"
    "20220718/dA.txt"
    "20220718/dB.txt"
    "20220718/dC.txt"
    "20220718/dD.txt"
    "20220718/dE.txt"
    "20220718/eA.txt"
    "20220718/fA.txt"
    ];

  files{3} = [        % New Mexico
    "20230721/JA.txt"
    "20230721/JB.txt"
    "20230721/JC.txt"
    "20230721/JD.txt"
    "20230721/JE.txt"
    "20230721/JF.txt"
    "20230721/JG.txt"
    "20230721/KA.txt"
    "20230721/KB.txt"
    "20230721/KC.txt"
    "20230721/cA.txt"
    "20230721/cB.txt"
    "20230721/dA.txt"
    "20230721/eA.txt"
    "20230721/fA.txt"
    "20230721/fB.txt"
    "20230721/fC.txt"
    "20230721/fD.txt"
    "20230721/fE.txt"
    "20230721/fF.txt"
    "20230721/fG.txt"
    "20230721/gA.txt"
    "20230721/gB.txt"
    "20230721/hA.txt"
    "20230721/iA.txt"
    "20230722/lA.txt"
    "20230722/lB.txt"
    "20230722/lC.txt"
    "20230722/lD.txt"
    "20230722/lE.txt"
    "20230722/mA.txt"
    "20230722/mB.txt"
    "20230722/mC.txt"
    "20230722/mD.txt"
    "20230722/nA.txt"
    "20230722/nB.txt"
    "20230722/nC.txt"
    "20230722/nD.txt"
    "20230722/nE.txt"
    "20230722/nF.txt"
    "20230722/nG.txt"
    "20230722/oA.txt"
    "20230722/oB.txt"
    "20230722/oC.txt"
    "20230722/oD.txt"
    "20230722/oE.txt"
    "20230722/oF.txt"
    "20230722/pA.txt"
    "20230722/pB.txt"
    "20230722/pC.txt"
    "20230724/qA.txt"
    "20230724/rA.txt"
    "20230725/SA.txt"
    "20230725/SB.txt"
    "20230725/SC.txt"
    "20230725/SD.txt"
    "20230725/SE.txt"
    "20230725/SF.txt"
    "20230725/SG.txt"
    "20230725/TA.txt"
    "20230725/TB.txt"
    "20230725/TC.txt"
    "20230725/TD.txt"
    "20230725/TE.txt"
    "20230725/TF.txt"
    "20230725/TG.txt"
    "20230725/UA.txt"
    "20230725/UB.txt"
    "20230725/UC.txt"
    "20230725/UD.txt"
    "20230725/UE.txt"
    "20230725/UF.txt"
    "20230725/UG.txt"
    "20230725/VA.txt"];

  files{4} = [
    "20230623/aA.txt"
    "20230623/bA.txt"
    "20230623/cA.txt"
    "20230623/dA.txt"
    "20230623/eA.txt"
    "20230623/fA.txt"
    "20230623/fB.txt"
    "20230717/MA.txt" 
    "20230717/aA.txt" 
    "20230717/bA.txt"
    "20230717/cA.txt"
    "20230717/cvA.txt"
    "20230717/dA.txt" 
    "20230717/eA.txt"
    "20230717/fA.txt" 
    "20230717/gA.txt"
    "20230717/hA.txt" 
    "20230717/iA.txt" 
    "20230717/jA.txt" 
    "20230717/kA.txt"
    "20230717/lA.txt"
    "20230717/nA.txt"
    "20230717/oA.txt"   
    "20220718/bA.txt"
    "20220718/cB.txt"
    "20220718/dB.txt"
    "20220718/dC.txt"
    "20220718/dD.txt"
    "20220718/dE.txt"  
    "20230720/bA.txt"
    "20230720/cA.txt"
    "20230720/dA.txt"
    "20230720/eA.txt"
    "20230720/fA.txt"
    "20230720/gA.txt"
    "20230720/gaA.txt"
    "20230720/gbA.txt"
    "20230720/ghA.txt"
    "20230720/gvA.txt"
    "20230720/hA.txt"
    "20230720/iA.txt"
    "20230720/jA.txt"
    "20230720/kA.txt"
    "20230720/lA.txt"
    "20230720/mA.txt"
    "20230720/nA.txt"
    "20230720/oA.txt"
    "20230801/abA.txt"
    "20230801/uA.txt"
    "20230801/vA.txt"
    "20230801/wA.txt"
    "20230801/xA.txt"
    "20230801/yA.txt"
    "20230801/zA.txt"     
    "20230821/cA.txt"
    "20230821/eB.txt"
    "20230913/aA.txt"
    "20230913/bA.txt"
    "20230913/cA.txt"
    "20230914/aA.txt"
    "20230914/bA.txt"
    "20230914/cA.txt"
    "20230914/dA.txt"
    "20230914/eA.txt"
    "20230914/fA.txt"
    "20230914/gA.txt"
    "20230914/hA.txt"
    "20230914/jA.txt"
    "20230914/kA.txt"
    "20230914/lA.txt"
    "20230914/uA.txt"
    "20230915/aA.txt"
    "20230915/bA.txt"
    "20230915/cA.txt"
    "20230915/dA.txt"	
% BSAfiles
    "20231107/aA.txt"
    "20231107/bA.txt"
    "20231107/cA.txt"
    "20231107/dA.txt"
    "20231107/eA.txt"
    "20231107/fA.txt"
    "20231107/gA.txt"
    "20231107/hA.txt"
    "20231127/aA.txt"
    "20231127/bA.txt"
    "20231127/cA.txt"
    "20231127/dA.txt"
    "20231127/eA.txt"
    "20231127/fA.txt"
    "20231128/gA.txt"
    "20231128/gB.txt"
    "20231128/hA.txt"
    "20231128/iA.txt"
    "20231128/jA.txt"
    "20231128/jB.txt"
    "20231128/kA.txt"
    "20231128/lA.txt"
    "20231129/mA.txt"
    "20231129/mB.txt"
    "20231129/mC.txt"
    "20231129/nA.txt"
    "20231129/oA.txt"
    "20231129/pA.txt"
    "20231129/qA.txt"
    "20231129/rA.txt"
    "20231129/sA.txt"
    "20231129/tA.txt"
    "20231129/uA.txt"
    
    ];

  
  Tlist = [];
  for k = 1:numel(files)
    if ismember(file,files{k})
      Tlist = lists{k};
      break;
    end
  end
