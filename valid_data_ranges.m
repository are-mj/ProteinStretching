function r = valid_data_ranges(file)
% Start and end indices for the parts of input data to be analysed.
%   data index normally correspond to line number in the file minus 4.
%   r: n by 2 array
%
% Used by analyse_file.m
%
% To specify parts in a new file, use read_experiment_file.m 
%   and plot f.  Use the zoom and data tip tools to specify start
%   and end indices.  If no r is specified for a given file, the whole 
%   data set is used by analyse_file

  file = char(file);  % Convert if string
  file = strrep(file,'\','/');
  s = strfind(file,'/');
  if numel(s)>1
    file = file(s(end-1)+1:end);  % Short form of file name
  end

  r = [];

  switch file
    case '02032022/bA.txt'
      r = [   64   6011
            7011   8849
            8800  11201];
    case '02032022/eA.txt'
      r = [   60  49100];
    case '02112022/oA.txt'
      r = [15651 18550];
    case '02112022/rA.txt'
      r = [2280  49900];
    case '02142022/bA.txt'
      r = [   277 424401
           470887 531743];
    case '02152022/eA.txt'
      r = [25          103232
           104693      107866
          109125      112550
          114061      127006
          130532      141966];
    case '02172022/pA.txt'
      r =  [124   17792
            18150 45438];   
    case '04082022/jA.txt'
      r = [     1 153833
           168100 175073
           176200 250000];     
    case '04082022/uA.txt'
      r = [1300    60160
           61500  174500];      
    case '07012022/zA.txt'
      r = [
        39       74982
       75139      111905
      111996      389519
      390307      454518
      390106      454472
      455049      527258
      527488      573551      
      574113      806401
];
    case '07012022/zB.txt'
      r = [
         347       51691
       52076       75898
       76360      115138
      115215      117450 ];
    case '07012022/zAB.txt'
      r = [
           1       74982
       75349      111905
      111996      385701
      390000      455518
      455623      528000
      528025      574190
      574195      858801
      859186      883008
      883470      921003]; 
    case '07082022/eA.txt'
      r = [3659      635222];   
    case '07082022/cA.txt'
      r = [
         102      110312
      110272      141405
      141424      184163
      185531      193999
      194032      205718
      205771      320851
      321065      328015
      328313      338845
      338968      347779
      347890      367793 ];  
    case '07122022/cA.txt'
        r = [
         338      125469          
      125594      204008 
      206835      448216]; 
    case '07142022/aA.txt'
      r = [
      1586      197178
      197238      269031
      271332      353365];
    case '07142022/aAB.txt'
      r = [
      1         197100
      197200    268810
      270000    352000];
    case '07152022/bAB.txt'
      r = [
      2331      130822
      131065    283423
      283536    938686];
    case '07182022/bA.txt'
      r = [64       61427];
    case '07182022/cAB.txt'
      r = [222      234858];
    case '07182022/dAE.txt'
      r = [1        656143 
          656672   770084];
    case '07182022/eA.txt'
      r = [178     112233];
    case '07182022/fA.txt'
      r = [89      174225];
    case 'Molecula qA/250322qA.txt'
      r = [1    288361
        290791  382894
        387714  426644
        442180  521289
        522839  662561];
    case 'yB/250322yB.txt'
      r = [1    259339
        289217  293722
        295618  333449];
    case '20230721/JAG.txt'
      r = [1    205853
        206149  1052170
        1052360 1058510
        1059630 1378960];
    case '20230721/JF.txt'
      r = [1    13502
        13529   20000
        20101  206889];     
    case '20230721/KA.txt'
      r = [1276 39750
        40000   215700];
    case '20230721/fA.txt'
      r = [1    45242
        48383   48465
        48500   131472
        131500  140629
        140800  188980];
    case '20230722/oA.txt'
      r = [5353 23155
        23216   26955
        27000   196000];
    case '20230722/oC.txt'
      r = [1    42634
        58105   76794
        76900   82104
        82255   192000];
    case '20230725/SF.txt'
      r = [1    145876
        146000  150795
        150824  192598];
     case '20230725/TA.txt'
       r = [1   2473
         2918   200753];
      case '20230725/TG.txt'
       r = [1   57056
         57244  60646];
       case '20230725/UC.txt'
         r = [1 32480
           32840  191856];
end