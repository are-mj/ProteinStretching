function tbl = select(tbl,criteria)
% Select rows in Matlab table that satisfy one or more criteria 
% combining column names and relational operators
%   relational operators: '<', '>' or '='
%     NOTE that Matlab's standard equality operator '==' is replaced by '='
%     (out of covenience for the programmer!)
%   Single relations have the form <column><relation><value>
%   Combine lower and upper limits with the form: lowlim<column<hilim
% Two or more relations may be combined with '&'.
% Exanple: 
%  T2 = select(Tun,'Pullingspeed>200 & 9<Temperature<15');

  crits = regexp(criteria,'&','split');
  for ii = 1:numel(crits)
    crit = strrep(crits{ii},' ','');  % remove blanks
    pos = regexp(crit,'[<>=]');  % operator position
    switch length(pos)
      case 1
        column = crit(1:pos-1);
        operator = crit(pos);
        value = str2double(crit(pos+1:end));
        tbl = modify(tbl,column,operator,value,1);
      case 2
        column = crit(pos(1)+1:pos(2)-1);
        operator = crit(pos(1));
        value = str2double(crit(1:pos(1)-1));
        tbl = modify(tbl,column,operator,value,-1);  

        operator = crit(pos(2));
        value = str2double(crit(pos(2)+1:end));
        tbl = modify(tbl,column,operator,value,1);      
      otherwise
        error('Illegal criteria expression')
    end
  end
end
  
function newtbl = modify(tbl,column,operator,value,direction)
  if direction < 0
    if operator == '>'
      operator = '<';
    elseif operator == '<'
      operator = '>';
    end
  end
    switch operator
      case '='
        newtbl = tbl(tbl.(column)==value,:);
      case '<'
        newtbl = tbl(tbl.(column) < value,:);
      case '>'
        newtbl = tbl(tbl.(column) > value,:);
      otherwise
        error(['Unsupported operator in ',crit])
    end
end

  


  