function [red,blue,Fred,Fblu] = necklace(red0,blue0)
% Remove beads from a necklace so that so that no neighbors are same color.
% Keep only the rightmost bead from same color sequences.  
% The necklace should start and end with red beads.
% red0 and blue0 are are sortet arrays of initial distances from the left 
% end of the nexklace to the red and blue beads, respectively. 
% red and blue are the distances of the remaining beads.  All reamining
% beads keep their original position.

% If the input is the timing of peaks amd valleys in a zig-zag time series
% it will eliminate missing pairs where a peak or a valley is missing 

  red = red0(:)';
  blue = blue0(:)';
  % Eliminate red beads with no blue between them:
  [~,i] = sort([blue,red]);
% Display input:
%   for j = 1:numel(i)
%     if i(j)<=numel(blue)
%       fprintf(' b ');
%     else 
%       fprintf(' r ');
%     end
%   end
%   fprintf('\n');
  redpos = find(i>numel(blue));
  reddiff = diff(redpos);
  Fred = find(reddiff<2);
  if ~isempty(Fred)
    red(Fred) = [];
  end
  
  % Eliminate blue beads with no red between them
  blupos = find(i<=numel(blue));
  bludiff = diff(blupos);
  Fblu = find([bludiff<2,0]|blue<red(1)|blue>red(end));
  if ~isempty(Fblu)
    blue(Fblu) = [];
  end
% Display output:
%   [~,i] = sort([blue,red]);
%   for j = 1:numel(i)
%     if i(j)<=numel(blue)
%       fprintf(' b ');
%     else 
%       fprintf(' r ');
%     end
%   end
%   fprintf('\n');
end


