%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [runtime] = compute_WTM_map(T,compile,filename,AT_filename)
%
% input: 
%
%        T = the WTM threshold
%        compile: 1 = compile c++ code (required first time)
%                 0 = don't compile c++ code
%        filename = name of file that contains edge list
%        AT_filename = name  of file wherein the activation times will be saved   
%
% output: runtime = the time required for the code to run
% 
% DRT May 13, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [runtime] = compute_WTM_map(T,compile,filename,AT_filename)

%% Possibly compile the c++ code


%this code creates a makefile to comile c++ code   

%you will need to modify this code to fit your computer's c++ compiler
if compile
     
   fid = fopen('Makefile', 'w');           
   fprintf(fid,'all: compute_WTM_map_v3.cpp\n');
   fprintf(fid,'\tg++ -g -Wall -o compute_WTM_map_v3 compute_WTM_map_v3.cpp\n');
   fprintf(fid,'clean:\n');
   fprintf(fid,'\t$(RM) compute_WTM_map_v3');

   %compile the c++ code 
   %cd WTM_map
   unix('make');
   %cd ..
   
end


%% Run the c++ code and pass it the input and output file names

   %cd WTM_map%set current folder as WTM_map
   tic
   unix(['./compute_WTM_map_v3 -t ',num2str(T),...
      ' -i ',filename,...
      ' -o ',AT_filename]);
   runtime = toc;% calculate run time for creating the WTM map
   %cd ..%return out of the current folder, WTM_map

end









