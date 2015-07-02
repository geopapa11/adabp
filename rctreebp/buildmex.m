
% % // WARNING! WARNING! WARNING!  This MEX file is not for the faint of heart.
% % //  It attempts to circumvent Matlab's memory-handling functions in order
% % //  to give more efficient code, and uses undocumented functions to do so.
% % //  It comes with no guarantees whatsoever.  It may not work on other
% % //  versions of Matlab, or other platforms, or even at all.  It may result
% % //  in memory leaks, segmentation faults, destroy data, set your computer on
% % //  fire, or suck it into a black hole, for all I know.  If you do not accept
% % //  these possible risks, do not use this MEX file.


fprintf('Warning!!! These MEX files use undocumented functions to circumvent Matlabs\n');
fprintf('  memory management system.  It comes with absolutely no guarantees whatsoever.\n');
fprintf('  See the C++ files for more dire admonitions.  If you are uncomfortable with\n');
fprintf('  this approach, do not use these MEX files.\n');
fprintf('Otherwise, hit any key to continue.\n');
pause;

mex sUnion32.cpp
mex sMember32.cpp
mex sDiff32.cpp
mex modCell.cpp
mex modU32Array.cpp
