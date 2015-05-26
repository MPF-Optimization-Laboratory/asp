% Compile the OCT/MEX interfaces to the heap routines.
% These are required for l2maxline.
% $Id: as_setup.m 495 2009-09-29 00:48:00Z saunders $

cd private
try
   if exist('OCTAVE_VERSION','builtin')
      mkoctfile heapify.cc heap.c
      mkoctfile heapdel.cc heap.c
   else
      mex heapify.c heap.c
      mex heapdel.c heap.c
   end
   fprintf('\nInterface successfully compiled')
   fprintf('\nAvailable mex files are')
   dir heap*.m*
catch
   fprintf('\nInterface did not compile!')
end
cd ..
 