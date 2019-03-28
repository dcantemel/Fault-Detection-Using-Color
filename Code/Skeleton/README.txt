MATLAB Implementation of Skeletonization
------------------------------------------------------------

Written by Nicholas R. Howe
Packaged May 2006

This package comes with no warranty of any kind (see below).


Description
-----------

The files in this package comprise the Matlab implementation of a
method for skeletonizing binary images.  These are my own implementation
of a skeletonization method communicated to me by Alex Telea
(http://www.win.tue.nl/~alext/).  See his papers for more details.

This implementation is not optimized -- it is designed to be slow but
reliable.  Alex has a faster (non-Matlab) implementation, but I found
that it didn't give the results I wanted in certain situations.

Execute demo.m line by line while reading the comments to see what
is going on.


Note from the author: http://www.cs.smith.edu/~nhowe/research/code/

Note: I have no published paper on this technique. The idea came to me orally from Alex Telea, and I believe it may be contained in one of his papers, but I am not sure which one. However, I can express the concept fairly easily. A point on the skeleton sits at the center of a circle that touches the edge of the figure at multiple points. The intensity in the gray-level skeleton image is based on the shortest distance you would have to travel around the perimeter of the figure to connect the most distant two points. Thus spurs in the skeleton due to small edge perturbations will have low intensity even if they are very long. If the circle touches two disconnected edges (like the inside and the outside of an O shape) then the skeleton will be infinite at that point. To get a final skeleton, you should pick your threshold depending on the expected size of any noisy protrusions in the silhouette.

Copyright Notice
----------------

This software comes with no warranty, expressed or implied, including
but not limited to merchantability or fitness for any particular
purpose.

All files are copyright Nicholas R. Howe.  Permission is
granted to use the material for noncommercial and research purposes.
