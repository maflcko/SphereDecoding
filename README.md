Sphere Decoder
==============
An implementation of the eight Sphere Decoder algorithms by A. Ghasemmehdi and E. Agrell, 2011 ("Faster Recursions in Sphere Decoding," IEEE Trans. Inf. Theory, vol. 57, no. 6, pp. 3530-3536).

Usage
-----
All algorithms take a lower triangular generator matrix __G__ or __G__^-1 and a vector __y__ as input. They return __u__ such that __u__ * __G__ is the closest lattice point to __y__.

For instance:

    >> u = decode5([-1.2 2.8], eye(2))
    u =
        -1     3

Support
-------
Currently this works with Matlab; Octave support is on the TODO-list.

License
-------
You can use all files I am the author of. The only condition is attributing them to me (Marco Falke). Other files (another author is specified) may be restricted in usage.
