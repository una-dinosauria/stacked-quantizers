Additive Quantization
---

This is an implementation of Additive Quantization by [Julieta](http://www.cs.ubc.ca/~julm/) as described in

Artem Babenko and Victor Lempitsky: *Additive Quantization for Extreme Vector Compression*, CVPR 14.

### TODOs

* The code for encoding under `AQ_encoding.m` can surely benefit from a C/C++ implementation.
* The codeboook update under `AQ_update_codebooks.m` runs a least-squares uptimization on sparse matrix of code assignments. My implementation currently creates a full matrix and then makes it sparse, which is quite readable but not very good for speed/memory.
