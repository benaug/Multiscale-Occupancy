# Multiscale-Occupancy
Basic multiscale occupancy model that samples latent states from marginal distributions

This repository has multiscale occupancy models that sample the z_i and w_{i,j} latent indicator variables by computing z|w,y and sampling z' from that, then sampling w' from w|z',y. This is approach may have most or all the benefit of fitting the model by marginalizing z and w, but it is faster because you only do these computations once per iteration instead of once per parameter per iteration. This strategy could work for other occupancy models, or models with independent latent indicator variables.

There are 2 versions, one that uses 2D data (no replicate or occasion-level variation in p) and another that used 3D data. The former is faster if p does not vary across replicates/occasions, and more so the more there are.

I suspect that sampling z|w,y and w|z,y as nimble will do by default will sometimes not let you fully explore the posterior, but I haven't found a case of this, yet. If so, this approach would most likely prevent that (as would typical marginalization). I have seen this behavior in more complex multiscale occupancy models with latent species IDs.