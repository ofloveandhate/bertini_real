This example is for showing how to do optimization of a non-algebraic function over a 1-dimensional algebraic set using Bertini_real and Chebfun.

1. dependencies.
  * make sure that `brakelab` is installed, and added to your Matlab path.
  * make sure that chebfun is downloaded and on your Matlab path.
  * make sure that `bertini` is on your path for making system calls.  [Help here](https://de.mathworks.com/matlabcentral/answers/10451-how-do-i-add-to-the-path-of-system).  Probably you need to modify your `startup.m` to include `/usr/local/bin` or something.
2. decompose the curve.  `bertini`, then `bertini_real`.  no need to run the `sampler`.
3. in matlab, `gather_br_samples`
4. run `demoplot`.
