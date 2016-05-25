COMP3320 Assignment 2
Tim Crundall

In order to make the sequential optimised version run:
make kernel_main_opt

In order to make the SSE vectorised version run:
make kernel_main_sse

note that the SSE version does not run correctly. The values of the forces
and potential energies quickly blow up to large numbers. The incorrectness
has been overlooked since it was beleived that useful, relevant peformance data
could be gathered with this incorrect code.


The OpenMP part has been omitted due to lack of time.
