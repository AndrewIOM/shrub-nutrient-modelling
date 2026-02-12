### Dependent libraries

This project uses a slightly modified version of DiffSharp from the current (Jan 2026)
nuget package. This is because the modified version drastically speeds up casting
from a boolean tensor to a float, by changing the ordering of the discriminated union
cases in the tensor creation method. This subverts bottlenecks in bisection.

The modified DiffSharp (from AndrewIOM/DiffSharp) also has apple silicon support for
Torch, although testing found that this was slower than the reference backend.

