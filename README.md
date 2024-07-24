# topology-estimation-mcp
This repository includes simulation code for network topology estimation with $`\ell^1`$ norm and MCP

[Paper]  
F. Matsuzaki and T. Ikeda, "[Sparse topology estimation for consensus network systems via minimax concave penalty](https://ieeexplore.ieee.org/document/10542337)," _IEEE Control Systems Letters (L-CSS)_, vol. 8, pp. 1012–1017, 2024.


## Description
### DataGenerate.m
MATLAB code for generating data of consensus networks that include probabilistic noise.
To solve the stochastic differential equation, the Euler-Maruyama scheme is implemented.

The generated data is stored in `Data/<sample_duration>.mat`

### L1.m
MATLAB code for estimating the topology with $`\ell^1`$ penalty.
The estimation result is stored in `EstimateResult/L1_<sample_duration>/L1_<sample_duration> <date>.mat`.

### MCP.m
MATLAB code for estimating the topology with MCP.
The estimation result is stored in `EstimateResult/MCP_<sample_duration>/MCP_<sample_duration> <date>.mat`.
