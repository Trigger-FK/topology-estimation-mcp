# topology-estimation-mcp
Test code for network topology estimation with $`\ell^1`$ norm and MCP

## Description
### DataGenerate.m
MATLAB code for generating data of consensus networks that include probabilistic noise.
To solve the stochastic differential equation, the Euler-Maruyama sheme is implemented.

The generated data is stored in `Data/<sample_duration>.mat`

### L1.m
MATLAB code for estimating the topology with $`\ell^1`$ penalty.
The estimation result is stored in `EstimateResult/L1_<smaple_duration>/L1_<smaple_duration> <date>.mat`.

### MCP.m
MATLAB code for estimating the topology with MCP.
The estimation result is stored in `EstimateResult/MCP_<smaple_duration>/MCP_<smaple_duration> <date>.mat`.
