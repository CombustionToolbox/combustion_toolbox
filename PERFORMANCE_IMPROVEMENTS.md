# Performance Improvements for Combustion Toolbox

This document outlines the performance improvements made to the Combustion Toolbox and provides recommendations for future optimizations.

## Summary of Changes

This optimization effort focused on identifying and addressing performance bottlenecks in the codebase, particularly in hot paths that are executed repeatedly during equilibrium calculations.

### Changes Implemented

#### 1. Caching Optimizations (High Impact)

**Problem**: Species lookup in thermodynamic calculations used `find(strcmp())` which is O(n) complexity.

**Solution**: Replaced with `containers.Map` for O(1) lookup complexity.

**Files Modified**:
- `+combustiontoolbox/+utils/+thermo/getGibbsEnergyArray.m`
- `+combustiontoolbox/+utils/+thermo/getEnthalpyArray.m`
- `+combustiontoolbox/+core/@ChemicalSystem/updatePropertiesMatrixThermo.m`

**Impact**: Significant speedup for calculations involving repeated species lookups. For a typical equilibrium calculation with 50 species called 1000 times, this reduces lookup time from O(50,000) to O(1000) operations.

**Before**:
```matlab
index = find(strcmp(cachedSpecies, species), 1);
```

**After**:
```matlab
persistent cachedMap
if isempty(cachedMap)
    cachedMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
end

if isKey(cachedMap, species)
    index = cachedMap(species);
else
    % Cache new species
    cachedMap(species) = index;
end
```

#### 2. Array Preallocation (Medium Impact)

**Problem**: `set_prop_DB.m` was growing arrays backward without preallocation.

**Solution**: Preallocate array with known size before the loop.

**Files Modified**:
- `utils/databases/set_prop_DB.m`

**Impact**: Reduced memory reallocation overhead during database property extraction.

**Before**:
```matlab
for i = length(listSpecies):-1:1
    value(1, i) = DB.(species).(property);
end
```

**After**:
```matlab
numSpecies = length(listSpecies);
value = zeros(1, numSpecies);
for i = 1:numSpecies
    value(i) = DB.(species).(property);
end
```

#### 3. Vector Operation Optimization (Medium Impact)

**Problem**: Redundant vector multiplications in equilibrium solver inner loops.

**Solution**: Compute intermediate products once and reuse.

**Files Modified**:
- `+combustiontoolbox/+equilibrium/@EquilibriumSolver/equilibriumGibbs.m`
- `+combustiontoolbox/+equilibrium/@EquilibriumSolver/equilibriumHelmholtz.m`

**Impact**: Reduces computational operations in the most frequently executed function during equilibrium calculations.

**Before**:
```matlab
b1 = (NatomE - bi + sum(A0(indexGas, :) .* N(indexGas) .* muRT(indexGas)))';
b3 = NP + sum(N(indexGas) .* muRT(indexGas) - N(indexGas));
```

**After**:
```matlab
N_gas = N(indexGas);
muRT_gas = muRT(indexGas);
N_gas_muRT = N_gas .* muRT_gas;  % Compute once

b1 = (NatomE - bi + sum(A0(indexGas, :) .* N_gas_muRT))';
b3 = NP + sum(N_gas_muRT - N_gas);
```

#### 4. Code Clarity (Low Impact)

**Problem**: Minor inefficiency in `findIndex.m` logic.

**Solution**: Improved code organization and added clarifying comments.

**Files Modified**:
- `+combustiontoolbox/+utils/findIndex.m`

**Impact**: Marginal performance improvement with better code readability.

## Performance Testing

To measure the impact of these changes, run the validation suite with timing:

```matlab
% Run basic performance test
results = unitTest().run;

% Run computation time validation
run_computation_time('DET', 'C2H2_acetylene', LS, 9, 1, 3);
```

## Additional Optimization Opportunities

The following areas could benefit from future optimization work but were not modified in this iteration to maintain minimal changes:

### 1. Matrix Operations in Equilibrium Solver

**Location**: `equilibriumGibbs.m`, lines 383-411

**Opportunity**: The Jacobian matrix construction could potentially be optimized by:
- Caching matrix structure when species composition doesn't change significantly
- Using sparse matrix operations where appropriate
- Pre-computing constant submatrices

**Estimated Impact**: Low to Medium (depends on problem size)

### 2. Species Finding Algorithm

**Location**: `+core/@ChemicalSystem/findProducts.m`

**Opportunity**: The species search algorithm iterates through the entire database. Could be optimized with:
- Building index structures for common search patterns
- Caching results of common reactant combinations

**Estimated Impact**: Medium (primarily affects initialization time)

### 3. Vectorization Opportunities

**Location**: Multiple files using element-wise operations

**Opportunity**: Some loops could potentially be vectorized:
- Species property calculations across multiple temperatures
- Batch processing of multiple equilibrium calculations

**Estimated Impact**: Medium (depends on use case)

### 4. Parallel Processing

**Location**: Problem-solving loops in validation scripts

**Opportunity**: Multiple independent equilibrium calculations could be parallelized using MATLAB's Parallel Computing Toolbox:
- Parametric studies (e.g., varying equivalence ratios)
- Multiple initial conditions
- Rocket performance calculations

**Estimated Impact**: High (for parametric studies)

**Implementation Consideration**: Would require `parfor` loop modifications and ensuring thread-safety of caching mechanisms.

## Best Practices for Performance

When adding new code to the Combustion Toolbox, follow these guidelines:

### 1. Preallocate Arrays
```matlab
% Good
n = 1000;
result = zeros(n, 1);
for i = 1:n
    result(i) = compute(i);
end

% Avoid
for i = 1:n
    result(i) = compute(i);  % Growing array
end
```

### 2. Use Efficient Lookups
```matlab
% Good - O(1) lookup
speciesMap = containers.Map(speciesList, 1:length(speciesList));
index = speciesMap(speciesName);

% Avoid - O(n) lookup in hot paths
index = find(strcmp(speciesList, speciesName), 1);
```

### 3. Minimize Repeated Computations
```matlab
% Good
temp = expensive_computation();
result1 = temp * factor1;
result2 = temp * factor2;

% Avoid
result1 = expensive_computation() * factor1;
result2 = expensive_computation() * factor2;
```

### 4. Vectorize When Possible
```matlab
% Good
result = sin(x) .* cos(y);

% Avoid (when x, y are vectors)
for i = 1:length(x)
    result(i) = sin(x(i)) * cos(y(i));
end
```

### 5. Profile Before Optimizing
```matlab
% Use MATLAB profiler to identify bottlenecks
profile on
your_function();
profile viewer
```

## Benchmarking

To compare performance before and after optimizations:

```matlab
% Create benchmark
benchmark = combustiontoolbox.utils.Benchmark();

% Run test case
tic;
results = solve_equilibrium(parameters);
benchmark.add_result('equilibrium_solver', toc);

% Compare with baseline
benchmark.compare('baseline_results.mat');
```

## Conclusion

The optimizations implemented in this PR focus on high-impact, low-risk improvements:
- Replacing O(n) operations with O(1) lookups in hot paths
- Eliminating redundant calculations
- Proper memory preallocation

These changes maintain code correctness while improving performance, particularly for large-scale parametric studies and repeated equilibrium calculations.

## References

- MATLAB Performance Documentation: https://www.mathworks.com/help/matlab/performance.html
- Gordon & McBride (1994): NASA Reference Publication 1311
- Cuadra et al. (2024): Combustion Toolbox arXiv:2409.15086
