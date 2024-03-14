

# Submodular Subset Maximization

This repository contains code to exactly solve the **Cardinality-Constrained Submodular Monotone Subset Maximization** problem.

Given a universe $\mathcal{U}$ consisting of $n$ arbitrary items. Let $f : 2^{\mathcal{U}} \to \mathbb{R}$ be a set function.    
Let $\Delta(e \mid A) \coloneqq f(A \cup \{e\}) - f(A)$ be the **marginal gain** of $e$.        
Set function $f$ is **submodular** if for every $A \subseteq B \subseteq \mathcal{U}$ and $e \in \mathcal{U} \setminus B$ it holds that $\Delta(e \mid A) \geq \Delta(e \mid B)$.      
Set function $f$ is **monotone (increasing)** if for every $A \subseteq \mathcal{U}$ and every $e \in \mathcal{U}$ it holds that $\Delta(e \mid A) \geq 0$.

#### Cardinality-Constrained Submodular Monotone Subset Maximization

Given a submodular and monotone set function $f$ and an integer $k$ determine a set $S$ with $|S| = k$ and      
$$S = \text{arg}\max_{S' \subseteq \mathcal{U}, |S'| = k} f(S').$$

## Available Functions

Currently, four functions can be optimized.

- Group Farness : Given a graph $G=(V, E)$ determine a set $S \subseteq V$ of size $k$, such that $$f(S) \coloneqq \sum_{v \in V \setminus S} \min_{u \in S} \text{dist}(u, v)$$ is minimal. Since this originally is a minimization problem, we will use Negative Group Farness $-f(S)$ as the function to maximize.

- Partial Dominating Set: Given a graph $G=(V, E)$ determine a set $S \subseteq V$ of size $k$, such that $$f(S) \coloneqq |\bigcup_{v \in S} N[v]|$$ is maximized. $N[v] = \{v\} \cup \{u \mid \{u, v\} \in E\}$ is the closed neighborhood of $v$.

- $k$-Medoid Clustering: Given a set $X$ consisting of $n$ datapoints, each of dimensionality $d$, determine a set $S \subseteq X$ of size $k$ such that $$f(S) \coloneqq \sum_{i=1}^{n} \min_{x_j \in S} d(x_i, x_j)$$ is minimized. $$d(x, y) = \sqrt{ \sum_{i = 1}^{d} (x_i - y_i)^{2} }$$ is the euclidian distance. Like Group Farness, we will maximize $-f(S)$.

- Facility Location: Given a set of $n$ locations $N$ and a set of $m$ customers $M$. By $g_{ij} \geq 0$ we denote the benefit for customer $j$ when building a facility at location $i$. Determine a set $S \subseteq N$ of size $k$, such that $$f(S) = \sum_{j \in M} \max_{i \in S} g_{ij}$$ is maximized.

## Installation

### Requirements

The [Boost](https://www.boost.org/) library version 1.74.0 or greater is required.    
The Minimum CMake version required is 3.24.

Komogorv's [Blossom V](https://pub.ista.ac.at/~vnk/software.html#BLOSSOM5)

> Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum
> cost perfect matching algorithm." In Mathematical Programming
> Computation (MPC), July 2009, 1(1):43-67.

implementation is required. Its redistribution is prohibited, so we give a small description on how to include it in the project.
Extract the folder into `3rd_party_tools`. The folder structure should look like

```
SSM/
├─ 3rd_party_tools/
│  ├─ blossom5-v2.05.src/
│  │  ├─ GEOM/
│  │  ├─ MinCost/
│  │  ├─ block.h
│  │  ├─ ...
│  │  ├─ PerfectMatching.h
│  │  ├─ ...
│  │  ├─ USAGE.TXT
├─ ...

```

Afterwads modify `PerfectMatching.h` on line 40 from `//#define PERFECT_MATCHING_DOUBLE` to `#define PERFECT_MATCHING_DOUBLE` (effectively uncommenting the line).

### Build
To build the binary use
```
cmake -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -S . -B ./build
cmake --build build --target subsetoptimization
```

The binary `subsetoptimization` will be in the folder `build`.

## Usage

Use `./subsetoptimization -h [ --help  ]` to print a help message.

#### Needed Arguments

- `-t [ --type  ] {graph | k-medoid | facility}`: To read given input.
- `-i [ --input-file  ] arg`: Path to the file containing the graph, datapoints, etc.
- `-k [ --k  ] arg`: The desired size of the set.
- `-s [ --score-function ] arg`: The desired score function to optimize. Currently available are `negative-group-farness` or `partial-dominating-set` if `graph` was specified, `euclidian-distance` if `k-medoid` was specified and `benefits`, if `facility` was specified.
- `-o [ --output-file  ] arg`: Path to the file, which will hold the results.

Example: `./subsetoptimization -t graph -i input.edges -k 10 -s partial-dominating-set -o result.JSON` will optimize the Partial Dominating Set function on the graph in input.edges for the set size $k = 10$ and write the results to results.JSON.

#### Optional Arguments

- `--time-limit arg`: Time-limit in seconds (0 equals infinity). The algorithm will return the best-found solution if the time limit is reached.
- `--SUB {0, 1}`: Set to 1 to enable Simple-Upper-Bound (recommended), or 0 to disable it.
- `--CR {0, 1}`: Set to 1 to enable Candidate-Reduction (recommended), or 0 to disable it.
- `--Greedy-Local-Search {0, 1}`: Set to 1 to enable an additional greedy local search algorithm at initialization (not recommended), or 0 to disable it.
- `--BFThreshold n0 k0`: If the algorithm encounters a node in the SE Tree with less-equal than $n_0$ and less-equal $k_0$ elements still have to be chosen, the algorithm switches to Brute-Force.
- `--Plain {0, 1}`: Set to 1 to disable all heuristics.

#### Lazy Evaluation Arguments

`--LE Disabled` to disable Lazy Evaluation or  
`--LE Avg*y1 {or, and} {n*y2, k*y2, r_n*y2, r_k*y2}` for any $y_1, y_2 \geq 0$.

- The first argument initializes the **Average Update Scheme**
- The second argument determines, whether to use the **and** or **or Update Scheme**
- The third argument initializes the **Rank Update Scheme**.
- A recommended setting would be `--LE Avg*1 and n*1`.

#### Upper Bound 2D Arguments

`--UB2D Disabled` to disable Upper Bound 2D or  
`--UB2D {Sqrt[x]*y, x*y, y} l_max {Greedy, Matching, BForce, Dynamic} {A, B, AB} {0, 1} {0, 1} lazyStart lazyAdd lowDepthPercentage highDepthPercentage subBoundPercentage` for $x \in \{n, k, r\_n, r\_k\}$, $y \geq 0$ (float), $l_{\max} \geq 0$ (int), $0 \leq \text{lazyStart} \leq 1$ (float), $0 \leq \text{lazyAdd} \leq 1$ (float), $0 \leq \text{lowDepthPercentage} \leq \text{highDepthPercentage} \leq 1$ (float) and $0 \leq \text{subBoundPercentage} \leq 1$ (float).

- The 1st argument determines how to caculate $\ell$.
- The 2nd argument determines the maximum value for $\ell$.
- The 3rd argument determines which algorithm computes the matching.
- The 4th argument determines which Odd-Strategy to use.
- The 5th argument determines if Candidate Reduction should be used with the heuristic.
- The 6th argument determines if Safe-Skips should be used.
- The 7th and 8th argument determines how to initialize Lazy-Skips.
- The 9th and 10th argument determines at which depth the heuristic should be activated.
- The 11th argument determines the percentage to start the heuristic, based on SUB dynamically.
- In general, it is recommended to disable the setting.

#### Partial Brute Force Arguments

`--PBF Disabled` to disable Partial Brute Force or  
`--PBF {Sqrt[x1]*y1, x1*y1, y1} n_max {Sqrt[x2]*y2, x2*y2, y2} l_max {BForce, Dynamic} {0, 1} {0, 1} lazyStart lazyAdd lowDepthPercentage highDepthPercentage subBoundPercentage` for $x_1, x_2 \in \{n, k, r\_n, r\_k\}$, $y_1, y_2 \geq 0$ (float), $n_{\max} \geq 0$ (int), $l_{\max} \geq 0$ (int), $0 \leq \text{lazyStart} \leq 1$ (float), $0 \leq \text{lazyAdd} \leq 1$ (float), $0 \leq \text{lowDepthPercentage} \leq \text{highDepthPercentage} \leq 1$ (float) and $0 \leq \text{subBoundPercentage} \leq 1$ (float).

- The 1st argument determines how to calculate $n$ ($\eta$ in thesis).
- The 2nd argument determines a maximum value for $n$.
- The 3rd argument determines how to calculate $\ell$ ($\lambda$ in thesis).
- The 4th argument determines a maximum value for $\ell$.
- The 5th argument determines which algorithm to use.
- The arguments 6 to 12 are identical to arguments 5 to 11 from Upper Bound 2D.
- In general, it is recommended to disable the setting.

#### Divide-And-Conquer Arguments

`--DAC Disabled` to disable it or  
`--DAC {Sqrt[x]*y, x*y, y} {l_max} {0, 1} {0, 1} lazyStart lazyAdd lowDepthPercentage highDepthPercentage subBoundPercentage` for $x \in \{n, k, r\_n, r _k\}$, $y \geq 0$ (float), $\ell_{\max} \geq 0$ (int), $0 \leq \text{lazyStart} \leq 1$ (float), $0 \leq \text{lazyAdd} \leq 1$ (float), $0 \leq \text{lowDepthPercentage} \leq \text{highDepthPercentage} \leq 1$ (float) and $0 \leq \text{subBoundPercentage} \leq 1$ (float).

- The 1st argument determines how to calculate $\ell$ (size of the first partition).
- The 2nd argument detmerines a maximum value for $\ell$.
- The arguments 3 to 9 are identical to arguments 5 to 11 from Upper Bound 2D.
- In general, it is recommended to disable the setting.

#### Developer Arguments

- `--measure-oracle arg`: Set to 1, if we want to measure the oracle timings for each heuristic.

## File Format

#### Graph

The file format should be

``` e11 e12 e21 e22 ... en1 en2 ```  
with `ei1 ei2` denoting edge $i$ of the graph.

- Each $e_{ij}$ should be a positive integer $\geq 0$.
- The algorithm assumes that the vertex with the smallest ID is 0.
- Lines starting with a `%` are comments and will be ignored.

#### Datapoints

The file format should be

```  
x11 x12 ... x1d  
x21 x22 ... x2d  
...  
xn1 xn2 ... xnd  
```  

each line denotes one point of the dataset.

- Each $x_{ij}$ should be a double value.
- Lines starting with a `%` are comments and will be ignored.

#### Facility Location

The file format should be

```  
g11 g12 ... g1m  
g21 g22 ... g2m  
...  
gn1 gn2 ... gnm  
```  

each line denotes the benefit a location provides to the customers.

- Each $g_{ij}$ should be a double value $\geq 0$.
- Currently, the format does not allow comments.

## License

Apache License 2.0

## Project status

The project development has finished, and there are no plans to continue working on this repository.
