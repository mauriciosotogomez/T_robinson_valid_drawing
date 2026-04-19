# Caterpillar Metric

Tools to determine whether a given similarity/distance matrix admits a **caterpillar tree metric** (a.k.a. T-Robinsonian structure), and to compute the closest such metric via Linear Programming.

## Background

A symmetric matrix is *T-Robinsonian* if its entries are consistent with distances on a weighted caterpillar tree. This code formulates the problem as an LP (using [PuLP](https://coin-or.github.io/pulp/)) and finds either a valid embedding or certifies infeasibility.

## Files

| File | Description |
|------|-------------|
| `caterpillar_metric.py` | Main script: reads input, builds and solves the LP, checks and reports results |
| `utils.py` | Helper functions: center-matrix generation, Robinson matrix generation, solution checking, text/graph plotting |

## Requirements

- Python 3.x
- `numpy`, `pandas`, `matplotlib`, `networkx`
- `pulp` (LP solver wrapper; uses CBC by default)

Install dependencies:

```bash
pip install numpy pandas matplotlib networkx pulp
```

## Usage

```
python caterpillar_metric.py [options]
```

### Options

| Flag | Description |
|------|-------------|
| `-e N` | Use built-in example number N (0–5) |
| `-f FILE` | Read matrix from a CSV file |
| `-g SIZE` | Generate a random Robinson distance matrix of given size (requires R + `create_Robinson.r`) |
| `-m SIZE` | Generate a random center matrix of given size |
| `-d` | Treat input as a **distance** matrix (default: similarity) |
| `-c` | Treat input as a **center** matrix |
| `-p` | Restrict embedding to a **path** (no hanging legs) |
| `-v` | Verbose output: print variables, dual values, and draw the tree |

### Examples

Run on built-in example 0 (similarity matrix):

```bash
python caterpillar_metric.py -e 0
```

Run on a CSV distance matrix with verbose output:

```bash
python caterpillar_metric.py -f my_matrix.csv -d -v
```

Generate and test a random center matrix of size 6:

```bash
python caterpillar_metric.py -m 6
```

## Output

- **Console**: LP status, solution correctness (green = valid caterpillar metric, red = not), timing, and (with `-v`) the tree structure and variable values.
- **`caterpillar_metric_out.txt`**: The output distance matrix realised by the caterpillar tree.
- **`A.csv`, `B.csv`, `Bt.csv`, `K.csv`**: Intermediate constraint and kernel matrices (useful for analysis).
- **LP files** (`caterpillar_metric.lp`, `aux_caterpillar_metric.lp`): The LP instances in standard format.

## How it works

1. **Center matrix** (`max_closer`): For each pair (i, j), identifies the node on the backbone path closest to both — the "Robinson center".
2. **Auxiliary LP**: Solves an unconstrained (no non-negativity) version to extract the kernel structure of the feasible cone.
3. **Main LP**: Solves the non-negative version to find actual caterpillar distances `d[i]` (backbone positions) and leg lengths `l[i]`.
4. **Verification**: Checks that the output matrix preserves the same ordering as the input (correct Robinson structure).

## License

MIT
