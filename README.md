# heihgm

heihgm is a program to solve matching  in hypergraphs  in (semi-)streaming fashion.

## Compiling & Running one of computations

We use bazel (https://bazel.build) as build system. For one-off computations run:

```
bazel build -c opt app
```
Running
```
./bazel-bin/app/app --command_textproto 'command:"run" hypergraph {    file_path: "path/to/hgr"    format: "hgr" } config { algorithm_configs{ data_structure:"<data_structure>" algorithm_name:"<algorithm_name>" }   }' --seed 1234
```
Valid `<data_structure>` values are `from_mem_stream_hypergraph` (loading to memory and stream then),`from_disk_stream_hypergraph` (streaming from disk) and `simple_matching` (only `greedy_slim`).

Algorithms for streaming (replace `<algorithm_name>`):
- `naive` Naive baseline algorithm
- `greedy_set`  Swap Set algorithm. Params: `double_params{key:"alpha" value:0.5}`
- `best_evict` SwapSet algorithm with `alpha=best`, needs a second pass.
- `greedy` Guarantee stack based algorithm Params: `double_params{key:"eps" value:0.5}`
- `lenient`  Lenient update function. Params: `double_params{key:"eps" value:0.5}`

Params need to be put in the `algorithm_configs{}` field
### Gurobi Compile

If you want to use Gurobi
```
export GUROBI_HOME=/path/to/gurobi1100/linux64
bazel build -c opt app --define gurobi=enabled
```
Running:

```
./bazel-bin/app/app --command_textproto 'command:"run" hypergraph {    file_path: "path/to/hgr"    format: "hgr" } config { algorithm_configs{ data_structure:"simple_matching" algorithm_name:"gurobi_exact" double_params{key:"timeout" value: 3600}}   }' --seed 1234
```

## Running experiments

### Compile

```
bazel build -c opt runner:fork_runner # --define gurobi=enabled # optional for gurobi
```

### Running

```
./bazel-bin/runner/fork_runner --experiment_path=/path/to/folder/containing/experiment  --random_order 1 --max_alloc_memory <max_mem_to_consume>  --max_alloc_memory_per_process <max_memory_for_one_process>  --cycles_to_queue_new 1 --seed 1234 --root_path=/path/to/folder/containing/hypergraphs
```
Please choose sensible values for <max_mem_to_consume> and <max_memory_for_one_process> in MB. So something like 50000 and 100000.

### Plotting
Please install python3.12 via spack (install via spack.io):
```
cd code
spack env activate .
spack install
```

```
spack env activate .
bazel run -c opt tools/plot:plot_cc /absolute/path/to/results/folder /absolute/path/to/visualisation.textproto 
```
