# Circuit file

Each circuit file is organized as follows:

```
n <number of nodes>
m <number of branches>
<branch 1>
<branch 2>
...
<branch m>
```

where each branch is represented as follows:

```
<start node id> <end node id> <current source in A> <resistance in Ohms> <voltage in V>
```

and node IDs are within the range of `[0, n - 1]`.
If a branch is missing one of the components, simply put a `0`.

*Note: You need to make sure of your circuit's polarity.*
