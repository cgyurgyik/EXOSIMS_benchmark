# EXOSIMS_benchmark
Initial benchmarks on EXOSIMS C/C++ code using Google Benchmark.

### KeplerSTM
On my not-so-fast Macbook Air:
```
2020-01-26 18:57:54
Run on (4 X 1600 MHz CPU s)
CPU Caches:
  L1 Data 32 KiB (x2)
  L1 Instruction 32 KiB (x2)
  L2 Unified 256 KiB (x2)
  L3 Unified 3072 KiB (x1)
Load Average: 4.15, 3.13, 3.06
***WARNING*** Library was built as DEBUG. Timings may be affected.
--------------------------------------------------------
Benchmark              Time             CPU   Iterations
--------------------------------------------------------
Old_KeplerSTM       9422 ns         8759 ns        74960
New_KeplerSTM       7612 ns         7298 ns        84318

```

