#!/usr/bin/env python3
"""
Benchmark script comparing schubmult (C++ executable) vs schubmult_py (Python CLI).

Generates test cases with permutations in S_and measures execution time
for both implementations, displaying results in a comparison table.
"""

import subprocess
import time
import sys
from pathlib import Path

# Test permutation pairs (two perms separated by hyphen)
# Each is a tuple: (perm1, perm2, description)
TEST_CASES = [
    # S_cases
    ("3 8 1 9 7 5 6 4 2 10 11 14 13 12 - 6 1 5 3 2 4 7 8 9 10 11 14 13 12", "S_moderate"),
    ("10 9 8 5 2 6 3 7 1 4 11 13 12 - 6 11 2 9 8 4 7 1 3 5 10 13 12", "S_dense"),
    ("4 1 2 3 5 6 7 8 9 10 11 14 13 12 - 5 3 2 1 4 6 7 8 9 10 11 14 13 12", "S_pattern"),
    
    ("3 1 5 4 2 6 7 8 9 10 11 14 13 12 - 6 1 5 3 2 4 7 8 9 10 11 14 13 12", "S_extended"),
    ("7 2 4 1 3 5 6 8 9 10 11 14 13 12 - 2 5 1 3 7 4 6 8 9 10 11 14 13 12", "S_dense"),
    ("6 1 2 3 4 5 7 8 9 10 11 14 13 12 - 4 3 2 1 5 6 7 8 9 10 11 14 13 12", "S_pattern"),
    
    ("3 1 5 4 2 6 7 8 9 10 11 14 13 12 - 6 1 5 2 3 4 7 8 9 10 11 14 13 12", "S_extended"),
    ("8 2 4 1 3 5 6 7 9 10 11 14 13 12 - 2 5 1 3 7 4 6 8 9 10 11 14 13 12", "S_dense"),
    ("7 1 2 3 4 5 6 8 9 10 11 14 13 12 - 5 4 3 2 1 6 7 8 9 10 11 14 13 12", "S_pattern"),
    
    ("3 1 5 4 2 6 7 8 9 10 11 14 13 12 - 6 3 1 5 2 4 7 8 9 10 11 14 13 12", "S_extended"),
    ("9 2 4 1 3 5 6 7 8 10 11 14 13 12 - 2 5 1 3 7 4 6 8 9 10 11 14 13 12", "S_dense"),
    ("8 1 2 3 4 5 6 7 9 10 11 14 13 12 - 6 5 4 3 2 1 7 8 9 10 11 14 13 12", "S_pattern"),
]


def run_benchmark(executable: str, perm_input: str, num_runs: int = 1, use_conda: bool = False, conda_env: str = "schubmult_312") -> float:
    """
    Run executable with given input and return total time in seconds.
    
    Args:
        executable: Path/name of executable to run
        perm_input: Input string with permutations
        num_runs: Number of times to run (default 1)
        use_conda: Whether to run within conda environment
        conda_env: Name of conda environment to use
    
    Returns:
        Total execution time in seconds
    """
    total_time = 0.0
    import schubmult._scripts.schubmult_py as py
    for _ in range(num_runs):
        if use_conda:            
            cmd = [executable] + perm_input.split() + ["-np"]
        else:
            cmd = [executable] + perm_input.split()
        
        start = time.time()
        try:
            if use_conda:
                result = py.main(cmd) 
                result = 0 if result is not None else 1
            else:
                result = subprocess.run(
                    cmd,
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                    timeout=30,
                )
            elapsed = time.time() - start
            total_time += elapsed
            
            # if (use_conda and result != 0) or (not use_conda and result.returncode != 0):
            #     print(f"Error running {executable}: {result.stderr[:100]}", file=sys.stderr)
            #     return None
                
        except subprocess.TimeoutExpired:
            print(f"Timeout running {executable}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Exception running {executable}: {e}", file=sys.stderr)
            return None
    
    return total_time


def format_time(seconds: float) -> str:
    """Format time for display."""
    if seconds is None:
        return "ERROR"
    if seconds < 0.001:
        return f"{seconds*1e6:.1f} µs"
    elif seconds < 1:
        return f"{seconds*1000:.2f} ms"
    else:
        return f"{seconds:.3f} s"


def main():
    print("Benchmarking schubmult vs schubmult_py")
    print("=" * 90)
    
    # First pass: warmup and check if both work
    print("Running warmup... ", end="", flush=True)
    run_benchmark("schubmult", TEST_CASES[0][0])
    run_benchmark("schubmult_py", TEST_CASES[0][0], use_conda=True)
    print("done")
    print()
    
    # Table header
    print(f"{'Example':<40} {'schubmult':<15} {'schubmult_py':<15} {'Ratio':<8}")
    print("-" * 90)
    
    results = []
    
    for perm_input, description in TEST_CASES:
        # Time schubmult (C++)
        time_cpp = run_benchmark("schubmult", perm_input, num_runs=3)
        
        # Time schubmult_py (Python) - run in conda environment
        time_py = run_benchmark("schubmult_py", perm_input, num_runs=3, use_conda=True)
        
        if time_cpp is None or time_py is None:
            print(f"{description:<40} {'FAILED':<15} {'FAILED':<15} {'-':<8}")
            continue
        
        # Average times
        avg_cpp = time_cpp / 3
        avg_py = time_py / 3
        
        # Ratio (Python / C++)
        ratio = avg_cpp / avg_py if avg_py > 0 else float('inf')
        
        results.append((description, avg_cpp, avg_py, ratio))
        
        print(
            f"{description:<40} {format_time(avg_cpp):<15} {format_time(avg_py):<15} {ratio:>6.2f}x"
        )
    
    print("-" * 90)
    
    # Summary statistics
    if results:
        ratios = [r[3] for r in results]
        avg_ratio = sum(ratios) / len(ratios)
        min_ratio = min(ratios)
        max_ratio = max(ratios)
        
        total_cpp = sum(r[1] for r in results)
        total_py = sum(r[2] for r in results)
        
        print(f"{'Summary':<40} {format_time(total_cpp):<15} {format_time(total_py):<15} {avg_ratio:>6.2f}x")
        print()
        print(f"Average ratio: {avg_ratio:.2f}x (min: {min_ratio:.2f}x, max: {max_ratio:.2f}x)")
        print(f"Total C++ time: {format_time(total_cpp)}")
        print(f"Total Python time: {format_time(total_py)}")


if __name__ == "__main__":
    main()


# import subprocess
# import timeit

# # Define your executables and the arguments you want to test
# # Replace "programA" and "scriptB.py" with your actual targets
# PROGRAM_A = ["schubmult"]
# SCRIPT_B = ["schubmult_py"]

# # The set of arguments you want to pass to both programs
# ARGUMENTS = ["--verbose", "--threads=4", "input_file.dat"]

# # Number of times to run each program to calculate the average
# NUMBER_OF_RUNS = 5


# def run_target(cmd_base, args):
#     """Combines the executable base with the arguments and runs it."""
#     full_command = cmd_base + args
#     # stdout/stderr are captured and discarded to prevent console flooding
#     # during benchmarking
#     subprocess.run(
#         full_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
#     )


# def main():
#     print(f"Benchmarking across {NUMBER_OF_RUNS} runs...")
#     print(f"Arguments being tested: {ARGUMENTS}\n")

#     # 1. Time Program A
#     # We use lambda to pass arguments into the timeit environment easily
#     timer_a = timeit.timeit(
#         lambda: run_target(PROGRAM_A, ARGUMENTS), number=NUMBER_OF_RUNS
#     )
#     avg_a = timer_a / NUMBER_OF_RUNS
#     print(f"🏁 Program A Total Time: {timer_a:.4f} seconds")
#     print(f"   Program A Avg Time:   {avg_a:.4f} seconds/run\n")

#     # 2. Time Script B
#     timer_b = timeit.timeit(
#         lambda: run_target(SCRIPT_B, ARGUMENTS), number=NUMBER_OF_RUNS
#     )
#     avg_b = timer_b / NUMBER_OF_RUNS
#     print(f"🏁 Script B Total Time: {timer_b:.4f} seconds")
#     print(f"   Script B Avg Time:   {avg_b:.4f} seconds/run\n")

#     # 3. Quick Comparison Summary
#     if avg_a < avg_b:
#         diff = avg_b / avg_a
#         print(f"🚀 Program A is roughly {diff:.1f}x faster than Script B.")
#     else:
#         diff = avg_a / avg_b
#         print(f"🚀 Script B is roughly {diff:.1f}x faster than Program A.")


# if __name__ == "__main__":
#     main()
