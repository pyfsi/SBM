import tracemalloc as tm
from utils import os
from utils import SBM_OUTPUT

def write_memory_profiler_result(stats_comparison):
    # Save the comparison results to a text file
    mem_prof_out_path = os.path.join(SBM_OUTPUT, 'memory_profiler.log')
    with open(mem_prof_out_path, 'a') as mem_f:
        mem_f.write("[ Memory usage increase from snapshot 1 to snapshot 2 ]\n")
        for stat in stats_comparison[:10]:
            mem_f.write(f"{stat}\n")

        # Detailed traceback for the top memory consumers
        mem_f.write("\n[ Detailed traceback for the top memory consumers ]\n")
        for stat in stats_comparison[:1]:
            mem_f.write('\n'.join(stat.traceback.format()) + '\n')

def mem_profile(do_mem_prof=False):
    def mem_profile_inner(func):
        def wrapper(*args):
            snapshot1 = tm.take_snapshot()
            func(*args)
            snapshot2 = tm.take_snapshot()
            stats1_vs_2 = snapshot2.compare_to(snapshot1, 'lineno')
            write_memory_profiler_result(stats1_vs_2)

        # use wrapper if in memory profile mode
        if do_mem_prof:
            return wrapper
        else:
            return func

    return mem_profile_inner



