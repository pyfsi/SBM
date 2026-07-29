from utils import os, shutil
import pstats

def check_stats_existence(stats_path) -> bool:
    return os.path.exists(stats_path)

def print_prof_stats(prof_path, stats_path):
    out_path = os.path.join(stats_path, "stats_cumulative.txt")
    with open(out_path, "w") as f:
        ps = pstats.Stats(prof_path, stream=f)
        ps.strip_dirs()
        ps.sort_stats('cumulative')
        ps.print_stats(r'\.py')

    out_path = os.path.join(stats_path, "stats_time.txt")
    with open(out_path, "w") as f:
        ps = pstats.Stats(prof_path, stream=f)
        ps.strip_dirs()
        ps.sort_stats('time')
        ps.print_stats(r'\.py')

def print_sbm_callers(prof_path, stats_path):
    out_path = os.path.join(stats_path, "callers.txt")
    with open(out_path, "w") as f:
        ps = pstats.Stats(prof_path, stream=f)
        ps.strip_dirs()
        ps.sort_stats('cumulative')
        ps.print_callers()

if __name__ == "__main__":
    # paths
    cwd = os.getcwd()
    stats_sbm_path = os.path.join(cwd, "sbm_files", "stats_sbm")
    prof_path = os.path.join(cwd, 'sbm.prof')

    if (os.path.exists(stats_sbm_path)):
        print(f"File path {stats_sbm_path} already exists.")
        user_input = input("Do you want to overwrite this directory? [y/n]")
        if (user_input.lower() == "y"):
            shutil.rmtree(stats_sbm_path)
            os.mkdir(stats_sbm_path)
            print(f"Previous statistics files deleted.")
        else:
            raise UserWarning(f"Previous statistics files not deleted. Ending print_profiler_stats.py.")
    else:
        os.mkdir(stats_sbm_path)

    # print stats
    print_prof_stats(prof_path, stats_sbm_path)
    print_sbm_callers(prof_path, stats_sbm_path)
    print(f"Profile statistics printed to {stats_sbm_path}.")