import argparse
import subprocess
import sys
import os
from pathlib import Path


def run_command(cmd_list, **kwargs):
    try:
        subprocess.run(cmd_list, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        print(
            f"Error: Command '{' '.join(cmd_list)}' failed with exit code {e.returncode}"
        )
        sys.exit(e.returncode)


def main():
    parser = argparse.ArgumentParser(description="DiPPER2 Python wrapper")
    parser.add_argument("-f", "--fold", required=True, help="Results folder")
    parser.add_argument("-d", "--assem_f", required=True, help="Assemblies folder")
    parser.add_argument("-l", "--target", required=True, help="Target list")
    parser.add_argument(
        "-F", "--forward", required=True, type=str, help="The forward primer sequence"
    )
    parser.add_argument(
        "-R", "--reverse", type=str, required=True, help="The reverse primer sequence"
    )
    parser.add_argument("-o", "--out", help="Output prefix")
    parser.add_argument(
        "-c", "--del_concat", default=1, type=int, help="Delete concatenated files"
    )
    parser.add_argument("-r", "--ref", help="Reference for BED files")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.1.0')

    args = parser.parse_args()

    print("Welcome to DiPPER2 Python wrapper")

    # get folders
    script_dir = Path(__file__).resolve().parent
    parent_dir = script_dir.parent

    # check if environment exists, otherwise create
    if not Path("dipper2").exists():
        run_command(["python3", "-m", "venv", "dipper2"])

    run_command(["dipper2/bin/pip", "install", "--upgrade", "pip"])
    run_command(
        ["dipper2/bin/pip", "install", "-r", str(parent_dir / "requirements.txt")]
    )

    python_env = script_dir / "dipper2/bin/python3"

    # change to folder with assemblies
    os.chdir(args.assem_f)

    # run target_move_module_optimized
    run_command(
        [
            python_env,
            str(script_dir / "target_move_module_optimized.py"),
            "-t",
            args.target,
            "-f",
            args.fold,
        ]
    )

    # Test primers
    test_args = [
        python_env,
        str(script_dir / "Primer_Testing_module_standalone.py"),
        "-f",
        args.fold,
        "-F",
        args.forward,
        "-R",
        args.reverse,
    ]
    if args.out:
        test_args += ["-o", args.out]
    if args.ref:
        test_args += ["-r", args.ref]
    if args.del_concat == 0:
        test_args += ["-c", str(args.del_concat)]
    run_command(test_args)
    print("Finished in silico PCR and target definition")

    # Summarize results
    summarize_args = [
        python_env,
        str(script_dir / "Summarize_results_module_standalone.py"),
        "-f",
        args.fold,
    ]
    #run 
    run_command(summarize_args)
    print("Finished! All steps completed successfully.")


if __name__ == "__main__":
    main()
