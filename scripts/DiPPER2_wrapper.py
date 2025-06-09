import argparse
import subprocess
import sys
import os
from pathlib import Path


def run_command(cmd_list, **kwargs):
    try:
        subprocess.run(cmd_list, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command '{' '.join(cmd_list)}' failed with exit code {e.returncode}")
        sys.exit(e.returncode)



def main():
    parser = argparse.ArgumentParser(description="DiPPER2 Python wrapper")
    parser.add_argument("-f", "--fold", required=True, help="Results folder")
    parser.add_argument("-d", "--assem_f", required=True, help="Assemblies folder")
    parser.add_argument("-l", "--target", required=True, help="Target list")
    parser.add_argument("-o", "--out", help="Output prefix")
    parser.add_argument("-p", "--fur", help="FUR parameters")
    parser.add_argument("-t", "--p3", help="Primer3 parameters; default:  primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=100 prodMaxSize=200 Oligo=1")
    parser.add_argument("-q", "--qpcr", default="n", help="qPCR (y) or conventional PCR (n)")
    parser.add_argument("-c", "--del_concat", default=1, type=int, help="Delete concatenated files")
    parser.add_argument("-r", "--ref", help="Reference for BED files")
    parser.add_argument("-a", "--ref_fur", help="Assembly used as FUR reference")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.1.0')

    args = parser.parse_args()

    print("Welcome to DiPPER2 Python wrapper")

    #get folders
    script_dir = Path(__file__).resolve().parent
    parent_dir = script_dir.parent

    #check if environment exists, otherwise create
    if not Path("dipper2").exists():
        run_command(["python3", "-m", "venv", "dipper2"])

    run_command(["dipper2/bin/pip", "install", "--upgrade", "pip"])
    run_command(["dipper2/bin/pip", "install", "-r", str(parent_dir / "requirements.txt")])
    
    python_env= script_dir / "dipper2/bin/python3"

    #change to folder with assemblies
    os.chdir(args.assem_f)

    #run target_move_module_optimized
    run_command([python_env, str(script_dir / "target_move_module_optimized.py"), "-t", args.target, "-f", args.fold])

    #run FUR
    fur_args = [python_env, str(script_dir / "FUR_module_optimized.py"), "-f", args.fold]
    if args.out:
        fur_args += ["-o", args.out]
    if args.fur:
        fur_args += ["-p", args.fur]
    if args.ref_fur:
        fur_args += ["-r", args.ref_fur]
    run_command(fur_args)

    print("Finished running FUR")

    #run Primer3
    p3_args = [python_env, str(script_dir / "Primer3_module_optimized.py"), "-f", args.fold]
    if args.out:
        p3_args += ["-o", args.out]

    if args.p3:
        p3_args += ["-p", args.p3]
    elif args.qpcr == "n":
        conv_pcr = "primMinTm=58 primOptTm=60 primMaxTm=62 inMinTm=63 inOptTm=65 inMaxTm=67 prodMinSize=200 prodMaxSize=1000 Oligo=0"
        p3_args += ["-p", conv_pcr]
    
    if args.qpcr == "y":
        p3_args += ["-q", args.qpcr]

    run_command(p3_args)
    print("Finished primer picking")

    #Test primers
    test_args = [python_env, str(script_dir / "Primer_Testing_module_optimized.py"), "-f", args.fold]
    if args.out:
        test_args += ["-o", args.out]
    if args.ref:
        test_args += ["-r", args.ref]
    if args.del_concat == 0:
        test_args += ["-c", str(args.del_concat)]
    run_command(test_args)
    print("Finished in silico PCR and target definition")

    #Summarize args
    summarize_args = [python_env, str(script_dir / "Summarize_results_module_improved.py"), "-f", args.fold]
    if args.qpcr == "y":
        summarize_args += ["-q", args.qpcr]

    run_command(summarize_args)
    print("Finished! All steps completed successfully.")


if __name__ == "__main__":
    main()
