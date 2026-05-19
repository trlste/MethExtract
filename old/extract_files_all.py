#!/usr/bin/env python3

import argparse
import os
import subprocess
from multiprocessing import Pool
from pathlib import Path


def process_single_beta_file(args_tuple):
    """
    Wrapper function to process a single beta file by calling the original script.
    
    Args:
        args_tuple: Tuple of (beta_file, original_script, indir, outdir)
    
    Returns:
        Tuple of (success, beta_file)
    """
    beta_file, original_script, indir, outdir = args_tuple
    
    beta_file = str(beta_file)
    basename = os.path.basename(beta_file)
    
    print(f"Starting processing: {basename}")
    
    # Build command to call original script
    cmd = [
        "python3",
        original_script,
        "--beta", beta_file,
        "--indir", indir,
        "--outdir", outdir
    ]
    
    try:
        # Run the original script
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Print output from the original script
        if result.stdout:
            print(result.stdout)
        
        print(f"Successfully processed: {basename}\n")
        return True, basename
        
    except subprocess.CalledProcessError as e:
        print(f"Error processing {basename}:")
        print(f"  Return code: {e.returncode}")
        if e.stdout:
            print(f"  Stdout: {e.stdout}")
        if e.stderr:
            print(f"  Stderr: {e.stderr}")
        print()
        return False, basename
    
    except Exception as e:
        print(f"✗ Unexpected error processing {basename}: {str(e)}\n")
        return False, basename


def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script to process multiple beta files in parallel using the original converter script"
    )

    parser.add_argument(
        "--beta-dir",
        required=True,
        help="Input directory containing beta files"
    )

    parser.add_argument(
        "--script",
        required=True,
        help="Path to the original beta-to-bigwig conversion script"
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="Input directory containing CpG.bed.gz and CpG.chrome.size"
    )

    parser.add_argument(
        "--outdir",
        default=f'{os.getcwd()}/roadmap',
        help="Output directory (default: ./roadmap)"
    )

    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of processes to use (default: number of CPU cores)"
    )

    parser.add_argument(
        "--pattern",
        default="*.beta",
        help="File pattern to match beta files (default: *.beta)"
    )

    args = parser.parse_args()

    # Validate script exists
    if not os.path.isfile(args.script):
        print(f"Error: Script not found at {args.script}")
        return 1

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Find all beta files
    beta_dir = Path(args.beta_dir)
    beta_files = sorted(beta_dir.glob(args.pattern))
    
    if not beta_files:
        print(f"No beta files found in {args.beta_dir} matching pattern '{args.pattern}'")
        return 1
    
    print(f"Found {len(beta_files)} beta file(s) to process")
    print(f"Original script: {args.script}")
    print(f"Output directory: {args.outdir}")
    print(f"Reference directory: {args.indir}")
    print(f"Using {args.processes or 'default number of'} processes")
    print()
    print("=" * 80)
    print()

    # Prepare arguments for each process
    task_args = [
        (beta_file, args.script, args.indir, args.outdir)
        for beta_file in beta_files
    ]

    # Process files using multiprocessing
    with Pool(processes=args.processes) as pool:
        results = pool.map(process_single_beta_file, task_args)
    
    # Summary
    print()
    print("=" * 80)
    print("Processing Summary:")
    print(f"Total files: {len(results)}")
    print(f"Successfully processed: {sum(1 for success, _ in results if success)}")
    print(f"Failed: {sum(1 for success, _ in results if not success)}")
    
    failed_files = [name for success, name in results if not success]
    if failed_files:
        print("\nFailed files:")
        for f in failed_files:
            print(f"  - {f}")
        return 1
    
    return 0


exit(main())
