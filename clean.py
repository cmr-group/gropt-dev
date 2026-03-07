import shutil
from pathlib import Path

# --- Configuration ---
# Get the project's root directory (where this script is located)
# This makes the script runnable from anywhere, as long as it's found by python.
project_root = Path(__file__).parent.resolve()

# List of directories to remove completely.
# '.scikit-build' and 'build' contain the CMake cache — only remove for a full rebuild.
dirs_to_remove = [
    'dist',
]

# Removing these forces a full CMake recompile of all C++ files.
cmake_dirs_to_remove = [
    '.scikit-build',
    'build',
]

# Glob pattern for other directories to remove.
glob_patterns_to_remove = [
    '*.egg-info',
    '*.so',  # Add this if you want to clean up stray .so files in the root
    '*.pyd',  # Add this if you want to clean up stray .pyd files in the root
    '**/gropt_wrapper*.pyd',  # Remove all compiled wrapper binaries to ensure only one exists after rebuild
    '**/gropt_wrapper*.so',
]

# Glob patterns for directories searched recursively under the project root.
recursive_glob_patterns_to_remove = [
    'gropt/__pycache__',
]


# --- Main Execution ---
def main():
    """Find and removes build artifacts in the project root."""
    import sys
    full = '--full' in sys.argv
    if full:
        print('--- Starting FULL cleanup script (CMake cache will be removed) ---')
    else:
        print('--- Starting cleanup script ---')

    all_dirs = dirs_to_remove + (cmake_dirs_to_remove if full else [])

    # Remove specified directories
    for dir_name in all_dirs:
        path = project_root / dir_name
        if path.is_dir():
            print(f'Removing directory: {path.relative_to(project_root)}')
            shutil.rmtree(path)
        elif path.exists():
            print(f'Removing file: {path.relative_to(project_root)}')
            path.unlink()

    # Remove directories/files matching glob patterns
    for pattern in glob_patterns_to_remove:
        for path in project_root.glob(pattern):
            if path.is_dir():
                print(f'Removing glob-matched directory: {path.relative_to(project_root)}')
                shutil.rmtree(path)
            elif path.is_file():
                print(f'Removing glob-matched file: {path.relative_to(project_root)}')
                path.unlink()

    # Remove specific recursive paths
    for pattern in recursive_glob_patterns_to_remove:
        path = project_root / pattern
        if path.is_dir():
            print(f'Removing directory: {path.relative_to(project_root)}')
            shutil.rmtree(path)

    print('--- Cleanup finished successfully ---')


if __name__ == '__main__':
    main()
