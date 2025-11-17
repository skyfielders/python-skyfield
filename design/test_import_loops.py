#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path

def main():
    my_path = Path(sys.argv[0])
    design_dir = my_path.parent
    project_dir = design_dir.parent
    skyfield_dir = project_dir / 'skyfield'

    for path in sorted(skyfield_dir.glob('*.py')):
        module_name = path.stem
        command = ['python3', '-c', 'import skyfield.' + module_name]
        print(module_name)
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError:
            exit(1)

if __name__ == '__main__':
    main()
