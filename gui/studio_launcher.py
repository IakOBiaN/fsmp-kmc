"""PyInstaller entry point for FSMP-kMC Studio (see studio.spec).

Besides starting the GUI, the frozen executable doubles as the detached
run launcher: once frozen there is no python to hand `-c` to, so
fsmp_gui.runs.launch() re-enters this file with
    --run-waiter <engine> <parameter file>
which must stay first in argv and skip the GUI entirely."""

import sys

if len(sys.argv) >= 4 and sys.argv[1] == "--run-waiter":
    from fsmp_gui.runs import waiter_main

    waiter_main(sys.argv[2], sys.argv[3])
else:
    from fsmp_gui.app import main

    main()
