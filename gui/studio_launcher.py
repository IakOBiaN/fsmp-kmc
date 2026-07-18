"""PyInstaller entry point for FSMP-kMC Studio (see studio.spec).

Besides starting the GUI, the frozen executable doubles as two windowless
tools; either flag must stay first in argv and skips the GUI entirely:

    --run-waiter <engine> <parameter file>
        the detached run launcher (once frozen there is no python to hand
        `-c` to, so fsmp_gui.runs.launch() re-enters this file);

    --selftest [report file]
        the bundle health check the release workflow runs after every
        build (see fsmp_gui.selftest)."""

import sys

if len(sys.argv) >= 4 and sys.argv[1] == "--run-waiter":
    from fsmp_gui.runs import waiter_main

    waiter_main(sys.argv[2], sys.argv[3])
elif len(sys.argv) >= 2 and sys.argv[1] == "--selftest":
    from fsmp_gui.selftest import run

    sys.exit(run(sys.argv[2] if len(sys.argv) > 2 else None))
else:
    from fsmp_gui.app import main

    main()
