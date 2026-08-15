## What this changes

<!-- One or two sentences: what a user can do after this that they could not
     before, or what stops going wrong. -->

## How it was checked

<!-- Which suites you ran and on what platform. If you could not run one,
     say so; CI will run all of them. -->

- [ ] `make test` (or `build test`) passes, and the pinned energies did not move
- [ ] `QT_QPA_PLATFORM=offscreen python -m unittest discover -s gui/tests` passes, if the Studio was touched
- [ ] The warning gate is clean: `-Wall -Wextra -Werror` under GCC and Clang
- [ ] Documentation is updated if behaviour changed (README, `configs/` comments, CHANGELOG)

## Anything to be aware of

<!-- A pin that moved and why, a new dependency, a format change, a decision
     you were unsure about. -->
