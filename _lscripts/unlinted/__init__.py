# Minimal package for CLI scripts. Keep lightweight (no heavy imports) so package
# import during installation is cheap. Individual script modules should import
# heavy dependencies only inside main() so they don't run on import.
__all__ = []
