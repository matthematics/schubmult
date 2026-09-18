"""Thin wrappers around the standard ``logging`` module."""

import logging


def init_logging(debug=False):
    """Configure root logging at DEBUG (if ``debug``) or ERROR with a timestamped ``file:line`` format."""
    if debug:
        logging.basicConfig(format="%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s", datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    else:
        logging.basicConfig(format="%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s", datefmt="%Y-%m-%d:%H:%M:%S", level=logging.ERROR)


def get_logger(name):
    """``logging.getLogger(name)``."""
    return logging.getLogger(name)
