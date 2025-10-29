import logging


def init_logging(debug=False):
    if debug:
        logging.basicConfig(format="%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s", datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    else:
        logging.basicConfig(format="%(asctime)s,%(msecs)03d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s", datefmt="%Y-%m-%d:%H:%M:%S", level=logging.ERROR)


def get_logger(name):
    return logging.getLogger(name)
