import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

def setup_logging(module_name:str, main_folder:Path, verbosity:bool):
    '''Set up the logger for the script. Currently global.'''
    logger = logging.getLogger(module_name)
    # where to put the log files
    log_folder= main_folder / "logs"
    log_folder.mkdir(parents=True, exist_ok=True)

    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(name)s - %(message)s")

    file_handler = RotatingFileHandler(
        filename = log_folder / "dipper2_run.log",
        maxBytes=10487650,
        backupCount=8
    )

    # set level of file handler (the actual log file)
    file_handler.setLevel(logging.DEBUG if verbosity else logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # stream handler is what user sees on command line
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if verbosity else logging.WARNING)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger

    #Creating a handler
# def handle_unhandled_exception(exc_type, exc_value, exc_traceback):
#     if issubclass(exc_type, KeyboardInterrupt):
#                 #Will call default excepthook
#         sys.__excepthook__(exc_type, exc_value, exc_traceback)
#         return
#         #Create a critical level log message with info from the except hook.
#     self.logger.critical("Unhandled exception", exc_info=(exc_type, exc_value, exc_traceback))
#Assign the excepthook to the handler
# sys.excepthook = handle_unhandled_exception

