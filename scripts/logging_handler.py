import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

class LoggingHandler:

    def __init__(self):
        # create a logger 
        self.logger = logging.getLogger(__name__)
        
    def setup_logging(self, main_folder:Path, verbosity:bool):
        '''Set up the logger for the script. Currently global.'''

        # where to put the log files
        log_folder= main_folder / "logs"
        log_folder.mkdir(parents=True, exist_ok=True)

        self.logger.setLevel(logging.DEBUG)

        file_handler = RotatingFileHandler(
            filename = log_folder / "dipper2_run.log",
            maxBytes=10487650,
            backupCount=8
        )

        # set level of file handler (the actual log file)
        file_handler.setLevel(logging.DEBUG if verbosity else logging.INFO)
        self.logger.addHandler(file_handler)

        # stream handler is what user sees on command line
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.DEBUG if verbosity else logging.WARNING)
        self.logger.addHandler(stream_handler)

    
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

