import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

class Logger:
    """A class to set up and manage logging configurations."""

    def __init__(
        self, 
        module_name: str = "default_logger", 
        log_folder: Path = Path("./dipper2/logs"), 
        verbosity: bool = False
    ):
        """
        Initializes the Logger instance and sets up logging.

        Args:
            module_name (str): Name of the module using the logger. Default is 'default_logger'.
            log_folder (Path): Path to the folder where logs should be stored. Default is './logs'.
            verbosity (bool): Whether to enable verbose logging. Default is False.
        """
        self.module_name = module_name
        self.log_folder = log_folder
        self.verbosity = verbosity

        # Ensure the log folder exists
        self.log_folder.mkdir(parents=True, exist_ok=True)

        # Set up the logger
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        """
        Configures the logger.

        Returns:
            logging.Logger: Configured logger instance.
        """
        logger = logging.getLogger(self.module_name)
        logger.setLevel(logging.DEBUG)  # Base logger level

        # Formatter for log messages
        formatter = logging.Formatter(
            "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
        )

        # File handler for saving logs to a file
        file_handler = RotatingFileHandler(
            filename=self.log_folder / "application.log",
            maxBytes=10 * 1024 * 1024,  # 10 MB
            backupCount=8,
        )
        file_handler.setLevel(logging.DEBUG if self.verbosity else logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # Stream handler for console output
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.DEBUG if self.verbosity else logging.WARNING)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

        return logger

    def get_logger(self) -> logging.Logger:
        """
        Returns the configured logger instance.

        Returns:
            logging.Logger: Configured logger instance.
        """
        return self.logger

