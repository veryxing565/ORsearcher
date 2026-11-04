"""Short operational logs plus opt-in, retained diagnostic detail."""
import logging
from pathlib import Path


LOGO = r'''   ___  ____                      _
  / _ \|  _ \ ___  ___  __ _ _ __| |__   ___ _ __
 | | | | |_) / __|/ _ \/ _` | '__| '_ \ / _ \ '__|
 | |_| |  _ <\__ \  __/ (_| | |  | | | |  __/ |
  \___/|_| \_\___/\___|\__,_|_|  |_| |_|\___|_|'''


def log_start(logger, version, args, root, sample_count):
    logger.info('\n%s', LOGO)
    logger.info('ORsearcher %s | Preset %s | %d thread(s), %d shard job(s), %d local worker(s)',
                version, getattr(args, 'preset', 'G0'), args.threads, args.jobs, args.workers)
    logger.info('Output: %s | Samples: %d | Workflow: %s', root, sample_count, args.command)


class ConciseFormatter(logging.Formatter):
    def format(self, record):
        # Keep tracebacks in debug.log, not in the short user-facing log.
        info, text = record.exc_info, record.exc_text
        record.exc_info = None
        record.exc_text = None
        try:
            return super().format(record)
        finally:
            record.exc_info, record.exc_text = info, text


def configure(root, level='info', quiet=False):
    root = Path(root)
    logger = logging.getLogger()
    # Only remove handlers owned by this package (safe for embedded callers).
    for handler in list(logger.handlers):
        if getattr(handler, '_orsearcher_owned', False):
            logger.removeHandler(handler)
            handler.close()
    logger.setLevel(logging.DEBUG)
    formatter = ConciseFormatter('%(asctime)s %(levelname)s %(message)s', '%Y-%m-%d %H:%M:%S')
    handlers = [logging.StreamHandler(),
                logging.FileHandler(root / 'batch.log', encoding='utf-8'),
                logging.FileHandler(root / 'debug.log', encoding='utf-8')]
    for index, handler in enumerate(handlers):
        handler._orsearcher_owned = True
        handler.setLevel(logging.DEBUG if index == 2 else logging.INFO if index == 1
                         else logging.WARNING if quiet else getattr(logging, level.upper()))
        handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
                             if index == 2 else formatter)
        logger.addHandler(handler)
