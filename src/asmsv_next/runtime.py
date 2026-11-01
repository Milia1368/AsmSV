"""Process lifetime helpers; no calling/scoring rules."""
import ctypes
import os
import signal
import sys
import time


def arm_parent_death(parent_pid):
    if sys.platform.startswith('linux'):
        libc=ctypes.CDLL(None,use_errno=True)
        if libc.prctl(1,signal.SIGTERM,0,0,0)!=0:raise OSError(ctypes.get_errno(),'PR_SET_PDEATHSIG failed')
        if os.getppid()!=parent_pid:raise SystemExit(143)


def worker_lifetime(parent_pid):
    signal.signal(signal.SIGTERM,signal.SIG_DFL)
    signal.signal(signal.SIGINT,signal.SIG_IGN)
    arm_parent_death(parent_pid)


def abort_pool(executor,futures):
    for future in futures:future.cancel()
    # Python 3.8–3.10 have no public terminate_workers(). Capture the owned
    # multiprocessing.Process handles before shutdown clears the private map.
    processes=list((getattr(executor,'_processes',None) or {}).values())
    for process in processes:
        if process.is_alive():process.terminate()
    deadline=time.monotonic()+2
    for process in processes:process.join(max(0,deadline-time.monotonic()))
    for process in processes:
        if process.is_alive():process.kill()
    for process in processes:process.join()
    executor.shutdown(wait=True)


def raise_termination(signum,frame):
    raise SystemExit(128+signum)
